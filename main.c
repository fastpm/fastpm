/**
 * main file: Time stepping and firing-sequences
 * 
 * This code was initially modified by Jun Koda, 
 * from the original serial COLA code
 * by Svetlin Tassev.
 *
 * The time stepping code has been largely re-written
 * to accompodate PM simulations.
 *
 * The firing sequences of snapshots and FOF is preserved.
 * 
 * Adpative PM size is added by Man-Yat Chu.
 * 
 * Yu Feng <rainwoodman@gmail.com>
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "lpt.h"
#include "msg.h"
#include "power.h"
#include "comm.h"
#include "domain.h"
#include "pm.h"
#include "stepping.h"
#include "fof.h"
#include "write.h"
#include "timer.h"
#include "mem.h"
#include "subsample.h"
#include "coarse_grid.h"

extern int write_runpb_snapshot(Snapshot * snapshot,  
        char * filebase,
        void * scratch, size_t scratch_bytes);
extern int read_runpb_ic(Parameters * param, double a_init, Particles * particles, 
        void * scratch, size_t scratch_bytes);

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
           double a_x, double a_v,
		   Particles const * const particles, 
		   Snapshot * const snapshot, 
		   const char fof_filename[], 
		   const char subsample_filename[], 
		   const char cgrid_filename[], const int cgrid_nc,
		   void* const mem1, const size_t size1,
		   const int write_longid,
		   const double fof_linking_factor);

void write_slice(const char filebase[], Particle* p, const int np, const float dz, const float boxsize);
void write_snapshot_slice(const char filebase[], Snapshot const * const snapshot, const float boxsize);
static double polval(const double * pol, const int degrees, const double x) {
    double rt = 0;
    int i;
    for(i = 0; i < degrees; i ++) {
        rt += pow(x, degrees - i - 1) * pol[i];
    }
    return rt;
}

static double measured_power(double k, Parameters * param) {
    size_t nk;
    double * power = pm_compute_power_spectrum(&nk);
    double k0 = 2 * 3.1415926 / param->boxsize;
    int ki = floor(k / k0);
    double u = k / k0 - ki;

    if(ki >= nk / 2 - 1) return 0;
    if(ki < 0) return 0;
    return power[ki] * (1-u) + power[ki] * u; // other corrections outsize.
}

int main(int argc, char* argv[])
{
  const int multi_thread= mpi_init(&argc, &argv);
  msg_init();
                                                      timer_set_category(Init);

  //
  // Initialization / Memory allocation
  //						      
  Parameters param;
  read_parameters(argc, argv, &param);
  int nc_factor = -1;
  const float change_pm= param.change_pm;
  const double OmegaM= param.omega_m;
  const double OmegaLambda= 1.0 - OmegaM;
  const double Hubble= param.h;
  const double sigma8= param.sigma8;

  msg_set_loglevel(param.loglevel);

  fft_init(multi_thread);
  const int nsteps= param.ntimestep;
  const double a_final= param.a_final;

  double a_init = 0.1;

  /* one extra item in the end; to avoid an if conditioni in main loop */
  double * A_X = (double*) calloc(nsteps + 2, sizeof(double));
  double * A_V = (double*) calloc(nsteps + 2, sizeof(double));

  switch(param.time_step) {
      case TIME_STEP_LOGA:
          msg_printf(info, "timestep linear in loga\n");
          for(int i = 0; i < nsteps + 2; i ++) {
              double dloga = (log(a_final) - log(a_init)) / nsteps;
              A_X[i] = exp(log(a_init) + i * dloga);
              A_V[i] = exp(log(a_init) + i * dloga - 0.5 * dloga);
          }
      break;
      case TIME_STEP_A:
          msg_printf(info, "timestep linear in a\n");
          for(int i = 0; i < nsteps + 2; i ++) {
              double da = (a_final - a_init) / nsteps;
              A_X[i] = a_init + (i * da);
              A_V[i] = a_init + (i * da - 0.5 * da);
          }
      break;
      case TIME_STEP_GROWTH:
          msg_printf(info, "timestep linear in growth factor\n");
          for(int i = 0; i < nsteps + 2; i ++) {
              /* best fit coefficients for unnormalized wmap7 cosmology dplus
               * shoudn't matter as long as we are close */
              const double a2d[] = {0.35630564,-1.58170455, 0.34493492, 3.65626615, 0};
              const double d2a[] = {0.00657966,-0.01088151, 0.00861105, 0.27018028, 0};
              double dfinal = polval(a2d, 5, a_final);
              double dinit = polval(a2d, 5, a_init);
              double dd = (dfinal - dinit) / nsteps;
              A_X[i] = polval(d2a, 5, dinit + i * dd);
              A_V[i] = polval(d2a, 5, dinit + i * dd - 0.5 * dd);
          }
      break;
      default:
          abort();
    
  }
  /* fix up the IC time */
  A_V[0] = A_X[0] = a_init;
  A_X[nsteps] = a_final;
  msg_printf(normal, "Drift points: \n");
  for(int i = 0; i < nsteps + 2; i ++) {
      msg_printf(normal, "%g, ", A_X[i]);
  }
  msg_printf(normal, "\n");

  power_init(param.power_spectrum_filename, a_init, 
	     sigma8, OmegaM, OmegaLambda);
  

  Memory mem; 
  allocate_shared_memory(param.nc, param.pm_nc_factor2, param.np_alloc_factor, &mem); 
  lpt_init(param.nc, mem.mem1, mem.size1);
  const int local_nx= lpt_get_local_nx();

  Particles* particles= 
     allocate_particles(param.nc, local_nx, param.np_alloc_factor);
  Snapshot* snapshot= allocate_snapshot(param.nc, local_nx, 
				  particles->np_allocated, mem.mem2, mem.size2);
  snapshot->boxsize= param.boxsize;
  snapshot->omega_m= OmegaM;
  snapshot->h= Hubble;
  //strncpy(snapshot->filename, param.snapshot_filename, 64);
  snapshot->filename= param.snapshot_filename;
  
  nc_factor = param.pm_nc_factor1;
  pm_set_diff_order(param.diff_order);
  pm_init(nc_factor*param.nc, nc_factor, param.boxsize,
	  mem.mem1, mem.size1, mem.mem2, mem.size2, param.nrealization>1);
  comm_init(nc_factor*param.nc, param.nc, param.boxsize);

  fof_init(particles->np_allocated, param.nc, mem.mem1, mem.size1);
  subsample_init(param.subsample_factor, param.random_seed);

  const int nout= param.n_zout;
  double* aout= malloc(sizeof(double)*nout);
  for(int i=0; i<nout; i++) {
    aout[i] = (double)(1.0/(1 + param.zout[i]));
    msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
	       i, param.zout[i], aout[i]);
  }

  if (param.qpm) {
      stepping_set_subtract_lpt(0);
  } else {
      stepping_set_subtract_lpt(1);
  }
  if (param.stdda) {
      stepping_set_std_da(1);
  } else {
      stepping_set_std_da(0);
  }
  if (param.nopm) {
      stepping_set_no_pm(1);
  } else {
      stepping_set_no_pm(0);
  }
  //
  // Many realizations with different initial conditions
  //
  for(int irealization=0; irealization<param.nrealization; irealization++) {
      MPI_Barrier(MPI_COMM_WORLD);
      msg_printf(verbose, "\n%d/%d realization.\n", 
              irealization+1, param.nrealization);
      int seed= param.random_seed + irealization;

      int iout= 0;
      domain_init(param.nc * nc_factor, param.boxsize);

      timer_set_category(LPT);
      if(param.readic_filename) {
          read_runpb_ic(&param, a_init, particles, mem.mem2, mem.size2);
      } else {
          // Sets initial grid and 2LPT desplacement
          lpt_set_displacement(a_init, OmegaM, seed,
                  param.boxsize, particles);
      }
      snapshot->seed= seed;

      // always do this because it intializes the initial velocity
      // correctly.
      set_nonstepping_initial(a_init, particles, snapshot);

      if(param.init_filename) {
          // Writes initial condition to file for e.g. Gadget N-body simulation
          char filename[256]; // TODO: use variable length for filename
          sprintf(filename, "%s_%05d", param.init_filename, seed);
          write_snapshot(filename, snapshot, param.write_longid);
      }
      timer_set_category(STEPPING);
    
      //
      // Time evolution loop
      //
      //   TODO: allow nstep=1?
      if(nout > 0 && nsteps > 1 && a_final > a_init) {
          msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
                  a_init, a_final, nsteps);

          int chk_change = 0;
          for (int istep=0; istep<= nsteps; istep++) {
              double a_v, a_x, a_v1, a_x1;


              a_v = A_V[istep];
              a_x = A_X[istep];
              a_v1= A_V[istep + 1];
              a_x1= A_X[istep + 1];

              if(a_x1 >= change_pm && chk_change != 1 && param.pm_nc_factor2 != param.pm_nc_factor1){

                  msg_printf(normal, "Switching to new pm factor: %d->%d\n",
                        param.pm_nc_factor1,
                        param.pm_nc_factor2);

                  nc_factor = param.pm_nc_factor2;
                  pm_finalize();
                  comm_finalize();
                  domain_finalize();

                  domain_init(param.nc * nc_factor, param.boxsize);
                  comm_init(nc_factor*param.nc, param.nc, param.boxsize);//what are these?
                  pm_init(nc_factor*param.nc, nc_factor, param.boxsize, mem.mem1, mem.size1, mem.mem2, mem.size2, param.nrealization>1);
                  chk_change = 1;
              }

              msg_printf(normal, "Timestep %d/%d\n", istep + 1, nsteps);

              timer_start(comm);
              // move particles to other nodes
              domain_wrap(particles);
              domain_decompose(particles, mem.mem1, mem.size1);

              timer_stop(comm);

              pm_calculate_forces(particles, mem.mem2, mem.size2); 

              if(param.measure_power_spectrum_filename) {
                  size_t nk = 0;
                  double * powerspectrum = pm_compute_power_spectrum(&nk);
                  char fname[9999];
                  sprintf(fname, "%s-%0.4f.txt", param.measure_power_spectrum_filename, a_x);
                  int myrank;
                  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
                  if(myrank == 0) {
                      FILE * fp = fopen(fname, "w");
                      if (fp == NULL) {
                        msg_abort(0020, "Unable to write to %s\n", fname);
                      }
                      for(int i = 0; i < nk; i ++) {
                          fprintf(fp, "%g %g\n", 3.1416 * 2 / param.boxsize * i, 
                                  powerspectrum[i] / pow(nc_factor * param.nc, 6) * pow(param.boxsize, 3.0)
                            );
                      }
                      fclose(fp);
                  }
                  double sigma8 = TopHatSigma2(8.0, measured_power, &param);
                  sigma8 = sigma8 / pow(param.nc * nc_factor, 6.0) * pow(param.boxsize, 3.0);
                  sigma8 /= pow(3.1415926 * 2, 3);
                  sigma8 = sqrt(sigma8);
                  msg_printf(verbose, "sigma8 = %g expected = %g\n", sigma8, param.sigma8 * GrowthFactor(1.0, a_x));
              }

              while(iout < nout && a_v <= aout[iout] && aout[iout] <= a_x) {
                  // Time to write output
                  snapshot_time(aout[iout], iout, a_x, a_v, particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
                  iout++;
              }
              if(iout >= nout) break;

              // Leap-frog "kick" -- velocities updated
              stepping_kick(particles, OmegaM, a_v, a_v1, a_x);

              while(iout < nout && a_x < aout[iout] && aout[iout] <= a_v1) {
                  // Time to write output
                  snapshot_time(aout[iout], iout, a_x, a_v, particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
                  iout++;
              }
              if(iout >= nout) break;

              // Leap-frog "drift" -- positions updated
              stepping_drift(particles, OmegaM, a_x, a_x1, a_v1);
          }
      }
      timer_print();
  }

  //move_particles(particles);
  //write_slice("slice", particles->p, particles->np_local, param.boxsize/32, param.boxsize);
  

  MPI_Finalize();
  return 0;
}

int mpi_init(int* p_argc, char*** p_argv)
{
  // MPI+OpenMP paralellization: MPI_THREAD_FUNNELED
  // supported by mpich2 1.4.1, but now by openmpi 1.2.8

#ifdef _OPENMP
  int thread_level, hybrid_parallel;
  MPI_Init_thread(p_argc, p_argv, MPI_THREAD_FUNNELED, &thread_level);
  hybrid_parallel = (thread_level >= MPI_THREAD_FUNNELED);

  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank == 0) {
    if(hybrid_parallel)
      printf("MPI + multi thread supported (MPI_THREAD_FUNNELED).\n");
    else
      printf("Warning: MPI + multi thread not supported. 1 thread per node.\n");
  }
	
  return hybrid_parallel;
#else
  MPI_Init(p_argc, p_argv);
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0)
    printf("MPI only without OpenMP\n");
  return 0;
#endif

}

void fft_init(int threads_ok)
{
  // Initialize FFTW3

#ifdef _OPENMP
  if(threads_ok)
    threads_ok= fftwf_init_threads();
  if(!threads_ok)
    msg_printf(warn, "Multi-thread FFTW not supported.\n");
#endif
    
  fftwf_mpi_init();

#ifdef _OPENMP
  if(threads_ok) {
    int nthreads= omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    msg_printf(info, "Multi-threaded FFTW: %d threads\n", nthreads);
  }
#endif

}


void write_slice(const char filebase[], Particle* p, const int np, const float dz, const float boxsize)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char filename[128];
  sprintf(filename, "%s%d.txt", filebase, myrank);
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(0020, "Unable to write to %s\n", filename);

  
  for(int i=0; i<np; i++) {
    float x= fmod(p->x[0], boxsize);
    float y= fmod(p->x[1], boxsize);
    float z= fmod(p->x[2], boxsize);
    if(0.0f < z && z < dz) {
      fprintf(fp, "%e %e %e\n",  x, y, z);
    }
    p++;
  }
  fclose(fp);
}

void write_snapshot_slice(const char filebase[], Snapshot const * const snapshot, const float boxsize)
{
  msg_printf(normal, "Writing snapshot slices %s\n", filebase);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char filename[128];
  sprintf(filename, "%s_%d.txt", filebase, myrank);
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(0020, "Unable to write to %s\n", filename);

  ParticleMinimum const * p= snapshot->p;
  const int np= snapshot->np_local;
  
  for(int i=0; i<np; i++) {
    float x= fmod(p->x[0], boxsize);
    float y= fmod(p->x[1], boxsize);
    float z= fmod(p->x[2], boxsize);
    if(0.0f < z && z < boxsize/32) {
      fprintf(fp, "%e %e %e\n",  x, y, z);
    }
    p++;
  }
  fclose(fp);
}
    
void snapshot_time(const float aout, const int iout, 
           double a_x, double a_v,
		   Particles const * const particles, 
		   Snapshot * const snapshot,
		   const char fof_filename[], 
		   const char subsample_filename[], 
		   const char cgrid_filename[], const int cgrid_nc,
		   void* const mem1, const size_t size1,
		   const int write_longid,
		   const double fof_linking_factor)
{
  // Halo finding and snapshot outputs

  char filebase[1024];      // TODO: make 256 to variable number...?
  const int isnp= iout+1;
                                                       timer_set_category(Snp);
  stepping_set_snapshot(aout, a_x, a_v, particles, snapshot);

  const int nc= snapshot->nc; assert(nc > 0);
  const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);


                                                       timer_start(write);
  // Gadget snapshot for all particles
  // periodic wrapup not done, what about after fof? what about doing move_particle_min here?
  if(snapshot->filename) {
    //sprintf(filebase, "%s%05d%c", snapshot->filename, snapshot->seed, suffix);
    //write_snapshot(filebase, snapshot, write_longid);
    sprintf(filebase, "%s%05d_%0.04f", snapshot->filename, snapshot->seed, snapshot->a);
    write_runpb_snapshot(snapshot, filebase, mem1, size1);
  }
                                                       timer_stop(write);
						       timer_start(sub);
  // particle subsample (binary)
  if(subsample_filename) {
    sprintf(filebase, "%s%05d_%0.04f", subsample_filename, snapshot->seed, snapshot->a);
    //write_subsample(filebase, subsample_factor, snapshot, mem1, size1); // this is regular subsamle but not used. Random subampling is used.
    // TODO: periodic wrapup not done. What about after fof?
    write_random_sabsample(filebase, snapshot, mem1, size1);
  }

  // coarse mesh (binary)
  if(cgrid_filename) {
    sprintf(filebase, "%s%05d_%0.04f.b", cgrid_filename, snapshot->seed, snapshot->a);
    coarse_grid2(filebase, snapshot, cgrid_nc, mem1, size1);
  }

                                                       timer_stop(sub);
  // text file of a slice
  //sprintf(filebase, "slice%02d", isnp); temp
  //write_snapshot_slice(filebase, snapshot, boxsize); temp


  const float ll= (float)(fof_linking_factor*boxsize/nc); // FoF linking length
  fof_find_halos(snapshot, ll);

  // FoF halo catalogue (binary file)
  if(fof_filename) {
    // comment: move_particles done here
    /* the binary format is mass, x, x, x, v, v, v*/
    /* mass is in 1e10 Msun/h */
    const double rho_crit = 27.7455;
    const double M0 = snapshot->omega_m*rho_crit*pow(snapshot->boxsize / snapshot->nc, 3.0);
  
    sprintf(filebase, "%s%05d_%0.04f.bin", fof_filename, snapshot->seed, snapshot->a);
    fof_write_halos(filebase, 1, M0);
  }

  const double z_out= 1.0/aout - 1.0;
  msg_printf(normal, "snapshot %d written z=%4.2f a=%5.3f\n", 
	     isnp, z_out, aout);
                                                       timer_set_category(STEPPING);
}
