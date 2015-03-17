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
#include "pm.h"
#include "cola.h"
#include "fof.h"
#include "write.h"
#include "timer.h"
#include "mem.h"
#include "move.h"
#include "subsample.h"
#include "coarse_grid.h"

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
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
  const int nc_factor= param.pm_nc_factor;
  const double OmegaM= param.omega_m;
  const double OmegaLambda= 1.0 - OmegaM;
  const double Hubble= param.h;
  const double sigma8= param.sigma8;

  msg_set_loglevel(param.loglevel);

  fft_init(multi_thread);
  comm_init(param.pm_nc_factor*param.nc, param.nc, param.boxsize);

  const int nsteps= param.ntimestep;
  const double a_final= param.a_final;

  double a_init;
  if(param.loga_step) {
    a_init = 0.1;
  } else {
    double da = a_final / nsteps;
    a_init = da;
  }

  power_init(param.power_spectrum_filename, a_init, 
	     sigma8, OmegaM, OmegaLambda);


  Memory mem; 
  allocate_shared_memory(param.nc, nc_factor, param.np_alloc_factor, &mem); 
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
  

  pm_init(nc_factor*param.nc, nc_factor, param.boxsize,
	  mem.mem1, mem.size1, mem.mem2, mem.size2);
  fof_init(particles->np_allocated, param.nc, mem.mem1, mem.size1);
  subsample_init(param.subsample_factor, param.random_seed);

  const int nout= param.n_zout;
  float* aout= malloc(sizeof(float)*nout);
  for(int i=0; i<nout; i++) {
    aout[i] = (float)(1.0/(1 + param.zout[i]));
    msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
	       i, param.zout[i], aout[i]);
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

      // Sets initial grid and 2LPT desplacement
      timer_set_category(LPT);
      lpt_set_displacement(a_init, OmegaM, seed,
              param.boxsize, particles);
      snapshot->seed= seed;

      if(param.init_filename) {
          // Writes initial condition to file for e.g. Gadget N-body simulation
          char filename[256]; // TODO: use variable length for filename
          sprintf(filename, "%s_%05d", param.init_filename, seed);
          set_noncola_initial(particles, snapshot);
          write_snapshot(filename, snapshot, param.write_longid);
      }

      timer_set_category(COLA);

      if(param.loga_step) {
          const double dloga= (log(a_final) - log(a_init))/nsteps;
          particles->a_v= exp(log(a_init) - 0.5*dloga);
          msg_printf(info, "timestep linear in loga\n");
      } else {
          const double da = a_final / nsteps;
          //  particles->a_v= 0.5*da;
          //  I thought 0.5*da is the natural initial time for leap frog integration,
          //  but da gives much better matter power. I don't know why. (Feb 12, 2014)
          particles->a_v= da;
          msg_printf(info, "timestep linear in a\n");
      }

      //
      // Time evolution loop
      //
      //   TODO: allow nstep=1?
      if(nout > 0 && nsteps > 1 && a_final > a_init) {
          msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
                  a_init, a_final, nsteps);
          for (int istep=1; istep<=nsteps; istep++) {
              msg_printf(normal, "Timestep %d/%d\n", istep, nsteps);

              timer_start(comm);
              // move particles to other nodes
              move_particles2(particles, param.boxsize, mem.mem1, mem.size1 );

              timer_stop(comm);

              pm_calculate_forces(particles); 

              double avel0, apos0, avel1, apos1;
              if(param.loga_step) {
                  const double dloga= (log(a_final) - log(a_init))/nsteps;
                  avel0= exp(log(a_init) + (istep-0.5)*dloga);
                  apos0= exp(log(a_init) + istep*dloga);

                  avel1= exp(log(a_init) + (istep+0.5)*dloga);
                  apos1= exp(log(a_init) + (istep+1)*dloga);
              } else {
                  const double da = a_final / nsteps;
                  avel0= (istep-0.5)*da;
                  apos0=  istep*da;

                  avel1= (istep+0.5)*da;
                  apos1= (istep+1.0)*da;
              } 

              while(iout < nout && avel0 <= aout[iout] && aout[iout] <= apos0) {
                  // Time to write output
                  snapshot_time(aout[iout], iout, particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
                  iout++;
              }
              if(iout >= nout) break;

              // Leap-frog "kick" -- velocities updated
              cola_kick(particles, OmegaM, avel1);

              while(iout < nout && apos0 < aout[iout] && aout[iout] <= avel1) {
                  // Time to write output
                  snapshot_time(aout[iout], iout,  particles, snapshot, param.fof_filename, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid, param.fof_linking_factor);
                  iout++;
              }
              if(iout >= nout) break;

              // Leap-frog "drift" -- positions updated
              cola_drift(particles, OmegaM, apos1);
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

  char filebase[256];      // TODO: make 256 to variable number...?
  const int isnp= iout+1;
  char suffix= 'a' + iout; // TODO: not suited for large nout -- change to #
                                                       timer_set_category(Snp);
  cola_set_snapshot(aout, particles, snapshot);

  const int nc= snapshot->nc; assert(nc > 0);
  const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);


                                                       timer_start(write);
  // Gadget snapshot for all particles
  // periodic wrapup not done, what about after fof? what about doing move_particle_min here?
  if(snapshot->filename) {
    sprintf(filebase, "%s%05d%c", snapshot->filename, snapshot->seed, suffix);
    write_snapshot(filebase, snapshot, write_longid);
  }
                                                       timer_stop(write);
						       timer_start(sub);
  // particle subsample (binary)
  if(subsample_filename) {
    sprintf(filebase, "%s%05d%c.b", subsample_filename, snapshot->seed, suffix);
    //write_subsample(filebase, subsample_factor, snapshot, mem1, size1); // this is regular subsamle but not used. Random subampling is used.
    // TODO: periodic wrapup not done. What about after fof?
    write_random_sabsample(filebase, snapshot, mem1, size1);
  }

  // coarse mesh (binary)
  if(cgrid_filename) {
    sprintf(filebase, "%s%05d%c.b", cgrid_filename, snapshot->seed, suffix);
    coarse_grid2(filebase, snapshot, cgrid_nc, mem1, size1);
  }

                                                       timer_stop(sub);
  // text file of a slice
  //sprintf(filebase, "slice%02d", isnp); temp
  //write_snapshot_slice(filebase, snapshot, boxsize); temp


  const float ll= (float)(fof_linking_factor*boxsize/nc); // FoF linking length
  fof_find_halos(snapshot, ll);

  // FoF halo catalogue (text file)
  if(fof_filename) {
    // comment: move_particles done here
    sprintf(filebase, "%s%05d%c.txt", fof_filename, snapshot->seed, suffix);
    fof_write_halos(filebase);
  }

  const double z_out= 1.0/aout - 1.0;
  msg_printf(normal, "snapshot %d (%c) written z=%4.2f a=%5.3f\n", 
	     isnp, suffix, z_out, aout);
                                                       timer_set_category(COLA);
}
