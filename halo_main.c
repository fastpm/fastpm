//
// Main file for program "halo"
//   Read Gadget snapshot and does data analysis only, 
//   FoF halo finding/density field/particle subsampling
//

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
#include "comm.h"
#include "fof.h"
#include "write.h"
#include "timer.h"
#include "mem.h"
#include "move.h"
#include "subsample.h"
#include "coarse_grid.h"
#include "read.h"

// TODO: update this with parameters in param
const double OmegaLambda = 0.727;
const double OmegaM      = 0.273;
const double Hubble      = 0.705;
const double sigma8      = 0.812;

static const int nc_factor= 3;

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
		   Particles const * const particles, 
		   Snapshot * const snapshot, 
		   const char subsample_filename[], 
		   const char cgrid_filename[], const int cgrid_nc,
		   void* const mem1, const size_t size1,
		   const int write_longid);

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
  msg_set_loglevel(param.loglevel);

  fft_init(multi_thread);
  comm_init(nc_factor*param.nc, param.nc, param.boxsize);

  //const int nsteps= param.ntimestep;
  //const double a_final= param.a_final;

  //power_init("camb0_matterpower.dat", a_init, sigma8, OmegaM, OmegaLambda);
  //power_init(param.power_spectrum_file, a_init, sigma8, OmegaM, OmegaLambda);


  Memory mem; 
  allocate_shared_memory(param.nc, nc_factor, param.np_alloc_factor, &mem); 
  //lpt_init(param.nc, mem.mem1, mem.size1, param.ext, param.dump);
  //const int local_nx= lpt_get_local_nx();

  const int nc= param.nc;
  ptrdiff_t local_nx, local_x_start;
  fftwf_mpi_local_size_3d(nc, nc, nc, MPI_COMM_WORLD,
			  &local_nx, &local_x_start);


  Particles* particles= 
     allocate_particles(param.nc, local_nx, param.np_alloc_factor);
  Snapshot* snapshot= allocate_snapshot(param.nc, local_nx, 
				  particles->np_allocated, mem.mem2, mem.size2);
  snapshot->boxsize= param.boxsize;
  snapshot->omega_m= OmegaM;
  snapshot->h= Hubble;
  snapshot->seed= 0;//***
  strncpy(snapshot->filename, param.snapshot_filename, 64);

  fof_init(particles->np_allocated, param.nc, mem.mem1, mem.size1);
  subsample_init(param.subsample_factor, param.random_seed);

  const int this_node= comm_this_node();
  //

  // for each directories given
  char dirname[256];

  FILE* fp_dir= 0;
  if(this_node == 0) {
    fp_dir= fopen("dirs.txt", "r");
    if(fp_dir == 0)
      msg_abort(0, "Unable to open list of directories dirs.txt\n");
  }    
                                                       timer_set_category(Snp);
  while(1) {
    if(this_node == 0) {
      char* sget= fgets(dirname, 255, fp_dir);
      if(sget == NULL)
	dirname[0]= '\0';
      else {
	const int last= strlen(dirname)-1;
	if(dirname[last] == '\n') dirname[last]= '\0';
      }
    }

    MPI_Bcast(dirname, 255, MPI_BYTE, 0, MPI_COMM_WORLD);

    if(dirname[0] == '\0')
      break;

    msg_printf(verbose, "Entering directory %s\n", dirname);

    char filebase[256];
    for(int iout=0; iout<100000; iout++) { 
      //const int iout= 0;
      //char suffix= 'a' + iout;
      
      if(snapshot->filename)
	msg_abort(1, "filename snapshot is not set in the parameter file\n");

      sprintf(filebase, "%s/%s_%03d", dirname, param.snapshot_filename, iout);

      int file_found= 
	read_snapshot(filebase, snapshot, mem.mem1, mem.size1);

      if(file_found == 0) {
	if(iout == 0)
	  msg_abort(2, "Unable to read snapshot %s\n", filebase);
	break;
      }
						       
      const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);

						       timer_start(sub);
      // subsample
      if(param.subsample_filename) {
	sprintf(filebase, "%s/%s_%03d.b", dirname, param.subsample_filename, iout);
	//sprintf(filebase, "%s%05d%c.b", param.subsample_filename, snapshot->seed, suffix);
	//write_subsample(filebase, subsample_factor, snapshot, mem1, size1);
	write_random_sabsample(filebase, snapshot, mem.mem1, mem.size1);
      }

      // coarse mesh
      if(param.cgrid_filename) {
	sprintf(filebase, "%s/%s_%03d.b", dirname, param.cgrid_filename, iout);
	//sprintf(filebase, "%s%05d%c.b", param.cgrid_filename, snapshot->seed, suffix);
	coarse_grid2(filebase, snapshot, param.cgrid_nc, mem.mem1, mem.size1);
      }

                                                       timer_stop(sub);

      const float ll= 0.2*boxsize/nc; // FOF linking length
      fof_find_halos(snapshot, ll);

      // text file of fof halo
      sprintf(filebase, "%s/fof_%03d.txt", dirname, iout);
      //sprintf(filebase, "fof%05d%c.txt", snapshot->seed, suffix);
      fof_write_halos(filebase);

                                                       timer_set_category(COLA);



      timer_print();
    }
  }

  
  if(this_node == 0)
    fclose(fp_dir);
  
  msg_printf(normal, "Processed all directories. halo done.\n");
  // done
  /*
  const int nout= 4;
  float aout[nout];
  for(int i=0; i<nout; i++)
    aout[i] = 1.0/(1+z_out[i]);

  for(int irealization=0; irealization<param.nrealization; irealization++) {
    MPI_Barrier(MPI_COMM_WORLD);
    msg_printf(verbose, "\n%d/%d realization.\n", irealization+1, param.nrealization);
    int seed= param.random_seed + irealization;

    int iout= 0;

    // Set initial grid and 2LPT desplacement
    //timer_set_category(LPT);
    //lpt_set_displacement(a_init, OmegaM, seed,
    //param.boxsize, particles);
    //snapshot->seed= seed;
    //cola_set_snapshot(a_init, particles, snapshot);

    // write initial condition to file

    if(param.init_filename) {
      set_noncola_initial(particles, snapshot);
      write_snapshot(param.init_filename, snapshot, Hubble, 
		     param.write_longid);
    }

    timer_set_category(COLA);

#ifndef LOGTIMESTEP
    particles->a_v= 0.5*da;
    msg_printf(info, "timestep linear in a\n");
#else
    const double loga_init= log(a_init);
    const double dloga= (log(a_final) - log(a_init))/nsteps;
    particles->a_v= exp(log(a_init) - 0.5*dloga);
    msg_printf(info, "timestep linear in loga\n");
#endif


    // Time evolution loop
    if(nsteps > 1 && a_final > a_init) {
      msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
		 a_init, a_final, nsteps);
      for (int istep=1; istep<=nsteps; istep++) {
	msg_printf(normal, "Timestep %d/%d\n", istep, nsteps);
      
                                                            timer_start(comm);
        // move particles to other nodes
        move_particles2(particles, param.boxsize, mem.mem1, mem.size1 );

                                                            timer_stop(comm);

        pm_calculate_forces(particles); 

#ifndef LOGTIMESTEP
	double avel0= (istep-0.5)*da;
	double apos0=  istep*da;
	
	double avel1= (istep+0.5)*da;
	double apos1= (istep+1.0)*da;
#else
	float avel0= exp(loga_init + (istep-0.5)*dloga);
	float apos0= exp(loga_init + istep*dloga);
	
	float avel1= exp(loga_init + (istep+0.5)*dloga);
	float apos1= exp(loga_init + (istep+1)*dloga);
#endif
      
	while(iout < nout && avel0 <= aout[iout] && aout[iout] <= apos0) {
	  snapshot_time(aout[iout], iout, particles, snapshot, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid);
	  iout++;
	}
	if(iout >= nout) break;

	cola_kick(particles, OmegaM, avel1);

	while(iout < nout && apos0 < aout[iout] && aout[iout] <= avel1) {
	  snapshot_time(aout[iout], iout,  particles, snapshot, param.subsample_filename, param.cgrid_filename, param.cgrid_nc, mem.mem1, mem.size1, param.write_longid);
	  iout++;
	}
	if(iout >= nout) break;

	cola_drift(particles, OmegaM, apos1);
      }
    }

  }

  // write ****
  //write_snapshot("snp", particles, param.boxsize, OmegaM, Hubble);

  // FOF halo finding ****

  //move_particles(particles);
  //write_slice("slice", particles->p, particles->np_local, param.boxsize/32, param.boxsize);
  

  //const float ll= 0.2*param.boxsize/param.nc;
  //fof_find_halos(particles, param.boxsize, ll);


  //write_slice("slice_ev", particles->p, particles->np_local, param.boxsize/32);

  // finalize:
  */

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

    
