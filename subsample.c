//
// 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "msg.h"
#include "write.h"

static gsl_rng* Random_generator= 0;
static double SubsampleFactor= 0.0;
static int ThisTask;
static int NTask;

void subsample_init(const double subsample_factor, const unsigned int seed)
{
  Random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  SubsampleFactor= subsample_factor;
    
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  const int this_node= ThisTask;
  gsl_rng_set(Random_generator, 2*seed + 100*this_node);

  msg_printf(verbose, "Subsampling initialized. Factor= %le, seed= %u\n", 
	     subsample_factor, 2*seed);  
}

void subsample_finalize(void)
{
   gsl_rng_free(Random_generator);
}

void write_random_sabsample(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size)
{
  ParticleMinimum const * const p= snapshot->p;
  const int np= snapshot->np_local;
  const int nbuf = mem_size/sizeof(ParticleMinimum);
  
  if(SubsampleFactor*np + 5*sqrt(SubsampleFactor*np) > nbuf) { 
    msg_abort(9300, 
	      "Not enough memory for local subsampling ~ %.2lf particles\n",
	      SubsampleFactor*np);
  }

  ParticleMinimum* out= (ParticleMinimum*) mem;
  int nsub= 0; // number of subsampled particles in this node
  for(int i=0; i<np; i++) {
    if(gsl_rng_uniform(Random_generator) < SubsampleFactor)
      out[nsub++]= p[i];
  }

  //
  // Prepare to gather particles. Get numbers.
  // 
  const int this_node= ThisTask;
  const int nnode= NTask;

  int* const nsub_recv= malloc(sizeof(int)*nnode*2); assert(nsub_recv);
  int* disp= nsub_recv + nnode;
  int ret= 
    MPI_Gather(&nsub, 1, MPI_INT, nsub_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  ParticleMinimum* buf1= 0;
  int ns= 0; // Total number of subsample particles
  if(this_node == 0) {
    for(int i=0; i<nnode; i++) {
      disp[i] = ns*sizeof(ParticleMinimum);
      ns  += nsub_recv[i];
      nsub_recv[i] *= sizeof(ParticleMinimum);
    }
    if(nbuf < nsub+ns)
      msg_abort(9100, "Not enough space to gather subsample particles %d (local) + %d (global)\n", nsub, ns);

    buf1= out+nsub; // Buffer for all particles

    msg_printf(verbose, "Subsampled particles %d (average %.2lf)\n", 
	       ns, SubsampleFactor*snapshot->np_total);
  }

    
  // Gather subsample particles to node 0
  ret= MPI_Gatherv(out, nsub*sizeof(ParticleMinimum), MPI_BYTE,
		   buf1, nsub_recv, disp, MPI_BYTE,
		   0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
  
  if(this_node == 0) {
    Snapshot subsample;
    subsample.p= buf1;
    subsample.np_local= ns;
    subsample.np_allocated= ns;
    subsample.np_total= ns;
    subsample.np_average= ns;
    subsample.a= snapshot->a;
    subsample.nc= 0;
    subsample.boxsize= snapshot->boxsize;
    subsample.omega_m= snapshot->omega_m;
    subsample.h= snapshot->h;
    subsample.seed= snapshot->seed;

    //write_particles_binary(filename, &subsample);
    //write_snapshot1(filename, &subsample); // Gadget Format
  }
  free(nsub_recv);
}

