#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include "msg.h"
#include "move.h"
#include "comm.h"
#include "particle.h"

//
// Move particles ver 2
//

static int sendrecv_particles(const int dix, Particle* const buf, const int nbuf, Particles* const particles);

void move_particles2(Particles* const particles, const float BoxSize, void* const vbuf, const size_t vsize)
{
  // particles->p and particles->np_local will be modified
  Particle* buf= (Particle*) vbuf;
  const int nbuf= (int)(vsize/sizeof(Particle));

  Particle* const p= particles->p;  
  const int np_local= particles->np_local;

  const int ThisNode= comm_this_node();
  const int NNode= comm_nnode();
  const float x_left= comm_xleft(0);
  const float x_right= comm_xright(0);

  int nsend= 0, move_to= 0;

  msg_printf(verbose, "Inter-node particle transfer\n");

  for(int i=0; i< np_local; i++) {
    if(p[i].x[0] < 0.0f) p[i].x[0] += BoxSize;
    else if(p[i].x[0] >= BoxSize) p[i].x[0] -= BoxSize;

    if(p[i].x[1] < 0.0f) p[i].x[1] += BoxSize;
    else if(p[i].x[1] >= BoxSize) p[i].x[1] -= BoxSize;

    if(p[i].x[2] < 0.0f) p[i].x[2] += BoxSize;
    else if(p[i].x[2] >= BoxSize) p[i].x[2] -= BoxSize;

    if(p[i].x[0] < x_left || p[i].x[0] >= x_right) {
      if(nsend >= nbuf)
	msg_abort(6030, "Error: Not enough space to move particles: "
		        "%d particles moving\n", nsend);
      buf[nsend++]= p[i];
    }
    else
      p[move_to++]= p[i];
  }

  const int np_not_moved= move_to;
  particles->np_local= move_to;


  msg_printf(verbose, "%d particles needs to be moved (%.3f%%).\n", 
	     nsend, (float)nsend/np_local*100.0f);

  for(int dix=1; dix < NNode; dix++) {
    nsend= sendrecv_particles(dix, buf, nsend, particles);
    nsend= sendrecv_particles(-dix, buf, nsend, particles);

    int nsend_max= comm_share_int(nsend, MPI_MAX);
    if(nsend_max == 0) break;
    
    msg_printf(verbose, "%d particles remaining, continue step %d.\n", 
	       nsend_max, dix+1);
  }

  int np_max= comm_reduce_int(particles->np_local, MPI_MAX);
  msg_printf(info, "Particle load imbalance %.2f %d/%.1f\n", np_max/particles->np_average, np_max, particles->np_average);

  // post-condition check
  // particles are all in the volume of this node
  const int np_new= particles->np_local;
  for(int i=np_not_moved; i<np_new; i++) {
    assert(x_left - 1.0001f*BoxSize <= p[i].x[0] && 
	                               p[i].x[0] < x_right + 0.0001f*BoxSize);
    assert(-0.0001f*BoxSize <= p[i].x[1] && p[i].x[1] < 1.0001f*BoxSize);
    assert(-0.0001f*BoxSize <= p[i].x[2] && p[i].x[2] < 1.0001f*BoxSize);
  }

  // Sanity checks
  // Total number of particles checked
  long long np= particles->np_local, np_global;
  MPI_Reduce(&np, &np_global, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisNode == 0 && np_global != particles->np_total)
    msg_abort(6040, 
         "Error: Number of particles not conserved after transfer %lld %lld\n", 
	 np_global, particles->np_total);
  else
    msg_printf(debug, 
	       "%lld particles after exchange_particles. Total number OK.\n", 
	       np_global);

}

int sendrecv_particles(const int dix, Particle* const buf, const int nbuf, Particles* const particles)
{
  Particle* const p= particles->p;
  const int np_local= particles->np_local;
  const int np_alloc= particles->np_allocated;

  const float x_left= comm_xleft(dix);
  const float x_right= comm_xright(dix);

  //fprintf(stderr, "x_left x_right sendrecv particle %f %f\n", x_left, x_right);

  // particles in buffer are moved such that
  // index 0..nstay-1 are not in x_left <= x < x_right
  // index nstay..nbuf-1 are in that range -> send to node dix

  Particle tmp;
  int i=0,j=nbuf-1;
  while(1) {
    while(i<nbuf && !(x_left <= buf[i].x[0] && buf[i].x[0] < x_right))
      i++;
    while(j>=0 && x_left <= buf[j].x[0] && buf[j].x[0] < x_right)
      j--;

    if(i < j) {  // swap
      tmp= buf[i]; buf[i]= buf[j]; buf[j]= tmp;
      i++; j--;
    }
    else
      break;
  }
  //printf("%d %d\n", i, j);

  int nsend, nstay;
  if(j == -1) {
    nsend= nbuf;
    nstay= 0;
  }
  else if(i == nbuf) {
    nsend= 0;
    nstay= nbuf;
  }
  else {
    assert(j+1 == i);
    nstay= i;
    nsend= nbuf-i;
  }
  Particle* sendbuf= buf+nstay;


  int nrecv;
  int tag= 800+5*dix;
  MPI_Status status;

  const int node_to= comm_node(dix);
  const int node_from= comm_node(-dix);

  // *** debug!!!
  //fprintf(stderr, "n%d Sending %d particles to node %d (step %+d)\n", 
  //	  comm_this_node(), nsend, node_to, dix);

  //msg_printf(verbose, "Sending %d particles to node %d (step %+d)\n", 
  //nsend, node_to, dix);

  MPI_Barrier(MPI_COMM_WORLD); // debug!!!
  MPI_Sendrecv(&nsend, 1, MPI_INT, node_to, tag, 
	       &nrecv, 1, MPI_INT, node_from, tag, MPI_COMM_WORLD, &status); 
  tag++;

  if(np_local + nrecv > np_alloc) {
    msg_printf(fatal, "%d + %d > %d\n", np_local, nrecv, np_alloc);
    msg_abort(6041, "Error: Not enough space for particles (np_alloc)\n");
  }

  MPI_Sendrecv(sendbuf,    nsend*sizeof(Particle), MPI_BYTE, node_to,   tag,
               p+np_local, nrecv*sizeof(Particle), MPI_BYTE, node_from, tag,
	       MPI_COMM_WORLD, &status); tag++;

  msg_printf(verbose, "Received %d particles from node %d\n", nrecv, node_from);
  
  particles->np_local= np_local + nrecv;

  return nstay;
}
