//
// MPI Communication
//
// Written by Jun Koda. 
// Needed by the Friend of Friend finder; otherwise
// replaced by domain.c

// #include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3-mpi.h>

#include "msg.h"
#include "comm.h"


static int Ngrid;
static int ThisNode, NNode, LeftNode=-1, RightNode=-1;
static int Local_nx, Local_x_start;
static int Tag;
static int Node[2];

static float BoxSize;

typedef struct {
  int node, ix, nx;
} Slice;

static Slice* Slices;
static int ThisSlice;

int compare_slice(void const * const a, void const * const b)
{
  return ((Slice*)a)->ix - ((Slice*)b)->ix;
}

void comm_init(const int nc_pm, const int nc_p, const float boxsize)
{
  /*
  ptrdiff_t local_nx, local_x_start;
  fftwf_mpi_local_size_3d(nc_p, nc_p, nc_p, MPI_COMM_WORLD,
			  &local_nx, &local_x_start);

  Np_local_average= local_nx*nc_p*nc_p;
  */

  //** copied from pm.c code is redundunt now

  msg_printf(verbose, "comm initialization\n");

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisNode);
  MPI_Comm_size(MPI_COMM_WORLD, &NNode);

  ptrdiff_t local_nx, local_x_start, local_ny, local_y_start;
  int* local_x_table  = malloc(sizeof(int)*NNode*2); assert(local_x_table);
  int* local_nx_table = local_x_table + NNode;

  // Print LPT initial domein decomposition (just for information**)
  fftwf_mpi_local_size_3d(nc_p, nc_p, nc_p, MPI_COMM_WORLD,
			  &local_nx, &local_x_start);

  MPI_Allgather(&local_nx, 1, MPI_INT, local_nx_table, 1, MPI_INT, 
		MPI_COMM_WORLD);
  MPI_Allgather(&local_x_start, 1, MPI_INT, local_x_table, 1, MPI_INT, 
		MPI_COMM_WORLD);

  for(int i=0; i<NNode; i++)
    msg_printf(debug, "LPT Task=%d x=%d..%d\n", i, local_x_table[i],
                      local_x_table[i]+local_nx_table[i]-1);

  // PM domain decomposition
  Ngrid= nc_pm;
  BoxSize= boxsize;
  Tag= 600;

  fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid, MPI_COMM_WORLD,
			  &local_nx, &local_x_start, &local_ny, &local_y_start);

  Local_nx= local_nx;
  Local_x_start= local_x_start;


  MPI_Allgather(&Local_nx, 1, MPI_INT, local_nx_table, 1, MPI_INT, 
		MPI_COMM_WORLD);
  MPI_Allgather(&Local_x_start, 1, MPI_INT, local_x_table, 1, MPI_INT, 
		MPI_COMM_WORLD);

  for(int i=0; i<NNode; i++)
    msg_printf(debug, "Task=%d x=%d..%d\n", i, local_x_table[i],
                      local_x_table[i]+local_nx_table[i]-1);

  for(int i=0; i<NNode; i++) {
    if((local_x_table[i] + local_nx_table[i]) % Ngrid == Local_x_start)
      LeftNode= i;
    if((Local_x_start + Local_nx) % Ngrid == local_x_table[i])
      RightNode= i;
  }
  msg_printf(debug, "Node %d LeftNode= %d, RightNode= %d\n", 
	     ThisNode, LeftNode, RightNode);
  assert(LeftNode >= 0 && RightNode >= 0);
  
  Node[ToRight]= RightNode;
  Node[ToLeft]= LeftNode;

  // Slice information
  Slices= malloc(sizeof(Slice)*NNode); assert(Slices);
  for(int i=0; i<NNode; i++) {
    Slices[i].node= i;
    Slices[i].ix= local_x_table[i];
    Slices[i].nx= local_nx_table[i];
  }
  qsort(Slices, NNode, sizeof(Slice), compare_slice);

  for(int i=0; i<NNode; i++) {
    if(Slices[i].node == ThisNode) {
      ThisSlice= i; break;
    }
  }

  //printf("Node initialized Node[%d]= %d, Node[%d]= %d\n", ToRight, RightNode, ToLeft, LeftNode);

  free(local_x_table);
}

void comm_finalize(void){

  free(Slices);

}


float comm_xmin(void)
{
  return BoxSize*Local_x_start/Ngrid;
}

float comm_xmax(void)
{
  return BoxSize*(Local_nx + Local_x_start)/Ngrid;
}

int comm_right_edge(void)
{
  return Local_x_start + Local_nx == Ngrid;
}

int comm_get_nrecv(enum Direction direction, int nsend)
{
  int nrecv;
  MPI_Status status;

  //printf("Node %d %d\n", Node[0], Node[1]);
  //printf("%d: %d -> %d\n", ThisNode, Node[direction], Node[!direction]);
  fflush(stdout);

  assert(Ngrid);
  int ret= MPI_Sendrecv(&nsend, 1, MPI_INT, Node[direction], Tag, 
			&nrecv, 1, MPI_INT, Node[!direction], Tag, 
			MPI_COMM_WORLD, &status); 
  assert(ret == MPI_SUCCESS);

  Tag++;

  return nrecv;
}

void comm_sendrecv(enum Direction direction, void* send, int nsend, void* recv, int nrecv, MPI_Datatype datatype)
{
  MPI_Status status;

  MPI_Sendrecv(send, nsend, datatype, Node[direction], Tag,
	       recv, nrecv, datatype, Node[!direction], Tag,
	       MPI_COMM_WORLD, &status);
  Tag++;
}

int comm_get_total_int(int x)
{
  int x_global;
  MPI_Reduce(&x, &x_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  return x_global;
}

int comm_reduce_int(int x, MPI_Op op)
{
  int x_global;
  MPI_Reduce(&x, &x_global, 1, MPI_INT, op, 0, MPI_COMM_WORLD);

  return x_global;
}

int comm_share_int(int x, MPI_Op op)
{
  int x_global;
  MPI_Reduce(&x, &x_global, 1, MPI_INT, op, 0, MPI_COMM_WORLD);

  MPI_Bcast(&x_global, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  return x_global;
}


int comm_this_node(void)
{
  return ThisNode;
}

int comm_nnode(void)
{
  return NNode;
}

float comm_xleft(const int dix) {
  const int i= (ThisSlice + dix + NNode) % NNode;

  if(Slices[i].ix == 0)
    return -0.0001*BoxSize;

  return BoxSize*Slices[i].ix/Ngrid;
}

float comm_xright(const int dix) {
  const int i= (ThisSlice + dix + NNode) % NNode;
  if(Slices[i].ix + Slices[i].nx == Ngrid)
    return 1.0001*BoxSize;

  return BoxSize*(Slices[i].ix + Slices[i].nx)/Ngrid;
}

int comm_node(const int dix) {
  const int i= (ThisSlice + dix + NNode) % NNode;
  return Slices[i].node;
}

//int comm_local_nx(void)
//{
//  return Local_nx;
//}
