#ifndef COMM_H
#define COMM_H 1

#include <mpi.h>

enum Direction {ToRight=0, ToLeft=1};

void comm_init(const int nc_pm, const int nc_part, const float boxsize);

int comm_this_node(void);
int comm_nnode(void);
float comm_xmax(void);
float comm_xmin(void);

int comm_right_edge(void);
int comm_get_nrecv(enum Direction direction, int nsend);
void comm_sendrecv(enum Direction direction, void* send, int nsend, void* recv, int nrecv, MPI_Datatype datatype);
int comm_get_total_int(int x);
int comm_reduce_int(int x, MPI_Op op);
int comm_share_int(int x, MPI_Op op);

float comm_xleft(const int dix);
float comm_xright(const int dix);
int comm_node(const int dix);
//int comm_local_nx(void);

#endif
