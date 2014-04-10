//
// Reads GADGET snapshot and distribute to all nodes
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include "particle.h"
#include "comm.h"
#include "msg.h"
#include "gadget_file.h"

/*
typedef struct {
  ParticleMinimum* p;
  int np_local;
  int np_allocated;
  long long np_total;
  float np_average;
  float a; //, a_x, a_v;
  float boxsize;
  int nc;
  float omega_m, h;
  int seed;
  char filename[64];
} Snapshot;
*/

static int find_files(const char filename[]);

static void check_separator(FILE* fp, int expected_value)
{
  int n= -1;
  fread(&n, sizeof(int), 1, fp);
  
  if(n != expected_value) {
    msg_abort(15005, 
	      "Error: Unable to read snapshot correctly. Separator %d != %d\n",
	      n, expected_value);
  }  
}

int read_snapshot(const char filename[], Snapshot* snapshot, void* buf, size_t size)
{
  // returns 0 if file is not found
  const int this_node= comm_this_node();
  ParticleMinimum* const p= snapshot->p;
  float* const x= (float*) buf;

  int nread= 0;

  msg_printf(verbose, "Reading Gadget snapshot %s...\n", filename);

  int numfiles= 0;
  if(this_node == 0)
    numfiles= find_files(filename);
  
  MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(numfiles <= 0) 
    return 0;

  msg_printf(verbose, "%d files per snapshot.\n", numfiles);

  GadgetHeader header;
  float vfac_back; // Gadget convention back to km/s
      
  for(int num = 0; num<numfiles; num++) {
    int np_snapshot= 0;

    if(this_node == 0) {
      char filename_i[256];

      if(numfiles > 1)
	sprintf(filename_i, "%s.%d", filename, num);
      else
	sprintf(filename_i, "%s", filename);

      FILE* fp= fopen(filename_i, "r");
      if(fp == 0) {
	msg_abort(15000, "Error: can't open snapshot file %s (read.c)\n", 
		  filename_i);
      }
      
      check_separator(fp, 256);
      fread(&header, sizeof(GadgetHeader), 1, fp);
      check_separator(fp, 256);
      
      np_snapshot = header.np[1];   // This code only reads type 1 dark matter
      vfac_back= sqrt(header.time);

      msg_printf(verbose, "%d particles read from %s.\n", 
		 np_snapshot, filename_i);

      if(np_snapshot*sizeof(float)*6 > size) {
	msg_abort(15001, "Error: Not enough space to read snapshot. %d particles in %s, buffer size= %ld.", np_snapshot, filename_i, size); 
      }

      int ret;

      // position
      check_separator(fp, sizeof(float)*3*np_snapshot);
      ret= fread(x, sizeof(float), 3*np_snapshot, fp); 
      assert(ret == 3*np_snapshot);
      check_separator(fp, sizeof(float)*3*np_snapshot);

      const float boxsize= (float) header.boxsize;

      for(int i=0; i<3*np_snapshot; i++) {
	if(x[i] < 0.0f) x[i] += boxsize;
	if(x[i] >= boxsize) x[i] -= boxsize;
      }

      // velocity      
      check_separator(fp, sizeof(float)*3*np_snapshot);
      ret= fread(x + 3*np_snapshot, sizeof(float), 3*np_snapshot, fp);
      assert(ret == 3*np_snapshot);
      check_separator(fp, sizeof(float)*3*np_snapshot);

      fclose(fp);
    }
    
    MPI_Bcast(&np_snapshot, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&vfac_back, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, np_snapshot * 6, MPI_FLOAT, 0, MPI_COMM_WORLD);

    const float x_left= comm_xleft(0);
    const float x_right= comm_xright(0);


    float* v= x + 3*np_snapshot;

    for(int i=0; i<np_snapshot; ++i) {
      if(x_left <= x[3*i] && x[3*i] < x_right) {
	if(nread >= snapshot->np_allocated) {
	  msg_abort(15002, 
		   "Error: Not enough space (np_allocated) to read snapshot\n");
	}

	p[nread].x[0]= x[3*i    ];
	p[nread].x[1]= x[3*i + 1];
	p[nread].x[2]= x[3*i + 2];

	p[nread].v[0]= vfac_back*v[3*i    ];
	p[nread].v[1]= vfac_back*v[3*i + 1];
	p[nread].v[2]= vfac_back*v[3*i + 2];

	p[nread].id= 0;

	nread++;
      }
    }
  }

  MPI_Bcast(&header, sizeof(GadgetHeader), MPI_BYTE, 0, MPI_COMM_WORLD);

  snapshot->np_total= (((long long) header.np_total_highword[1]) << 32) +
                      (long long) header.np_total[1];
  snapshot->np_average= (float)(((double) snapshot->np_total)/comm_nnode());
  snapshot->a= (float) header.time;
  snapshot->boxsize= (float) header.boxsize;
  snapshot->omega_m= (float) header.omega0;
  snapshot->h= (float) header.hubble_param;
  snapshot->np_local= nread;

  /*
  // debug!!!
  for(int i=0; i<snapshot->np_local; i++) {
    if(!(comm_xleft(0) <= snapshot->p[i].x[0] &&
	 snapshot->p[i].x[0] < comm_xright(0))) {
      printf("assertion fail %e, node=%d\n", p[i].x[0], comm_this_node());
    }
  }
  printf("assert OK\n");
  */

  //printf("np_local %d\n", nread);
  long long np_local= nread, np_total;
  MPI_Reduce(&np_local, &np_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if(this_node == 0) {
    msg_printf(debug, "np_total %lld %lld\n", np_total, snapshot->np_total);
    assert(np_total == snapshot->np_total);
  }

  return 1;
}


int find_files(const char filename[])
{
  GadgetHeader header;

  FILE* fp= fopen(filename, "r");
  if(fp) {
    check_separator(fp, 256);
    fread(&header, sizeof(header), 1, fp);
    check_separator(fp, 256);
    fclose(fp);

    return header.num_files;
  }

  char filename_i[256];
  sprintf(filename_i, "%s.0", filename);

  fp= fopen(filename_i, "r");
  if(fp) {
    check_separator(fp, 256);
    fread(&header, sizeof(header), 1, fp);
    check_separator(fp, 256);
    fclose(fp);
    
    return header.num_files;
  }

  return 0;
}
