//
// Calculate density feild on coarse mesh nc and write as binary file
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "particle.h"
#include "msg.h"
#include "comm.h"
#include "coarse_grid.h"

static void cic_mass_assignment(Snapshot const * const snapshot, float* const grid, const int nc);
//static void write_grid(const char filename[], float* const grid, const int nc, const float boxsize);

/*
void coarse_grid(const char filename[], Snapshot const * const snapshot, const int nc, void* const mem1, const size_t size1)
{
  float* const grid= (float*) mem1;
  const int ngrid= nc*nc*nc;

  if(size1 < sizeof(float)*nc*nc*nc*2) {
    msg_abort(10010, 
	      "Not enough space for coarse grid. %ld Mbytes neseccarly > %ld\n",
	      sizeof(float)*nc*nc*nc*2/(1024*1024), size1/(1024*1024));
  }

  msg_printf(verbose, "Coarse density grid %d\n", nc);

  for(int i=0; i<ngrid; i++)
    grid[i]= 0.0f;

  cic_mass_assignment(snapshot, grid, nc);

  float* const global_grid= grid + ngrid;
  MPI_Reduce(grid, global_grid, ngrid, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  const int this_node= comm_this_node();
  if(this_node == 0) {
    double sum= 0.0, delta_sum= 0.0;

    // convert density to density constrast delta
    const float nbar_inv= (double) ngrid / snapshot->np_total;
    for(int i=0; i<ngrid; i++) {
      sum += global_grid[i];
      global_grid[i]= nbar_inv*global_grid[i] - 1.0f;
      delta_sum += global_grid[i];
    }
    
    //const double acc= 1.0e-1*snapshot->np_total/(nc*nc*nc);
    msg_printf(debug, "Coarse grid check sum %.2lf %lld\n", sum, snapshot->np_total);
    msg_printf(debug, "Coarse grid delta sum %.2lf\n", delta_sum);
    assert(fabs(sum - (double)snapshot->np_total) < 0.1*sum);

    write_grid(filename, global_grid, nc, snapshot->boxsize);
  }
}
*/

void coarse_grid2(const char filename[], Snapshot const * const snapshot, const int nc, void* const mem1, const size_t size1)
{
  // When coarse_grid() used too much memory via MPI_Reduce,
  // This version reduces little by little (nreduce_slices)
  float* grid= (float*) mem1;
  const int ngrid= nc*nc*nc;
  const int nreduce_slices=16;

  msg_printf(verbose, "Coase grid nc=%d, nreduce_slice=%d\n", nc, nreduce_slices);

  if(size1 < sizeof(float)*nc*nc*(nc + nreduce_slices)) {
    msg_abort(10010, 
	      "Not enough space for coarse grid. %ld Mbytes neseccarly > %ld\n",
	      sizeof(float)*nc*nc*nc*2/(1024*1024), size1/(1024*1024));
  }

  for(int i=0; i<ngrid; i++)
    grid[i]= 0.0f;
 

  cic_mass_assignment(snapshot, grid, nc);


  // Reduce grid data to node 0 and save
  const int this_node= comm_this_node();
  FILE* fp= 0; int ret;
  if(this_node == 0) {
    fp= fopen(filename, "w");
    if(fp == 0)
      msg_abort(10021, "Unable to write grid to file: %s\n", filename);

    ret= fwrite(&snapshot->boxsize, sizeof(float), 1, fp);   assert(ret == 1);
    ret= fwrite(&nc, sizeof(int), 1, fp);                    assert(ret == 1);
  }

  float* const grid_recv= grid + ngrid;
  const float nbar_inv= (double) ngrid / snapshot->np_total;
  int nrecv_total= 0;
  double sum= 0.0, delta_sum= 0.0;


  // Reduce grid data for each nreduce_slices to save buffer memory used by MPI
  for(int ix=0; ix<nc; ix+=nreduce_slices) {
    int nx= nc - ix < nreduce_slices ? nc - ix : nreduce_slices;
    int nrecv= nx*nc*nc;
    MPI_Reduce(grid, grid_recv, nrecv, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(this_node == 0) {
      // convert density to density constrast delta
      for(int i=0; i<nrecv; i++) {
	sum += grid_recv[i];
	grid_recv[i]= nbar_inv*grid_recv[i] - 1.0f;
	delta_sum += grid_recv[i];
      }
      ret= fwrite(grid_recv, sizeof(float), nrecv, fp); assert(ret == nrecv);

      nrecv_total += nrecv;
    }
    grid += nrecv;
  }

  msg_printf(verbose, "Coarse grid file written %s, nc=%d\n", filename, nc);

  if(this_node == 0) {
    assert(nrecv_total == nc*nc*nc);
  
    msg_printf(debug, "Coarse grid check sum %.2lf %%, %.2lf %lld\n", 
	       (sum - snapshot->np_total)/snapshot->np_total*100.0, 
	       sum, snapshot->np_total);
    msg_printf(debug, "Coarse grid delta sum %.2lf\n", delta_sum);
    assert(fabs(sum - (double)snapshot->np_total) < 0.1*sum);

    ret= fwrite(&nc, sizeof(int), 1, fp);                    assert(ret == 1);
    ret= fclose(fp);                                         assert(ret == 0);
  }
}

void cic_mass_assignment(Snapshot const * const snapshot, float* const grid, const int nc)
{
  const int np= snapshot->np_local;
  ParticleMinimum* const p= snapshot->p;
  const float dx_inv= nc/snapshot->boxsize;

  for(int i=0; i<np; i++) {
    int ix[3], ix0[3], ix1[3];
    float w[3];

    for(int j=0; j<3; ++j) {
      ix[j]= (int) floor(p[i].x[j]*dx_inv);
      w[j]= 1.0f - (p[i].x[j]*dx_inv - ix[j]);  // CIC weight for left point
      ix0[j]= (ix[j] + nc) % nc;              // left grid (periodic)
      ix1[j]= (ix[j] + 1 + nc) % nc;          // right grid (periodic)
    }
    
    grid[(ix0[0]*nc + ix0[1])*nc + ix0[2]] += w[0]*w[1]*w[2];
    grid[(ix0[0]*nc + ix1[1])*nc + ix0[2]] += w[0]*(1-w[1])*w[2];
    grid[(ix0[0]*nc + ix0[1])*nc + ix1[2]] += w[0]*w[1]*(1-w[2]);
    grid[(ix0[0]*nc + ix1[1])*nc + ix1[2]] += w[0]*(1-w[1])*(1-w[2]);

    grid[(ix1[0]*nc + ix0[1])*nc + ix0[2]] += (1-w[0])*w[1]*w[2];
    grid[(ix1[0]*nc + ix1[1])*nc + ix0[2]] += (1-w[0])*(1-w[1])*w[2];
    grid[(ix1[0]*nc + ix0[1])*nc + ix1[2]] += (1-w[0])*w[1]*(1-w[2]);
    grid[(ix1[0]*nc + ix1[1])*nc + ix1[2]] += (1-w[0])*(1-w[1])*(1-w[2]);
  }
}



/*
void write_grid(const char filename[], float* const grid, const int nc, const float boxsize)
{
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(10020, "Unable to write grid to file: %s\n", filename);

  const int ngrid= nc*nc*nc;

  int ret;
  ret= fwrite(&boxsize, sizeof(float), 1, fp);             assert(ret == 1);
  ret= fwrite(&nc, sizeof(int), 1, fp);                    assert(ret == 1);
  ret= fwrite(grid, sizeof(float), nc*nc*nc, fp);          assert(ret == ngrid);
  ret= fwrite(&nc, sizeof(int), 1, fp);                    assert(ret == 1);
  ret= fclose(fp);                                         assert(ret == 0);

  msg_printf(verbose, "Coarse grid file written %s, nc=%d\n", filename, nc);
}
*/
