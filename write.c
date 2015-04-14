#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "msg.h"
#include "comm.h"
#include "write.h"
#include "gadget_file.h"

//#ifdef LONGIDS
//typedef unsigned long long snpid_t;
//#else
//typedef unsigned int snpid_t;
//#endif

// ** particle mass
// ** 64-bit id option
void write_snapshot(const char filebase[], Snapshot const * const snapshot,
		    int use_long_id)
{
  int ThisTask;
  int NTask;

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  const double h= snapshot->h;
  char filename[256];
  sprintf(filename, "%s.%d", filebase, ThisTask);


  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const double boxsize= snapshot->boxsize;
  const double omega_m= snapshot->omega_m;

  /*
#ifdef LONGIDS
  msg_printf(normal, "LONGIDS used for snapshot. %d-byte.\n", sizeof(snpid_t));
#else
  msg_printf(normal, "ID is %d-byte unsigned int\n", sizeof(snpid_t));
#endif
  */

  if(use_long_id)
    msg_printf(normal, "Longid is used for GADGET snapshot. %d-byte.\n", 
	       sizeof(unsigned long long));
  else
    msg_printf(normal, "ID is %d-byte unsigned int\n", sizeof(unsigned int));


  long long np_send= np, np_total;
  MPI_Reduce(&np_send, &np_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&np_total, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

  GadgetHeader header; assert(sizeof(GadgetHeader) == 256);
  memset(&header, 0, sizeof(GadgetHeader));

  const double rho_crit = 27.7455;
  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np_total;
  
  header.np[1]= np;
  header.mass[1]= m;
  header.time= snapshot->a;
  header.redshift= 1.0/header.time - 1;
  header.np_total[1]= (unsigned int) np_total;
  header.np_total_highword[1]= (unsigned int) (np_total >> 32);
  //header.np_total[2]= (int) (np_total >> 32);
  header.num_files= NTask;
  header.boxsize= boxsize;
  header.omega0= omega_m;
  header.omega_lambda= 1.0 - omega_m;
  header.hubble_param= h;


  int blklen= sizeof(GadgetHeader);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  fwrite(&header, sizeof(GadgetHeader), 1, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // position
  blklen= np*sizeof(float)*3;
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++)
    fwrite(p[i].x, sizeof(float), 3, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // velocity
  const float vfac= 1.0/sqrt(snapshot->a); // Gadget convention

  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    float vout[]= {vfac*p[i].v[0], vfac*p[i].v[1], vfac*p[i].v[2]};
    fwrite(vout, sizeof(float), 3, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // id
  if(use_long_id) {
    blklen= np*sizeof(unsigned long long);
    fwrite(&blklen, sizeof(blklen), 1, fp);
    for(int i=0; i<np; i++) {
      unsigned long long id_out= p[i].id;
      fwrite(&id_out, sizeof(unsigned long long), 1, fp); 
    }
  }
  else {
    blklen= np*sizeof(unsigned int);
    fwrite(&blklen, sizeof(blklen), 1, fp);
    for(int i=0; i<np; i++) {
      unsigned int id_out= p[i].id;
      fwrite(&id_out, sizeof(unsigned int), 1, fp); 
    }
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);


  
  fclose(fp);  

  msg_printf(normal, "snapshot %s written\n", filebase);
}

// Writing Gadget file for subsampled particles
//  *number of file = 1
//  *ID is int
void write_snapshot1(const char filename[], Snapshot const * const snapshot)
{
  int ThisTask;
  int NTask;

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const double boxsize= snapshot->boxsize;
  const double omega_m= snapshot->omega_m;

  GadgetHeader header; assert(sizeof(GadgetHeader) == 256);
  memset(&header, 0, sizeof(GadgetHeader));

  const double rho_crit = 27.7455;
  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np;
  
  header.np[1]= np;
  header.mass[1]= m;
  header.time= snapshot->a;
  header.redshift= 1.0/header.time - 1;
  header.np_total[1]= (unsigned int) np;
  header.np_total_highword[1]= 0;
  header.num_files= NTask;
  header.boxsize= boxsize;
  header.omega0= omega_m;
  header.omega_lambda= 1.0 - omega_m;
  header.hubble_param= snapshot->h;


  int blklen= sizeof(GadgetHeader);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  fwrite(&header, sizeof(GadgetHeader), 1, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // position
  blklen= np*sizeof(float)*3;
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++)
    fwrite(p[i].x, sizeof(float), 3, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // velocity
  const float vfac= 1.0/sqrt(snapshot->a); // Gadget convention

  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    float vout[]= {vfac*p[i].v[0], vfac*p[i].v[1], vfac*p[i].v[2]};
    fwrite(vout, sizeof(float), 3, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // id
  blklen= np*sizeof(int);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    int id_out= p[i].id;
    fwrite(&id_out, sizeof(int), 1, fp); 
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);
  
  fclose(fp);  

  msg_printf(normal, "subsample %s written\n", filename);
}


//
// Write particle force
//

void write_force(const char filebase[], Particles const * const particles)
{
  int ThisTask;
  int NTask;

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  char filename[256];
  int inode= ThisTask;

  sprintf(filename, "%s.%d", filebase, inode);
  FILE* fp= fopen(filename, "w");
  if(fp == 0) {
    msg_abort(9010, "Unable to write force to %s\n", filename);
  }

  Particle const * const p= particles->p;
  float3 * const f= particles->force;
  const int np= particles->np_local;

  fwrite(&np, sizeof(int), 1, fp);
  for(int i=0; i<np; i++) {
    int id= (int) p[i].id;
    fwrite(&id, sizeof(int), 1, fp);
    fwrite(p[i].x, sizeof(float), 3, fp);
    fwrite(f[i], sizeof(float), 3, fp);
  }
  fwrite(&np, sizeof(int), 1, fp);

  fclose(fp);
}

// Writing binary file for subsampled particles
//  *number of file = 1
//  *ID is int
void write_particles_binary(const char filename[], Snapshot const * const snapshot)
{
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const float boxsize= snapshot->boxsize; 
  const float omega_m= snapshot->omega_m;
  const double rho_crit = 27.7455;
  const float m= omega_m*rho_crit*pow(boxsize, 3.0)/np;  
  const float redshift= 1.0/snapshot->a - 1.0;

  // Header 6 floats
  fwrite(&snapshot->boxsize, sizeof(float), 1, fp);
  fwrite(&m, sizeof(float), 1, fp);
  fwrite(&snapshot->omega_m, sizeof(float), 1, fp);
  fwrite(&snapshot->h, sizeof(float), 1, fp);
  fwrite(&snapshot->a, sizeof(float), 1, fp);
  fwrite(&redshift, sizeof(float), 1, fp);

  fwrite(&np, sizeof(int), 1, fp);

  // positions, velocities
  for(int i=0; i<np; i++) {
    fwrite(p[i].x, sizeof(float), 3, fp);
    fwrite(p[i].v, sizeof(float), 3, fp);
  }

  fwrite(&np, sizeof(int), 1, fp);
  
  int ret= fclose(fp); assert(ret == 0);

  msg_printf(normal, "subsample binary %s written\n", filename);
}
