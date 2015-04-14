//
// Particle Mesh gravitational force calculation
//
// This code is a modification to the original serial COLA code
// by Svetlin Tassev. See below.
//

/*
    Copyright (c) 2011-2013       Svetlin Tassev
                           Harvard University, Princeton University
 
    This file is part of COLAcode.

    COLAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    COLAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with COLAcode.  If not, see <http://www.gnu.org/licenses/>.
*/


// This file contains some standard functions for a PM code. 
// Nothing COLA-specific.

#include <stdlib.h>
#include <alloca.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "pm.h"
#include "msg.h"
#include "comm.h"
#include "timer.h"
#include "domain.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static int ThisNode, NNode, LeftNode=-1, RightNode=-1;

static int Ngrid;
static size_t NgridL; // for index to avoid 4byte integer overflow
static int PM_factor;
static int Local_nx, Local_x_start;
static int Local_ny_td, Local_y_start_td;  // transposed
static float BoxSize;
//static float NP_Average;

static fftwf_complex* fftdata;
static fftwf_plan p0, p11;

static double * PowerSpectrumVariable;

//static fftwf_complex *P3D;
//static float *density;
//static fftwf_plan p0;
//static fftwf_plan p11; //,p12,p13;

//static fftwf_complex *FN11;
//static float* N11;

//static void* Buf;
//static Particle* ParticleBuffer;
//static int NParticleBuffer;
static long long NParticleTotal;

/* the transfer functions for force in fourier space applied to potential */
/* super lanzcos in CH6 P 122 Digital Filters by Richard W. Hamming */
static double super_lanzcos_diff_kernel_3(double w) {
/* order N = 3*/
    return 1. / 594 * 
       (126 * sin(w) + 193 * sin(2 * w) + 142 * sin (3 * w) - 86 * sin(4 * w));
}
static double super_lanzcos_diff_kernel_2(double w) {
/* order N = 2*/
    return 1 / 126. * (58 * sin(w) + 67 * sin (2 * w) - 22 * sin(3 * w));
}
static double super_lanzcos_diff_kernel_1(double w) {
/* order N = 1 */
/* 
 * This is the same as GADGET-2 but in fourier space: 
 * see gadget-2 paper and Hamming's book.
 * c1 = 2 / 3, c2 = 1 / 12
 * */
    return 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}
static double super_lanzcos_diff_kernel_0(double w) {
    return w;
}

static double (*diff_kernel)(double w);

static inline void WRtPlus(float * const d, 
		    const int i, const int j, const int k, const float f)
{
#ifdef _OPENMP
  #pragma omp atomic
#endif
  d[k + 2*(NgridL/2 + 1)*(j + NgridL*i)] += f;
}

static inline float REd(float const * const d, const int i, const int j, const int k)
{
  return d[k + 2*(NgridL/2 + 1)*(j + NgridL * i)];
}
void pm_set_diff_order(int order) {
    double (*kernels[4])(double w) = {
          super_lanzcos_diff_kernel_0,
          super_lanzcos_diff_kernel_1,
          super_lanzcos_diff_kernel_2,
          super_lanzcos_diff_kernel_3};
    if(order > 4) {
        msg_abort(8888, "differentiation of order %d unsupported\n", order);
    } 
    if(order < 0) order = 0;
    diff_kernel = kernels[order];
}
void pm_init(const int nc_pm, const int nc_pm_factor, const float boxsize,
	     void* const mem1, const size_t size1,
	     void* const mem2, const size_t size2, int many)
{
  // Assume FFTW3 is initialized in lpt.c

  Ngrid= nc_pm; NgridL= nc_pm;
  PowerSpectrumVariable = malloc(sizeof(double) * Ngrid / 2);
  BoxSize= boxsize;
  PM_factor= nc_pm_factor;

  ptrdiff_t local_nx, local_x_start, local_ny, local_y_start;
  ptrdiff_t total_size= 
    fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid, MPI_COMM_WORLD,
	                 &local_nx, &local_x_start, &local_ny, &local_y_start);

  Local_nx= local_nx;
  Local_x_start= local_x_start;
  Local_ny_td= local_ny;
  Local_y_start_td= local_y_start;

  size_t bytes= sizeof(fftwf_complex)*total_size;

  if(mem1 == 0) {
    fftdata=  fftwf_alloc_complex(total_size);
    msg_printf(info, 
	       "%d Mbytes allocated for particle mesh (density)\n", 
	       (int)(bytes/(1024*1024)));
  }
  else {
    assert(size1 >= total_size*sizeof(fftwf_complex));
    fftdata= (fftwf_complex*) mem1;
  }
  //density= (float*) P3D;

  /*
  if(mem2 == 0) {
    FN11= fftwf_alloc_complex(total_size);
    msg_printf(info, 
	       "%d Mbytes allocated for particle mesh (force)\n", 
	       (int)(bytes/(1024*1024)));
  }
  else {
    assert(size2 >= total_size*sizeof(fftwf_complex));
    FN11= (fftwf_complex*) mem2;
  }
  */

  if(fftdata == 0)
    msg_abort(4001, "Unable to allocate memory for particle mesh: %d Mbytes\n",
	      (int)(2*bytes/(1024*1024)));

  double pm_index_max= (double) Ngrid*Ngrid*(Ngrid/2+1);
  if(pm_index_max > pow(2.0, 31))
    msg_printf(info, 
	       "Local number of PM Mesh exceeds 4-byte integer limit %ld\n", 
	       2*NgridL*NgridL*(NgridL/2+1));

  //FN11= P3D;
  //N11= (float*) FN11;

  msg_printf(verbose, "Setting up FFTW3 plans\n");

  p0=  fftwf_mpi_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, 
				 (float*) fftdata, fftdata,
		    MPI_COMM_WORLD, (many?FFTW_MEASURE:FFTW_ESTIMATE) | FFTW_MPI_TRANSPOSED_OUT);
  
  // inverse FFT
  p11= fftwf_mpi_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, fftdata, (float*)fftdata,
		     MPI_COMM_WORLD, (many?FFTW_MEASURE:FFTW_ESTIMATE) | FFTW_MPI_TRANSPOSED_IN);

  // FFTW_MEASURE is probably better for multiple realization **

  // Find Adjacent Nodes
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisNode);
  MPI_Comm_size(MPI_COMM_WORLD, &NNode);

  int* local_x_table  = malloc(sizeof(int)*NNode*2);
  int* local_nx_table = local_x_table + NNode;

  MPI_Allgather(&Local_nx, 1, MPI_INT, local_nx_table, 1, MPI_INT, 
		MPI_COMM_WORLD);
  MPI_Allgather(&Local_x_start, 1, MPI_INT, local_x_table, 1, MPI_INT, 
		MPI_COMM_WORLD);

  //for(int i=0; i<NNode; i++)
  //  msg_printf(debug, "Task=%d x=%d..%d\n", i, local_x_table[i],
  //                    local_x_table[i]+local_nx_table[i]-1);

  for(int i=0; i<NNode; i++) {
    if((local_x_table[i] + local_nx_table[i]) % Ngrid == Local_x_start)
      LeftNode= i;
    if((Local_x_start + Local_nx) % Ngrid == local_x_table[i])
      RightNode= i;
  }
  //msg_printf(debug, "Node %d LeftNode= %d, RightNode= %d\n", 
  //     ThisNode, LeftNode, RightNode);
  assert(LeftNode >= 0 && RightNode >= 0);
  
  free(local_x_table);
  //
  // buffer for particle transfer
  const int ncp= nc_pm/nc_pm_factor;
  NParticleTotal= (long long) ncp*ncp*ncp;


}

void pm_finalize(void)
{
  //fftwf_free(fftdata);
  //fftwf_destroy_plan(plan);
  //fftwf_free(P3D);
  //fftwf_free(FN11);

  fftwf_destroy_plan(p0);
  fftwf_destroy_plan(p11);
}


void PtoMesh(const Particle P[], const int np, float* const density)
{
  // ** precondition
  //   particles are assumed to be periodiclly wraped up in y,z direction
  //float* const density= (float*) fftdata;

  msg_printf(verbose, "Calculating PtoMesh\n");

  const float scaleBox=((float) Ngrid)/((float) BoxSize);
  const float WPAR= pow(PM_factor, 3);

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Ngrid; j++)
      for(int k = 0; k < Ngrid; k++)
	density[(i*NgridL + j)*2*(NgridL/2 + 1) + k] = 0.0f;

    int c1 = 0;
    int c2 = 0;
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0; i<np; i++) {
    float X=P[i].x[0]*scaleBox;
    float Y=P[i].x[1]*scaleBox;
    float Z=P[i].x[2]*scaleBox;

    int iI=(int) floorf(X); // without floor, -1 < X < 0 is mapped to iI=0
    int J=(int) floorf(Y);          // Assumes Y,Z are positive
    int K=(int) floorf(Z);
    float D1=X-((float) iI);
    float D2=Y-((float) J);
    float D3=Z-((float) K);
    float T1=1.-D1;
    float T2=1.-D2;
    float T3=1.-D3;

    float T2W =T2*WPAR;
    float D2W =D2*WPAR;

#ifdef CHECK
    assert(Y >= 0.0f && Z >= 0.0f);
#endif
            
    // Do periodic wrapup in all directions. 
    // Buffer particles are copied from adjacent nodes
    while(iI >= Ngrid) iI -= Ngrid;
    while(J >= Ngrid) J -= Ngrid;
    while(K >= Ngrid) K -= Ngrid;
            
    int I1=iI+1; if(I1 >= Ngrid) I1=0;
    int J1=J+1; if(J1 >= Ngrid) J1=0; // assumes y,z < BoxSize
    int K1=K+1; if(K1 >= Ngrid) K1=0;

    iI -= Local_x_start;
    I1 -= Local_x_start;

    if(0 <= iI && iI < Local_nx) {
      WRtPlus(density, iI, J,  K,  T3*T1*T2W);
      WRtPlus(density, iI, J,  K1, D3*T1*T2W);
      WRtPlus(density, iI, J1, K,  T3*T1*D2W);
      WRtPlus(density, iI, J1, K1, D3*T1*D2W);
      c1 += 1;
    }

    if(0 <= I1 && I1 < Local_nx) {
      WRtPlus(density, I1, J,  K,  T3*D1*T2W);
      WRtPlus(density, I1, J,  K1, D3*D1*T2W);
      WRtPlus(density, I1, J1, K,  T3*D1*D2W);
      WRtPlus(density, I1, J1, K1, D3*D1*D2W);
      c2 += 1;
    }
  }
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Ngrid; j++)
      for(int k = 0; k < Ngrid; k++)
	density[(i*NgridL + j)*2*(NgridL/2 + 1) + k] += -1.0f;


  msg_printf(verbose, "CIC density assignment finished.\n");
}

static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}
// FFT density mesh and copy it to density_k
void compute_density_k(fftwf_complex * density_k)
{
    // FFT density(x) mesh -> density(k)
    fftwf_mpi_execute_dft_r2c(p0, (float*) fftdata, fftdata);

    // copy density(k) in fftdata to density_k

    /* the CIC deconvolution kernel is
     *
     * sinc_unnormed(k_x L / 2 All.Nmesh) ** 2
     *
     * k_x = kpos * 2pi / L
     *
     * */

    /* 
     * we deconvolve CIC twice here.
     * which means, in powerspectrum we need to convolve once.
     * */
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int Jl=0; Jl<Local_ny_td; Jl++) {
        int J = Jl + Local_y_start_td;
        int J0 = J <= (Ngrid/2) ? J : J - Ngrid;
        for(int iI=0; iI<Ngrid; iI++) {
            int I0 = iI <= (Ngrid/2) ? iI : iI - Ngrid;
            for(int K=0; K<Ngrid/2+1; K++){
                int K0 = K;
                size_t index = K + (NgridL/2+1)*(iI + NgridL*Jl);
                density_k[index][0] = fftdata[index][0];
                density_k[index][1] = fftdata[index][1];
            }
        }
    }
}
double * pm_compute_power_spectrum(size_t * nk) {
    double * power = PowerSpectrumVariable;
    nk[0] = Ngrid / 2;
    return power;
}
void compute_power_spectrum(fftwf_complex * density_k) {
    msg_printf(verbose, "Calculating power spectrum...\n");
    double * power = PowerSpectrumVariable;
    double * count = (double*) alloca(sizeof(double) * (Ngrid / 2));
    for(int i = 0; i < Ngrid / 2; i ++) {
        count[i] = 0;
        power[i] = 0;
    }
    for(int Jl=0; Jl<Local_ny_td; Jl++) {
        int J = Jl + Local_y_start_td;
        int J0 = J <= (Ngrid/2) ? J : J - Ngrid;
        const double fy = sinc_unnormed(J0 * M_PI / Ngrid);
        for(int iI=0; iI<Ngrid; iI++) {
            int I0 = iI <= (Ngrid/2) ? iI : iI - Ngrid;
            const double fx = sinc_unnormed(I0 * M_PI / Ngrid);
            for(int K=0; K<Ngrid/2; K++){
                int K0 = K;
                int i = sqrt(1.0 * I0 * I0 + 1.0 * J0 * J0 + 1.0 * K0 * K0);
                if (i >= Ngrid / 2) continue;
                const double fz = 1/sinc_unnormed(K0 * M_PI / Ngrid);
                double fxyz = fx * fy * fz;
                fxyz *= fxyz;
                /* we revert one extra CIC */
                size_t index = K + (NgridL/2+1)*(iI + NgridL*Jl);
                double a = density_k[index][0] / fxyz;
                double b = density_k[index][1] / fxyz;
                double p = a * a + b * b;
                power[i] += p;                
                count[i] += 1;
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, power, Ngrid / 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, count, Ngrid / 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int i = 0; i < Ngrid / 2; i ++) {
        power[i] /= (count[i]);
    }
}
// Calculate one component of force mesh from precalculated density(k)
void compute_force_mesh(const int axes, fftwf_complex * const P3D)
{
                                                       timer_start(force_mesh);

  fftwf_complex* const FN11= fftdata;

  //k=0 zero mode force is zero
  FN11[0][0]= 0.0f;
  FN11[0][1]= 0.0f;

  const float scale=2.0*M_PI/BoxSize;
  //const float dens_fac= 1.0/(pow(Ngrid, 3));

  const float f1= -1.0/pow(Ngrid, 3.0)/scale;

  //complex float di[3];
  //#shared(FN11, P3D, Local_ny_td, Local_y_start_td, Ngrid, NgridL)

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for(int Jl=0; Jl<Local_ny_td; Jl++) {
      float di[3];
      float di2_sinc[3];

      int J= Jl + Local_y_start_td;
      int J0= J <= (Ngrid/2) ? J : J - Ngrid;
      //di[1]= (complex float) J0;
      const double fy = sinc_unnormed(J0 * M_PI / Ngrid);

      di[1]= diff_kernel(J0 * M_PI * 2.0 / Ngrid) * Ngrid / (M_PI * 2.0);
      di2_sinc[1]= (J0 * fy) * (J0 * fy);

      for(int iI=0; iI<Ngrid; iI++) {
          int I0= iI <= (Ngrid/2) ? iI : iI - Ngrid;
          const double fx = sinc_unnormed(I0 * M_PI / Ngrid);
          //di[0]= (complex float) I0;

          di[0]= diff_kernel(I0 * M_PI * 2.0 / Ngrid) * Ngrid / (M_PI * 2.0);
          di2_sinc[0]= (I0 * fx) * (I0 * fx);

          int KMIN= (iI==0 && J==0); // skip (0,0,0) because FN=0 for k=(0,0,0)

          for(int K=KMIN; K<Ngrid/2+1; K++){
              const double fz = sinc_unnormed(K * M_PI / Ngrid);
              //di[2]= (complex float) K;
              di[2]= diff_kernel(K * M_PI * 2.0 / Ngrid) * Ngrid / (M_PI * 2.0);

              di2_sinc[2]= (K * fz) * (K * fz);
              //float RK =(float) (K*K + I0*I0 + J0*J0);
              const float di2 = di2_sinc[0] + di2_sinc[1] + di2_sinc[2];
              float f2= f1/di2 * di[axes];

              size_t index= K + (NgridL/2+1)*(iI + NgridL*Jl);
              FN11[index][0]= -f2*P3D[index][1];
              FN11[index][1]=  f2*P3D[index][0];
              //complex float dens = -P3D[K + (NgridL/2+1)*(iI + NgridL*Jl)];
              //dens *= dens_fac/RK;

              //#ifdef FILTER
              //	int KK= RK*Scale*Scale; // what is this factor?
              //#else
              //int KK= 1;// ** factors around here can be cleaned/optimized
              //#endif
              //FN11[K + (NgridL/2+1)*(iI + NgridL*Jl)]= dens*I*di[axes]/Scale*KK;

              //FN11[K + (NgridL/2+1)*(iI + NgridL*Jl)]= dens*I*di[axes]/scale;
              //FN11[K + (NgridL/2+1)*(iI + NgridL*Jl)][0]= 0.0f; //dens[0];
          }
      }
  }
  timer_stop(force_mesh);


  //msg_printf("Inverse FFT to Force, axes= %d ...", axes);
  timer_start(fft);
  //fftwf_mpi_execute_dft_c2r(p11, FN11, N11);
  // Force(k) -> Force(x)
  fftwf_mpi_execute_dft_c2r(p11, fftdata, (float*) fftdata);
  timer_stop(fft);
  //msg_printf("done\n");
}

// Does 3-linear interpolation
// particles= Values of mesh at particle positions P.x
void force_at_particle_locations(const Particle P[], const int np, 
				 const int axes, 
		                 const float fmesh[], float3 f[])
{
  const float scaleBox=((float) Ngrid)/((float) BoxSize);
     
  msg_printf(verbose, "Calculating MP...\n");

#ifdef _OPENMP
  #pragma omp parallel for default(shared)     
#endif
  for(int i=0; i<np; i++) {
    float X=P[i].x[0]*scaleBox;
    float Y=P[i].x[1]*scaleBox;
    float Z=P[i].x[2]*scaleBox;
            
    int iI= (int) floorf(X);
    int J=  (int) floorf(Y);
    int K=  (int) floorf(Z);
    float D1= X - (float) iI;
    float D2= Y - (float) J;
    float D3= Z - (float) K;
    float T1= 1.0f - D1;
    float T2= 1.0f - D2;
    float T3= 1.0f - D3;

    if(J >= Ngrid) J=0;
    if(K >= Ngrid) K=0;
            
    int I1=iI+1; if(I1 >= Ngrid) I1=0;
    int J1=J+1; if(J1 >= Ngrid) J1=0;
    int K1=K+1; if(K1 >= Ngrid) K1=0;

    iI -= Local_x_start;
    I1 -= Local_x_start;

    f[i][axes]= 0.0f;

    if(0 <= iI && iI < Local_nx) {
      f[i][axes] += 
	REd(fmesh, iI, J,  K )*T3*T1*T2 +
	REd(fmesh, iI, J,  K1)*D3*T1*T2 +
	REd(fmesh, iI, J1, K )*T3*T1*D2 +
	REd(fmesh, iI, J1, K1)*D3*T1*D2;
    }
    if(0 <= I1 && I1 < Local_nx) {
      f[i][axes] += 
	REd(fmesh, I1, J,  K )*T3*D1*T2 +
	REd(fmesh, I1, J,  K1)*D3*D1*T2 +
	REd(fmesh, I1, J1, K )*T3*D1*D2 +
	REd(fmesh, I1, J1, K1)*D3*D1*D2;
    }
  }
}


void check_total_density(float const * const density)
{
  double sum= 0.0;

  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Ngrid; j++)
      for(int k = 0; k < Ngrid; k++)
	sum += density[(i*NgridL + j)*2*(NgridL/2 + 1) + k];

  double sum_global;
  MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

 
  if(ThisNode == 0) {
    double tol= 1.0e-7*pow(PM_factor, 3)*NParticleTotal;

    if(fabs(sum_global) > tol)
      msg_abort(6040, "Error: total CIC density error is large: %le > %le\n", 
		sum_global, tol);

    msg_printf(debug, 
	      "Total CIC density OK within machine precision: %lf (< %lf).\n",
	       sum_global, tol);


  }
  MPI_Barrier(MPI_COMM_WORLD);
}
  
// Calculates force on particles, particles->f, using particle mesh method
void pm_calculate_forces(Particles* particles, void * mem2, size_t size2)
{

    fftwf_complex * density_k;
    if(mem2 == 0) {
        density_k= fftwf_alloc_complex((NgridL/2+1)*NgridL*Local_ny_td);
    }
    else {
        assert(size2 >= sizeof(fftwf_complex)*(NgridL/2+1)*NgridL*Local_ny_td);
        density_k= (fftwf_complex*) mem2;
    }

    int nghosts = domain_create_ghosts(particles, BoxSize/Ngrid, mem2, size2);
    int np_plus_buffer= particles->np_local + nghosts;

                                                            timer_start(assign);
  // x_i -> density(x) = fftdata
    PtoMesh(particles->p, np_plus_buffer, (float*) fftdata);
                                                            timer_stop(assign);
                                                            timer_start(check);
    check_total_density((float*) fftdata);
                                                            timer_stop(check);


                                                            timer_start(fft);
  // density(x) -> density(k)
    compute_density_k(density_k);
  //fftwf_mpi_execute_dft_r2c(p0, density, P3D);
                                                            timer_stop(fft);

    compute_power_spectrum(density_k);

    for(int axes=0; axes<3; axes++) {
    // density(k) -> f(x_i) [fftdata]
        compute_force_mesh(axes, density_k);

							    timer_start(pforce);
        force_at_particle_locations(particles->p, np_plus_buffer, axes,
				(float*) fftdata, particles->force);
                                                            timer_stop(pforce);
    }

                                                            timer_start(comm);
    domain_annihilate_ghosts(particles, nghosts, particles->force, mem2, size2);
         
                                                            timer_stop(comm);
}

