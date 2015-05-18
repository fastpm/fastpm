/**
 * Particle Mesh gravitational force calculation
 *
 * This code was initially modified by Jun Koda, 
 * from the original serial COLA code
 * by Svetlin Tassev.
 *
 * The CIC kernel remains the orignal form.
 *
 * Yu Feng <rainwoodman@gmail.com>
 ***/

// This file contains some standard functions for a PM code. 
// Nothing COLA-specific.

#include <stdlib.h>
#include <alloca.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "pm.h"
#include "power.h"
#include "msg.h"
#include "timer.h"
#include "domain.h"
#include "heap.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static int Ngrid;
static size_t NgridL; // for index to avoid 4byte integer overflow
static int PM_factor;
static int Local_nx, Local_x_start;
static int Local_ny_td, Local_y_start_td;  // transposed
static size_t total_size;
static float BoxSize;

static fftwf_plan p0, p11;

static double * PowerSpectrumVariable;

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

static inline float WRtPlus(float * const d, 
        const int i, const int j, const int k, const float f)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
    d[k + 2*(NgridL/2 + 1)*(j + NgridL*i)] += f;
    return f;
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

void pm_init(const int nc_pm, const int nc_pm_factor, const float boxsize, int many)
{
    // Assume FFTW3 is initialized in lpt.c

    Ngrid= nc_pm; NgridL= nc_pm;
    PowerSpectrumVariable = malloc(sizeof(double) * Ngrid / 2);
    BoxSize= boxsize;
    PM_factor= nc_pm_factor;

    ptrdiff_t local_nx, local_x_start, local_ny, local_y_start;
    total_size= 
        fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid / 2 + 1, MPI_COMM_WORLD,
                &local_nx, &local_x_start, &local_ny, &local_y_start);

    Local_nx= local_nx;
    Local_x_start= local_x_start;
    Local_ny_td= local_ny;
    Local_y_start_td= local_y_start;


    fftwf_complex * fftdata = heap_allocate(sizeof(fftwf_complex) * total_size);

    msg_printf(verbose, "Setting up FFTW3 plans\n");

    p0=  fftwf_mpi_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, 
            (float*) fftdata, fftdata,
            MPI_COMM_WORLD, (many?FFTW_MEASURE:FFTW_ESTIMATE) | FFTW_MPI_TRANSPOSED_OUT);

    // inverse FFT
    p11= fftwf_mpi_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, fftdata, (float*)fftdata,
            MPI_COMM_WORLD, (many?FFTW_MEASURE:FFTW_ESTIMATE) | FFTW_MPI_TRANSPOSED_IN);

    const int ncp= nc_pm/nc_pm_factor;
    NParticleTotal= (long long) ncp*ncp*ncp;
    heap_return(fftdata);
}

void pm_finalize(void)
{
    fftwf_destroy_plan(p0);
    fftwf_destroy_plan(p11);
}


void PtoMesh(float (*Px)[3], int np, float* const density)
{
    // ** precondition
    //   particles are assumed to be periodiclly wraped up in y,z direction

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

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0; i<np; i++) {
        float X=Px[i][0]*scaleBox;
        float Y=Px[i][1]*scaleBox;
        float Z=Px[i][2]*scaleBox;

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

        int I1=iI+1; while(I1 >= Ngrid) I1-=Ngrid;
        int J1=J+1; while(J1 >= Ngrid) J1-=Ngrid; // assumes y,z < BoxSize
        int K1=K+1; while(K1 >= Ngrid) K1-=Ngrid;

        iI -= Local_x_start;
        I1 -= Local_x_start;

        if(0 <= iI && iI < Local_nx) {
            WRtPlus(density, iI, J,  K,  T3*T1*T2W);
            WRtPlus(density, iI, J,  K1, D3*T1*T2W);
            WRtPlus(density, iI, J1, K,  T3*T1*D2W);
            WRtPlus(density, iI, J1, K1, D3*T1*D2W);
        }

        if(0 <= I1 && I1 < Local_nx) {
            WRtPlus(density, I1, J,  K,  T3*D1*T2W);
            WRtPlus(density, I1, J,  K1, D3*D1*T2W);
            WRtPlus(density, I1, J1, K,  T3*D1*D2W);
            WRtPlus(density, I1, J1, K1, D3*D1*D2W);
        }
    }
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
void compute_density_k(float * density, fftwf_complex * density_k)
{
    fftwf_complex * fftdata = (fftwf_complex*) density;
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
        for(int iI=0; iI<Ngrid; iI++) {
            for(int K=0; K<Ngrid/2+1; K++){
                size_t index = K + (NgridL/2+1)*(iI + NgridL*Jl);
                density_k[index][0] = fftdata[index][0];
                density_k[index][1] = fftdata[index][1];
            }
        }
    }
    /* we check the consistency of CIC here. 
     * delta_k[0] is summed in more accurately in FFTW
     * than simple linear summing.
     * */
    double Norm = 0.0;
    /* remove the mean  */
    if(Local_x_start == 0) {
        Norm = density_k[0][0];
        density_k[0][0] = 0;
        density_k[0][1] = 0;
    }
    MPI_Bcast(&Norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    msg_printf(verbose, "delta_k(0) = %g ; CIC err is = %g\n", 
            Norm, pow(Ngrid, 3), Norm / pow(Ngrid, 3) - 1);
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
    float * ff = (float*) alloca(sizeof(float) * (Ngrid));
    float * i2 = (float*) alloca(sizeof(float) * (Ngrid));
    for(int J = 0; J < Ngrid; J ++) {
        int J0 = J <= (Ngrid/2) ? J : J - Ngrid;
        double tmp = sin(M_PI * 0.5 * J0 / ( 0.5 * Ngrid));
        ff[J] = 1 - 2. / 3 * tmp * tmp;
        i2[J] = (float)J0 * J0;
    }
    for(int Jl=0; Jl<Local_ny_td; Jl++) {
        int J = Jl + Local_y_start_td;
        for(int iI=0; iI<Ngrid; iI++) {
            for(int K=0; K<Ngrid/2; K++){
                int i = sqrt(i2[J] + i2[K] + i2[iI]);
                if (i >= Ngrid / 2) continue;
                const float fx2 = ff[iI];
                const float fy2 = ff[J];
                const float fz2 = ff[K];
                const float fxyz2 = fx2 * fy2 * fz2;

                size_t index = K + (NgridL/2+1)*(iI + NgridL*Jl);
                float a = density_k[index][0];
                float b = density_k[index][1];
                float p = a * a + b * b;
                power[i] += p / (fxyz2);                
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
void compute_force_mesh(const int axes, fftwf_complex * fftdata, fftwf_complex * P3D)
{
    timer_start(force_mesh);

    fftwf_complex * FN11 = fftdata;
    //k=0 zero mode force is zero
    FN11[0][0]= 0.0f;
    FN11[0][1]= 0.0f;

    const float scale=2.0*M_PI/BoxSize;
    //const float dens_fac= 1.0/(pow(Ngrid, 3));

    const float f1= -1.0/pow(Ngrid, 3.0)/scale;

    float * diff = alloca(sizeof(float) * Ngrid);
    float * di2 = alloca(sizeof(float) * Ngrid);
    float * ff = alloca(sizeof(float) * Ngrid);
    for(int J = 0 ; J < Ngrid; J ++) {
        int J0= J <= (Ngrid/2) ? J : J - Ngrid;
        diff[J] = diff_kernel(J0 * M_PI * 2.0 / Ngrid) * Ngrid / (M_PI * 2.0);
        ff[J] = sinc_unnormed(J0 * M_PI / Ngrid);
        di2[J] = J0 * ff[J] * J0 * ff[J];
    }
    //complex float di[3];
    //#shared(FN11, P3D, Local_ny_td, Local_y_start_td, Ngrid, NgridL)

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int Jl=0; Jl<Local_ny_td; Jl++) {
        int J= Jl + Local_y_start_td;
        for(int iI=0; iI<Ngrid; iI++) {
            int KMIN= (iI==0 && J==0); // skip (0,0,0) because FN=0 for k=(0,0,0)

            for(int K=KMIN; K<Ngrid/2+1; K++){
                //float RK =(float) (K*K + I0*I0 + J0*J0);
                const float tmp = di2[J] + di2[iI] + di2[K];
                const int ind[] = {iI, J, K};
                float f2= f1/tmp * diff[ind[axes]];

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
#if 0
    This will dump force into files
        static int step = 0;
    {
        int ThisTask;
        MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
        char buf[1024];
        sprintf(buf, "densitydump-%d-%d.%d", step, axes, ThisTask);
        FILE * fp = fopen(buf, "w");
        fwrite(fftdata, sizeof(float), Local_nx * NgridL * (NgridL /2 + 1) * 2, fp);
        fclose(fp);
        if (axes == 2) step ++;
    }
#endif
}

// Does 3-linear interpolation
// particles= Values of mesh at particle positions P.x
void force_at_particle_locations(float (*Px)[3], const int np, 
        const int axes, 
        const float fmesh[], float (*f)[3])
{
    const float scaleBox=((float) Ngrid)/((float) BoxSize);

    msg_printf(verbose, "Calculating MP...\n");

#ifdef _OPENMP
#pragma omp parallel for default(shared)     
#endif
    for(int i=0; i<np; i++) {
        float X=Px[i][0]*scaleBox;
        float Y=Px[i][1]*scaleBox;
        float Z=Px[i][2]*scaleBox;

        int iI= (int) floorf(X);
        int J=  (int) floorf(Y);
        int K=  (int) floorf(Z);
        float D1= X - (float) iI;
        float D2= Y - (float) J;
        float D3= Z - (float) K;
        float T1= 1.0f - D1;
        float T2= 1.0f - D2;
        float T3= 1.0f - D3;

        while(iI >= Ngrid) iI -= Ngrid;
        while(J >= Ngrid) J -= Ngrid;
        while(K >= Ngrid) K -= Ngrid;

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


// Calculates force on particles, particles->f, using particle mesh method
void pm_calculate_forces(Particles* particles)
{

    int nghosts = domain_create_ghosts(particles, BoxSize/Ngrid);
    int np_plus_buffer= particles->np_local + nghosts;

    fftwf_complex * density_k = heap_allocate(total_size * sizeof(fftwf_complex));
    fftwf_complex * fftdata = heap_allocate(total_size * sizeof(fftwf_complex));

    timer_start(assign);
    // x_i -> density(x) = fftdata
    PtoMesh(particles->x, np_plus_buffer, (float*) fftdata);
    timer_stop(assign);

    timer_start(fft);
    // density(x) -> density(k)
    compute_density_k((float *) fftdata, density_k);

    timer_stop(fft);

    timer_start(powerspectrum);
    compute_power_spectrum(density_k);
    timer_stop(powerspectrum);

        for(int axes=0; axes<3; axes++) {
            // density(k) -> f(x_i) [fftdata]
            compute_force_mesh(axes, fftdata, density_k);

            timer_start(pforce);
            force_at_particle_locations(particles->x, np_plus_buffer, axes,
                    (float*) fftdata, particles->force);
            timer_stop(pforce);
        }

    heap_return(fftdata);
    heap_return(density_k);
    timer_start(comm);
    domain_annihilate_ghosts(particles, nghosts, particles->force);

    timer_stop(comm);
}

