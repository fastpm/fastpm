/**
 * main file: Time stepping and firing-sequences
 * 
 * This code was initially modified by Jun Koda, 
 * from the original serial COLA code
 * by Svetlin Tassev.
 *
 * The time stepping code has been largely re-written
 * to accompodate PM simulations.
 *
 * The firing sequences of snapshots and FOF is preserved.
 * 
 * Adpative PM size is added by Man-Yat Chu.
 * 
 * Yu Feng <rainwoodman@gmail.com>
 *
 */
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
#include "power.h"
#include "domain.h"
#include "pm.h"
#include "stepping.h"
#include "timer.h"
#include "heap.h"

extern int write_runpb_snapshot(Particles * snapshot,  
        char * filebase);
extern int read_runpb_ic(Parameters * param, double a_init, Particles * particles);

Particles* allocate_snapshot(const int nc, const int np_alloc);
Particles* allocate_particles(const int nc, const int nx, const double np_alloc_factor);

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
        double a_x, double a_v,
        Particles * particles, 
        Particles * snapshot
        );

static double polval(const double * pol, const int degrees, const double x) {
    double rt = 0;
    int i;
    for(i = 0; i < degrees; i ++) {
        rt += pow(x, degrees - i - 1) * pol[i];
    }
    return rt;
}

static double measured_power(double k, Parameters * param) {
    size_t nk;
    double * power = pm_compute_power_spectrum(&nk);
    double k0 = 2 * 3.1415926 / param->boxsize;
    int ki = floor(k / k0);
    double u = k / k0 - ki;

    if(ki >= nk / 2 - 1) return 0;
    if(ki < 0) return 0;
    return power[ki] * (1-u) + power[ki] * u; // other corrections outsize.
}

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
    int nc_factor = -1;
    const float change_pm= param.change_pm;
    const double OmegaM= param.omega_m;
    const double OmegaLambda= 1.0 - OmegaM;
    const double Hubble= param.h;
    const double sigma8= param.sigma8;

    msg_set_loglevel(param.loglevel);

    fft_init(multi_thread);
    const int nsteps= param.ntimestep;
    const double a_final= param.a_final;

    double a_init = 0.1;

    /* one extra item in the end; to avoid an if conditioni in main loop */
    double * A_X = (double*) calloc(nsteps + 2, sizeof(double));
    double * A_V = (double*) calloc(nsteps + 2, sizeof(double));

    switch(param.time_step) {
        case TIME_STEP_LOGA:
            msg_printf(info, "timestep linear in loga\n");
            for(int i = 0; i < nsteps + 2; i ++) {
                double dloga = (log(a_final) - log(a_init)) / nsteps;
                A_X[i] = exp(log(a_init) + i * dloga);
                A_V[i] = exp(log(a_init) + i * dloga - 0.5 * dloga);
            }
            break;
        case TIME_STEP_A:
            msg_printf(info, "timestep linear in a\n");
            for(int i = 0; i < nsteps + 2; i ++) {
                double da = (a_final - a_init) / nsteps;
                A_X[i] = a_init + (i * da);
                A_V[i] = a_init + (i * da - 0.5 * da);
            }
            break;
        case TIME_STEP_GROWTH:
            msg_printf(info, "timestep linear in growth factor\n");
            for(int i = 0; i < nsteps + 2; i ++) {
                /* best fit coefficients for unnormalized wmap7 cosmology dplus
                 * shoudn't matter as long as we are close */
                const double a2d[] = {0.35630564,-1.58170455, 0.34493492, 3.65626615, 0};
                const double d2a[] = {0.00657966,-0.01088151, 0.00861105, 0.27018028, 0};
                double dfinal = polval(a2d, 5, a_final);
                double dinit = polval(a2d, 5, a_init);
                double dd = (dfinal - dinit) / nsteps;
                A_X[i] = polval(d2a, 5, dinit + i * dd);
                A_V[i] = polval(d2a, 5, dinit + i * dd - 0.5 * dd);
            }
            break;
        default:
            abort();

    }
    /* fix up the IC time */
    A_V[0] = A_X[0] = a_init;
    A_X[nsteps] = a_final;
    msg_printf(normal, "Drift points: \n");
    for(int i = 0; i < nsteps + 2; i ++) {
        msg_printf(normal, "%g, ", A_X[i]);
    }
    msg_printf(normal, "\n");

    power_init(param.power_spectrum_filename, a_init, 
            sigma8, OmegaM, OmegaLambda);


    heap_init(0);
    lpt_init(param.nc);
    const int local_nx= lpt_get_local_nx();

    Particles* particles= allocate_particles(param.nc, local_nx, param.np_alloc_factor);
    Particles* snapshot= allocate_snapshot(param.nc, particles->np_allocated);

    snapshot->boxsize= param.boxsize;
    snapshot->omega_m= OmegaM;
    snapshot->h= Hubble;
    //strncpy(snapshot->filename, param.snapshot_filename, 64);
    snapshot->filename= param.snapshot_filename;

    nc_factor = param.pm_nc_factor1;
    pm_set_diff_order(param.diff_order);
    pm_init(nc_factor*param.nc, nc_factor, param.boxsize, param.nrealization>1);

    const int nout= param.n_zout;
    double* aout= malloc(sizeof(double)*nout);
    for(int i=0; i<nout; i++) {
        aout[i] = (double)(1.0/(1 + param.zout[i]));
        msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
                i, param.zout[i], aout[i]);
    }

    if (param.qpm) {
        stepping_set_subtract_lpt(0);
    } else {
        stepping_set_subtract_lpt(1);
    }
    if (param.stdda) {
        stepping_set_std_da(1);
    } else {
        stepping_set_std_da(0);
    }
    if (param.nopm) {
        stepping_set_no_pm(1);
    } else {
        stepping_set_no_pm(0);
    }
    //
    // Many realizations with different initial conditions
    //
    for(int irealization=0; irealization<param.nrealization; irealization++) {
        MPI_Barrier(MPI_COMM_WORLD);
        msg_printf(verbose, "\n%d/%d realization.\n", 
                irealization+1, param.nrealization);
        int seed= param.random_seed + irealization;

        int iout= 0;
        domain_init(param.nc * nc_factor, param.boxsize);

        timer_set_category(LPT);
        if(param.readic_filename) {
            read_runpb_ic(&param, a_init, particles);
        } else {
            // Sets initial grid and 2LPT desplacement
            lpt_set_displacement(a_init, OmegaM, seed,
                    param.boxsize, particles);
        }
        snapshot->seed= seed;

        // always do this because it intializes the initial velocity
        // correctly.
        stepping_set_initial(a_init, OmegaM, particles);

        timer_set_category(STEPPING);

        //
        // Time evolution loop
        //
        //   TODO: allow nstep=1?
        if(nout > 0 && nsteps > 1 && a_final > a_init) {
            msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
                    a_init, a_final, nsteps);

            int chk_change = 0;
            for (int istep=0; istep<= nsteps; istep++) {
                double a_v, a_x, a_v1, a_x1;


                a_v = A_V[istep];
                a_x = A_X[istep];
                a_v1= A_V[istep + 1];
                a_x1= A_X[istep + 1];

                if(a_x1 >= change_pm && chk_change != 1 && param.pm_nc_factor2 != param.pm_nc_factor1){

                    msg_printf(normal, "Switching to new pm factor: %d->%d\n",
                            param.pm_nc_factor1,
                            param.pm_nc_factor2);

                    nc_factor = param.pm_nc_factor2;
                    pm_finalize();
                    domain_finalize();

                    domain_init(param.nc * nc_factor, param.boxsize);
                    pm_init(nc_factor*param.nc, nc_factor, param.boxsize, param.nrealization>1);
                    chk_change = 1;
                }

                msg_printf(normal, "Timestep %d/%d\n", istep + 1, nsteps);

                timer_start(comm);
                // move particles to other nodes
                domain_wrap(particles);
                domain_decompose(particles);

                timer_stop(comm);

                pm_calculate_forces(particles); 

                if(param.measure_power_spectrum_filename) {
                    size_t nk = 0;
                    double * powerspectrum = pm_compute_power_spectrum(&nk);
                    char fname[9999];
                    sprintf(fname, "%s-%0.4f.txt", param.measure_power_spectrum_filename, a_x);
                    int myrank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
                    if(myrank == 0) {
                        FILE * fp = fopen(fname, "w");
                        if (fp == NULL) {
                            msg_abort(0020, "Unable to write to %s\n", fname);
                        }
                        for(int i = 0; i < nk; i ++) {
                            fprintf(fp, "%g %g\n", 3.1416 * 2 / param.boxsize * (i + 0.5), 
                                    powerspectrum[i] / pow(nc_factor * param.nc, 6) * pow(param.boxsize, 3.0)
                                   );
                        }
                        fclose(fp);
                    }
                    double sigma8 = TopHatSigma2(8.0, (void*) measured_power, &param);
                    sigma8 = sigma8 / pow(param.nc * nc_factor, 6.0) * pow(param.boxsize, 3.0);
                    sigma8 /= pow(3.1415926 * 2, 3);
                    sigma8 = sqrt(sigma8);
                    msg_printf(verbose, "sigma8 = %g expected = %g\n", sigma8, param.sigma8 * GrowthFactor(1.0, a_x));
                }

                while(iout < nout && a_v <= aout[iout] && aout[iout] <= a_x) {
                    // Time to write output
                    snapshot_time(aout[iout], iout, a_x, a_v, particles, snapshot);
                    iout++;
                }
                if(iout >= nout) break;

                // Leap-frog "kick" -- velocities updated
                stepping_kick(particles, OmegaM, a_v, a_v1, a_x);

                while(iout < nout && a_x < aout[iout] && aout[iout] <= a_v1) {
                    // Time to write output
                    snapshot_time(aout[iout], iout, a_x, a_v, particles, snapshot);
                    iout++;
                }
                if(iout >= nout) break;

                // Leap-frog "drift" -- positions updated
                stepping_drift(particles, OmegaM, a_x, a_x1, a_v1);
                msg_printf(verbose, "Max memory = %td bytes\n", heap_get_max_usage());
            }
        }
        timer_print();
    }

    //move_particles(particles);


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


void snapshot_time(const float aout, const int iout, 
        double a_x, double a_v,
        Particles * particles, 
        Particles * snapshot)
{
    // Halo finding and snapshot outputs

    char filebase[1024];      // TODO: make 256 to variable number...?
    const int isnp= iout+1;

    msg_printf(verbose, "Taking a snapshot...\n");

    snapshot->x = heap_allocate(sizeof(float) * 3 * particles->np_allocated);
    snapshot->v = heap_allocate(sizeof(float) * 3 * particles->np_allocated);
    snapshot->id = heap_allocate(sizeof(int64_t) * particles->np_allocated);

    timer_set_category(Snp);
    stepping_set_snapshot(aout, a_x, a_v, particles, snapshot);

    const int nc= snapshot->nc; assert(nc > 0);
    const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);

    timer_start(comm);
    domain_wrap(snapshot);
    domain_decompose(snapshot);
    timer_stop(comm);

    timer_start(write);
    // Gadget snapshot for all particles
    // periodic wrapup not done, what about after fof? what about doing move_particle_min here?
    if(snapshot->filename) {
        sprintf(filebase, "%s%05d_%0.04f.bin", snapshot->filename, snapshot->seed, snapshot->a);
        write_runpb_snapshot(snapshot, filebase);
    }
    timer_stop(write);
    // text file of a slice
    //sprintf(filebase, "slice%02d", isnp); temp
    //write_snapshot_slice(filebase, snapshot, boxsize); temp

    const double rho_crit = 27.7455;
    const double M0 = snapshot->omega_m*rho_crit*pow(snapshot->boxsize / snapshot->nc, 3.0);
    msg_printf(verbose, "mass of a particle is %g 1e10 Msun/h\n", M0); 

    const double z_out= 1.0/aout - 1.0;
    heap_return(snapshot->id);
    heap_return(snapshot->v);
    heap_return(snapshot->x);

    msg_printf(normal, "snapshot %d written z=%4.2f a=%5.3f\n", 
            isnp, z_out, aout);
    timer_set_category(STEPPING);
}

Particles* allocate_particles(const int nc, const int nx, const double np_alloc_factor)
{
    Particles* particles= malloc(sizeof(Particles));

    const int np_alloc= (int)(np_alloc_factor*nc*nc*(nx));

    particles->x= heap_allocate(sizeof(float) * 3 *np_alloc);
    particles->v= heap_allocate(sizeof(float) * 3 *np_alloc);
    particles->force= heap_allocate(sizeof(float) * 3 *np_alloc);
    particles->dx1= heap_allocate(sizeof(float) * 3 *np_alloc);
    particles->dx2= heap_allocate(sizeof(float) * 3 *np_alloc);
    particles->id= heap_allocate(sizeof(int64_t) * 3 *np_alloc);

    particles->force= heap_allocate(sizeof(float)*3*np_alloc);

    particles->np_allocated= np_alloc;

    particles->np_total= (long long) nc*nc*nc;

    int NTask; 
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    particles->np_average= (float)(pow((double) nc, 3) / NTask);

    return particles;
}

Particles * allocate_snapshot(const int nc, const int np_alloc)
{
    Particles * snapshot= malloc(sizeof(Particles));
    memset(snapshot, 0, sizeof(Particles)); 
    snapshot->np_allocated= np_alloc;
    long long nc_long= nc;
    snapshot->np_total= nc_long*nc_long*nc_long;
    snapshot->nc= nc;
    snapshot->a= 0.0f; //snapshot->a_v= 0.0f; snapshot->a_x= 0.0f;

    return snapshot;
}
