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
#include "mond.h"
#include "version.h"

extern int write_runpb_snapshot(Parameters * param, Particles * snapshot,  
        char * filebase);
extern int read_runpb_ic(Parameters * param, double a_init, Particles * particles);

Particles* allocate_particles(Parameters * param, int allocate_memory);
static size_t calculate_heap_size(Parameters * param);
static void write_powerspectrum(Parameters * param, double a_x, int nc_factor);

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(Parameters * param, float aout, int iout, 
        double a_x, double a_v,
        Particles * particles, 
        Particles * snapshot
        );

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
    const double sigma8= param.sigma8;

    msg_set_loglevel(param.loglevel);

    msg_printf(verbose, "This fPM. Version: %s\n", FPM_VERSION);

    fft_init(multi_thread);

    power_init(param.power_spectrum_filename, param.time_step[0], 
            sigma8, OmegaM, OmegaLambda);

    heap_init(calculate_heap_size(&param));
    Particles * particles = allocate_particles(&param, 1);
    Particles * snapshot = allocate_particles(&param, 0);

    nc_factor = param.pm_nc_factor1;
    pm_set_diff_order(param.diff_order);
    pm_set_poisson_order(param.poisson_order);
    if(param.pm_mond_mode != PM_MOND_NONE) {
        mond_init(&param);
        pm_set_mond(mond_get_modulator(), &param);
    }
    pm_init(param.boxsize, param.nc);
    domain_init(param.boxsize, param.nc);
    lpt_init(param.nc);

    stepping_init(&param);
    int nsteps = stepping_get_nsteps();
    int nout = param.n_zout;

    double* aout = malloc(sizeof(double)*nout);
    for(int i=0; i<nout; i++) {
        aout[i] = (double)(1.0/(1 + param.zout[i]));
        msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
                i, param.zout[i], aout[i]);
    }

    //
    // Many realizations with different initial conditions
    //
    for(int irealization=0; irealization<param.nrealization; irealization++) {
        MPI_Barrier(MPI_COMM_WORLD);
        msg_printf(verbose, "\n%d/%d realization.\n", 
                irealization+1, param.nrealization);

        msg_printf(normal, "Time integration a= %g -> %g, %d steps\n", 
                param.time_step[0], param.time_step[param.n_time_step-1], param.n_time_step);

        int seed= param.random_seed + irealization;

        int iout= 0;

        pm_set_size(nc_factor);
        domain_set_size(nc_factor);

        timer_set_category(LPT);

        particles->dx1 = heap_allocate(sizeof(float) * 3 * particles->np_allocated);
        particles->dx2 = heap_allocate(sizeof(float) * 3 * particles->np_allocated);

        if(param.readic_filename) {
            read_runpb_ic(&param, param.time_step[0], particles);
        } else {
            // Sets initial grid and 2LPT desplacement
            lpt_set_displacement(param.time_step[0], OmegaM, seed,
                    param.boxsize, particles);
        }
        snapshot->seed= seed;

        // always do this because it intializes the initial velocity
        // correctly.
        stepping_set_initial(param.time_step[0], particles);

        if(param.force_mode == FORCE_MODE_PM) {
            heap_return(particles->dx2);
            heap_return(particles->dx1);
        }
        timer_set_category(STEPPING);

        double broadband_init = 0;
        double a_init = 0; 
        //
        // Time evolution loop
        //
        //   TODO: allow nstep=1?
        int chk_change = 0;
        for (int istep = 0; istep <= nsteps; istep++) {
            double a_v, a_x, a_v1, a_x1;

            stepping_get_times(istep,
                &a_x, &a_x1, &a_v, &a_v1);

            if(a_x1 >= change_pm && chk_change != 1 && param.pm_nc_factor2 != param.pm_nc_factor1){

                msg_printf(normal, "Switching to new pm factor: %d->%d\n",
                        param.pm_nc_factor1,
                        param.pm_nc_factor2);

                nc_factor = param.pm_nc_factor2;
                domain_set_size(nc_factor);
                pm_set_size(nc_factor);
                chk_change = 1;
            }

            msg_printf(normal, "Timestep %d/%d a_x:%g a_x1:%g a_v:%g a_v1:%g\n", 
                    istep + 1, nsteps, a_x, a_x1, a_v, a_v1);

            timer_start(comm);
            // move particles to other nodes
            domain_wrap(particles);
            domain_decompose(particles);

            timer_stop(comm);

            if(param.force_mode & FORCE_MODE_PM) {
                mond_set_time(a_x);
                if(param.enforce_broadband && istep > 0) {
                    double growth = GrowthFactor(a_init, a_x);
                    pm_enforce_broadband(broadband_init * growth * growth);
                }
                pm_calculate_forces(particles); 
                if(istep == 0) {
                    broadband_init = pm_get_broadband();
                    a_init = a_x;
                }
                write_powerspectrum(&param, a_x, nc_factor);
                
            }

            while(iout < nout && a_v <= aout[iout] && aout[iout] <= a_x) {
                // Time to write output
                snapshot_time(&param, aout[iout], iout, a_x, a_v, particles, snapshot);
                iout++;
            }

            msg_printf(verbose, "Memory utilization = %td / %td bytes at root rank\n", 
                heap_get_max_usage(),
                heap_get_total_bytes());

            if(iout >= nout) break;

            // Leap-frog "kick" -- velocities updated
            stepping_kick(particles, a_v, a_v1, a_x);

            while(iout < nout && a_x < aout[iout] && aout[iout] <= a_v1) {
                // Time to write output
                snapshot_time(&param, aout[iout], iout, a_x, a_v, particles, snapshot);
                iout++;
            }
            if(iout >= nout) break;

            // Leap-frog "drift" -- positions updated
            stepping_drift(particles, a_x, a_x1, a_v1);
        }
        timer_print();
    }


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


void snapshot_time(Parameters * param, float aout, int iout, 
        double a_x, double a_v,
        Particles * particles, 
        Particles * snapshot)
{
    char filebase[1024];      // TODO: make 256 to variable number...?
    const int isnp= iout+1;

    msg_printf(verbose, "Taking a snapshot...\n");

    snapshot->np_allocated = particles->np_allocated;
    snapshot->x = heap_allocate(sizeof(float) * 3 * particles->np_allocated);
    snapshot->v = heap_allocate(sizeof(float) * 3 * particles->np_allocated);
    snapshot->id = heap_allocate(sizeof(int64_t) * particles->np_allocated);

    timer_set_category(Snp);
    stepping_set_snapshot(aout, a_x, a_v, particles, snapshot);

    const int nc= snapshot->nc; assert(nc > 0);
    const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);

    timer_start(comm);
    domain_wrap(snapshot);

    timer_stop(comm);

    timer_start(write);
    // periodic wrapup not done, what about after fof? what about doing move_particle_min here?
    if(snapshot->filename) {
        sprintf(filebase, "%s%05d_%0.04f.bin", snapshot->filename, snapshot->seed, snapshot->a);
        write_runpb_snapshot(param, snapshot, filebase);
    }
    timer_stop(write);

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

Particles* allocate_particles(Parameters * param, int allocate_memory)
{
    Particles * particles= malloc(sizeof(Particles));
    memset(particles, 0, sizeof(Particles));
    particles->boxsize = param->boxsize;
    particles->omega_m = param->omega_m;
    particles->h = param->h;
    particles->filename = param->snapshot_filename;

    int nc = param->nc;
    particles->nc = nc;
    particles->np_total= ((size_t) nc)*nc*nc;

    int NTask;
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    int np_alloc = (int) (particles->np_total * param->np_alloc_factor / NTask);

    if(allocate_memory) {
        particles->x = heap_allocate(sizeof(float) * 3 *np_alloc);
        particles->v = heap_allocate(sizeof(float) * 3 *np_alloc);
        particles->force = heap_allocate(sizeof(float) * 3 *np_alloc);
        particles->id = heap_allocate(sizeof(int64_t) * np_alloc);

        particles->np_allocated = np_alloc;
    } else {
        particles->np_allocated = 0;
    }


    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    return particles;
}

static size_t calculate_heap_size(Parameters * param) {
    size_t nc = param->nc;
    size_t nm = param->pm_nc_factor2 * nc;
    int NTask;
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    
    size_t np_alloc = param->np_alloc_factor * param->nc * param->nc * param->nc / NTask;
    size_t psize;
    if (param->force_mode == FORCE_MODE_PM) {
        /* dx1 and dx2 are not needed */
        psize = np_alloc * (sizeof(float) * 3 * 3 + 8);
    } else {
        psize = np_alloc * (sizeof(float) * 3 * 5 + 8);
    }
    size_t fftsize = nm * nm * (nm / 2 + 1) * sizeof(float) * 2 * 2 / NTask; 
    size_t snapshotsize = (np_alloc) * ( sizeof(float) * 3 * 2 + 8 + sizeof(int));
    if (snapshotsize > fftsize) {
        fftsize = snapshotsize;
    }
    size_t total = psize + fftsize + 4096 * 128; 
    msg_printf(verbose, "Particles: %td bytes\nFFT: %td bytes for FFT\nTotal: %td\n",
        psize, fftsize, total);
    return total;
}
static void write_powerspectrum(Parameters * param, double a_x, int nc_factor) {
    if(param->measure_power_spectrum_filename) {
        size_t nk = 0;
        double * powerspectrum = pm_compute_power_spectrum(&nk);
        char fname[9999];
        sprintf(fname, "%s-%0.4f.txt", param->measure_power_spectrum_filename, a_x);
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if(myrank == 0) {
            FILE * fp = fopen(fname, "w");
            if (fp == NULL) {
                msg_abort(0020, "Unable to write to %s\n", fname);
            }
            for(int i = 0; i < nk; i ++) {
                fprintf(fp, "%g %g\n", 3.1416 * 2 / param->boxsize * (i + 0.5), 
                        powerspectrum[i] / pow(nc_factor * param->nc, 6) * pow(param->boxsize, 3.0)
                       );
            }
            fclose(fp);
        }
        double sigma8 = TopHatSigma2(8.0, (void*) measured_power, param);
        sigma8 = sigma8 / pow(param->nc * nc_factor, 6.0) * pow(param->boxsize, 3.0);
        sigma8 /= pow(3.1415926 * 2, 3);
        sigma8 = sqrt(sigma8);
        msg_printf(verbose, "Non-Linear sigma8 = %g Linear sigma8 = %g\n", sigma8, param->sigma8 * GrowthFactor(1.0, a_x));
    }
}
