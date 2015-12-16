#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <omp.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

#include "pmpfft.h"
#include "vpm.h"
#include "pmic.h"
#include "pm2lpt.h"
#include "parameters.h"
#include "pmsteps.h"
#include "msg.h"
#include "power.h"
#include "walltime.h"

int read_runpb_ic(Parameters * param, double a_init, PMStore * p, MPI_Comm comm);
int write_runpb_snapshot(double boxsize, double omega_m, PMStore * p, double aa,
        char * filebase, MPI_Comm comm);

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void 
ensure_dir(char * path);

static void
fix_linear_growth(PMStore * p, double correction);

/* Snapshots */
typedef struct {
    MPI_Comm comm;
    PMStore * p;
    char template[1024];    
    int nout;
    double * aout;
    int iout;
    double boxsize;
    double omega_m;
} SNPS;

static int 
snps_interp(SNPS * snps, PMStore * p, double a_x, double a_v, PMStepper * stepper);
static void 
snps_init(SNPS * snps, Parameters * prr, MPI_Comm comm);
static void 
snps_start(SNPS * snps);



/* Useful stuff */
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);

int fastpm(Parameters * prr, MPI_Comm comm) {

    int NTask; 
    int ThisTask;
    MPI_Comm_size(comm, &NTask);
    MPI_Comm_rank(comm, &ThisTask);
    struct ClockTable CT;
    msg_init(comm);
    walltime_init(&CT);

    msg_set_loglevel(verbose);

    const double rho_crit = 27.7455;
    const double M0 = prr->omega_m*rho_crit*pow(prr->boxsize / prr->nc, 3.0);
    msg_printf(verbose, "mass of a particle is %g 1e10 Msun/h\n", M0); 

    PMStore pdata;
    VPM * vpm_list;
    PMInit baseinit = {
            .Nmesh = prr->nc,
            .BoxSize = prr->boxsize, .NprocY = prr->NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = prr->UseFFTW,
        };

    power_init(prr->power_spectrum_filename, 
            prr->time_step[0], 
            prr->sigma8, 
            prr->omega_m, 
            1 - prr->omega_m, comm);

    pm_store_init(&pdata);

    double alloc_factor;
    alloc_factor = prr->np_alloc_factor;
    msg_printf(info, "Using alloc factor of %g\n", alloc_factor);

    pm_store_alloc_evenly(&pdata, pow(prr->nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC, 
        alloc_factor, comm);

    walltime_measure("/Init/Misc");

    vpm_list = vpm_create(prr->n_pm_nc_factor, 
                           prr->pm_nc_factor, 
                           prr->change_pm,
                           &baseinit, &pdata.iface, comm);

    walltime_measure("/Init/Plan");

    if(prr->readic_filename) {
        read_runpb_ic(prr, prr->time_step[0], &pdata, comm);
        walltime_measure("/Init/ReadIC");
    } else {
        double shift[3] = {
            prr->boxsize / prr->nc * 0.5,
            prr->boxsize / prr->nc * 0.5,
            prr->boxsize / prr->nc * 0.5,
            };

        PM pm;

        pm_init_simple(&pm, &pdata, prr->nc, prr->boxsize, comm);

        pm_store_set_lagrangian_position(&pdata, &pm, shift);

        FastPMFloat * delta_k = pm_alloc(&pm);

        if(prr->readnoise_filename) {
            pmic_read_gaussian(&pm, delta_k, prr->readnoise_filename, PowerSpecWithData, NULL);
        } else {
            pmic_fill_gaussian_gadget(&pm, delta_k, prr->random_seed, PowerSpecWithData, NULL);
        }

        /* read out values at locations with an inverted shift */
        pm_2lpt_solve(&pm, delta_k, &pdata, shift);

        pm_free(&pm, delta_k);
        pm_destroy(&pm);

        walltime_measure("/Init/2LPT");
    }

    pm_2lpt_evolve(prr->time_step[0], &pdata, prr->omega_m);

    walltime_measure("/Init/Drift");

    SNPS snps;

    snps_init(&snps, prr, comm);

    snps_start(&snps);

    PMStepper stepper;
    stepping_init(&stepper, prr->omega_m, prr->force_mode, prr->cola_stdda);

    int istep;
    int nsteps = prr->n_time_step;

    snps_interp(&snps, &pdata, prr->time_step[0], prr->time_step[0], &stepper);

    walltime_measure("/Init/Start");

    double Plin0 = 0;
    /* The last step is the 'terminal' step */

    for (istep = 0; istep < nsteps; istep++) {
        double a_v, a_x, a_v1, a_x1;

        /* begining and ending of drift(x) and kick(v)*/
        pm_get_times(istep, prr->time_step, prr->n_time_step,
            &a_x, &a_x1, &a_v, &a_v1);

        /* Find the Particle Mesh to use for this time step */
        VPM * vpm = vpm_find(vpm_list, a_x);
        PM * pm = &vpm->pm;

        /* watch out: boost the density since mesh is finer than grid */
        double density_factor =  pow(vpm->pm_nc_factor, 3); 

        msg_printf(normal, "==== Step %d a_x = %6.4f a_x1 = %6.4f a_v = %6.4f a_v1 = %6.4f Nmesh = %d ====\n", 
                    istep, a_x, a_x1, a_v, a_v1, pm->init.Nmesh);

        PowerSpectrum ps;

        walltime_measure("/Stepping/Start");

        /* apply periodic boundary and move particles to the correct rank */
        pm_store_wrap(&pdata, pm->BoxSize);
        walltime_measure("/Stepping/Periodic");

        pm_store_decompose(&pdata, to_rank, pm, comm);

        size_t np_max;
        size_t np_min;
        double np_mean = pow(prr->nc, 3) / NTask;
        MPI_Allreduce(&pdata.np, &np_max, 1, MPI_LONG, MPI_MAX, comm);
        MPI_Allreduce(&pdata.np, &np_min, 1, MPI_LONG, MPI_MIN, comm);

        msg_printf(info, "Load imbalance is - %g / + %g\n",
            np_min / np_mean, np_max / np_mean);

        walltime_measure("/Stepping/Decompose");

        /* Calculate PM forces. */
        FastPMFloat * delta_k = pm_alloc(pm);

        pm_calculate_forces(&pdata, pm, delta_k, density_factor);

        /* calculate the power spectrum */
        power_spectrum_init(&ps, pm->Nmesh[0] / 2);

        pm_calculate_powerspectrum(pm, delta_k, &ps);

        double Plin = pm_calculate_linear_power(pm, delta_k, prr->enforce_broadband_kmax);
        Plin /= pow(stepper_get_growth_factor(&stepper, a_x), 2.0);
        if(istep == 0) {
            Plin0 = Plin;
        }
        double correction = sqrt(Plin0 / Plin);
        if(!prr->enforce_broadband) correction = 1.0;
        msg_printf(info, "<P(k<%g)> = %g Linear Theory = %g, correction=%g\n", 
            prr->enforce_broadband_kmax,
            Plin,
            Plin0,
            correction
        ); 
        fix_linear_growth(&pdata, correction);

        walltime_measure("/PowerSpectrum/Measure");
        
        if(prr->measure_power_spectrum_filename) {
            if(pm->ThisTask == 0) {
                ensure_dir(prr->measure_power_spectrum_filename);
                power_spectrum_write(&ps, pm, ((double)prr->nc * prr->nc * prr->nc), 
                    prr->measure_power_spectrum_filename, prr->random_seed, a_x);
            }
        }

        MPI_Barrier(comm);
        walltime_measure("/PowerSpectrum/Write");
        power_spectrum_destroy(&ps);

        pm_free(pm, delta_k);

        /* take snapshots if needed, before the kick */
        snps_interp(&snps, &pdata, a_x, a_v, &stepper);

        /* never go beyond 1.0 */
        if(a_x >= 1.0) break; 
        
        // Leap-frog "kick" -- velocities updated

        stepping_kick(&stepper, &pdata, &pdata, a_v, a_v1, a_x);
        walltime_measure("/Stepping/kick");

        /* take snapshots if needed, before the drift */
        snps_interp(&snps, &pdata, a_x, a_v1, &stepper);
        
        // Leap-frog "drift" -- positions updated
        stepping_drift(&stepper, &pdata, &pdata, a_x, a_x1, a_v1);
        walltime_measure("/Stepping/drift");

        /* no need to check for snapshots here, it will be checked next loop.  */
    }

    pm_store_destroy(&pdata);

    msg_printf(info, "Total Time\n");
    walltime_summary(0, comm);
    walltime_report(stdout, 0, comm);

    pfft_cleanup();
    return 0;
}

static int 
to_rank(void * pdata, ptrdiff_t i, void * data) 
{
    PMStore * p = (PMStore *) pdata;
    PM * pm = (PM*) data;
    double pos[3];
    p->iface.get_position(p, i, pos);
    return pm_pos_to_rank(pm, pos);
}


#if 0
static void 
smooth_density(PM * pm, double r_s) 
{
    /* 
     *  This function smooth density by scale r_s. There could be a factor of sqrt(2)
     *  It is not used. */

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);
    {
        /* fill in the extra 'smoothing kernels' we will take the product */
        ptrdiff_t ind;
        int d;
        for(d = 0; d < 3; d++)
        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            fac[d][ind].extra = exp(-0.5 * fac[d][ind].kk * r_s * r_s);
        }
    }

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double smth = 1.;
            double kk = 0.;
            for(d = 0; d < 3; d++) {
                smth *= fac[d][i[d] + pm->ORegion.start[d]].extra;
                kk += fac[d][i[d] + pm->ORegion.start[d]].kk;
            }
            /* - i k[d] / k2 */
            if(LIKELY(kk> 0)) {
                pm->workspace[ind + 0] = pm->canvas[ind + 0] * smth;
                pm->workspace[ind + 1] = pm->canvas[ind + 1] * smth;
            } else {
                pm->workspace[ind + 0] = 0;
                pm->workspace[ind + 1] = 0;
            }
            pm_inc_o_index(pm, i);
        }
    }

    pm_destroy_k_factors(pm, fac);
}
#endif

static int
snps_interp(SNPS * snps, PMStore * p, double a_x, double a_v, PMStepper * stepper)
{
    /* interpolate and write snapshots, assuming p 
     * is at time a_x and a_v. */
    PMStore snapshot;

    while(snps->iout < snps->nout && (
        /* after a kick */
        (a_x < snps->aout[snps->iout] && snps->aout[snps->iout] < a_v)
        ||
        /* after a drift */
        (a_x >= snps->aout[snps->iout] && snps->aout[snps->iout] >= a_v)
        )) {

        pm_store_init(&snapshot);

        pm_store_alloc(&snapshot, p->np_upper, PACK_ID | PACK_POS | PACK_VEL);

        msg_printf(verbose, "Taking a snapshot...\n");

        double aout = snps->aout[snps->iout];
        int isnp= snps->iout+1;

        stepping_set_snapshot(stepper, p, &snapshot, aout, a_x, a_v);
        walltime_measure("/Snapshot/KickDrift");

        char filebase[1024];
        sprintf(filebase, snps->template, aout);
        ensure_dir(filebase);
        write_runpb_snapshot(snps->boxsize, stepper->omega_m, &snapshot, aout, filebase, snps->comm);
        
        walltime_measure("/Snapshot/IO");

        MPI_Barrier(snps->comm);
        walltime_measure("/Snapshot/Wait");

        const double z_out= 1.0/aout - 1.0;

        msg_printf(normal, "snapshot %d written z = %6.4f a = %6.4f\n", 
                isnp, z_out, aout);

        snps->iout ++;
        pm_store_destroy(&snapshot);
    }
    return (snps->iout == snps->nout);
}

static void 
snps_init(SNPS * snps, Parameters * prr, MPI_Comm comm) 
{
    snps->iout = 0;
    snps->nout = prr->n_zout;
    snps->comm = comm;

    snps->boxsize = prr->boxsize;
    snps->omega_m = prr->omega_m;

    snps->aout = malloc(sizeof(double)*snps->nout);

    sprintf(snps->template, "%s%05d_%%0.04f.bin", prr->snapshot_filename, prr->random_seed);

    int i;
    for(i=0; i<snps->nout; i++) {
        snps->aout[i] = (double)(1.0/(1 + prr->zout[i]));
        msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
                i, prr->zout[i], snps->aout[i]);
    }
}

static void 
snps_start(SNPS * snps) 
{
    snps->iout = 0;
}


static void 
_mkdir(const char *dir) 
{
        char * tmp = strdup(dir);
        char *p = NULL;
        size_t len;

        len = strlen(tmp);
        if(tmp[len - 1] == '/')
                tmp[len - 1] = 0;
        for(p = tmp + 1; *p; p++)
                if(*p == '/') {
                        *p = 0;
                        mkdir(tmp, S_IRWXU | S_IRWXG | S_IRWXO);
                        *p = '/';
                }
        mkdir(tmp, S_IRWXU | S_IRWXG | S_IRWXO);
        free(tmp);
}
static void 
ensure_dir(char * path) 
{
    int i = strlen(path);
    char * dup = strdup(path);
    char * p;
    for(p = i + dup; p >= dup && *p != '/'; p --) {
        continue;
    }
    /* plain file name in current directory */
    if(p < dup) return;
    
    /* p == '/', so set it to NULL, dup is the dirname */
    *p = 0;
    _mkdir(dup);
    free(dup);
}

static void
fix_linear_growth(PMStore * p, double correction)
{
    ptrdiff_t i;
    int d;
    for(d = 0; d < 3; d ++) {
        for(i = 0; i < p->np; i ++) {
            p->acc[d][i] *= correction;
        }
    }
}
