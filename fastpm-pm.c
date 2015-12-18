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

#include "parameters.h"
#include "fastpm-internal.h"
#include "fastpm-pm.h"

void 
fastpm_init(FastPM * fastpm, Parameters * prr, MPI_Comm comm);

void 
fastpm_destroy(FastPM * fastpm);

/*
static void DUMP(PM * pm, char * filename, FastPMFloat *data) {
    char fn2[1024];
    if(pm->NTask > 1) {
        sprintf(fn2, "%s.%03d", filename, pm->ThisTask);
        printf("%d: %td %td %td\n", pm->ThisTask,
                    pm->IRegion.strides[0],
                    pm->IRegion.strides[1],
                    pm->IRegion.strides[2]);
    } else {
        sprintf(fn2, "%s", filename);
    }
    FILE * fp = fopen(fn2, "w");
    fwrite(data, sizeof(FastPMFloat), pm->allocsize, fp);
    fclose(fp);
}
*/

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, void * template);

static void 
ensure_dir(char * path);

static void
fix_linear_growth(PMStore * p, double correction);

/* Useful stuff */
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);

void fastpm_init(FastPM * fastpm, Parameters * prr, MPI_Comm comm) {
    PMInit baseinit = {
            .Nmesh = prr->nc,
            .BoxSize = prr->boxsize, .NprocY = prr->NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = prr->UseFFTW,
        };
    double alloc_factor = prr->np_alloc_factor;

    fastpm->comm = comm;
    MPI_Comm_rank(comm, &fastpm->ThisTask);
    MPI_Comm_size(comm, &fastpm->NTask);

    fastpm->nc = prr->nc;
    fastpm->boxsize = prr->boxsize;
    fastpm->omega_m = prr->omega_m;

    fastpm->p= malloc(sizeof(PMStore));
    fastpm->pm_2lpt = malloc(sizeof(PM));
    fastpm->time_step = prr->time_step;
    fastpm->n_time_step = prr->n_time_step;
    fastpm->USE_COLA = prr->force_mode == FORCE_MODE_COLA;
    fastpm->USE_NONSTDDA = !prr->cola_stdda;
    fastpm->USE_LINEAR_THEORY = prr->enforce_broadband;
    fastpm->nLPT = -2.5f;
    fastpm->K_LINEAR = prr->enforce_broadband_kmax;

    pm_store_init(fastpm->p);

    pm_store_alloc_evenly(fastpm->p, pow(1.0 * fastpm->nc, 3), 
        PACK_POS | PACK_VEL | PACK_ID | PACK_DX1 | PACK_DX2 | PACK_ACC, 
        alloc_factor, comm);

    fastpm->vpm_list = vpm_create(prr->n_pm_nc_factor, 
                           prr->pm_nc_factor, 
                           prr->change_pm,
                           &baseinit, &fastpm->p->iface, comm);

    pm_init_simple(fastpm->pm_2lpt, fastpm->p, fastpm->nc, fastpm->boxsize, comm);

    walltime_measure("/Init/Plan");
}

void 
fastpm_prepare_ic(FastPM * fastpm, FastPMFloat * delta_k) 
{
    double shift[3] = {
        fastpm->boxsize / fastpm->nc * 0.5,
        fastpm->boxsize / fastpm->nc * 0.5,
        fastpm->boxsize / fastpm->nc * 0.5,
        };

    pm_store_set_lagrangian_position(fastpm->p, fastpm->pm_2lpt, shift);

    /* read out values at locations with an inverted shift */
    pm_2lpt_solve(fastpm->pm_2lpt, delta_k, fastpm->p, shift);
}

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

    int istep;
    int nsteps = prr->n_time_step;

    FastPM * fastpm = alloca(sizeof(FastPM));

    fastpm_init(fastpm, prr, comm);

    walltime_measure("/Init/Misc");

    if(prr->readic_filename) {
        fastpm_read_runpb_ic(fastpm, fastpm->p, prr->readic_filename);

        walltime_measure("/Init/ReadIC");
    } else {
        FastPMFloat * delta_k = pm_alloc(fastpm->pm_2lpt);

        power_init(prr->power_spectrum_filename, 
                prr->time_step[0], 
                prr->sigma8, 
                prr->omega_m, 
                1 - prr->omega_m, comm);

        if(prr->readnoise_filename) {
            pmic_read_gaussian(fastpm->pm_2lpt, delta_k, prr->readnoise_filename, PowerSpecWithData, NULL);
        } else {
            pmic_fill_gaussian_gadget(fastpm->pm_2lpt, delta_k, prr->random_seed, PowerSpecWithData, NULL);
        }
        fastpm_prepare_ic(fastpm, delta_k);

        pm_free(fastpm->pm_2lpt, delta_k);

        walltime_measure("/Init/2LPT");
    }

    pm_2lpt_evolve(fastpm->time_step[0], fastpm->p, fastpm->omega_m);

    walltime_measure("/Init/Evolve");

    char TEMPLATE[1024];
    sprintf(TEMPLATE, "%s%05d_%%0.04f.bin", prr->snapshot_filename, prr->random_seed);
    fastpm_interp(fastpm, prr->aout, prr->n_aout, take_a_snapshot, TEMPLATE);

    walltime_measure("/Init/Start");

    double Plin0 = 0;
    /* The last step is the 'terminal' step */

    for (istep = 0; istep < nsteps; istep++) {
        double a_v, a_x, a_v1, a_x1;

        /* begining and ending of drift(x) and kick(v)*/
        fastpm_set_time(fastpm, istep, 
                    &a_x, &a_x1, &a_v, &a_v1);

        msg_printf(normal, "==== Step %d a_x = %6.4f a_x1 = %6.4f a_v = %6.4f a_v1 = %6.4f Nmesh = %d ====\n", 
                    istep, a_x, a_x1, a_v, a_v1, fastpm->pm->init.Nmesh);

        walltime_measure("/Stepping/Start");

        fastpm_decompose(fastpm);

        walltime_measure("/Stepping/Decompose");

        /* Calculate PM forces. */
        FastPMFloat * delta_k = pm_alloc(fastpm->pm);

        /* watch out: boost the density since mesh is finer than grid */
        double density_factor = fastpm->pm->Norm / pow(1.0 * fastpm->nc, 3);

        pm_calculate_forces(fastpm->p, fastpm->pm, delta_k, density_factor);

        PowerSpectrum ps;
        /* calculate the power spectrum */
        power_spectrum_init(&ps, fastpm->pm->Nmesh[0] / 2);

        pm_calculate_powerspectrum(fastpm->pm, delta_k, &ps);

        walltime_measure("/PowerSpectrum/Measure");

        if(prr->measure_power_spectrum_filename) {
            if(fastpm->ThisTask == 0) {
                ensure_dir(prr->measure_power_spectrum_filename);
                power_spectrum_write(&ps, fastpm->pm, ((double)prr->nc * prr->nc * prr->nc), 
                    prr->measure_power_spectrum_filename, prr->random_seed, a_x);
            }
        }
        power_spectrum_destroy(&ps);
        MPI_Barrier(comm);
        walltime_measure("/PowerSpectrum/Write");

        double Plin = pm_calculate_linear_power(fastpm->pm, delta_k, fastpm->K_LINEAR);

        Plin /= pow(fastpm_growth_factor(fastpm, a_x), 2.0);
        if(istep == 0) {
            Plin0 = Plin;
        }

        double correction = sqrt(Plin0 / Plin);

        if(!fastpm->USE_LINEAR_THEORY) correction = 1.0;
        msg_printf(info, "<P(k<%g)> = %g Linear Theory = %g, correction=%g\n", 
                          fastpm->K_LINEAR, Plin, Plin0, correction); 
        fix_linear_growth(fastpm->p, correction);

        pm_free(fastpm->pm, delta_k);

        /* take snapshots if needed, before the kick */
        fastpm_interp(fastpm, prr->aout, prr->n_aout, take_a_snapshot, TEMPLATE);

        /* never go beyond 1.0 */
        if(a_x >= 1.0) break; 
        
        // Leap-frog "kick" -- velocities updated

        fastpm_kick(fastpm, fastpm->p, fastpm->p, a_v1);
        walltime_measure("/Stepping/kick");

        /* take snapshots if needed, before the drift */
        fastpm_interp(fastpm, prr->aout, prr->n_aout, take_a_snapshot, TEMPLATE);
        //
        // Leap-frog "drift" -- positions updated
        fastpm_drift(fastpm, fastpm->p, fastpm->p, a_x1);
        walltime_measure("/Stepping/drift");

        /* no need to check for snapshots here, it will be checked next loop.  */
    }
    fastpm_destroy(fastpm);

    msg_printf(info, "Total Time\n");
    walltime_summary(0, comm);
    walltime_report(stdout, 0, comm);

    pfft_cleanup();
    return 0;
}
void 
fastpm_destroy(FastPM * fastpm) 
{
    pm_store_destroy(fastpm->p);
    /* FIXME: free VPM and stuff. */
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

void
fastpm_decompose(FastPM * fastpm) {
    /* apply periodic boundary and move particles to the correct rank */
    pm_store_wrap(fastpm->p, fastpm->pm->BoxSize);
    pm_store_decompose(fastpm->p, to_rank, fastpm->pm, fastpm->comm);
    size_t np_max;
    size_t np_min;
    /* FIXME move NTask to somewhere else. */
    double np_mean = pow(fastpm->nc, 3) / fastpm->pm->NTask;
    MPI_Allreduce(&fastpm->p->np, &np_max, 1, MPI_LONG, MPI_MAX, fastpm->comm);
    MPI_Allreduce(&fastpm->p->np, &np_min, 1, MPI_LONG, MPI_MIN, fastpm->comm);

    msg_printf(info, "Load imbalance is - %g / + %g\n",
        np_min / np_mean, np_max / np_mean);

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

void 
fastpm_interp(FastPM * fastpm, double * aout, int nout, 
            fastpm_interp_action action, void * userdata) 
{
    /* interpolate and write snapshots, assuming p 
     * is at time a_x and a_v. */
    double a_x = fastpm->p->a_x;
    double a_v = fastpm->p->a_v;
    int iout;
    for(iout = 0; iout < nout; iout ++) {
        if(
        ! /* after a kick */
        (a_x < aout[iout] && aout[iout] < a_v)
        &&
        ! /* after a drift */
        (a_x >= aout[iout] && aout[iout] >= a_v)
        ) continue;

        PMStore * snapshot = alloca(sizeof(PMStore));
        pm_store_init(snapshot);
        pm_store_alloc(snapshot, fastpm->p->np_upper, PACK_ID | PACK_POS | PACK_VEL);

        msg_printf(verbose, "Taking a snapshot...\n");

        fastpm_set_snapshot(fastpm, fastpm->p, snapshot, aout[iout]);
        walltime_measure("/Snapshot/KickDrift");

        action(fastpm, snapshot, userdata);

        pm_store_destroy(snapshot);

    }
}

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, void * template) 
{
    char filebase[1024];
    double aout = snapshot->a_x;
    double z_out= 1.0/aout - 1.0;

    sprintf(filebase, template, aout);
    ensure_dir(filebase);
    fastpm_write_runpb_snapshot(fastpm, snapshot, filebase);
    
    walltime_measure("/Snapshot/IO");
    MPI_Barrier(fastpm->comm);
    walltime_measure("/Snapshot/Wait");

    msg_printf(normal, "snapshot %s written z = %6.4f a = %6.4f\n", 
            filebase, z_out, aout);
    return 0;
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

void 
fastpm_set_time(FastPM * fastpm, 
    int istep,
    double * a_x,
    double * a_x1,
    double * a_v,
    double * a_v1) 
{
    double * time_step = fastpm->time_step;
    int nstep = fastpm->n_time_step;

    /* The last step is the terminal step. */
    *a_x = time_step[(istep >= nstep)?(nstep - 1):istep];
    *a_x1 = time_step[(istep + 1 >= nstep)?(nstep - 1):(istep + 1)];

    double a_xm1 = time_step[(istep > 0)?(istep - 1):0];
    *a_v = sqrt(a_xm1 * *(a_x));
    *a_v1 = sqrt(*a_x * *a_x1);
    fastpm->istep = istep;
    VPM * vpm = vpm_find(fastpm->vpm_list, *a_x);
    fastpm->pm = &vpm->pm;
}

/*
void 
fastpm_paint (FastPM * fastpm, FastPMFloat * delta_x, double density_factor) {
    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL); 
    pm_paint(pm, delta_x, p, p->np + pgd->nghosts, density_factor);
    pm_ghosts_free(pgd);
}
*/
