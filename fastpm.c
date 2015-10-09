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

#include "pmpfft.h"
#include "vpm.h"
#include "pmic.h"
#include "pm2lpt.h"
#include "parameters.h"
#include "pmsteps.h"
#include "msg.h"
#include "power.h"
#include "pmtimer.h"

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void 
rungdb(const char* fmt, ...);

static void 
ensure_dir(char * path);

/* Snapshots */
typedef struct {
    Parameters * param;
    MPI_Comm comm;
    PMStore * p;
    int nout;
    double * aout;
    int iout;
} SNPS;

static int 
snps_interp(SNPS * snps, double a_x, double a_v);
static void 
snps_init(SNPS * snps, Parameters * prr, PMStore * p, MPI_Comm comm);
static void 
snps_start(SNPS * snps);



typedef struct {
    size_t size;
    double *k;
    double *p;
    double *N;
} PowerSpectrum;

static void 
power_spectrum_init(PowerSpectrum * ps, size_t size);

static void 
power_spectrum_destroy(PowerSpectrum * ps);

static void
write_power_spectrum(PowerSpectrum * ps, PM * pm, double ntotal, char * basename, int random_seed, double aout);

/* Useful stuff */
static void 
do_pm(PMStore * p, VPM * vpm, PowerSpectrum * ps);
static int 
to_rank(void * pdata, ptrdiff_t i, void * data);
static double 
sinc_unnormed(double x);
static double 
estimate_alloc_factor(double Volume, int NTask, double failure_rate);

int fastpm(Parameters * prr, MPI_Comm comm) {

    int NTask; 
    MPI_Comm_size(comm, &NTask);

    msg_init(comm);
    msg_set_loglevel(verbose);

    PMStore pdata;
    PM pm;
    VPM * vpm_list;
    PMInit baseinit = {
            .Nmesh = prr->nc,
            .BoxSize = prr->boxsize,
            .NprocY = prr->NprocY, /* 0 for auto, 1 for slabs */
            .transposed = 1,
            .use_fftw = prr->UseFFTW,
        };

    timer_set_category(INIT);

    stepping_init(prr);

    power_init(prr->power_spectrum_filename, 
            prr->time_step[0], 
            prr->sigma8, 
            prr->omega_m, 
            1 - prr->omega_m, comm);

    pm_store_init(&pdata);

    double alloc_factor;
    alloc_factor = prr->np_alloc_factor;
    msg_printf(info, "Using alloc factor of %g\n", alloc_factor);

    pm_store_alloc_evenly(&pdata, pow(prr->nc, 3), alloc_factor, comm);

    vpm_list = vpm_create(prr->n_pm_nc_factor, 
                           prr->pm_nc_factor, 
                           prr->change_pm,
                           &baseinit, &pdata.iface, comm);


    timer_set_category(LPT);

    if(prr->readic_filename) {
        read_runpb_ic(prr, prr->time_step[0], &pdata, comm);
    } else {
        PM pm;

        pm_2lpt_init(&pm, &pdata, prr->nc, prr->boxsize, comm);

        pm_start(&pm);

        pm_ic_fill_gaussian_gadget(&pm, prr->random_seed, PowerSpecWithData, NULL);

        pm_2lpt_main(&pm, &pdata, comm);

        pm_destroy(&pm);
    }

    double shift[3] = {
        prr->boxsize / prr->nc * 0.5,
        prr->boxsize / prr->nc * 0.5,
        prr->boxsize / prr->nc * 0.5,
        };


    pm_2lpt_set_initial(prr->time_step[0], &pdata, shift, prr->omega_m);

    if(prr->force_mode != FORCE_MODE_PM) {
        /* for COLA and COLA1, v cancels out such that the initial is zero */
        memset(pdata.v, 0, sizeof(pdata.v[0]) * pdata.np);
    }

    SNPS snps;

    snps_init(&snps, prr, &pdata, comm);

    snps_start(&snps);

    timer_set_category(STEPPING);

    int istep;
    int nsteps = stepping_get_nsteps();

    snps_interp(&snps, prr->time_step[0], prr->time_step[0]);

    for (istep = 0; istep <= nsteps; istep++) {
        double a_v, a_x, a_v1, a_x1;

        /* begining and ending of drift(x) and kick(v)*/
        stepping_get_times(istep,
            &a_x, &a_x1, &a_v, &a_v1);

        /* Find the Particle Mesh to use for this time step */
        VPM * vpm = vpm_find(vpm_list, a_x);
        PM * pm = &vpm->pm;
        msg_printf(debug, "Using PM of size %td\n", pm->init.Nmesh);

        PowerSpectrum ps;

        /* apply periodic boundary and move particles to the correct rank */
        timer_start("comm");
        pm_store_wrap(&pdata, pm->BoxSize);
        pm_store_decompose(&pdata, to_rank, pm, comm);
        timer_stop("comm");

        size_t np_max;
        size_t np_min;
        double np_mean = pow(prr->nc, 3) / NTask;
        MPI_Allreduce(&pdata.np, &np_max, 1, MPI_LONG, MPI_MAX, comm);
        MPI_Allreduce(&pdata.np, &np_min, 1, MPI_LONG, MPI_MIN, comm);

        msg_printf(info, "Load imbalance is - %g / + %g\n",
            np_min / np_mean, np_max / np_mean);

        /* Calculate PM forces, only if needed. */
        power_spectrum_init(&ps, pm->Nmesh[0] / 2);

        if(prr->force_mode & FORCE_MODE_PM) {
            /* watch out: boost the density since mesh is finer than grid */
            do_pm(&pdata, vpm, &ps);
        }
        if(prr->measure_power_spectrum_filename) {
            if(pm->ThisTask == 0)
                ensure_dir(prr->measure_power_spectrum_filename);
                write_power_spectrum(&ps, pm, ((double)prr->nc * prr->nc * prr->nc), 
                    prr->measure_power_spectrum_filename, prr->random_seed, a_x);
        }
        power_spectrum_destroy(&ps);
#if 0
        fwrite(pdata.x, sizeof(pdata.x[0]), pdata.np, fopen("x.f8x3", "w"));
        fwrite(pdata.v, sizeof(pdata.v[0]), pdata.np, fopen("v.f4x3", "w"));
        fwrite(pdata.id, sizeof(pdata.id[0]), pdata.np, fopen("id.i8", "w"));
        fwrite(pdata.acc, sizeof(pdata.acc[0]), pdata.np, fopen("acc.f4x3", "w"));
#endif
        /* take snapshots if needed, before the kick */
        if(snps_interp(&snps, a_x, a_v)) break;

        // Leap-frog "kick" -- velocities updated

        timer_start("evolve");  
        stepping_kick(&pdata, &pdata, a_v, a_v1, a_x);
        timer_stop("evolve");  

        /* take snapshots if needed, before the drift */
        if(snps_interp(&snps, a_x, a_v1)) break;
        
        // Leap-frog "drift" -- positions updated
        timer_start("evolve");  
        stepping_drift(&pdata, &pdata, a_x, a_x1, a_v1);
        timer_stop("evolve");  

        /* no need to check for snapshots here, it will be checked next loop.  */
    }

    pm_store_destroy(&pdata);
    timer_print();
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

static double 
diff_kernel(double w) 
{
    /* order N = 1 super lanzcos kernel */
    /* 
     * This is the same as GADGET-2 but in fourier space: 
     * see gadget-2 paper and Hamming's book.
     * c1 = 2 / 3, c2 = 1 / 12
     * */
    return 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}


typedef struct {
    float k_finite; /* i k, finite */
    float kk_finite; /* k ** 2, on a mesh */
    float kk;  /* k ** 2 */
    float cic;  /* 1 - 2 / 3 sin^2 ( 0.5 k L / N)*/
    float extra;  /* any temporary variable that can be useful. */
} KFactors;

static void 
create_k_factors(PM * pm, KFactors * fac[3]) 
{ 
    /* This function populates fac with precalculated values that
     * are useful for force calculation. 
     * e.g. k**2 and the finite differentiation kernels. 
     * precalculating them means in the true kernel we only need a 
     * table look up. watch out for the offset ORegion.start
     * */
    int d;
    ptrdiff_t ind;
    for(d = 0; d < 3; d++) {
        fac[d] = malloc(sizeof(fac[0][0]) * pm->Nmesh[d]);
        double CellSize = pm->BoxSize[d] / pm->Nmesh[d];
        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            float k = pm->MeshtoK[d][ind];
            float w = k * CellSize;
            float ff = sinc_unnormed(0.5 * w);

            fac[d][ind].k_finite = 1 / CellSize * diff_kernel(w);
            fac[d][ind].kk_finite = k * k * ff * ff;
            fac[d][ind].kk = k * k;
            double tmp = sin(0.5 * k * CellSize);
            fac[d][ind].cic = 1 - 2. / 3 * tmp * tmp;
        }
    } 
}

static void 
destroy_k_factors(PM * pm, KFactors * fac[3]) 
{
    int d;
    for(d = 0; d < 3; d ++) {
        free(fac[d]);
    }
}

static void 
prepare_omp_loop(PM * pm, ptrdiff_t * start, ptrdiff_t * end, ptrdiff_t i[3]) 
{ 
    /* static schedule the openmp loops. start, end is in units of 'real' numbers.
     *
     * i is in units of complex numbers.
     *
     * We call pm_unravel_o_index to set the initial i[] for each threads,
     * then rely on pm_inc_o_index to increment i, because the former is 
     * much slower than pm_inc_o_index and would eliminate threading advantage.
     *
     * */
    int nth = omp_get_num_threads();
    int ith = omp_get_thread_num();

    *start = ith * pm->ORegion.total / nth * 2;
    *end = (ith + 1) * pm->ORegion.total / nth * 2;

    /* do not unravel if we are not looping at all. 
     * This fixes a FPE when
     * the rank has ORegion.total == 0 
     * -- with PFFT the last transposed dimension
     * on some ranks will be 0 */
    if(*end > *start) 
        pm_unravel_o_index(pm, *start / 2, i);

#if 0
        msg_aprintf(info, "ith %d nth %d start %td end %td pm->ORegion.strides = %td %td %td\n", ith, nth,
            *start, *end,
            pm->ORegion.strides[0],
            pm->ORegion.strides[1],
            pm->ORegion.strides[2]
            );
#endif

}

static void 
apply_force_kernel(PM * pm, int dir) 
{
    /* This is the force in fourier space. - i k[dir] / k2 */

    KFactors * fac[3];

    create_k_factors(pm, fac);

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k_finite;
            double kk_finite = 0;
            double kk = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
            }
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                pm->workspace[ind + 0] =   pm->canvas[ind + 1] * (k_finite / kk_finite);
                pm->workspace[ind + 1] = - pm->canvas[ind + 0] * (k_finite / kk_finite);
            } else {
                pm->workspace[ind + 0] = 0;
                pm->workspace[ind + 1] = 0;
            }
    //        pm->workspace[ind + 0] = pm->canvas[ind + 0];
     //       pm->workspace[ind + 1] = pm->canvas[ind + 1];
            pm_inc_o_index(pm, i);
        }
    }
    destroy_k_factors(pm, fac);
}

static void 
smooth_density(PM * pm, double r_s) 
{
    /* 
     *  This function smooth density by scale r_s. There could be a factor of sqrt(2)
     *  It is not used. */

    KFactors * fac[3];

    create_k_factors(pm, fac);
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

        prepare_omp_loop(pm, &start, &end, i);

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

    destroy_k_factors(pm, fac);
}

static void calculate_powerspectrum(PM * pm, PowerSpectrum * ps, double density_factor) {
    KFactors * fac[3];

    create_k_factors(pm, fac);

    memset(ps->p, 0, sizeof(ps->p[0]) * ps->size);
    memset(ps->k, 0, sizeof(ps->k[0]) * ps->size);
    memset(ps->N, 0, sizeof(ps->N[0]) * ps->size);

    double k0 = 2 * M_PI / pm->BoxSize[0];

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double kk = 0.;
            double cic = 1.0;
            for(d = 0; d < 3; d++) {
                kk += fac[d][i[d] + pm->ORegion.start[d]].kk;
                cic *= fac[d][i[d] + pm->ORegion.start[d]].cic;
            }

            double real = pm->canvas[ind + 0];
            double imag = pm->canvas[ind + 1];
            double value = real * real + imag * imag;
            double k = sqrt(kk);
            ptrdiff_t bin = floor(k / k0);
            if(bin >= 0 && bin < ps->size) {
                int w = 2;
                if(i[2] == 0) w = 1;
                ps->N[bin] += w;
                ps->p[bin] += w * value; /// cic;
                ps->k[bin] += w * k;
            }
            pm_inc_o_index(pm, i);
        }
    }


    MPI_Allreduce(MPI_IN_PLACE, ps->p, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->N, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->k, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);

    ptrdiff_t ind;
    for(ind = 0; ind < ps->size; ind++) {
        ps->k[ind] /= ps->N[ind];
        ps->p[ind] /= ps->N[ind];
        ps->p[ind] *= pm->Volume / (pm->Norm * pm->Norm) * (density_factor * density_factor);
    }

    destroy_k_factors(pm, fac);
}

static void 
do_pm(PMStore * p, VPM * vpm, PowerSpectrum * ps)
{
    PM * pm = &vpm->pm;
    double density_factor =  pow(vpm->pm_nc_factor, 3); 

    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
    };
    pm_start(pm);

    timer_start("ghosts1");
    pm_append_ghosts(&pgd);
    timer_stop("ghosts1");

    timer_start("paint");    

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * we compensate this later at readout by density_factor.
     * */
    pm_paint(pm, p, p->np + pgd.nghosts);
    
    timer_stop("paint");    

#if 0
    fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen("density.f4", "w"));
#endif
    timer_start("fft");
    pm_r2c(pm);
    timer_stop("fft");

    /* calculate the power spectrum */

    timer_start("power");
    calculate_powerspectrum(pm, ps, density_factor);
    timer_stop("power");
    
#if 0
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("density-k.f4", "w"));
#endif

    /* calculate the forces save them to p->acc */

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Y};
    for(d = 0; d < 3; d ++) {
        timer_start("transfer");
        apply_force_kernel(pm, d);
        timer_stop("transfer");

#if 0
        char * fname[] = { "acc-0.f4", "acc-1.f4", "acc-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname[d], "w"));
#endif
        timer_start("fft");
        pm_c2r(pm);
        timer_stop("fft");

#if 0
        char * fname2[] = { "accr-0.f4", "accr-1.f4", "accr-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname2[d], "w"));
#endif


        timer_start("readout");

#pragma omp parallel for
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            /* compensate the density is less than the true density */
            p->acc[i][d] = pm_readout_one(pm, p, i) * (density_factor / pm->Norm);
        }
        timer_stop("readout");

        timer_start("ghosts2");
        pm_reduce_ghosts(&pgd, ACC[d]); 
        timer_stop("ghosts2");
    }
    pm_destroy_ghosts(&pgd);
    pm_stop(pm);
}    

static void rungdb(const char* fmt, ...){
    /* dumpstack(void) Got this routine from http://www.whitefang.com/unix/faq_toc.html
 *     ** Section 6.5. Modified to redirect to file to prevent clutter
 *         */
    /* This needs to be changed... */
    char dbx[160];
    char cmd[160];
    char * tmpfilename;
    extern const char *__progname;
    va_list va;
    va_start(va, fmt);
    
    vsprintf(cmd, fmt, va);
    va_end(va);

    tmpfilename = tempnam(NULL, NULL);

    sprintf(dbx, "echo '%s\n' > %s", cmd, tmpfilename);
    system(dbx);

    sprintf(dbx, "echo 'where\ndetach' | gdb -batch --command=%s %s %d", tmpfilename, __progname, getpid() );
    system(dbx);
    unlink(tmpfilename);
    free(tmpfilename);

    return;
}
static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}

static int
snps_interp(SNPS * snps, double a_x, double a_v)
{
    /* interpolate and write snapshots, assuming snps->p 
     * is at time a_x and a_v. */
    char filebase[1024];    
    PMStore * p = snps->p;
    Parameters * param = snps->param;
    PMStore snapshot;
    double BoxSize[3] = {param->boxsize, param->boxsize, param->boxsize};

    timer_set_category(SNP);

    while(snps->iout < snps->nout && (
        /* after a kick */
        (a_x < snps->aout[snps->iout] && snps->aout[snps->iout] < a_v)
        ||
        /* after a drift */
        (a_x >= snps->aout[snps->iout] && snps->aout[snps->iout] >= a_v)
        )) {

        pm_store_init(&snapshot);

        pm_store_alloc_bare(&snapshot, p->np_upper);

        msg_printf(verbose, "Taking a snapshot...\n");

        double aout = snps->aout[snps->iout];
        int isnp= snps->iout+1;

        stepping_set_snapshot(p, &snapshot, aout, a_x, a_v);

        timer_start("comm");
        pm_store_wrap(&snapshot, BoxSize);
        timer_stop("comm");

        MPI_Barrier(snps->comm);
        timer_start("write");

        if(param->snapshot_filename) {
            ensure_dir(param->snapshot_filename);
            sprintf(filebase, "%s%05d_%0.04f.bin", param->snapshot_filename, param->random_seed, aout);
            write_runpb_snapshot(param, &snapshot, aout, filebase, snps->comm);
        }
        MPI_Barrier(snps->comm);
        timer_stop("write");

        const double rho_crit = 27.7455;
        const double M0 = param->omega_m*rho_crit*pow(param->boxsize / param->nc, 3.0);
        msg_printf(verbose, "mass of a particle is %g 1e10 Msun/h\n", M0); 

        const double z_out= 1.0/aout - 1.0;

        msg_printf(normal, "snapshot %d written z=%4.2f a=%5.3f\n", 
                isnp, z_out, aout);

        snps->iout ++;
        pm_store_destroy(&snapshot);
    }
    timer_set_category(STEPPING);
    return (snps->iout == snps->nout);
}

static void 
snps_init(SNPS * snps, Parameters * prr, PMStore * p, MPI_Comm comm) 
{
    snps->iout = 0;
    snps->nout = prr->n_zout;
    snps->param = prr;
    snps->p = p;
    snps->comm = comm;

    snps->aout = malloc(sizeof(double)*snps->nout);
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
power_spectrum_init(PowerSpectrum * ps, size_t size) 
{
    ps->size = size;
    ps->k = malloc(sizeof(ps->k[0]) * size);
    ps->p = malloc(sizeof(ps->p[0]) * size);
    ps->N = malloc(sizeof(ps->N[0]) * size);
}
static void 
power_spectrum_destroy(PowerSpectrum * ps) {
    free(ps->N);
    free(ps->p);
    free(ps->k);
}

static void
write_power_spectrum(PowerSpectrum * ps, PM * pm, double ntotal, char * basename, int random_seed, double aout) 
{
    char buf[1024];
    sprintf(buf, "%s%05d_%0.04f.txt", basename, random_seed, aout);
    FILE * fp = fopen(buf, "w");
    int i;
    fprintf(fp, "# k p N \n");
    for(i = 0; i < ps->size; i ++) {
        fprintf(fp, "%g %g %g\n", ps->k[i], ps->p[i], ps->N[i]);
    }
    fprintf(fp, "# metadata 7\n");
    fprintf(fp, "# volume %g float64\n", pm->Volume);
    fprintf(fp, "# shotnoise %g float64\n", pm->Volume / ntotal);
    fprintf(fp, "# N1 %g int\n", ntotal);
    fprintf(fp, "# N2 %g int\n", ntotal);
    fprintf(fp, "# Lz %g float64\n", pm->BoxSize[2]);
    fprintf(fp, "# Lx %g float64\n", pm->BoxSize[0]);
    fprintf(fp, "# Ly %g float64\n", pm->BoxSize[1]);
    fclose(fp);
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

static double 
estimate_alloc_factor(double Volume, int NTask, double failure_rate)
{

    double R = pow(4 * M_PI / 3 * Volume, 0.3333);
    double sigma = sqrt(Sigma2(R));
    double factor = 1.001;
    double probfail;
    while (1) {
        double x = (factor - 1.0);
        probfail = erfc(x / (1.414) / sigma);
        /* 1 - (1 - p) ** k == - [ exp(log(1-p) * k) - 1];
         * Because p is small, at large cores directly calculating (1-p)**k
         * can be very inaccurate. */

        probfail = - expm1(log1p(-probfail) * NTask);

        if (probfail < failure_rate) break;
        factor *= 1.01;
    }
    msg_printf(info, "Sigma %g AllocFactor %g Overrun probability %g\n", sigma, factor, probfail); 

    return factor;
}
