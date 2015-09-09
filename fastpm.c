#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "pmpfft.h"
#include "pm2lpt.h"
#include "parameters.h"
#include "pmsteps.h"
#include "msg.h"
#include "power.h"
#include "pmtimer.h"

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void rungdb(const char* fmt, ...);
static void calculate_forces(PMStore * p, PM * pm, double density_factor);
static int to_rank(PMStore * pdata, ptrdiff_t i, void * data);
static double sinc_unnormed(double x);

typedef struct {
    Parameters * param;
    PMStore * p;
    int nout;
    double * aout;
    int iout;
} SNPS;

static int snps_interp(SNPS * snps, double a_x, double a_v);
static void snps_init(SNPS * snps, Parameters * prr, PMStore * p);
static void snps_start(SNPS * snps);

typedef struct {
    PM pm;
    double a_start;
    int pm_nc_factor;
} VPM;

static VPM vpm[2];
    
static void vpm_init(Parameters * prr, PMIFace * iface, MPI_Comm comm) {
    /* plan for the variable PMs */
    int i;
    for (i = 0; i < 2; i ++) {
        if (i == 0) {
            vpm[i].pm_nc_factor = prr->pm_nc_factor1;
            vpm[i].a_start = 0;
        } else {
            vpm[i].pm_nc_factor = prr->pm_nc_factor2;
            vpm[i].a_start = prr->change_pm;
        }
        PMInit pminit = {
            .Nmesh = (int)(prr->nc * vpm[i].pm_nc_factor),
            .BoxSize = prr->boxsize,
            .NprocX = 0, /* 0 for auto, 1 for slabs */
            .transposed = 1,
        };
        pm_pfft_init(&vpm[i].pm, &pminit, iface, comm);
    }
}

static VPM vpm_find(double a) {
    /* find the PM object for force calculation at time a*/
    int i;
    for (i = 0; i < 2; i ++) {
        if(vpm[i].a_start > a) break;
    }
    return vpm[i-1];
}

int main(int argc, char ** argv) {

    MPI_Init(&argc, &argv);

    msg_init();
    msg_set_loglevel(verbose);

    Parameters prr;
    PMStore pdata;
    PM pm;

    timer_set_category(INIT);

    read_parameters(argc, argv, &prr);

    stepping_init(&prr);

    MPI_Comm comm = MPI_COMM_WORLD;

    int NTask;
    MPI_Comm_size(comm, &NTask);
    power_init(prr.power_spectrum_filename, 
            prr.time_step[0], 
            prr.sigma8, 
            prr.omega_m, 
            1 - prr.omega_m);

    
    pm_store_init(&pdata);

    pm_store_alloc(&pdata, 1.0 * prr.nc * prr.nc * prr.nc / NTask * prr.np_alloc_factor);

    vpm_init(&prr, &pdata.iface, comm);

    timer_set_category(LPT);

    pm_2lpt_main(&pdata, prr.nc, prr.boxsize, PowerSpecWithData, prr.random_seed, NULL, comm);

    double shift[3] = {
        prr.boxsize / prr.nc * 0.5,
        prr.boxsize / prr.nc * 0.5,
        prr.boxsize / prr.nc * 0.5,
        };

    stepping_set_initial(prr.time_step[0], &pdata, shift);

    SNPS snps;

    snps_init(&snps, &prr, &pdata);

    snps_start(&snps);

    timer_set_category(STEPPING);

    int istep;
    int nsteps = stepping_get_nsteps();

    snps_interp(&snps, prr.time_step[0], prr.time_step[0]);

    for (istep = 0; istep <= nsteps; istep++) {
        double a_v, a_x, a_v1, a_x1;

        stepping_get_times(istep,
            &a_x, &a_x1, &a_v, &a_v1);

        VPM vpm = vpm_find(a_x);
        PM * pm = &vpm.pm;

        timer_start("comm");
        pm_store_wrap(&pdata, pm->BoxSize);
        pm_store_decompose(&pdata, to_rank, pm, comm);
        timer_stop("comm");

        if(prr.force_mode & FORCE_MODE_PM) {
            calculate_forces(&pdata, pm, pow(vpm.pm_nc_factor, 3)); 
        }
#if 0
        fwrite(pdata.x, sizeof(pdata.x[0]), pdata.np, fopen("x.f8x3", "w"));
        fwrite(pdata.v, sizeof(pdata.v[0]), pdata.np, fopen("v.f4x3", "w"));
        fwrite(pdata.id, sizeof(pdata.id[0]), pdata.np, fopen("id.i8", "w"));
        fwrite(pdata.acc, sizeof(pdata.acc[0]), pdata.np, fopen("acc.f4x3", "w"));
#endif
        if(snps_interp(&snps, a_x, a_v)) break;

        // Leap-frog "kick" -- velocities updated

        timer_start("evolve");  
        stepping_kick(&pdata, a_v, a_v1, a_x);
        timer_stop("evolve");  

        if(snps_interp(&snps, a_x, a_v1)) break;
        
        // Leap-frog "drift" -- positions updated
        timer_start("evolve");  
        stepping_drift(&pdata, a_x, a_x1, a_v1);
        timer_stop("evolve");  

    }
    pm_store_destroy(&pdata);
    timer_print();
    pfft_cleanup();
    MPI_Finalize();
}

static int to_rank(PMStore * pdata, ptrdiff_t i, void * data) {
    PM * pm = (PM*) data;
    double pos[3];
    pdata->iface.get_position(pdata, i, pos);
    return pm_pos_to_rank(pm, pos);
}
static double diff_kernel(double w) {
    /* order N = 1 */
    /* 
     * This is the same as GADGET-2 but in fourier space: 
     * see gadget-2 paper and Hamming's book.
     * c1 = 2 / 3, c2 = 1 / 12
     * */
    return 1 / 6.0 * (8 * sin (w) - sin (2 * w));
}
static void apply_force_kernel(PM * pm, int dir) {
    ptrdiff_t ind;
    int d;

    struct {
        float k_finite; /* i k, finite */
        float kk_finite; /* k ** 2, on a mesh */
        float kk;  /* k ** 2 */
    } * fac[3];
    
    for(d = 0; d < 3; d++) {
        fac[d] = alloca(sizeof(fac[0][0]) * pm->Nmesh[d]);
        double CellSize = pm->BoxSize[d] / pm->Nmesh[d];
        for(ind = 0; ind < pm->Nmesh[d]; ind ++) {
            float k = pm->MeshtoK[d][ind];
            float w = k * CellSize;
            float ff = sinc_unnormed(0.5 * w);

            fac[d][ind].k_finite = 1 / CellSize * diff_kernel(w);
            fac[d][ind].kk_finite = k * k * ff * ff;
            fac[d][ind].kk = k * k;
        }
    } 
    
    ptrdiff_t i[3] = {0};
    for(ind = 0; ind < pm->ORegion.total * 2; 
        ind += 2, pm_inc_o_index(pm, i)
    ) {
        double k_finite = fac[dir][i[dir] + pm->ORegion.start[dir]].k_finite;
        double kk_finite = 0;
        double kk = 0;
        for(d = 0; d < 3; d++) {
            kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
//            kk += fac[d][i[d]].kk;
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
    }
}

static void calculate_forces(PMStore * p, PM * pm, double density_factor) {

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
     * it is less than the density (a cell is smaller than the mean seperation between particles. */
    pm_paint(pm, p, p->np + pgd.nghosts);
    
    timer_stop("paint");    

#if 0
    fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen("density.f4", "w"));
#endif
    timer_start("fft");
    pm_r2c(pm);
    timer_stop("fft");

#if 0
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("density-k.f4", "w"));
#endif

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Y};
    for(d = 0; d < 3; d ++) {
        timer_start("transfer");
        apply_force_kernel(pm, d);
        timer_stop("transfer");

#if 0
        char * fname[] = {
                "acc-0.f4",
                "acc-1.f4",
                "acc-2.f4",
            };

        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname[d], "w"));
#endif
        timer_start("fft");
        pm_c2r(pm);
        timer_stop("fft");

#if 0
        char * fname2[] = {
                "accr-0.f4",
                "accr-1.f4",
                "accr-2.f4",
            };

        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname2[d], "w"));
#endif


        timer_start("readout");

        for(i = 0; i < p->np + pgd.nghosts; i ++) {
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
    char filebase[1024];    
    PMStore * p = snps->p;
    Parameters * param = snps->param;
    PMStore snapshot;
    double BoxSize[3] = {param->boxsize, param->boxsize, param->boxsize};

    timer_set_category(SNP);
    while(snps->iout < snps->nout && (
        /* after a kick */
        (a_x < snps->aout[snps->iout] && snps->aout[snps->iout] <= a_v)
        ||
        /* after a drift */
        (a_x >= snps->aout[snps->iout] && snps->aout[snps->iout] >= a_v)
        )) {

        pm_store_init(&snapshot);

        pm_store_alloc_bare(&snapshot, p->np_upper);

        msg_printf(verbose, "Taking a snapshot...\n");

        double aout = snps->aout[snps->iout];
        int isnp= snps->iout+1;

        stepping_set_snapshot(aout, a_x, a_v, p, &snapshot);

        timer_start("comm");
        pm_store_wrap(&snapshot, BoxSize);
        timer_stop("comm");

        timer_start("write");

        if(param->snapshot_filename) {
            sprintf(filebase, "%s%05d_%0.04f.bin", param->snapshot_filename, param->random_seed, aout);
            write_runpb_snapshot(param, &snapshot, aout, filebase);
        }
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
static void snps_init(SNPS * snps, Parameters * prr, PMStore * p) {
    snps->iout = 0;
    snps->nout = prr->n_zout;
    snps->param = prr;
    snps->p = p;

    snps->aout = malloc(sizeof(double)*snps->nout);
    int i;
    for(i=0; i<snps->nout; i++) {
        snps->aout[i] = (double)(1.0/(1 + prr->zout[i]));
        msg_printf(verbose, "zout[%d]= %lf, aout= %f\n", 
                i, prr->zout[i], snps->aout[i]);
    }
}

static void snps_start(SNPS * snps) {
    snps->iout = 0;
}
