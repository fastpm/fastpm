#include "pmpfft.h"
#include "msg.h"
#include "walltime.h"

static void 
apply_force_kernel(PM * pm, int dir) 
{
    /* This is the force in fourier space. - i k[dir] / k2 */

    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

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
//            pm->workspace[ind + 0] = pm->canvas[ind + 0];
//            pm->workspace[ind + 1] = pm->canvas[ind + 1];
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);
}

void 
pm_calculate_forces(PMStore * p, PM * pm, double density_factor)
{
    PMGhostData pgd = {
        .pm = pm,
        .pdata = p,
        .np = p->np,
        .np_upper = p->np_upper,
        .attributes = PACK_POS,
        .nghosts = 0,
        .get_position = p->iface.get_position,
    };
    walltime_measure("/Force/Init");

    pm_append_ghosts(&pgd);
    walltime_measure("/Force/AppendGhosts");

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * we compensate this later at readout by density_factor.
     * */
    pm_paint(pm, p, p->np + pgd.nghosts);
    walltime_measure("/Force/Paint");
    
    pm_r2c(pm);
    walltime_measure("/Force/FFT");

#if 0
    fwrite(pm->canvas, sizeof(pm->canvas[0]), pm->allocsize, fopen("density-k.f4", "w"));
#endif

    /* calculate the forces save them to p->acc */

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        apply_force_kernel(pm, d);
        walltime_measure("/Force/Transfer");

#if 0
        char * fname[] = { "acc-0.f4", "acc-1.f4", "acc-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname[d], "w"));
#endif
        pm_c2r(pm);
        walltime_measure("/Force/FFT");

#if 0
{
    char buf[1000];
    sprintf(buf, "accr-%d.f4-rank-%d", d, pm->ThisTask);
    fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(buf, "w"));
}
#endif

#if 0
        char * fname2[] = { "accr-0.f4", "accr-1.f4", "accr-2.f4", };
        fwrite(pm->workspace, sizeof(pm->workspace[0]), pm->allocsize, fopen(fname2[d], "w"));
#endif


#pragma omp parallel for
        for(i = 0; i < p->np + pgd.nghosts; i ++) {
            /* compensate the density is less than the true density */
            p->acc[i][d] = pm_readout_one(pm, p, i) * (density_factor / pm->Norm);
        }
        walltime_measure("/Force/Readout");

        pm_reduce_ghosts(&pgd, ACC[d]); 
        walltime_measure("/Force/ReduceGhosts");
    }
    pm_destroy_ghosts(&pgd);
    walltime_measure("/Force/Finish");

    MPI_Barrier(pm->Comm2D);
    walltime_measure("/Force/Wait");
}    

