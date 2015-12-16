#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "cosmology.h"
#include "pmpfft.h"
#include "pm2lpt.h"
#include "msg.h"

static void 
apply_za_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir) 
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
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
            }
            /* - i k[d] / k2 */
            if(LIKELY(kk_finite > 0)) {
                to[ind + 0] = - from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = + from[ind + 0] * (k_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);
}


static void 
apply_2lpt_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir1, int dir2) 
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
            double k1_finite = fac[dir1][i[dir1] + pm->ORegion.start[dir1]].k_finite;
            double k2_finite = fac[dir2][i[dir2] + pm->ORegion.start[dir2]].k_finite;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += fac[d][i[d] + pm->ORegion.start[d]].kk_finite;
            }
            if(kk_finite > 0) {
                to[ind + 0] = from[ind + 0] * (-k1_finite * k2_finite / kk_finite);
                to[ind + 1] = from[ind + 1] * (-k1_finite * k2_finite / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
            pm_inc_o_index(pm, i);
        }
    }
    pm_destroy_k_factors(pm, fac);
}


void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, PMStore * p, double shift[3]) 
{
/* calculate dx1, dx2, for initial fluctuation delta_k.
 * shift: martin has shift = 0.5, 0.5, 0.5.
 * Use shift of 0, 0, 0 if in doublt. 
 *   */

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL);

    ptrdiff_t i;
    int d;

    FastPMFloat * workspace = pm_alloc(pm);
    FastPMFloat * source =  pm_alloc(pm);
    memset(source, 0, sizeof(source[0]) * pm->allocsize);

    FastPMFloat * field[3];

    for(d = 0; d < 3; d++ )
        field[d] = pm_alloc(pm);

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] -= shift[d];
        }
    }
     
    int DX1[] = {PACK_DX1_X, PACK_DX1_Y, PACK_DX1_Z};
    int DX2[] = {PACK_DX2_X, PACK_DX2_Y, PACK_DX2_Z};
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};

    for(d = 0; d < 3; d++) {
        msg_printf(info, "Solving for DX1 axis = %d\n", d);

        apply_za_transfer(pm, delta_k, workspace, d);

        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {        
            p->dx1[i][d] = pm_readout_one(pm, workspace, p, i);
        }
        pm_ghosts_reduce(pgd, DX1[d]);
    } 

    for(d = 0; d< 3; d++) {
        msg_printf(info, "Solving for 2LPT axes = %d %d .\n", d, d);
        apply_2lpt_transfer(pm, delta_k, field[d], d, d);
        pm_c2r(pm, field[d]);
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }    
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
        msg_printf(info, "Solving for 2LPT axes = %d %d .\n", d1, d2);
        apply_2lpt_transfer(pm, delta_k, workspace, d1, d2);
        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] -= workspace[i] * workspace[i];
        }
    } 
    pm_r2c(pm, source, workspace);
    pm_assign(pm, workspace, source);

    for(d = 0; d < 3; d++) {
        msg_printf(info, "Solving for DX2 axis = %d .\n", d);
        /* 
         * We absorb some the negative factor in za transfer to below;
         *
         * */
        apply_za_transfer(pm, source, workspace, d);
        pm_c2r(pm, workspace);

#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {        
            /* this ensures x = x0 + dx1(t) + 3/ 7 dx2(t) */
            p->dx2[i][d] = pm_readout_one(pm, workspace, p, i) / pm->Norm ;
        }
        pm_ghosts_reduce(pgd, DX2[d]);
    }
    double dx1disp[3] = {0};
    double dx2disp[3] = {0};

#pragma omp parallel for
    for(i = 0; i < p->np; i ++) {
        int d;
        for(d =0; d < 3; d++) {
#pragma omp atomic
            dx1disp[d] += p->dx1[i][d] * p->dx1[i][d];
#pragma omp atomic
            dx2disp[d] += p->dx2[i][d] * p->dx2[i][d];
        } 
    }
    uint64_t Ntot = p->np;
    MPI_Allreduce(MPI_IN_PLACE, dx1disp, 3, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, dx2disp, 3, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, &Ntot,   1, MPI_LONG,  MPI_SUM, pm->Comm2D);
    for(d =0; d < 3; d++) {
        dx1disp[d] /= Ntot;
        dx1disp[d] = sqrt(dx1disp[d]);
        dx2disp[d] /= Ntot;
        dx2disp[d] = sqrt(dx2disp[d]);
    }
    msg_printf(info, "dx1 disp : %g %g %g %g\n", 
            dx1disp[0], dx1disp[1], dx1disp[2],
            (dx1disp[0] + dx1disp[1] + dx1disp[2]) / 3.0);
    msg_printf(info, "dx2 disp : %g %g %g %g\n", 
            dx2disp[0], dx2disp[1], dx2disp[2],
            (dx2disp[0] + dx2disp[1] + dx2disp[2]) / 3.0);

#ifdef PM_2LPT_DUMP
    fwrite(p->dx1, sizeof(p->dx1[0]), p->np, fopen("dx1.f4x3", "w"));
    fwrite(p->dx2, sizeof(p->dx2[0]), p->np, fopen("dx2.f4x3", "w"));
#endif

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += shift[d];
        }
    }

    for(d = 0; d < 3; d ++) {
        pm_free(pm, field[2-d]);
    }
    pm_free(pm, source);
    pm_free(pm, workspace);

    pm_ghosts_free(pgd);

}

// Interpolate position and velocity for snapshot at a=aout
void 
pm_2lpt_evolve(double aout, PMStore * p, double Omega)
{
    int np= p->np;

    msg_printf(verbose, "Evolveing by 2lpt to a= %6.4f (z=%6.4f).\n", aout, 1.0f/aout-1);
    Cosmology c = {
            .OmegaM = Omega,
            .OmegaLambda = 1 - Omega,
        };

    const float Dplus = GrowthFactor(aout, c);

    const double omega=OmegaA(aout, c);
    const double D2 = Dplus*Dplus*pow(omega/Omega, -1.0/143.0);
    const double D20 = pow(Omega, -1.0/143.0);
    

    float Dv=DprimeQ(aout, 1.0, c); // dD_{za}/dy
    float Dv2=GrowthFactor2v(aout, c);   // dD_{2lpt}/dy

    int i;
#pragma omp parallel for 
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            /* Use the more accurate 2LPT dx2 term */
            p->dx2[i][d] *= 3.0 / 7.0 * D20;

            p->x[i][d] += Dplus * p->dx1[i][d] + D2 * p->dx2[i][d];

            p->v[i][d] = (p->dx1[i][d]*Dv + p->dx2[i][d]*Dv2);
        }
    }
}

