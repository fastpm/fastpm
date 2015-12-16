#include <string.h>
#include "pmpfft.h"
#include "msg.h"
#include "walltime.h"

static void 
apply_force_kernel(PM * pm, float_t * from, float_t * to, int dir) 
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
                to[ind + 0] =   from[ind + 1] * (k_finite / kk_finite);
                to[ind + 1] = - from[ind + 0] * (k_finite / kk_finite);
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
pm_calculate_forces(PMStore * p, PM * pm, float_t * delta_k, double density_factor)
{
    walltime_measure("/Force/Init");

    PMGhostData * pgd = pm_ghosts_create(pm, p, PACK_POS, NULL); 

    walltime_measure("/Force/AppendGhosts");

    float_t * canvas = pm_alloc(pm);

    /* Watch out: this paints number of particles per cell. when pm_nc_factor is not 1, 
     * it is less than the density (a cell is smaller than the mean seperation between particles. 
     * We thus have to boost the density by density_factor.
     * */
    pm_paint(pm, canvas, p, p->np + pgd->nghosts, density_factor);
    walltime_measure("/Force/Paint");
    
    pm_r2c(pm, canvas, delta_k);
    walltime_measure("/Force/FFT");

    /* calculate the forces save them to p->acc */

    int d;
    ptrdiff_t i;
    int ACC[] = {PACK_ACC_X, PACK_ACC_Y, PACK_ACC_Z};
    for(d = 0; d < 3; d ++) {
        apply_force_kernel(pm, delta_k, canvas, d);
        walltime_measure("/Force/Transfer");

        pm_c2r(pm, canvas);
        walltime_measure("/Force/FFT");

#pragma omp parallel for
        for(i = 0; i < p->np + pgd->nghosts; i ++) {
            p->acc[i][d] = pm_readout_one(pm, canvas, p, i) / pm->Norm;
        }
        walltime_measure("/Force/Readout");

        pm_ghosts_reduce(pgd, ACC[d]); 
        walltime_measure("/Force/ReduceGhosts");
    }
    pm_free(pm, canvas);

    pm_ghosts_free(pgd);
    walltime_measure("/Force/Finish");

    MPI_Barrier(pm->Comm2D);
    walltime_measure("/Force/Wait");
}    

void 
pm_calculate_powerspectrum(PM * pm, float_t * delta_k, PowerSpectrum * ps) 
{
    PMKFactors * fac[3];

    pm_create_k_factors(pm, fac);

    memset(ps->p, 0, sizeof(ps->p[0]) * ps->size);
    memset(ps->k, 0, sizeof(ps->k[0]) * ps->size);
    memset(ps->N, 0, sizeof(ps->N[0]) * ps->size);

    double k0 = 2 * M_PI / pm->BoxSize[0];

#pragma omp parallel 
    {
        ptrdiff_t ind;
        ptrdiff_t start, end;
        ptrdiff_t i[3];

        pm_prepare_omp_loop(pm, &start, &end, i);

        for(ind = start; ind < end; ind += 2) {
            int d;
            double kk = 0.;
            double cic = 1.0;
            for(d = 0; d < 3; d++) {
                kk += fac[d][i[d] + pm->ORegion.start[d]].kk;
                cic *= fac[d][i[d] + pm->ORegion.start[d]].cic;
            }

            double real = delta_k[ind + 0];
            double imag = delta_k[ind + 1];
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
        ps->p[ind] *= pm->Volume / (pm->Norm * pm->Norm);
    }

    pm_destroy_k_factors(pm, fac);
}

void 
power_spectrum_init(PowerSpectrum * ps, size_t size) 
{
    ps->size = size;
    ps->k = malloc(sizeof(ps->k[0]) * size);
    ps->p = malloc(sizeof(ps->p[0]) * size);
    ps->N = malloc(sizeof(ps->N[0]) * size);
}

void 
power_spectrum_destroy(PowerSpectrum * ps) {
    free(ps->N);
    free(ps->p);
    free(ps->k);
}

void
power_spectrum_write(PowerSpectrum * ps, PM * pm, double ntotal, char * basename, int random_seed, double aout) 
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


