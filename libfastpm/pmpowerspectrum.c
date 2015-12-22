#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"

void 
pm_calculate_powerspectrum(PM * pm, FastPMFloat * delta_k, PowerSpectrum * ps) 
{
    double Volume = 1.0;
    double Norm = 1.0;
    int d;
    for (d = 0; d < 3; d ++) {
        Volume *= pm_boxsize(pm)[d];
        Norm *= pm_nmesh(pm)[d];
    }
    memset(ps->p, 0, sizeof(ps->p[0]) * ps->size);
    memset(ps->k, 0, sizeof(ps->k[0]) * ps->size);
    memset(ps->N, 0, sizeof(ps->N[0]) * ps->size);

    double k0 = 2 * M_PI / pm_boxsize(pm)[0];

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk = 0.;
            for(d = 0; d < 3; d++) {
                kk += kiter.fac[d][kiter.iabs[d]].kk;
            }

            ptrdiff_t ind = kiter.ind;

            double real = delta_k[ind + 0];
            double imag = delta_k[ind + 1];
            double value = real * real + imag * imag;
            double k = sqrt(kk);
            ptrdiff_t bin = floor(k / k0);
            if(bin >= 0 && bin < ps->size) {
                int w = 2;
                /* fixme: older version of code has this bug. */
                if(kiter.i[2] == 0) w = 1;
                #pragma omp atomic
                ps->N[bin] += w;
                #pragma omp atomic
                ps->p[bin] += w * value; /// cic;
                #pragma omp atomic
                ps->k[bin] += w * k;
            }
        }
    }


    MPI_Allreduce(MPI_IN_PLACE, ps->p, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->N, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->k, ps->size, MPI_DOUBLE, MPI_SUM, pm->Comm2D);

    ptrdiff_t ind;
    for(ind = 0; ind < ps->size; ind++) {
        ps->k[ind] /= ps->N[ind];
        ps->p[ind] /= ps->N[ind];
        ps->p[ind] *= Volume / (Norm * Norm);
    }
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
    double * BoxSize = pm_boxsize(pm);
    double Volume = BoxSize[0] * BoxSize[1] * BoxSize[1];
    fprintf(fp, "# metadata 7\n");
    fprintf(fp, "# volume %g float64\n", Volume);
    fprintf(fp, "# shotnoise %g float64\n", Volume / ntotal);
    fprintf(fp, "# N1 %g int\n", ntotal);
    fprintf(fp, "# N2 %g int\n", ntotal);
    fprintf(fp, "# Lz %g float64\n", BoxSize[2]);
    fprintf(fp, "# Lx %g float64\n", BoxSize[0]);
    fprintf(fp, "# Ly %g float64\n", BoxSize[1]);
    fclose(fp);
}

