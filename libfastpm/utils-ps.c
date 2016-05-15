#include <string.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include "pmpfft.h"

void
fastpm_powerspectrum_init(FastPMPowerSpectrum * ps, PM * pm, FastPMFloat * delta1_k, FastPMFloat * delta2_k)
{
    /* N is used to store metadata -- the shot-noise level. */
    ps->size = pm_nmesh(pm)[0] / 2;
    ps->pm = pm;
    double Volume = 1.0;
    int d;
    for (d = 0; d < 3; d ++) {
        Volume *= pm_boxsize(pm)[d];
    }
    ps->Volume = Volume;
    ps->k = malloc(sizeof(ps->k[0]) * ps->size);
    ps->p = malloc(sizeof(ps->p[0]) * ps->size);
    ps->Nmodes = malloc(sizeof(ps->Nmodes[0]) * ps->size);
    ps->k0 = 6.28 / pm_boxsize(pm)[0];

    memset(ps->p, 0, sizeof(ps->p[0]) * ps->size);
    memset(ps->k, 0, sizeof(ps->k[0]) * ps->size);
    memset(ps->Nmodes, 0, sizeof(ps->Nmodes[0]) * ps->size);

    double k0 = 2 * M_PI / pm_boxsize(ps->pm)[0];

#pragma omp parallel 
    {
        PMKIter kiter;
        for(pm_kiter_init(ps->pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk = 0.;
            for(d = 0; d < 3; d++) {
                kk += kiter.kk[d][kiter.iabs[d]];
            }

            ptrdiff_t ind = kiter.ind;

            double k = sqrt(kk);
            ptrdiff_t bin = floor(k / k0);
            if(bin >= 0 && bin < ps->size) {
                double real1 = delta1_k[ind + 0];
                double imag1 = delta1_k[ind + 1];
                double real2 = delta2_k[ind + 0];
                double imag2 = delta2_k[ind + 1];
                double value = real1 * real2 + imag1 * imag2;
                int w = 2;
                /* fixme: older version of code has this bug. */
                if(kiter.i[2] == 0) w = 1;
                #pragma omp atomic
                ps->Nmodes[bin] += w;
                #pragma omp atomic
                ps->p[bin] += w * value; /// cic;
                #pragma omp atomic
                ps->k[bin] += w * k;
            }
        }
    }


    MPI_Allreduce(MPI_IN_PLACE, ps->p, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->Nmodes, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->k, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);

    ptrdiff_t ind;
    for(ind = 0; ind < ps->size; ind++) {
        ps->k[ind] /= ps->Nmodes[ind];
        ps->p[ind] /= ps->Nmodes[ind];
        ps->p[ind] *= ps->Volume;
    }
}

void
fastpm_powerspectrum_destroy(FastPMPowerSpectrum * ps) {
    free(ps->Nmodes);
    free(ps->p);
    free(ps->k);
}

void
fastpm_powerspectrum_write(FastPMPowerSpectrum * ps, char * filename, double N)
{
    FILE * fp = fopen(filename, "w");
    int i;
    fprintf(fp, "# k p N \n");
    for(i = 0; i < ps->size; i ++) {
        fprintf(fp, "%g %g %g\n", ps->k[i], ps->p[i], ps->Nmodes[i]);
    }
    double * BoxSize = pm_boxsize(ps->pm);
    fprintf(fp, "# metadata 7\n");
    fprintf(fp, "# volume %g float64\n", ps->Volume);
    fprintf(fp, "# shotnoise %g float64\n", ps->Volume / N);
    fprintf(fp, "# N1 %g int\n", N);
    fprintf(fp, "# N2 %g int\n", N);
    fprintf(fp, "# Lz %g float64\n", BoxSize[2]);
    fprintf(fp, "# Lx %g float64\n", BoxSize[0]);
    fprintf(fp, "# Ly %g float64\n", BoxSize[1]);
    fclose(fp);
}

double
fastpm_powerspectrum_large_scale(FastPMPowerSpectrum * ps, int Nmax)
{
    double kmax = Nmax * ps->k0;
    ptrdiff_t i;
    double Plin = 0;
    double Nmodes = 0;
    /* ignore zero mode ! */
    for(i = 1; (i == 1) || (i < ps->size && ps->k[i] <= kmax); i ++) {
        Plin += ps->p[i] * ps->Nmodes[i];
        Nmodes += ps->Nmodes[i];
    }
    Plin /= Nmodes;
    return Plin;
}
