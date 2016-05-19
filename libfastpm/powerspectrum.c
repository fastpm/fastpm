#include <string.h>
#include <mpi.h>
#include <math.h>
#include <alloca.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>

#include "pmpfft.h"

void
fastpm_powerspectrum_init(FastPMPowerSpectrum * ps, size_t size)
{
    ps->size = size;
    ps->pm = NULL;
    ps->k = malloc(sizeof(ps->k[0]) * ps->size);
    ps->p = malloc(sizeof(ps->p[0]) * ps->size);
    ps->Nmodes = malloc(sizeof(ps->Nmodes[0]) * ps->size);
}

void
fastpm_powerspectrum_init_from_delta(FastPMPowerSpectrum * ps, PM * pm, FastPMFloat * delta1_k, FastPMFloat * delta2_k)
{
    /* N is used to store metadata -- the shot-noise level. */
    fastpm_powerspectrum_init(ps, pm_nmesh(pm)[0] / 2);

    ps->pm = pm;
    double Volume = 1.0;
    int d;
    for (d = 0; d < 3; d ++) {
        Volume *= pm_boxsize(pm)[d];
    }
    ps->Volume = Volume;
    ps->k0 = 6.28 / pm_boxsize(pm)[0];

    double Norm = 0;

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
                if(kiter.iabs[0] == 0 &&
                   kiter.iabs[1] == 0 &&
                   kiter.iabs[2] == 0)
                    Norm = value;

                #pragma omp atomic
                ps->Nmodes[bin] += w;
                #pragma omp atomic
                ps->p[bin] += w * value; /// cic;
                #pragma omp atomic
                ps->k[bin] += w * k;
            }
        }
    }


    MPI_Allreduce(MPI_IN_PLACE, &Norm, 1, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->p, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->Nmodes, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->k, ps->size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);

    ptrdiff_t ind;
    for(ind = 0; ind < ps->size; ind++) {
        ps->k[ind] /= ps->Nmodes[ind];
        ps->p[ind] /= ps->Nmodes[ind];
        ps->p[ind] *= ps->Volume;
        if(Norm != 0)
            ps->p[ind] /= Norm;
    }
}

void
fastpm_transferfunction_init(FastPMPowerSpectrum * ps, PM * pm, FastPMFloat * src_k, FastPMFloat * dest_k)
{
    FastPMPowerSpectrum * ps2 = alloca(sizeof(*ps2));

    fastpm_powerspectrum_init_from_delta(ps, pm, src_k, src_k);
    fastpm_powerspectrum_init_from_delta(ps2, pm, dest_k, dest_k);

    ptrdiff_t i;
    for(i = 0; i < ps->size; i ++) {
        ps->p[i] = sqrt(ps2->p[i] / ps->p[i]);
    }

    fastpm_powerspectrum_destroy(ps2);
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

/* inverted calling signature for callbacks */
double
fastpm_powerspectrum_eval2(double k, FastPMPowerSpectrum * ps)
{
    return fastpm_powerspectrum_eval(ps, k);
}

double
fastpm_powerspectrum_eval(FastPMPowerSpectrum * ps, double k)
{
    int i;
    if(k == 0) return 0;

    /* ignore the 0 mode */

    int l = 1;
    int r = ps->size - 1;

    while(r - l > 1) {
        int m = (r + l) / 2;
        if(k < ps->k[m])
            r = m;
        else
            l = m;
    }
    double k2 = log(ps->k[r]),
           k1 = log(ps->k[l]);
    double p2 = log(ps->p[r]),
           p1 = log(ps->p[l]);

    k = log(k);

    double p = (k - k1) * p2 + (k2 - k) * p1;
    p /= (k2 - k1);

//    fastpm_info("Evaluating P(%g) = %g, l=%d, r=%d\n", exp(k), exp(p), l, r);
    return exp(p);
}

double
fastpm_powerspectrum_get(FastPMPowerSpectrum * ps, double k)
{
    int i;
    if(k == 0) return 0;

    /* ignore the 0 mode */

    int l = 1;
    int r = ps->size - 1;

    while(r - l > 1) {
        int m = (r + l) / 2;
        if(k < ps->k[m])
            r = m;
        else
            l = m;
    }

    return ps->p[l];
}

double
fastpm_powerspectrum_get2(double k, FastPMPowerSpectrum * ps)
{
    return fastpm_powerspectrum_get(ps, k);
}

struct sigma2_int {
    FastPMPowerSpectrum * ps;
    double R;
};
static
double sigma2_int(double k, struct sigma2_int * param)
{
    double kr, kr3, kr2, w, x;
    double r_tophat = param->R;

    kr = r_tophat * k;
    kr2 = kr * kr;
    kr3 = kr2 * kr;

    if(kr < 1e-8)
        return 0;

    w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
    x = 4 * M_PI * k * k * w * w * fastpm_powerspectrum_eval(param->ps, k);

    return x / pow(2 * M_PI, 3);
}

double
fastpm_powerspectrum_sigma(FastPMPowerSpectrum * ps, double R)
{
    const int WORKSIZE = 81920;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    struct sigma2_int param;
    param.ps = ps;
    param.R = R;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = (void*) &sigma2_int;
    F.params = &param;

    void * handler = gsl_set_error_handler_off ();
    //
    // note: 500/R is here chosen as (effectively) infinity integration boundary
    gsl_integration_qag(&F, 0, 500.0 * 1 / R,
            0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

    // high precision caused error
    gsl_integration_workspace_free(workspace);
    gsl_set_error_handler (handler);

    return sqrt(result);
}

void
fastpm_powerspectrum_scale(FastPMPowerSpectrum * ps, double factor)
{
    ptrdiff_t i;
    /* neglect the zero mode */
    for(i = 1; i < ps->size; i ++) {
        ps->p[i] *= factor;
    }
}

void
fastpm_powerspectrum_rebin(FastPMPowerSpectrum * ps, size_t newsize)
{
    /* this doesn't work */
    double * k1 = malloc(newsize * sizeof(k1[0]));
    double * p1 = malloc(newsize * sizeof(p1[0]));
    double * Nmodes1 = malloc(newsize * sizeof(Nmodes1[0]));

    ptrdiff_t i;

    k1[0] = 0;
    p1[0] = ps->p[0]; /* shall be 1 */
    Nmodes1[0] = ps->Nmodes[0]; /* shall be 1 */

    for(i = 1; i < newsize; i ++) {
        ptrdiff_t j1 = 1 + (i - 1) * (ps->size - 1) / (newsize - 1);
        ptrdiff_t j2 = 1 + (i) * (ps->size - 1) / (newsize - 1);
        ptrdiff_t j;
        k1[i] = 0;
        p1[i] = 0;
        Nmodes1[i] = 0;
        for(j = j1; j == j1 || j < j2; j ++) {
            k1[i] += ps->k[j] * ps->Nmodes[j];
            p1[i] += ps->p[j] * ps->Nmodes[j];
            Nmodes1[i] += ps->Nmodes[j];
        }
        k1[i] /= Nmodes1[j];
        p1[i] /= Nmodes1[j];
    }
    free(ps->k);
    free(ps->p);
    free(ps->Nmodes);
    ps->k = k1;
    ps->p = p1;
    ps->Nmodes = Nmodes1;
    ps->size = newsize;
}

int
fastpm_powerspectrum_init_from_string(FastPMPowerSpectrum * ps, const char * string)
{
    char ** list = fastpm_strsplit(string, "\n");
    char ** line;
    ptrdiff_t i;
    int pass = 0;
    /* two pass parsing, first pass for counting */
    /* second pass for assignment */
    while(pass < 2) {
        i = 0;
        for (line = list; *line; line++) {
            double k, p;
            if(2 == sscanf(*line, "%lg %lg", &k, &p)) {
                if(pass == 1) {
                    ps->k[i + 1] = k;
                    ps->p[i + 1] = p;
                    ps->Nmodes[i + 1] = 1;
                }
                i ++;
            }
        }

        if(pass == 0) {
            fastpm_powerspectrum_init(ps, i + 1);
            ps->k[0] = 0;
            ps->p[0] = 1;
            ps->Nmodes[0] = 1;
        }
        pass ++;
    }

    free(list);
    if(ps->size == 0) {
        /* Nothing is savagable in the file.*/
        return 0;
    }
    return 0;
}
