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
fastpm_powerspectrum_init(FastPMPowerSpectrum * ps, const size_t size)
{
    fastpm_funck_init(&ps->base, size);
    ps->pm = NULL;
    ps->edges = malloc(sizeof(ps->edges[0]) * (ps->base.size + 1));
    ps->Nmodes = malloc(sizeof(ps->Nmodes[0]) * ps->base.size);
}

void
fastpm_powerspectrum_init_from(FastPMPowerSpectrum * ps, const FastPMPowerSpectrum * other)
{
    fastpm_powerspectrum_init(ps, other->base.size);
    memcpy(ps->base.k, other->base.k, sizeof(ps->base.k[0]) * ps->base.size);
    memcpy(ps->base.f, other->base.f, sizeof(ps->base.f[0]) * ps->base.size);
    memcpy(ps->edges, other->edges, sizeof(ps->edges[0]) * (ps->base.size + 1));
    memcpy(ps->Nmodes, other->Nmodes, sizeof(ps->Nmodes[0]) * ps->base.size);
}

void
fastpm_powerspectrum_init_from_delta(FastPMPowerSpectrum * ps, PM * pm, const FastPMFloat * delta1_k, const FastPMFloat * delta2_k)
{
    /* This function measures powerspectrum from two overdensity or 1+overdensity fields */
    /* normalize them with fastpm_apply_normalize_transfer if needed before using this function.*/
    /* N is used to store metadata -- the shot-noise level. */
    fastpm_powerspectrum_init(ps, pm_nmesh(pm)[0] / 2);

    ps->pm = pm;
    double Volume = 1.0;
    int d;
    double k0 = 2 * M_PI / pm_boxsize(ps->pm)[0];
    for (d = 0; d < 3; d ++) {
        Volume *= pm_boxsize(pm)[d];
    }
    ps->Volume = Volume;
    ps->k0 = k0;

    memset(ps->base.f, 0, sizeof(ps->base.f[0]) * ps->base.size);
    memset(ps->base.k, 0, sizeof(ps->base.k[0]) * ps->base.size);
    memset(ps->edges, 0, sizeof(ps->edges[0]) * (ps->base.size + 1));
    memset(ps->Nmodes, 0, sizeof(ps->Nmodes[0]) * ps->base.size);

    int i;
    for(i = 0; i < ps->base.size + 1; i ++) {
        ps->edges[i] = i * k0;
    }

#pragma omp parallel
    {
        PMKIter kiter;
        for(pm_kiter_init(ps->pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            ptrdiff_t kk = 0.;
            for(d = 0; d < 3; d++) {
                double ik = kiter.iabs[d];
                if(ik > pm->Nmesh[d] / 2) ik -= pm->Nmesh[d];
                kk += ik * ik;
            }

            ptrdiff_t ind = kiter.ind;

            ptrdiff_t bin = ((ptrdiff_t)floor(sqrt(kk))) - 2;
            if(bin < 0) bin = 0;
            while((bin + 1) * (bin + 1) <= kk) {
                bin ++;
            }

            double k = sqrt(kk) * k0;

            if(bin >= 0 && bin < ps->base.size) {
                double real1 = delta1_k[ind + 0];
                double imag1 = delta1_k[ind + 1];
                double real2 = delta2_k[ind + 0];
                double imag2 = delta2_k[ind + 1];
                double value = real1 * real2 + imag1 * imag2;
                int w = 2;

                if(kiter.i[2] == 0 || kiter.i[2] == pm->Nmesh[2] / 2) w = 1;

                if(kiter.iabs[0] == 0 &&
                   kiter.iabs[1] == 0 &&
                   kiter.iabs[2] == 0) {

                } else {
                    #pragma omp atomic
                    ps->Nmodes[bin] += w;
                    #pragma omp atomic
                    ps->base.f[bin] += w * value; /// cic;
                    #pragma omp atomic
                    ps->base.k[bin] += w * k;
                }
            }
        }
    }


    MPI_Allreduce(MPI_IN_PLACE, ps->base.f, ps->base.size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->Nmodes, ps->base.size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);
    MPI_Allreduce(MPI_IN_PLACE, ps->base.k, ps->base.size, MPI_DOUBLE, MPI_SUM, ps->pm->Comm2D);

    ptrdiff_t ind;
    for(ind = 0; ind < ps->base.size; ind++) {
        if(ps->Nmodes[ind] == 0) continue;
        ps->base.k[ind] /= ps->Nmodes[ind];
        ps->base.f[ind] /= ps->Nmodes[ind];
        ps->base.f[ind] *= ps->Volume;
    }
}

void
fastpm_transferfunction_init(FastPMPowerSpectrum * ps, PM * pm, FastPMFloat * src_k, FastPMFloat * dest_k)
{
    FastPMPowerSpectrum * ps2 = alloca(sizeof(*ps2));

    fastpm_powerspectrum_init_from_delta(ps, pm, src_k, src_k);
    fastpm_powerspectrum_init_from_delta(ps2, pm, dest_k, dest_k);

    ptrdiff_t i;
    for(i = 0; i < ps->base.size; i ++) {
        ps->base.f[i] = sqrt(ps2->base.f[i] / ps->base.f[i]);
    }

    fastpm_powerspectrum_destroy(ps2);
}

void
fastpm_powerspectrum_destroy(FastPMPowerSpectrum * ps) {
    free(ps->edges);
    free(ps->Nmodes);
    fastpm_funck_destroy(&ps->base);
}

void
fastpm_powerspectrum_write(FastPMPowerSpectrum * ps, char * filename, double N)
{
    FILE * fp = fopen(filename, "w");
    int i;
    fprintf(fp, "# k p N \n");
    for(i = 0; i < ps->base.size; i ++) {
        fprintf(fp, "%g %g %g\n", ps->base.k[i], ps->base.f[i], ps->Nmodes[i]);
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
    for(i = 0; (i == 0) || (i < ps->base.size && ps->base.k[i] <= kmax); i ++) {
        Plin += ps->base.f[i] * ps->Nmodes[i];
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
    return fastpm_funck_eval(&ps->base, k);
}

double
fastpm_powerspectrum_get(FastPMPowerSpectrum * ps, double k)
{
    if(k == 0) return 1;

    /* ignore the 0 mode */

    int l = 0;
    int r = ps->base.size;

    while(r - l > 1) {
        int m = (r + l) / 2;
        /* if we are on the exact k, return value from there.*/
        if(k <= ps->edges[m])
            r = m;
        else
            l = m;
    }

    return ps->base.f[l];
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
    for(i = 1; i < ps->base.size; i ++) {
        ps->base.f[i] *= factor;
    }
}

void
fastpm_powerspectrum_rebin(FastPMPowerSpectrum * ps, size_t newsize)
{
    /* this doesn't work */
    double * k1 = malloc(newsize * sizeof(k1[0]));
    double * p1 = malloc(newsize * sizeof(p1[0]));
    double * Nmodes1 = malloc(newsize * sizeof(Nmodes1[0]));
    double * edges1 = malloc((newsize + 1) * sizeof(edges1[0]));

    ptrdiff_t i;

    for(i = 0; i < newsize; i ++) {
        ptrdiff_t j1 = (i    ) * ps->base.size / newsize;
        ptrdiff_t j2 = (i + 1) * ps->base.size / newsize;
        ptrdiff_t j;
        k1[i] = 0;
        p1[i] = 0;
        Nmodes1[i] = 0;
        edges1[i] = ps->edges[j1];
        edges1[i + 1] = ps->edges[j2];
        for(j = j1; j < j2; j ++) {
            k1[i] += ps->base.k[j] * ps->Nmodes[j];
            p1[i] += ps->base.f[j] * ps->Nmodes[j];
            Nmodes1[i] += ps->Nmodes[j];
        }
        if(Nmodes1[i] > 0) {
            k1[i] /= Nmodes1[i];
            p1[i] /= Nmodes1[i];
        }
    }
    free(ps->base.k);
    free(ps->base.f);
    free(ps->Nmodes);
    free(ps->edges);
    ps->base.k = k1;
    ps->base.f = p1;
    ps->Nmodes = Nmodes1;
    ps->edges = edges1;
    ps->base.size = newsize;
}

int
fastpm_powerspectrum_init_from_string(FastPMPowerSpectrum * ps, const char * string)
{
    int r = fastpm_funck_init_from_string(&ps->base, string);
    ps->edges = malloc(sizeof(ps->edges[0]) * (ps->base.size + 1));
    ps->Nmodes = malloc(sizeof(ps->Nmodes[0]) * ps->base.size);
    return r;
}

void
fastpm_funck_init(FastPMFuncK * fk, const size_t size)
{
    fk->size = size;
    fk->k = malloc(sizeof(fk->k[0]) * fk->size);
    fk->f = malloc(sizeof(fk->f[0]) * fk->size);
}

int
fastpm_funck_init_from_string(FastPMFuncK * fk, const char * string)
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
            double k, f;
            if(2 == sscanf(*line, "%lg\t%lg", &k, &f)) {   // note tab. Could make more general?
                if(pass == 1) {
                    fk->k[i] = k;
                    fk->f[i] = f;
                }
                i ++;
            }
        }
        if(pass == 0) {
            fastpm_funck_init(fk, i);
        }
        pass ++;
    }

    free(list);
    if(fk->size == 0) {
        /* Nothing is savagable in the file.*/
        return 0;
    }
    return 0;
}

/* inverted calling signature for callbacks */
double
fastpm_funck_eval2(double k, FastPMFuncK * fk)
{
    return fastpm_funck_eval(fk, k);
}

double
fastpm_funck_eval(FastPMFuncK * fk, double k)
{

    /* ignore the 0 mode */
    if(k == 0) return 1;

    int l = 0;
    int r = fk->size - 1;

    /* determine the right and left values around k */
    while(r - l > 1) {
        int m = (r + l) / 2;
        if(k < fk->k[m])
            r = m;
        else
            l = m;
    }
    double k2 = fk->k[r],
           k1 = fk->k[l];
    double f2 = fk->f[r],
           f1 = fk->f[l];

    if(l == r) {
        return fk->f[l];
    }

    if(f1 <= 0 || f2 <= 0 || k1 == 0 || k2 == 0) {
        /* if any of the f is zero, or negative, use linear interpolation */
        double f = (k - k1) * f2 + (k2 - k) * f1;
        f /= (k2 - k1);
        return f;
    } else {
        k = log(k);
        f1 = log(f1);
        f2 = log(f2);
        k1 = log(k1);
        k2 = log(k2);
        double f = (k - k1) * f2 + (k2 - k) * f1;
        f /= (k2 - k1);
        return exp(f);
    }
}

void
fastpm_funck_destroy(FastPMFuncK * fk)
{
    free(fk->f);
    free(fk->k);
}
