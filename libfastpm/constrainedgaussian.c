#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/constrainedgaussian.h>

#include "pmpfft.h"
#include <gsl/gsl_linalg.h>

double fastpm_2pcf_eval(FastPM2PCF* self, double r)
{
    double rMax = self->size * self->step_size;
    if(r > rMax)
        return 0.0;
    if(r == rMax)
        return self->xi[self->size];
    int i = (int)floor(r / self->step_size);
    double prev = self->xi[i], next = self->xi[i + 1];
    double deltaR = r - i * self->step_size;
    return prev + (next - prev) * deltaR / self->step_size;
}

/*
void
fastpm_generate_covariance_matrix(PM *pm, fastpm_fkfunc pkfunc, void * data, FastPMFloat *cov_x)
{
}
*/
void
fastpm_2pcf_from_powerspectrum(FastPM2PCF *self, fastpm_fkfunc pkfunc, void * data, double r_max, int steps)
{
    self->size = steps;
    self->step_size = r_max / steps;
    self->xi = malloc((steps + 1) * sizeof(double));

    double logKMin = -10, logKMax = 5;
    double logKSteps = 10000;
    double logKStepSize = (logKMax - logKMin) / logKSteps;
    double pi = 3.141593;
    int i, j;
    for(i = 0; i <= steps; ++i)
    {
        double r = i * self->step_size;
        double res = 0;
        double prev = 0;
        for(j = 1; j <= logKSteps; ++j)
        {
            double logK = logKMin + j * logKStepSize;
            double k = exp(logK);
            double kr = k * r;
            double func = 1;
            if(kr > 0)
                func = sin(kr) / kr;
            func *= pkfunc(k, data) * k * k * k;
            res += (prev + func) / 2;
            prev = func;
        }
        res *= logKStepSize / (2 * pi * pi);
        self->xi[i] = res;
    }
}

static void
_solve(int size, double * Cij, double * dfi, double * x)
{
    gsl_matrix_view m = gsl_matrix_view_array(Cij, size, size);
    gsl_vector_view b = gsl_vector_view_array(dfi, size);
    gsl_vector_view xx = gsl_vector_view_array(x, size);
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, &xx.vector);
}

static void
_readout(FastPMConstraint * constraints, int size, PM * pm, FastPMFloat * delta_x, double * dfi)
{
    int i;

    for(i = 0; i < size; ++i)
    {
        int ii[3];
        int d;
        int inBox = 1;
        int index = 0;
        for(d = 0; d < 3; ++d)
        {
            ii[d] = constraints[i].x[d] * pm->InvCellSize[d] - pm->IRegion.start[d];
            if(ii[d] < 0 || ii[d] > pm->IRegion.size[d])
                inBox = 0;
            index += ii[d] * pm->IRegion.strides[d];
        }

        if(inBox) {
            dfi[i] = delta_x[index];
        } else {
            dfi[i] = 0;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, dfi, size, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
}
static double
_sigma(PM * pm, FastPMFloat * delta_x)
{
    double d2 = 0.0;
    PMXIter xiter;
    for(pm_xiter_init(pm, &xiter);
       !pm_xiter_stop(&xiter);
        pm_xiter_next(&xiter))
    {
        double od = delta_x[xiter.ind] - 1;
        d2 += od * od;
    }
    MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    /* unbiased estimator of the variance. the mean is 1. */
    d2 /= (pm_norm(pm) - 1);
    return sqrt(d2);
}

void
fastpm_cg_apply_constraints(FastPMConstrainedGaussian *cg, PM * pm, FastPM2PCF *xi, FastPMFloat * delta_k)
{
    FastPMConstraint * constraints = cg->constraints;

    int size;
    for(size = 0; constraints[size].x[0] >= 0; size ++)
        continue;

    int i;
    fastpm_info("Constrained Gaussian with %d constraints\n", size);
    double dfi[size];
    double e[size];
    double Cij[size * size];
    double sigma = 0;
    FastPMFloat * delta_x = pm_alloc(pm);

    pm_assign(pm, delta_k, delta_x);
    pm_c2r(pm, delta_x);

    sigma = _sigma(pm, delta_x);
    fastpm_info("Measured sigma on the grid = %g\n", sigma);

    _readout(constraints, size, pm, delta_x, dfi);

    for(i = 0; i < size; ++i)
    {
        dfi[i] = (1 + constraints[i].c * sigma) - dfi[i];

        int j;
        for(j = i; j < size; ++j)
        {
            int d;
            double r = 0;
            for(d = 0; d < 3; ++d) {
                double dx = constraints[i].x[d] - constraints[j].x[d];
                
                if(dx > 0.5*pm->BoxSize[d]){
                    dx -= pm->BoxSize[d];
                }
                else if(dx < -0.5*pm->BoxSize[d]){
                    dx += pm->BoxSize[d];
                }
                
                r += dx * dx;
            }
            r = sqrt(r);
            double v = fastpm_2pcf_eval(xi, r);
            Cij[i * size + j] = v;
            Cij[j * size + i] = v;
        }
    }

    _solve(size, Cij, dfi, e);

    PMXIter xiter;
    for(pm_xiter_init(pm, &xiter);
       !pm_xiter_stop(&xiter);
        pm_xiter_next(&xiter))
    {
        double v = 0;
        for(i = 0; i < size; ++i)
        {
            int d;

            double r = 0;
            for(d = 0; d < 3; d ++) {
                double dx = xiter.iabs[d] * pm->CellSize[d] - constraints[i].x[d];
                
                if(dx > 0.5*pm->BoxSize[d]){
                    dx -= pm->BoxSize[d];
                }
                else if(dx < -0.5*pm->BoxSize[d]){
                    dx += pm->BoxSize[d];
                }
                
                r += dx * dx;
            }
            r = sqrt(r);

            v += e[i] * fastpm_2pcf_eval(xi, r);
        }
        delta_x[xiter.ind] += v;
    }

    _readout(constraints, size, pm, delta_x, dfi);
    for(i = 0; i < size; i ++) {
        fastpm_info("After constraints, Realization x[] = %g %g %g overdensity = %g, peak-sigma= %g\n",
                constraints[i].x[0],
                constraints[i].x[1],
                constraints[i].x[2],
                (dfi[i] - 1.0), (dfi[i] - 1.0) / sigma);
    }
    pm_r2c(pm, delta_x, delta_k);
    pm_free(pm, delta_x);
}
