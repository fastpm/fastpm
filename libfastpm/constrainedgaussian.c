#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/constrainedgaussian.h>

#include "pmpfft.h"
#include <gsl/gsl_linalg.h>

double fastpm_2pcf_eval(FastPM2PCF* self, double r)
{
    return 1.0;
}

/*
void
fastpm_generate_covariance_matrix(PM *pm, fastpm_fkfunc pkfunc, void * data, FastPMFloat *cov_x)
{
}
*/
void
fastpm_2pcf_from_powerspectrum(FastPM2PCF *self, fastpm_fkfunc pkfunc, void * data)
{
}

static void
_solve(int size, double * Cij, double * dfi, double * x)
{


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

void
fastpm_cg_induce_correlation(FastPMConstrainedGaussian *cg, PM * pm, FastPM2PCF *xi, FastPMFloat * delta_k)
{
    FastPMConstraint * constraints = cg->constraints;

    int size;
    for(size = 0; constraints[size].x[0] >= 0; size ++)
        continue;

    int i;
    fastpm_info("Constrained Gaussian with %d constraints\n", size);
    for(i = 0; i < size; i ++) {
        fastpm_info("x[] = %g %g %g ; c = %g\n",
                constraints[i].x[0],
                constraints[i].x[1],
                constraints[i].x[2],
                constraints[i].c);
    }

    double dfi[size];
    double e[size];
    double Cij[size * size];

    FastPMFloat * delta_x = pm_alloc(pm);

    pm_assign(pm, delta_k, delta_x);
    pm_c2r(pm, delta_x);

    _readout(constraints, size, pm, delta_x, dfi);

    for(i = 0; i < size; ++i)
    {
        dfi[i] = constraints[i].c - dfi[i];

        int j;
        for(j = i; j < size; ++j)
        {
            int d;
            double r = 0;
            for(d = 0; d < 3; ++d) {
                double dx = constraints[i].x[d] - constraints[j].x[d];
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
                r += dx * dx;
            }
            r = sqrt(r);

            v += e[i] * fastpm_2pcf_eval(xi, r);
        }
        delta_x[xiter.ind] += v;
    }

    _readout(constraints, size, pm, delta_x, dfi);
    for(i = 0; i < size; i ++) {
        fastpm_info("Realization x[] = %g %g %g v = %g\n",
                constraints[i].x[0],
                constraints[i].x[1],
                constraints[i].x[2],
                dfi[i]);
    }
    pm_r2c(pm, delta_x, delta_k);
    pm_free(pm, delta_x);
}
