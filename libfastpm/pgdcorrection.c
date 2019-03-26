#include <math.h>
#include <stdlib.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/prof.h>
#include "pmpfft.h"
#include "pmghosts.h"

double
fastpm_pgdc_get_alpha(FastPMPGDCorrection * pgdc, double a)
{
    return pgdc->alpha0*pow(10, pgdc->A*a*a-pgdc->B*a); 
}

double
fastpm_pgdc_get_ks(FastPMPGDCorrection * pgdc, double a)
{
    return pgdc->ks;
}

double
fastpm_pgdc_get_kl(FastPMPGDCorrection * pgdc, double a)
{
    return pgdc->kl;
}

static void
apply_pgdpot_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double alpha, double kl, double ks)
{
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
	double kl2 = kl*kl;
	double ks4 = ks*ks*ks*ks;
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk = 0;
            for(d = 0; d < 3; d++) {
                kk += kiter.kk[d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* - 1 / k2 */
            if(LIKELY(kk > 0)) {
		double fac = alpha * exp(-kl2/kk-kk*kk/ks4) / kk;
                to[ind + 0] = fac * from[ind + 0];
                to[ind + 1] = fac * from[ind + 1];
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

void
fastpm_pgdc_calculate(FastPMPGDCorrection * pgdc,
    PM * pm,
    FastPMStore * p,
    FastPMFloat * delta_k, double a, double fac)
{
    FastPMPainter reader[1];

    fastpm_painter_init(reader, pm, pgdc->PainterType, pgdc->PainterSupport);

    /* watch out: boost the density since mesh is finer than grid */
    long long np = p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG_LONG, MPI_SUM, pm_comm(pm));

    CLOCK(ghosts);
    PMGhostData * pgd = pm_ghosts_create(pm, p, p->attributes, reader->support);
    pm_ghosts_send(pgd, COLUMN_POS);
    LEAVE(ghosts);

    int d;

    FastPMFieldDescr C[] = {
         {COLUMN_PGDC, 0},
         {COLUMN_PGDC, 1},
         {COLUMN_PGDC, 2},
    };

    CLOCK(transfer);
    CLOCK(c2r);
    CLOCK(readout);
    CLOCK(reduce);

    FastPMFloat * canvas = pm_alloc(pm);

    double kl = fastpm_pgdc_get_kl(pgdc, a);
    double ks = fastpm_pgdc_get_ks(pgdc, a);
    double alpha = fastpm_pgdc_get_alpha(pgdc, a)*fac;

    for(d = 0; d < 3; d ++) {

        ENTER(transfer);
        /* apply the transfers */
        apply_pgdpot_transfer(pm, delta_k, canvas, alpha, kl, ks);
        fastpm_apply_diff_transfer(pm, canvas, canvas, C[d].memb);
        /* result saved to canvas. */
        LEAVE(transfer);

        ENTER(c2r);
        pm_c2r(pm, canvas);
        LEAVE(c2r);

        ENTER(readout);
        fastpm_readout_local(reader, canvas, p, p->np, C[d]);
        fastpm_readout_local(reader, canvas, pgd->p, pgd->p->np, C[d]);
        LEAVE(readout);

    }

    ENTER(reduce);
    pm_ghosts_reduce(pgd, COLUMN_PGDC, FastPMReduceAddFloat, NULL);
    LEAVE(reduce);

    pm_free(pm, canvas);

    pm_ghosts_free(pgd);
}

