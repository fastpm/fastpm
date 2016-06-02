/* libfastpm: */
#include <mpi.h>
#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/transfer.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "pm2lpt.h"
#include "pmghosts.h"

int 
fastpm_2lpt_init(FastPM2LPTSolver * solver, int nmesh, int nc, double boxsize, double alloc_factor, MPI_Comm comm)
{
    FastPMSolverBase * base = &solver->base;
    base->p = malloc(sizeof(PMStore));
    base->pm = malloc(sizeof(PM));

    pm_store_init(base->p);

    pm_store_alloc_evenly(base->p, pow(nc, 3),
        PACK_POS | PACK_VEL | PACK_ID | PACK_ACC | PACK_DX1 | PACK_DX2 | PACK_Q,
        alloc_factor, comm);

    solver->nc = nc;
    solver->boxsize = boxsize;

    PMInit pminit = {
        .Nmesh = nmesh,
        .BoxSize = boxsize,
        .NprocY = 0,
        .transposed = 1,
        .use_fftw = 0,
    };

    pm_init(base->pm, &pminit, comm);
    base->comm = comm;
    MPI_Comm_size(comm, &base->NTask);
    MPI_Comm_rank(comm, &base->ThisTask);

    return 0;
}

void
fastpm_2lpt_destroy(FastPM2LPTSolver * solver)
{
    FastPMSolverBase * base = &solver->base;

    pm_destroy(base->pm);
    pm_store_destroy(base->p);
    free(base->pm);
    free(base->p);
}

void
fastpm_2lpt_evolve(FastPM2LPTSolver * solver,
        FastPMFloat * delta_k_i, double aout, double omega_m)
{
    FastPMSolverBase * base = &solver->base;
    /* evolve particles by 2lpt to time a. pm->canvas contains rho(x, a) */
    double shift0;
    if(solver->USE_SHIFT) {
        shift0 = solver->boxsize / solver->nc * 0.5;
    } else {
        shift0 = 0;
    }

    double shift[3] = {shift0, shift0, shift0};
    int nc[3] = {solver->nc, solver->nc, solver->nc};

    pm_store_set_lagrangian_position(base->p, base->pm, shift, nc);

    pm_2lpt_solve(base->pm, delta_k_i, base->p, shift);

    if(solver->USE_DX1_ONLY) {
        ptrdiff_t i;
        for(i = 0; i < base->p->np; i ++) {
            int d;
            for(d = 0; d < 3; d ++) {
                base->p->dx2[i][d] = 0;
            }
        }
    }
    /* pdata->dx1 and pdata->dx2 are s1 and s2 terms 
     * S = D * dx1 + D2 * 3 / 7 * D20 * dx2; 
     *
     * See pmsteps.c 
     * */

    /* now shift particles to the correct locations. */

    /* predict particle positions by 2lpt */
    pm_2lpt_evolve(aout, base->p, omega_m, solver->USE_DX1_ONLY);
}
