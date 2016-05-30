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
    solver->p = malloc(sizeof(PMStore));
    solver->pm = malloc(sizeof(PM));
    pm_store_init(solver->p);

    pm_store_alloc_evenly(solver->p, pow(nc, 3),
        PACK_POS | PACK_VEL | PACK_ID | PACK_ACC | PACK_DX1 | PACK_DX2 | PACK_Q,
        alloc_factor, comm);

    solver->nc = nc;
    solver->nmesh = nmesh;

    PMInit pminit = {
        .Nmesh = nmesh,
        .BoxSize = boxsize,
        .NprocY = 0,
        .transposed = 1,
        .use_fftw = 0,
    };

    pm_init(solver->pm, &pminit, comm);
    solver->boxsize = boxsize;
    solver->comm = comm;
    return 0;
}

void
fastpm_2lpt_destroy(FastPM2LPTSolver * solver)
{
    pm_destroy(solver->pm);
    pm_store_destroy(solver->p);
    free(solver->pm);
    free(solver->p);
}

void
fastpm_2lpt_evolve(FastPM2LPTSolver * solver,
        FastPMFloat * delta_k_i, double aout, double omega_m)
{
    /* evolve particles by 2lpt to time a. pm->canvas contains rho(x, a) */
    double shift0;
    if(solver->USE_SHIFT) {
        shift0 = solver->boxsize / solver->nc * 0.5;
    } else {
        shift0 = 0;
    }

    double shift[3] = {shift0, shift0, shift0};
    int nc[3] = {solver->nc, solver->nc, solver->nc};

    pm_store_set_lagrangian_position(solver->p, solver->pm, shift, nc);

    pm_2lpt_solve(solver->pm, delta_k_i, solver->p, shift);

    if(solver->USE_DX1_ONLY) {
        ptrdiff_t i;
        for(i = 0; i < solver->p->np; i ++) {
            int d;
            for(d = 0; d < 3; d ++) {
                solver->p->dx2[i][d] = 0;
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
    pm_2lpt_evolve(aout, solver->p, omega_m, solver->USE_DX1_ONLY);
}
