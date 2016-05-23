#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>

#include <fastpm/prof.h>
#include <fastpm/logging.h>

#include "pmpfft.h"
#include "pmstore.h"
#include "pmghosts.h"
#include "pm2lpt.h"
#include "solver-pm-internal.h"

typedef struct {
    double * NonLinearGrowthRate;
    double * NonLinearGrowthRateTime;
    double * steps;
    int nsteps;
    int istep;
} FastPMModelPMPriv;

static void _before_kick(FastPM * solver, FastPMKick * kick, void * userdata)
{
    FastPMModel * model = userdata;
    FastPMModelPMPriv * priv = model->priv;
    double P = fastpm_model_measure_large_scale_power(model, solver->p);
    priv->NonLinearGrowthRateTime[priv->istep] = solver->p->a_x;
    priv->NonLinearGrowthRate[priv->istep] = P;
    priv->istep ++;
};

static void fastpm_model_pm_build(FastPMModel * model, PMStore * p, double ainit, double afinal)
{
    FastPMModelPMPriv * priv = model->priv;
    PMStore * psub = malloc(sizeof(PMStore));

    fastpm_model_create_subsample(model, psub, PACK_POS | PACK_DX1 | PACK_DX2);

    FastPM * solver = &(FastPM) {
        .FORCE_TYPE = FASTPM_FORCE_PM,
        .USE_NONSTDDA = 0,
        .USE_MODEL = FASTPM_MODEL_NONE, /* this does not use any model */
        .K_LINEAR = model->fastpm->K_LINEAR,
        .nc = model->fastpm->nc / model->factor,
        .vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = model->pm->Nmesh[0] / (model->fastpm->nc / model->factor)},
            {.a_start = -1, .pm_nc_factor = 0}},
        .omega_m = model->fastpm->omega_m,
        .alloc_factor = model->fastpm->alloc_factor,
        .boxsize = model->fastpm->boxsize
    };
    double stepsize = 1.02;
    double a = ainit;
    int i;

    /* big enough to hold all */
    int nsteps_guess = (int) (log(afinal / ainit) / log(stepsize) + 10);

    priv->steps = malloc(sizeof(double) * nsteps_guess);
    /* We actually just save the mean power spectrum at large scale here */
    priv->NonLinearGrowthRate = malloc(sizeof(double) * nsteps_guess);
    priv->NonLinearGrowthRateTime = malloc(sizeof(double) * nsteps_guess);

    for(a = ainit, i = 0; a <= afinal * stepsize && a < 1.0 * stepsize; i ++) {
        if(a > 1.0) a = 1.0;
        priv->steps[i] = a;
        a *= stepsize;
    }
    priv->nsteps = i;

    fastpm_info("Calibrating Non-linear growth with a low resolution simulation.\n");
    /* From this point we run a small PM simulation internally, so mute the logging. */
    fastpm_push_msg_handler(fastpm_void_msg_handler, model->fastpm->comm, NULL);

    fastpm_init(solver, 0, 0, model->fastpm->comm);

    if(psub->np > solver->p->np_upper) {
        fastpm_raise(-1, "Too many subsample particles. This is due to an internal limitation of memory. Easy fix.\n");
    }

    memcpy(solver->p->x, psub->x, sizeof(psub->x[0]) * psub->np);
    memcpy(solver->p->dx1, psub->dx1, sizeof(psub->dx1[0]) * psub->np);
    memcpy(solver->p->dx2, psub->dx2, sizeof(psub->dx2[0]) * psub->np);

    solver->p->np = psub->np;
    fastpm_add_extension(solver, FASTPM_EXT_BEFORE_KICK, _before_kick, model);

    /* evolve and calculate expected large scale power */
    priv->istep = 0;
    fastpm_evolve(solver, priv->steps, priv->nsteps);

    fastpm_pop_msg_handler();

    fastpm_destroy(solver);
    pm_store_destroy(psub);
    model->Pexpect = priv->NonLinearGrowthRate[0];
}

static void fastpm_model_pm_evolve(FastPMModel * model, double af)
{
    int i;
    FastPMModelPMPriv * priv = model->priv;
    for(i = 0; i < priv->nsteps - 1; i ++) {
        if(af > priv->NonLinearGrowthRateTime[i] &&
          af <= priv->NonLinearGrowthRateTime[i + 1]) {
            break;
        }
    }

    if(i == priv->nsteps - 1) {
        model->Pexpect = priv->NonLinearGrowthRate[i];
    } else {
        double u = (af - priv->NonLinearGrowthRateTime[i]);
        u /= (priv->NonLinearGrowthRateTime[i + 1] - priv->NonLinearGrowthRateTime[i]);
        double v = 1 - u;
        double x = sqrt(priv->NonLinearGrowthRate[i]);
        double y = sqrt(priv->NonLinearGrowthRate[i + 1]);
        model->Pexpect = pow(u * y + v * x, 2);
    }
}
static void fastpm_model_pm_destroy(FastPMModel * model)
{
    FastPMModelPMPriv * priv = model->priv;
    free(priv->NonLinearGrowthRateTime);
    free(priv->NonLinearGrowthRate);
    free(priv->steps);
    free(priv);
}
void fastpm_model_pm_init(FastPMModel * model)
{
    model->build = fastpm_model_pm_build;
    model->evolve = fastpm_model_pm_evolve;
    model->destroy = fastpm_model_pm_destroy;
    model->priv = malloc(sizeof(FastPMModelPMPriv));
}
