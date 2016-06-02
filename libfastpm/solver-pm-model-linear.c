#include <alloca.h>
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
    double Plin;
} FastPMModelLinearPriv;

static void fastpm_model_linear_build(FastPMModel * model, double ainit, double afinal)
{
    FastPMModelLinearPriv * priv = model->priv;
    PMStore * psub = alloca(sizeof(PMStore));
    fastpm_model_create_subsample(model, psub);
    pm_2lpt_evolve(ainit, psub, model->fastpm->omega_m, 0);
    pm_store_wrap(psub, model->pm->BoxSize);
    priv->Plin = fastpm_model_measure_large_scale_power(model, psub);
    priv->Plin /= pow(fastpm_growth_factor(model->fastpm, ainit), 2);
    pm_store_destroy(psub);
}

static void fastpm_model_linear_evolve(FastPMModel * model, double af)
{
    FastPMModelLinearPriv * priv = model->priv;
    model->Pexpect = priv->Plin * pow(fastpm_growth_factor(model->fastpm, af), 2);
}

static void fastpm_model_linear_destroy(FastPMModel * model)
{
    free(model->priv);
}
void fastpm_model_linear_init(FastPMModel * model)
{
    model->build = fastpm_model_linear_build;
    model->evolve = fastpm_model_linear_evolve;
    model->destroy = fastpm_model_linear_destroy;
    model->priv = malloc(sizeof(FastPMModelLinearPriv));
}
