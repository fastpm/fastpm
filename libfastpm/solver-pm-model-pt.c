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
    PMStore * psub;
} FastPMModelPTPriv;

static void fastpm_model_pt_build(FastPMModel * model, double ainit, double afinal)
{
    FastPMModelPTPriv * priv = model->priv;
    PMStore * psub = malloc(sizeof(PMStore));
    PMStore * p = model->fastpm->p;

    pm_store_init(psub);
    pm_store_alloc(psub, 1.0 * p->np_upper / model->factor, p->attributes);
    pm_store_create_subsample(psub, p, model->factor, model->fastpm->nc);
    pm_2lpt_evolve(ainit, psub, model->fastpm->omega_m, 0);
    pm_store_wrap(psub, model->pm->BoxSize);

    priv->psub = psub;
}

static void fastpm_model_pt_evolve(FastPMModel * model, double af)
{
    FastPMModelPTPriv * priv = model->priv;
    FastPMDrift drift;
    PMStore * psub = priv->psub;
    fastpm_drift_init(&drift, model->fastpm, psub, af);
    int i;
    for(i = 0; i < psub->np; i ++) {
        int d;
        for(d = 0; d < 3; d ++) {
            psub->x[i][d] += psub->dx1[i][d] * drift.da1;
            if(model->type == FASTPM_MODEL_2LPT)
                psub->x[i][d] += psub->dx2[i][d] * drift.da2;
        }
    }
    psub->a_x = af;
    psub->a_v = af;
    model->Pexpect = fastpm_model_measure_large_scale_power(model, psub);
}

static void fastpm_model_pt_destroy(FastPMModel * model)
{
    FastPMModelPTPriv * priv = model->priv;
    pm_store_destroy(priv->psub);
    free(priv->psub);
    free(model->priv);
}

void fastpm_model_pt_init(FastPMModel * model)
{
    model->build = fastpm_model_pt_build;
    model->evolve = fastpm_model_pt_evolve;
    model->destroy = fastpm_model_pt_destroy;
    model->priv = malloc(sizeof(FastPMModelPTPriv));
}
