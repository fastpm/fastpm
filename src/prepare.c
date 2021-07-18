#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <signal.h>
#include <getopt.h>
#include <limits.h>

#include <fastpm/libfastpm.h>

#include "lua-config.h"
#include "param.h"
#include "prepare.h"

void
prepare_cosmology(FastPMCosmology * c, LUAParameters * lua) {
    c->h = CONF(lua, h);
    c->Omega_m = CONF(lua, Omega_m);
    c->T_cmb = CONF(lua, T_cmb);
    c->Omega_k = CONF(lua, Omega_k);
    c->w0 = CONF(lua, w0);
    c->wa = CONF(lua, wa);
    c->N_eff = CONF(lua, N_eff);
    c->N_nu = CONF(lua, N_nu);
    c->N_ncdm = CONF(lua, n_m_ncdm);
    c->ncdm_matterlike = CONF(lua, ncdm_matterlike);
    c->ncdm_freestreaming = CONF(lua, ncdm_freestreaming);
    c->growth_mode = CONF(lua, growth_mode);

    int i;
    for(i = 0; i < CONF(lua, n_m_ncdm); i ++) {
        c->m_ncdm[i] = CONF(lua, m_ncdm)[i];
    }

}
