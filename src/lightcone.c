typedef struct {
    double LightSpeedFactor;

    /* storage of the particles on light cone */
    FastPMStore * p;
    /* need a table for drift factors */
    
} FastPMLightCone;

void fastpm_lc_init(FastPMLightCone * lc, size_t np_upper)
{
    /* allocation */
}

void fastpm_lc_destroy(FastPMLightCone * lc)
{

    /* free */
}

void fastpm_lc_write(FastPMLightCone * lc, const char * filename)
{

}

double fastpm_lc_horizon(FastPMLightCone * lc, double a)
{

}

double _unknown(double a, void * params)
{

    double xo[3];
    double horizon;

    horizon = fastpm_lc_horizon(lc, a);

    fastpm_drift_one(drift, p, i, xo, a);

    return horizon - xo[2];
}

int fastpm_lc_intersect(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p, int i, double * solution)
{

    /* XXX: culling */

    gsl_fsolve(solution);

    if( *solution >= drift->a_i
     && *solution <= drift->a_f)
    {
        return 1;
    }
    return 0;
}

int fastpm_lc_scan(FastPMLightCone * lc, FastPMDrift * drift, FastPMStore * p)
{
    for(i in p) {
        status = fastpm_lc_intersect(lc, drift, p, i, solution)
        if(dint is good)
            add(i, lc);
    }
}
