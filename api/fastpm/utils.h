FASTPM_BEGIN_DECLS

struct fastpm_powerspec_eh_params {
    double hubble_param;
    double omegam;
    double omegab;
    double Norm;
};

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param); /* Eisenstein & Hu */

double 
fastpm_utils_powerspec_white(double k, double * amplitude); /* white noise. */

void
fastpm_utils_smooth(PM * pm, FastPMFloat * delta_x, FastPMFloat * delta_smooth, double sml);

void 
fastpm_utils_dump(PM * pm , const char * filename, FastPMFloat *data);

void 
fastpm_utils_load(PM * pm , const char * filename, FastPMFloat *data);

double 
fastpm_utils_get_random(uint64_t id);

void
fastpm_utils_healpix_ra_dec (
                int nside,
                double **ra,
                double **dec,
                size_t * n,
                double fov,
                MPI_Comm comm
                );

FASTPM_END_DECLS

