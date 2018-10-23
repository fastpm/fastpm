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
                uint64_t **pix,
                size_t * n,
                FastPMLightCone * lc,
                MPI_Comm comm
                );

void
fastpm_gldot(double glmatrix[4][4], double xi[4], double xo[4]);

void
fastpm_gldotf(double glmatrix[4][4], float vi[4], float vo[4]);

int
MPIU_Any(MPI_Comm comm, int value);

int
MPIU_All(MPI_Comm comm, int value);

int
MPIU_stats(MPI_Comm comm, double value, const char * fmt, ...);

char *
MPIU_Bcast_string(MPI_Comm comm, char * src, int root, void (*free)(void * ptr));

FASTPM_END_DECLS

