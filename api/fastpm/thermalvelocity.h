typedef enum {
    FASTPM_NCDM_SPHERE_HEALPIX = 0,
    FASTPM_NCDM_SPHERE_FIBONACCI = 1,
} FastPMncdmSphereScheme;

typedef struct FastPMncdmInitData{
    double BoxSize;
    double Omega_ncdm;
    double m_ncdm[3];
    int n_ncdm;        /* number of ncdm species, each has a m_ncdm */
    double m_ncdm_sum; /* total eV mass of all ncdm */
    double z;          /* initialization redshift of ncdm species */

    int n_shells;    /* number of shells for splitting velocity magnitude distribution */
    FastPMncdmSphereScheme ncdm_sphere_scheme;
    int n_side;      /* n_side for healpix, also n_fib for fibonacci */
    size_t n_sphere; /* number of tiles on the sphere */
    size_t n_split;  /* total number of phase space splits */
    int lvk;         /* bool: switch on low velocity kernel? */
    
    /* a table for look up of the coherent thermal velocity and mass of split particles */
    double (* vel)[3];
    double * mass;
} FastPMncdmInitData;

FastPMncdmInitData *
fastpm_ncdm_init_create(double BoxSize, FastPMCosmology * c, double z, int n_shells, int n_side,
                        int lvk, FastPMncdmSphereScheme ncdm_sphere_scheme);

void
fastpm_ncdm_init_free(FastPMncdmInitData* nid);

void
fastpm_split_ncdm(FastPMncdmInitData* nid, FastPMStore * src, FastPMStore * dest, MPI_Comm comm);
