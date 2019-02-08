enum FastPMDivideSphere {     //FIX should implement this at some point
    FASTPM_DIVIDE_SPHERE_HEALPIX = 0,
    FASTPM_DIVIDE_SPHERE_FIBONACCI = 1,
};

typedef struct FastPMncdmInitData{
    double BoxSize;
    double Omega_ncdm;
    double m_ncdm[3];
    int n_ncdm; /* number of ncdm species, each has a m_ncdm */
    double m_ncdm_sum; /* total ev mass of all ncdm */
    double z;   /* initialization redshift of ncdm species */

    int n_shells;
    int lvk;    /* bool: switch on low velocity kernel? */

    union {
        int n_side; /* for healpix splits */
        int n_fib;  /* for fib splits */
    };

    size_t n_split; /* number of phase space splits of each initial position */
    /* a table for quick look up of the coherent thermal velocity and mass of split particles. */
    double (* vel)[3];
    double * mass;
} FastPMncdmInitData;

FastPMncdmInitData *
fastpm_ncdm_init_create(double BoxSize, double m_ncdm[3], int n_ncdm, double z, int n_shells, int n_side, int lvk);

void
fastpm_ncdm_init_free(FastPMncdmInitData* nid);

void
fastpm_split_ncdm(FastPMncdmInitData* nid, FastPMStore * src, FastPMStore * dest, int f_subsample_1d, MPI_Comm comm);

//unsigned int isqrt(int number);   //c
//void pix2vec (int pix, double *vec, int n_side);  //c do i need ot put in c?

//double fermi_dirac_kernel(double x);
//double low_vel_kernel(double x);
//double fermi_dirac_dispersion(double x);

//void divide_fd(double *vel_table, double *mass, int n_shells);

//void divide_sphere_healpix(double *vec_table, int n_side);
//void divide_sphere_fibonacci(double *vec_table, int n_fibonacci);

////double *create_nu_table_healpix(double m_nu, double redshift, int n_shells, int n_side);
////double *create_nu_table_fibonacci(double m_nu, double Redshift, int n_shells, int n_fibonacci);
