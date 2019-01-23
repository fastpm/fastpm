typedef struct FastPMncdmInitData{
    double m_ncdm;
    double z;
    int n_shells;
    int n_side;
    double (* vel)[3];
    double * mass;
} FastPMncdmInitData;

FastPMncdmInitData* 
fastpm_ncdm_init_create(double m_ncdm, double z, int n_shells, int n_side);

void fastpm_ncdm_init_free(FastPMncdmInitData* nid);

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