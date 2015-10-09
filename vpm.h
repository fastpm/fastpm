
/* Variable Particle Mesh */

typedef struct {
    PM pm;
    double a_start;
    int pm_nc_factor;
    int end;
} VPM;

VPM * 
vpm_create (int size, int * pm_nc_factors, double * change_pm, 
        PMInit * baseinit, PMIFace * iface, MPI_Comm comm);

VPM *
vpm_find(VPM * vpm, double a);

double
vpm_estimate_alloc_factor(double failure_rate);

