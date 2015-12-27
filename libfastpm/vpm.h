
/* Variable Particle Mesh */

typedef struct VPM {
    PM pm;
    double a_start;
    int pm_nc_factor;
    int end;
} VPM;

VPM * 
vpm_create (VPMInit * vpminit, PMInit * baseinit, PMIFace * iface, MPI_Comm comm);

VPM *
vpm_find(VPM * vpm, double a);

double
vpm_estimate_alloc_factor(double failure_rate);

void vpm_free (VPM * vpm);
