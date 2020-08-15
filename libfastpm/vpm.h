
/* Variable Particle Mesh */

struct VPM {
    PM pm;
    double a_start;
    double pm_nc_factor;
    int end;
};

VPM * 
vpm_create (VPMInit * vpminit, PMInit * baseinit, MPI_Comm comm);

VPM *
vpm_find(VPM * vpm, double a);

double
vpm_estimate_alloc_factor(double failure_rate);

void vpm_free (VPM * vpm);
