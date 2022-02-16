#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

// function to calculate subsampling rate
double calc_subrate(double ell_lim, double z, double res_box, FastPMCosmology * c)
{
    double res_lim = VolumeDensityFromEll(ell_lim, z, c);

    double rate_sub = res_lim / res_box;
    /* FIXME: in principle we can replicate particles to achieve a rate > 1.
     * probably want to move this clipping to the caller side.*/
    if (rate_sub > 1) {
        rate_sub = 1; // for low redshift, no subsample
    }
    return rate_sub;
}

int main(int argc, char * argv[]){
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    
    // Init the FD interpolation object
    FastPMFDInterp FDi;
    fastpm_fd_interp_init(&FDi);
    
    FastPMCosmology c[1] = {{
        .h = 0.6774,
        .Omega_m = 0.309,
        .T_cmb = 0,
        .N_eff = 3.04,
        .Omega_Lambda = 0.691,   // FIXME Note this test cosmology isn't properly closed   0.682407
        .Omega_cdm = 0.309,
        .wa = 0,
        .w0 = -1,
    }};

    fastpm_cosmology_init(c);
    FILE *f;
    f = fopen("subs_rate.txt","w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    //print box size, total particle number, resolution of the system
    double box_size_Mpch = 5000;
    int Npart = pow(2,13);
    double res_box = pow(Npart / box_size_Mpch, 3);
    double scale_f = 1 / (1 + 4.0);
    double comoving_l = ComovingDistance(scale_f, c) * HubbleDistance;

    //printf("%f: comoving distance\n",comoving_l);
    double Npart_tot = 4 * M_PI / 3 * pow(comoving_l, 3) * res_box;
    printf("%f: total particle number\n", Npart_tot);

    double output1 = Npart_tot / pow(10,12);
    double output2 = Npart_tot / pow(Npart, 3);
    printf("%f: particle/(Mpc/h)^3 resolution\n",res_box);
    printf("box size: %g (Mpc/h)\n",box_size_Mpch);
    printf("%f x 10^12 particles, %f x box size\n",output1,output2);

    //read ell_lim from command line
    double ell_lim;
    if (argc == 2){
        ell_lim = atof(argv[1]); 
    }
    else
        ell_lim = 2000;
    
    //calculate subsampling based on subsampling.c
    for (double z=0.1 ; z<=4 ; z+=0.01){
        double subsamplingrate = calc_subrate(ell_lim,z,res_box,c);
        fprintf(f, "%7.6lf, %10.7lf\n",z, subsamplingrate);
    }
    
    fastpm_cosmology_destroy(c);
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;

}
