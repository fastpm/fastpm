#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

int main(int argc, char * argv[]){
    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    
    // Init the FD interpolation object
    FastPMFDInterp FDi;
    fastpm_fd_interp_init(&FDi);

    FastPMCosmology c[1] = {{
        .h = 0.6711,
        .Omega_m = 0.3175,
        .T_cmb = 2.7255,
        .N_eff = 3.046,
        .N_nu = 3,
        .N_ncdm = 3,                    //0,
        .m_ncdm = {0.0666667, 0.0666667, 0.0666667}, //{},
        .Omega_Lambda = 0.682445,   // FIXME Note this test cosmology isn't properly closed   0.682407
        .Omega_cdm = 0.312725,
        .FDinterp = &FDi,
    }};

    FILE * pFile;
    pFile = fopen ("cosmology_test_mncdm20.txt","w");
    
    for (double loga=-6 ; loga<=0 ; loga+=0.001){
        double a = pow(10, loga);
        double z = 1./a - 1;
        double E = HubbleEa(a, c);
        double Oc = Omega_cdm_a(a, c);
        double Og = Omega_g(c) / (a*a*a*a) / (E*E);
        double Our = Omega_ur(c) / (a*a*a*a) / (E*E);
        double On = Omega_ncdmTimesHubbleEaSq(a, c) / (E*E);
        double Om = Omega_m(a, c);
        double OL = c->Omega_Lambda / (E*E);
        //double w = w_ncdm_i(a,0,c);   //look at w of species 0 as example.
        
        fprintf (pFile, "%g %g %g %g %g %g %g %g\n", z, E, Oc, Og, Our, On, Om, OL);
    }
    fclose (pFile);
    
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;

}