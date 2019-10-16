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
    
    FastPMCosmology c[1] = {{
        .h = 0.7,
        .Omega_cdm = 0.3,
        .T_cmb = 2.73,
        .N_eff = 3.046,
        .N_nu = 3,
        .m_ncdm = {0.05,0.05,0.05},
        .N_ncdm = 3
    }};

    FILE * pFile;
    pFile = fopen ("ncdm_test.txt","w");
    
    for (double loga=-4 ; loga<=0 ; loga+=0.01){
        double a = pow(10, loga);
        double On = Omega_ncdm(a, c);
        double Onm = Omega_ncdm_m(a, c);
        double Oc = Omega_cdm_a(a, c);
        double Om = Omega_m(a, c);
        double w = w_ncdm_i(a,0,c);   //look at w of species 0 as example.
        
        fprintf (pFile, "%g %g %g %g %g %g\n", a, On, Onm, Oc, Om, w);
    }
    fclose (pFile);
    
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;

}