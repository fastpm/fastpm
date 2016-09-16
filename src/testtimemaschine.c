#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/cosmology.h>
    

int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    //~ entry template[] = {
    //~ {0, 0, 1}, /* Kick */
    //~ {0, 1, 1}, /* Drift */
    //~ {0, 2, 1}, /* Drift */
    //~ {2, 2, 1}, /* Force */
    //~ {2, 2, 2}, /* Kick */
    //~ {-1, -1, -1} /* End of table */
    //~ };
    //~ entry template1[] = {
    //~ {0, 1, 0}, /* Drift */
    //~ {0, 1, 1}, /* Kick */
    //~ {1, 1, 1}, /* Force */
    //~ {-1, -1, -1} /* End of table */
    //~ };

    //~ entry template[] = {
    //~ {0, 0, 1}, /* Kick*/
    //~ {0, 1, 1}, /* Drift*/
    //~ {1, 1, 1}, /* Force */
    //~ {-1, -1, -1} /* End of table */
    //~ };
    
    FastPMTEEntry template[] = {
    {0, 0, 1}, /* Kick */
    {0, 1, 1}, /* Drift */
    {0, 2, 1}, /* Drift */
    {2, 2, 1}, /* Force */
    {2, 2, 2}, /* Kick */
    {-1, -1, -1} /* End of table */
    };

    double timesteps[5] = {0.2, 0.4, 0.6, 0.8, 1.0};
    FastPMTEStates *states = malloc(sizeof(FastPMTEStates));
    
    fastpm_tevo_generate_states(states, 5, template, timesteps);

    fastpm_tevo_print_states(states);
    fastpm_tevo_destroy_states(states);

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    libfastpm_cleanup();

    MPI_Finalize();
    
    return 0;
}
