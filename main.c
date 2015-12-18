#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <mpi.h>
#include <math.h>
#include <signal.h>
#include <getopt.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

#include "parameters.h"
#include "readparams.h"
#include "libfastpm.h"

#include "walltime.h"
#include "power.h"
#include "msg.h"
#include "pmic.h"


/* command-line arguments */
static char * ParamFileName;

static void 
parse_args(int argc, char ** argv, Parameters * prr);

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, double aout, void * template);

static void 
_mkdir(const char *dir);
static void 
ensure_dir(char * path);

int 
fastpm_read_runpb_ic(FastPM * fastpm, PMStore * p, char * filename);

int 
fastpm_write_runpb_snapshot(FastPM * fastpm, PMStore * p, char * filebase);

int fastpm(Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD; 

    Parameters prr;

    parse_args(argc, argv, &prr);

    struct ClockTable CT;
    msg_init(comm);
    walltime_init(&CT);

    msg_set_loglevel(verbose);

    read_parameters(ParamFileName, &prr, comm);

    fastpm(&prr, comm);

    walltime_summary(0, comm);

    msg_printf(info, "Total Time\n");
    walltime_report(stdout, 0, comm);

    MPI_Finalize();
    return 0;
}

static int default_step_action(FastPM * fastpm, Parameters * prr);
static int 
measure_power_action(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr);

int fastpm(Parameters * prr, MPI_Comm comm) {

    const double rho_crit = 27.7455;
    const double M0 = prr->omega_m*rho_crit*pow(prr->boxsize / prr->nc, 3.0);
    msg_printf(verbose, "mass of a particle is %g 1e10 Msun/h\n", M0); 

    FastPM * fastpm = & (FastPM) {
        .time_step = prr->time_step,
        .n_time_step = prr->n_time_step,
        .nc = prr->nc,
        .boxsize = prr->boxsize,
        .omega_m = prr->omega_m,
        .USE_COLA = prr->force_mode == FORCE_MODE_COLA,
        .USE_NONSTDDA = !prr->cola_stdda,
        .USE_LINEAR_THEORY = prr->enforce_broadband,
        .nLPT = -2.5f,
        .K_LINEAR = prr->enforce_broadband_kmax,
    };

    fastpm_init(fastpm, 
        prr->np_alloc_factor, prr->NprocY, prr->UseFFTW, 
        prr->pm_nc_factor, prr->n_pm_nc_factor, prr->change_pm, 
        comm);

    walltime_measure("/Init/Misc");

    if(prr->readic_filename) {
        fastpm_read_runpb_ic(fastpm, fastpm->p, prr->readic_filename);

        walltime_measure("/Init/ReadIC");
    } else {
        FastPMFloat * delta_k = pm_alloc(fastpm->pm_2lpt);

        power_init(prr->power_spectrum_filename, 
                prr->time_step[0], 
                prr->sigma8, 
                prr->omega_m, 
                1 - prr->omega_m, comm);

        if(prr->readnoise_filename) {
            pmic_read_gaussian(fastpm->pm_2lpt, delta_k, prr->readnoise_filename, PowerSpecWithData, NULL);
        } else {
            pmic_fill_gaussian_gadget(fastpm->pm_2lpt, delta_k, prr->random_seed, PowerSpecWithData, NULL);
        }
        fastpm_prepare_ic(fastpm, delta_k);

        pm_free(fastpm->pm_2lpt, delta_k);

        walltime_measure("/Init/2LPT");
    }

    fastpm_evolve(fastpm,
            (fastpm_after_force_action) measure_power_action, 
            (fastpm_after_drift_action) default_step_action, 
            (fastpm_after_kick_action) default_step_action,
            prr
            );
    fastpm_destroy(fastpm);

    return 0;
}

static int default_step_action(FastPM * fastpm, Parameters * prr) {
    char TEMPLATE[1024];
    sprintf(TEMPLATE, "%s%05d_%%0.04f.bin", prr->snapshot_filename, prr->random_seed);
    fastpm_interp(fastpm, prr->aout, prr->n_aout, take_a_snapshot, TEMPLATE);
    return 0;
}

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, double aout, void * template) 
{
    char filebase[1024];
    double z_out= 1.0/aout - 1.0;

    sprintf(filebase, template, aout);
    ensure_dir(filebase);
    fastpm_write_runpb_snapshot(fastpm, snapshot, filebase);
    
    walltime_measure("/Snapshot/IO");
    MPI_Barrier(fastpm->comm);
    walltime_measure("/Snapshot/Wait");

    msg_printf(normal, "snapshot %s written z = %6.4f a = %6.4f\n", 
            filebase, z_out, aout);
    return 0;
}

static int 
measure_power_action(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr) 
{
    PowerSpectrum ps;
    /* calculate the power spectrum */
    power_spectrum_init(&ps, pm_nmesh(fastpm->pm)[0] / 2);

    pm_calculate_powerspectrum(fastpm->pm, delta_k, &ps);

    walltime_measure("/PowerSpectrum/Measure");

    if(prr->measure_power_spectrum_filename) {
        if(fastpm->ThisTask == 0) {
            ensure_dir(prr->measure_power_spectrum_filename);
            power_spectrum_write(&ps, fastpm->pm, pow(fastpm->nc, 3.0),
                prr->measure_power_spectrum_filename, prr->random_seed, a_x);
        }
    }
    power_spectrum_destroy(&ps);
    MPI_Barrier(fastpm->comm);
    walltime_measure("/PowerSpectrum/Write");
    return 0;
}

static void 
parse_args(int argc, char ** argv, Parameters * prr) 
{
    char opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    ParamFileName = NULL;
    prr->NprocY = 0;    
    while ((opt = getopt(argc, argv, "h?y:f")) != -1) {
        switch(opt) {
            case 'y':
                prr->NprocY = atoi(optarg);
            break;
            case 'f':
                prr->UseFFTW = 1;
            break;
            case 'h':
            case '?':
            default:
                goto usage;
            break;
        }
    }
    if(optind >= argc) {
        goto usage;
    }

    ParamFileName = argv[optind];
    return;

usage:
    printf("Usage: fastpm [-f] [-y NprocY] paramfile\n"
    "-f Use FFTW \n"
    "-y Set the number of processes in the 2D mesh\n"
);
    MPI_Finalize();
    exit(1);
}

static void 
ensure_dir(char * path) 
{
    int i = strlen(path);
    char * dup = strdup(path);
    char * p;
    for(p = i + dup; p >= dup && *p != '/'; p --) {
        continue;
    }
    /* plain file name in current directory */
    if(p < dup) return;
    
    /* p == '/', so set it to NULL, dup is the dirname */
    *p = 0;
    _mkdir(dup);
    free(dup);
}

static void 
_mkdir(const char *dir) 
{
        char * tmp = strdup(dir);
        char *p = NULL;
        size_t len;

        len = strlen(tmp);
        if(tmp[len - 1] == '/')
                tmp[len - 1] = 0;
        for(p = tmp + 1; *p; p++)
                if(*p == '/') {
                        *p = 0;
                        mkdir(tmp, S_IRWXU | S_IRWXG | S_IRWXO);
                        *p = '/';
                }
        mkdir(tmp, S_IRWXU | S_IRWXG | S_IRWXO);
        free(tmp);
}
