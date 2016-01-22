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

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>

#include "parameters.h"
#include "power.h"


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
read_runpb_ic(FastPM * fastpm, PMStore * p, char * filename);

void 
read_grafic_gaussian(PM * pm, FastPMFloat * g_x, char * filename);

int 
write_runpb_snapshot(FastPM * fastpm, PMStore * p, char * filebase);

int 
read_parameters(char * filename, Parameters * param, MPI_Comm comm);

int run_fastpm(FastPM * fastpm, Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    MPI_Init(&argc, &argv);

    Parameters prr;

    parse_args(argc, argv, &prr);

    MPI_Comm comm = MPI_COMM_WORLD; 

    libfastpm_init();

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    read_parameters(ParamFileName, &prr, comm);


    /* convert parameter files pm_nc_factor into VPMInit */
    VPMInit * vpminit = alloca(sizeof(VPMInit) * (prr.n_pm_nc_factor + 1));
    int i;
    for(i = 0; i < prr.n_pm_nc_factor; i ++) {
        vpminit[i].a_start = prr.change_pm[i];
        vpminit[i].pm_nc_factor = prr.pm_nc_factor[i];
    }
    /* mark the end */
    vpminit[i].pm_nc_factor = 0;

    FastPM * fastpm = & (FastPM) {
        .nc = prr.nc,
        .alloc_factor = prr.np_alloc_factor, 
        .vpminit = vpminit,
        .boxsize = prr.boxsize,
        .omega_m = prr.omega_m,
        .USE_COLA = prr.force_mode == FORCE_MODE_COLA,
        .USE_NONSTDDA = !prr.cola_stdda,
        .USE_LINEAR_THEORY = prr.enforce_broadband,
        .USE_DX1_ONLY = prr.use_dx1_only,
        .nLPT = -2.5f,
        .K_LINEAR = prr.enforce_broadband_kmax,
    };

    run_fastpm(fastpm, &prr, comm);

    libfastpm_cleanup();

    MPI_Finalize();
    return 0;
}

static int 
check_snapshots(FastPM * fastpm, Parameters * prr);

static int 
measure_powerspectrum(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr);

int run_fastpm(FastPM * fastpm, Parameters * prr, MPI_Comm comm) {
    CLOCK(init);
    CLOCK(ic);
    CLOCK(evolve);

    const double rho_crit = 27.7455;
    const double M0 = prr->omega_m*rho_crit*pow(prr->boxsize / prr->nc, 3.0);
    fastpm_info("mass of a particle is %g 1e10 Msun/h\n", M0); 

    MPI_Barrier(comm);
    ENTER(init);

    fastpm_init(fastpm, 
        prr->NprocY, prr->UseFFTW, 
        comm);

    LEAVE(init);

    fastpm_add_extension(fastpm,
        FASTPM_EXT_AFTER_FORCE,
        measure_powerspectrum,
        prr);

    fastpm_add_extension(fastpm,
        FASTPM_EXT_AFTER_KICK,
        check_snapshots,
        prr);

    fastpm_add_extension(fastpm,
        FASTPM_EXT_AFTER_DRIFT,
        check_snapshots,
        prr);

    MPI_Barrier(comm);
    ENTER(ic);

    if(prr->readic_filename) {
        read_runpb_ic(fastpm, fastpm->p, prr->readic_filename);
        fastpm_setup_ic(fastpm, NULL, prr->time_step[0]);
    } else {
        FastPMFloat * delta_k = pm_alloc(fastpm->pm_2lpt);

        fastpm_info("Powerspecectrum file: %s\n", prr->power_spectrum_filename);

        power_init(prr->power_spectrum_filename, 
                prr->time_step[0], 
                prr->sigma8, 
                prr->omega_m, 
                1 - prr->omega_m, comm);

        if(prr->readnoise_filename) {
            fastpm_info("Reading grafic white noise file from '%s'.\n", prr->readnoise_filename);
            fastpm_info("GrafIC noise is Fortran ordering. FastPM is in C ordering.\n");
            fastpm_info("The simulation will be transformed x->z y->y z->x.\n");

            FastPMFloat * g_x = pm_alloc(fastpm->pm_2lpt);

            read_grafic_gaussian(fastpm->pm_2lpt, g_x, prr->readnoise_filename);

            fastpm_utils_induce_correlation(fastpm->pm_2lpt, g_x, delta_k, PowerSpecWithData, NULL);
            pm_free(fastpm->pm_2lpt, g_x);
        } else {

            fastpm_utils_fill_deltak(fastpm->pm_2lpt, delta_k, prr->random_seed, PowerSpecWithData, NULL, FASTPM_DELTAK_GADGET);
        }
        
        if(prr->writenoisek_filename) {
            fastpm_utils_dump(fastpm->pm_2lpt, prr->writenoisek_filename, delta_k);
        }

        fastpm_setup_ic(fastpm, delta_k, prr->time_step[0]);

        pm_free(fastpm->pm_2lpt, delta_k);

    }
    LEAVE(ic);

    MPI_Barrier(comm);
    ENTER(evolve);

    fastpm_evolve(fastpm, prr->time_step, prr->n_time_step);

    LEAVE(evolve);

    fastpm_destroy(fastpm);

    fastpm_clock_stat(comm);
    return 0;
}

static int check_snapshots(FastPM * fastpm, Parameters * prr) {
    char TEMPLATE[1024];
    sprintf(TEMPLATE, "%s%05d_%%0.04f.bin", prr->snapshot_filename, prr->random_seed);
    fastpm_interp(fastpm, prr->aout, prr->n_aout, take_a_snapshot, TEMPLATE);
    return 0;
}

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, double aout, void * template) 
{
    CLOCK(io);
    CLOCK(meta);

    char filebase[1024];
    double z_out= 1.0/aout - 1.0;

    sprintf(filebase, template, aout);

    ENTER(meta);
    ensure_dir(filebase);
    LEAVE(meta);

    MPI_Barrier(fastpm->comm);
    ENTER(io);

    write_runpb_snapshot(fastpm, snapshot, filebase);

    LEAVE(io);

    fastpm_info("snapshot %s written z = %6.4f a = %6.4f\n", 
            filebase, z_out, aout);
    return 0;
}

static int 
measure_powerspectrum(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr) 
{
    CLOCK(compute);
    CLOCK(io);

    FastPMPowerSpectrum ps;
    /* calculate the power spectrum */
    fastpm_power_spectrum_init(&ps, pm_nmesh(fastpm->pm)[0] / 2);

    MPI_Barrier(fastpm->comm);
    ENTER(compute);

    fastpm_utils_calculate_powerspectrum(fastpm->pm, delta_k, &ps, pow(fastpm->nc, 3.0));

    LEAVE(compute);

    MPI_Barrier(fastpm->comm);

    ENTER(io);
    if(prr->measure_power_spectrum_filename) {
        if(fastpm->ThisTask == 0) {
            ensure_dir(prr->measure_power_spectrum_filename);
            char buf[1024];
            sprintf(buf, "%s%05d_%0.04f.txt", prr->measure_power_spectrum_filename, prr->random_seed, a_x);
            fastpm_power_spectrum_write(&ps, fastpm->pm, buf);
        }
    }
    LEAVE(io);

    fastpm_power_spectrum_destroy(&ps);

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
