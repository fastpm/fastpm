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
parse_args(int * argc, char *** argv, Parameters * prr);

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, double aout, Parameters * prr);

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
write_snapshot(FastPM * fastpm, PMStore * p, char * filebase, char * parameters, int Nwriters);

int 
read_snapshot(FastPM * fastpm, PMStore * p, char * filebase);

int 
read_parameters(char * filename, Parameters * param, int argc, char ** argv, MPI_Comm comm);

int run_fastpm(FastPM * fastpm, Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    MPI_Init(&argc, &argv);

    Parameters prr = {0};

    parse_args(&argc, &argv, &prr);

    MPI_Comm comm = MPI_COMM_WORLD; 

    libfastpm_init();

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    read_parameters(ParamFileName, &prr, argc, argv, comm);

    /* convert parameter files pm_nc_factor into VPMInit */
    VPMInit * vpminit = alloca(sizeof(VPMInit) * (prr.n_pm_nc_factor + 1));
    int i;
    for(i = 0; i < prr.n_pm_nc_factor; i ++) {
        vpminit[i].a_start = prr.change_pm[i];
        vpminit[i].pm_nc_factor = prr.pm_nc_factor[i];
    }
    /* mark the end */
    vpminit[i].pm_nc_factor = 0;

    fastpm_info("np_alloc_factor = %g\n", prr.np_alloc_factor);
    FastPM * fastpm = & (FastPM) {
        .nc = prr.nc,
        .alloc_factor = prr.np_alloc_factor, 
        .vpminit = vpminit,
        .boxsize = prr.boxsize,
        .omega_m = prr.omega_m,
        .USE_COLA = prr.force_mode == FORCE_MODE_COLA,
        .USE_NONSTDDA = !prr.cola_stdda,
        .USE_DX1_ONLY = prr.use_dx1_only,
        .nLPT = -2.5f,
        .K_LINEAR = prr.enforce_broadband_kmax,
    };

    if(prr.enforce_broadband_mode == MODEL_LINEAR) {
        fastpm->USE_MODEL = FASTPM_MODEL_LINEAR;
    } else
    if(prr.enforce_broadband_mode == MODEL_2LPT) {
        fastpm->USE_MODEL = FASTPM_MODEL_2LPT;
    } else 
    if(prr.enforce_broadband_mode == MODEL_ZA) {
        fastpm->USE_MODEL = FASTPM_MODEL_ZA;
    } else 
    if(prr.enforce_broadband_mode == MODEL_NONE) {
        fastpm->USE_MODEL = FASTPM_MODEL_NONE;
    } else {
        fastpm_raise(-1, "wrong model type!\n");
    }

    run_fastpm(fastpm, &prr, comm);

    libfastpm_cleanup();

    MPI_Finalize();
    return 0;
}

static int 
check_snapshots(FastPM * fastpm, void * unused, Parameters * prr);

static int 
measure_powerspectrum(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr);

static void 
prepare_ic(FastPM * fastpm, Parameters * prr, MPI_Comm comm);

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
        FASTPM_EXT_BEFORE_KICK,
        check_snapshots,
        prr);

    fastpm_add_extension(fastpm,
        FASTPM_EXT_BEFORE_DRIFT,
        check_snapshots,
        prr);

    MPI_Barrier(comm);
    ENTER(ic);
    prepare_ic(fastpm, prr, comm);
    LEAVE(ic);

    MPI_Barrier(comm);
    ENTER(evolve);

    fastpm_evolve(fastpm, prr->time_step, prr->n_time_step);

    LEAVE(evolve);

    fastpm_destroy(fastpm);

    fastpm_clock_stat(comm);
    return 0;
}

static void 
prepare_ic(FastPM * fastpm, Parameters * prr, MPI_Comm comm) 
{
    /* we may need a read gadget ic here too */
    if(prr->read_runpbic) {
        read_runpb_ic(fastpm, fastpm->p, prr->read_runpbic);
        fastpm_setup_ic(fastpm, NULL);
        return;
    } 

    /* at this point generating the ic involves delta_k */
    FastPMFloat * delta_k = pm_alloc(fastpm->pm_2lpt);

    if(prr->read_noisek) {
        fastpm_info("Reading Fourier space noise from %s\n", prr->read_noisek);
        fastpm_utils_load(fastpm->pm_2lpt, prr->read_noisek, delta_k);
        goto finish;
    } 

    if(prr->read_noise) {
        fastpm_info("Reading Real space noise from %s\n", prr->read_noise);

        FastPMFloat * g_x = pm_alloc(fastpm->pm_2lpt);
        fastpm_utils_load(fastpm->pm_2lpt, prr->read_noise, g_x);
        pm_r2c(fastpm->pm_2lpt, g_x, delta_k);
        pm_free(fastpm->pm_2lpt, g_x);
        goto finish;
    }

    /* at this power we need a powerspectrum file to convolve the guassian */
    if(!prr->read_powerspectrum) {
        fastpm_raise(-1, "Need a power spectrum to start the simulation.\n");
    }

    fastpm_info("Powerspecectrum file: %s\n", prr->read_powerspectrum);

    power_init(prr->read_powerspectrum, 
            prr->time_step[0], 
            prr->sigma8, 
            prr->omega_m, 
            1 - prr->omega_m, comm);

    if(prr->read_grafic) {
        fastpm_info("Reading grafic white noise file from '%s'.\n", prr->read_grafic);
        fastpm_info("GrafIC noise is Fortran ordering. FastPM is in C ordering.\n");
        fastpm_info("The simulation will be transformed x->z y->y z->x.\n");

        FastPMFloat * g_x = pm_alloc(fastpm->pm_2lpt);

        read_grafic_gaussian(fastpm->pm_2lpt, g_x, prr->read_grafic);

        fastpm_utils_induce_correlation(fastpm->pm_2lpt, g_x, delta_k, PowerSpecWithData, NULL);
        pm_free(fastpm->pm_2lpt, g_x);
        goto finish;
    } 

    /* Nothing to read from, just generate a gadget IC with the seed. */

    fastpm_utils_fill_deltak(fastpm->pm_2lpt, delta_k, prr->random_seed, PowerSpecWithData, NULL, FASTPM_DELTAK_GADGET);

    /* our write out and clean up stuff.*/
finish:
    if(prr->inverted_ic) {
        ptrdiff_t i;
        for(i = 0; i < pm_size(fastpm->pm_2lpt); i ++) {
            delta_k[i] *= -1;
        }
    }
    if(prr->remove_cosmic_variance) {
        fastpm_utils_remove_cosmic_variance(fastpm->pm_2lpt, delta_k, PowerSpecWithData, NULL);
    }

    if(prr->write_noisek) {
        fastpm_info("Writing fourier space noise to %s\n", prr->write_noisek);
        ensure_dir(prr->write_noisek);
        fastpm_utils_dump(fastpm->pm_2lpt, prr->write_noisek, delta_k);
    }

    if(prr->write_noise) {
        FastPMFloat * g_x = pm_alloc(fastpm->pm_2lpt);
        pm_assign(fastpm->pm_2lpt, delta_k, g_x);
        pm_c2r(fastpm->pm_2lpt, g_x);

        fastpm_info("Writing real space noise to %s\n", prr->write_noise);
        ensure_dir(prr->write_noise);
        fastpm_utils_dump(fastpm->pm_2lpt, prr->write_noise, g_x);
        pm_free(fastpm->pm_2lpt, g_x);
    }
    fastpm_setup_ic(fastpm, delta_k);

    pm_free(fastpm->pm_2lpt, delta_k);
}

static int check_snapshots(FastPM * fastpm, void * unused, Parameters * prr) {
    fastpm_interp(fastpm, prr->aout, prr->n_aout, (fastpm_interp_action)take_a_snapshot, prr);
    return 0;
}

static int 
take_a_snapshot(FastPM * fastpm, PMStore * snapshot, double aout, Parameters * prr) 
{
    CLOCK(io);
    CLOCK(meta);

    if(prr->write_snapshot) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;
        int Nwriters = prr->Nwriters;
        if(Nwriters == 0) {
            MPI_Comm_size(fastpm->comm, &Nwriters);
        }
        sprintf(filebase, "%s_%0.04f", prr->write_snapshot, aout);

        fastpm_info("Writing snapshot %s at z = %6.4f a = %6.4f with %d writers\n", 
                filebase, z_out, aout, Nwriters);

        ENTER(meta);
        ensure_dir(filebase);
        LEAVE(meta);

        MPI_Barrier(fastpm->comm);
        ENTER(io);
        write_snapshot(fastpm, snapshot, filebase, prr->string, Nwriters);
        LEAVE(io);

        fastpm_info("snapshot %s written\n", filebase);
    }
    if(prr->write_runpb_snapshot) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;

        sprintf(filebase, "%s_%0.04f.bin", prr->write_runpb_snapshot, aout);
        ENTER(meta);
        ensure_dir(filebase);
        LEAVE(meta);

        MPI_Barrier(fastpm->comm);
        ENTER(io);
        write_runpb_snapshot(fastpm, snapshot, filebase);

        LEAVE(io);

        fastpm_info("snapshot %s written z = %6.4f a = %6.4f\n", 
                filebase, z_out, aout);

    }
    return 0;
}

static int 
measure_powerspectrum(FastPM * fastpm, FastPMFloat * delta_k, double a_x, Parameters * prr) 
{
    CLOCK(compute);
    CLOCK(io);

    fastpm_report_memory(fastpm->comm);

    FastPMPowerSpectrum ps;
    /* calculate the power spectrum */
    fastpm_power_spectrum_init(&ps, pm_nmesh(fastpm->pm)[0] / 2);

    MPI_Barrier(fastpm->comm);
    ENTER(compute);

    fastpm_utils_calculate_powerspectrum(fastpm->pm, delta_k, &ps, pow(fastpm->nc, 3.0));

    LEAVE(compute);

    MPI_Barrier(fastpm->comm);

    ENTER(io);
    if(prr->write_powerspectrum) {
        if(fastpm->ThisTask == 0) {
            ensure_dir(prr->write_powerspectrum);
            char buf[1024];
            sprintf(buf, "%s_%0.04f.txt", prr->write_powerspectrum, a_x);
            fastpm_power_spectrum_write(&ps, fastpm->pm, buf);
        }
    }
    LEAVE(io);

    fastpm_power_spectrum_destroy(&ps);

    return 0;
}

static void 
parse_args(int * argc, char *** argv, Parameters * prr) 
{
    char opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    ParamFileName = NULL;
    prr->NprocY = 0;    
    prr->Nwriters = 0;
    while ((opt = getopt(*argc, *argv, "h?y:fW:")) != -1) {
        switch(opt) {
            case 'y':
                prr->NprocY = atoi(optarg);
            break;
            case 'f':
                prr->UseFFTW = 1;
            break;
            case 'W':
                prr->Nwriters = atoi(optarg);
            break;
            case 'h':
            case '?':
            default:
                goto usage;
            break;
        }
    }
    if(optind >= *argc) {
        goto usage;
    }

    ParamFileName = (*argv)[optind];
    *argv += optind;
    *argc -= optind; 
    return;

usage:
    printf("Usage: fastpm [-W Nwriters] [-f] [-y NprocY] paramfile\n"
    "-f Use FFTW \n"
    "-y Set the number of processes in the 2D mesh\n"
    "-n Throttle IO (bigfile only) \n"
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

