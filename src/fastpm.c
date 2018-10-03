#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <signal.h>
#include <getopt.h>
#include <limits.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>
#include <fastpm/lightcone.h>
#include <fastpm/constrainedgaussian.h>
#include <fastpm/io.h>
#include <fastpm/fof.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lua-config.h"

/* c99 has no pi. */
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

typedef struct {
    int UseFFTW;
    int NprocY;
    int Nwriters;
    size_t MemoryPerRank;
    LuaConfig * config;
    char * string;
} Parameters;

extern char * 
lua_config_parse(char * entrypoint, char * filename, int argc, char ** argv, char ** error);

extern void
init_stacktrace();

#define CONF(prr, name) lua_config_get_ ## name (prr->config)
#define HAS(prr, name) lua_config_has_ ## name (prr->config)

/* command-line arguments */
static char * ParamFileName;

static void
_memory_peak_handler(FastPMMemory * mem, void * userdata);

static void 
parse_args(int * argc, char *** argv, Parameters * prr);

static void
free_parameters(Parameters * prr);

static int 
take_a_snapshot(FastPMSolver * fastpm, FastPMStore * snapshot, double aout, Parameters * prr);

static void
smesh_force_handler(FastPMSolver * solver, FastPMForceEvent * event, FastPMSMesh * smesh);

struct smesh_ready_handler_data {
    FastPMSolver * fastpm;
    Parameters * prr;
    int64_t * hist;
    double * aedges;
    int Nedges;
    int Nslices;
    FastPMSMeshSlice * slices;
};

static void
smesh_ready_handler(FastPMSMesh * mesh, FastPMLCEvent * lcevent, struct smesh_ready_handler_data * userdata);

struct usmesh_ready_handler_data {
    FastPMSolver * fastpm;
    Parameters * prr;
    FastPMStore tail[1];
    int64_t * hist;
    double * aedges;
    int Nedges;
};

static void
usmesh_ready_handler(FastPMUSMesh * mesh, FastPMLCEvent * lcevent, struct usmesh_ready_handler_data * userdata);

int 
read_runpb_ic(FastPMSolver * fastpm, FastPMStore * p, const char * filename);

void 
read_grafic_gaussian(PM * pm, FastPMFloat * g_x, const char * filename);

int
write_runpb_snapshot(FastPMSolver * fastpm, FastPMStore * p, const char * filebase);

int
read_complex(PM * pm, FastPMFloat * data, const char * filename, const char * blockname, int Nwriters);

int
read_parameters(char * filename, Parameters * param, int argc, char ** argv, MPI_Comm comm);

int
read_powerspectrum(FastPMPowerSpectrum *ps, const char filename[], const double sigma8, MPI_Comm comm);

static void
write_fof(FastPMSolver * fastpm, FastPMStore * snapshot, char * filebase, Parameters * prr);

static void
write_usmesh_fof(FastPMSolver * fastpm,
        FastPMStore * snapshot,
        char * filebase, Parameters * prr,
        FastPMLCEvent * lcevent,
        FastPMStore * tail,
        FastPMLightCone * lc);

int run_fastpm(FastPMConfig * config, Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    init_stacktrace();

    MPI_Init(&argc, &argv);

    Parameters * prr = alloca(sizeof(prr[0]));

    parse_args(&argc, &argv, prr);

    MPI_Comm comm = MPI_COMM_WORLD; 

    int Nwriters = prr->Nwriters;
    if(Nwriters == 0) {
        MPI_Comm_size(comm, &Nwriters);
        prr->Nwriters = Nwriters;
    }
    libfastpm_init();

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    fastpm_info("This is FastPM, with libfastpm version %s.\n", LIBFASTPM_VERSION);

    libfastpm_set_memory_bound(prr->MemoryPerRank * 1024 * 1024);
    fastpm_memory_set_handlers(_libfastpm_get_gmem(), NULL, _memory_peak_handler, &comm);

    read_parameters(ParamFileName, prr, argc, argv, comm);

    /* convert parameter files pm_nc_factor into VPMInit */
    VPMInit * vpminit = NULL;
    if(CONF(prr, ndim_pm_nc_factor) == 0) {
        vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = CONF(prr, pm_nc_factor)[0] },
            {.a_start = 1, .pm_nc_factor = 0},
            };
    } else
    if(CONF(prr, ndim_pm_nc_factor) == 2) {
        vpminit = alloca(sizeof(VPMInit) * (CONF(prr, shape_pm_nc_factor)[0] + 1));
        int i;
        for(i = 0; i < CONF(prr, n_pm_nc_factor); i ++) {
            vpminit[i].a_start = CONF(prr, pm_nc_factor)[2 * i];
            vpminit[i].pm_nc_factor = CONF(prr, pm_nc_factor)[2 * i + 1];
        }
        /* mark the end */
        vpminit[i].pm_nc_factor = 0;
    } else {
        fastpm_raise(-1, "Unknown format of pm_nc_factor, either a scalar or a 2d array. ");
    }

    fastpm_info("np_alloc_factor = %g\n", CONF(prr, np_alloc_factor));

    FastPMConfig * config = & (FastPMConfig) {
        .nc = CONF(prr, nc),
        .alloc_factor = CONF(prr, np_alloc_factor),
        .vpminit = vpminit,
        .boxsize = CONF(prr, boxsize),
        .omega_m = CONF(prr, omega_m),
        .hubble_param = CONF(prr, h),
        .USE_DX1_ONLY = CONF(prr, za),
        .nLPT = -2.5f,
        .USE_SHIFT = CONF(prr, shift),
        .FORCE_TYPE = CONF(prr, force_mode),
        .KERNEL_TYPE = CONF(prr, kernel_type),
        .DEALIASING_TYPE = CONF(prr, dealiasing_type),
        .PAINTER_TYPE = CONF(prr, painter_type),
        .painter_support = CONF(prr, painter_support),
        .NprocY = prr->NprocY,
        .UseFFTW = prr->UseFFTW,
        .ExtraAttributes = 0,
    };

    if(CONF(prr, compute_potential)) {
        config->ExtraAttributes |= PACK_POTENTIAL;
    }

    run_fastpm(config, prr, comm);

    libfastpm_cleanup();

    free_parameters(prr);

    MPI_Finalize();

    return 0;
}

static int 
check_snapshots(FastPMSolver * fastpm, FastPMInterpolationEvent * event, Parameters * prr);

static int 
check_lightcone(FastPMSolver * fastpm, FastPMInterpolationEvent * event, FastPMUSMesh * lc);

static int 
write_powerspectrum(FastPMSolver * fastpm, FastPMForceEvent * event, Parameters * prr);

static void 
prepare_ic(FastPMSolver * fastpm, Parameters * prr, MPI_Comm comm);

static void
prepare_lc(FastPMSolver * fastpm, Parameters * prr,
        FastPMLightCone * lc, FastPMUSMesh ** usmesh,
        FastPMSMesh ** smesh);

static int 
print_transition(FastPMSolver * fastpm, FastPMTransitionEvent * event, Parameters * prr);

int run_fastpm(FastPMConfig * config, Parameters * prr, MPI_Comm comm) {
    FastPMSolver fastpm[1];

    CLOCK(init);
    CLOCK(ic);
    CLOCK(evolve);
    CLOCK(io);
    CLOCK(sort);
    CLOCK(indexing);

    const double rho_crit = 27.7455;
    const double M0 = CONF(prr, omega_m) * rho_crit
                    * pow(CONF(prr, boxsize) / CONF(prr, nc), 3.0);
    fastpm_info("mass of a particle is %g 1e10 Msun/h\n", M0); 

    MPI_Barrier(comm);
    ENTER(init);

    fastpm_solver_init(fastpm, config, comm);

    fastpm_info("BaseProcMesh : %d x %d\n",
            pm_nproc(fastpm->basepm)[0], pm_nproc(fastpm->basepm)[1]);
#ifdef _OPENMP
    fastpm_info("%d Threads\n", omp_get_max_threads());
#endif

    LEAVE(init);

    fastpm_add_event_handler(&fastpm->event_handlers,
        FASTPM_EVENT_FORCE,
        FASTPM_EVENT_STAGE_AFTER,
        (FastPMEventHandlerFunction) write_powerspectrum,
        prr);

    fastpm_add_event_handler(&fastpm->event_handlers,
        FASTPM_EVENT_INTERPOLATION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) check_snapshots,
        prr);

    fastpm_add_event_handler(&fastpm->event_handlers,
        FASTPM_EVENT_TRANSITION,
        FASTPM_EVENT_STAGE_BEFORE,
        (FastPMEventHandlerFunction) print_transition,
        prr);

    /* initialize the lightcone */
    FastPMLightCone lc[1] = {{
        .speedfactor = CONF(prr, dh_factor),
        .cosmology = fastpm->cosmology,
        .fov = CONF(prr, lc_fov),
        .octants = {0, 0, 0, 0, 0, 0, 0, 0},
        .tol = 2. / CONF(prr, nc) / CONF(prr, lc_smesh_fraction),
    }};

    FastPMUSMesh * usmesh = NULL;
    FastPMSMesh * smesh = NULL;

    prepare_lc(fastpm, prr, lc, &usmesh, &smesh);

    MPI_Barrier(comm);

    ENTER(ic);
    prepare_ic(fastpm, prr, comm);

    fastpm_store_fill_subsample_mask(fastpm->p, CONF(prr, particle_fraction), fastpm->p->mask, comm);

    fastpm_info("dx1  : %g %g %g %g\n", 
            fastpm->info.dx1[0], fastpm->info.dx1[1], fastpm->info.dx1[2],
            (fastpm->info.dx1[0] + fastpm->info.dx1[1] + fastpm->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n", 
            fastpm->info.dx2[0], fastpm->info.dx2[1], fastpm->info.dx2[2],
            (fastpm->info.dx2[0] + fastpm->info.dx2[1] + fastpm->info.dx2[2]) / 3.0);

    LEAVE(ic);

    MPI_Barrier(comm);
    ENTER(evolve);
    fastpm_solver_evolve(fastpm, CONF(prr, time_step), CONF(prr, n_time_step));
    LEAVE(evolve);

    if(smesh)
        fastpm_smesh_destroy(smesh);

    if(usmesh)
        fastpm_usmesh_destroy(usmesh);

    free(smesh);
    free(usmesh);

    fastpm_lc_destroy(lc);

    fastpm_solver_destroy(fastpm);

    {
        FastPMMemory * g = _libfastpm_get_gmem();
        double max_used_bytes;
        double min_used_bytes;

        MPIU_stats(comm, g->peak_bytes, "<>",
                &min_used_bytes,
                &max_used_bytes);

        fastpm_log(INFO, "Peak memory usage max: %g MB min : %g MB\n",
                max_used_bytes / 1024. / 1024,
                min_used_bytes / 1024. / 1024
            );
    }

    fastpm_clock_stat(comm);

    return 0;
}

static void 
prepare_ic(FastPMSolver * fastpm, Parameters * prr, MPI_Comm comm) 
{
    /* we may need a read gadget ic here too */
    if(CONF(prr, read_runpbic)) {
        read_runpb_ic(fastpm, fastpm->p, CONF(prr, read_runpbic));
        fastpm_solver_setup_ic(fastpm, NULL, CONF(prr, time_step)[0]);
        return;
    } 

    /* at this point generating the ic involves delta_k */
    FastPMFloat * delta_k = pm_alloc(fastpm->basepm);

    if(CONF(prr, read_lineark)) {
        fastpm_info("Reading Fourier space linear overdensity from %s\n", CONF(prr, read_lineark));
        read_complex(fastpm->basepm, delta_k, CONF(prr, read_lineark), "LinearDensityK", prr->Nwriters);

        if(CONF(prr, inverted_ic)) {
            fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, -1);
        }
        goto produce;
    }

    /* at this power we need a powerspectrum file to convolve the guassian */
    if(!CONF(prr, read_powerspectrum)) {
        fastpm_raise(-1, "Need a power spectrum to start the simulation.\n");
    }

    FastPMPowerSpectrum linear_powerspectrum;

    read_powerspectrum(&linear_powerspectrum, CONF(prr, read_powerspectrum), CONF(prr, sigma8), comm);

    if(CONF(prr, read_grafic)) {
        fastpm_info("Reading grafic white noise file from '%s'.\n", CONF(prr, read_grafic));
        fastpm_info("GrafIC noise is Fortran ordering. FastPMSolver is in C ordering.\n");
        fastpm_info("The simulation will be transformed x->z y->y z->x.\n");

        FastPMFloat * g_x = pm_alloc(fastpm->basepm);

        read_grafic_gaussian(fastpm->basepm, g_x, CONF(prr, read_grafic));

        /* r2c will reduce the variance. Compensate here.*/
        fastpm_apply_multiply_transfer(fastpm->basepm, g_x, g_x, sqrt(pm_norm(fastpm->basepm)));
        pm_r2c(fastpm->basepm, g_x, delta_k);

        pm_free(fastpm->basepm, g_x);

        goto induce;
    }

    if(CONF(prr, read_whitenoisek)) {
        fastpm_info("Reading Fourier white noise file from '%s'.\n", CONF(prr, read_whitenoisek));

        read_complex(fastpm->basepm, delta_k, CONF(prr, read_whitenoisek), "WhiteNoiseK", prr->Nwriters);
        goto induce;
    }

    /* Nothing to read from, just generate a gadget IC with the seed. */
    fastpm_ic_fill_gaussiank(fastpm->basepm, delta_k, CONF(prr, random_seed), FASTPM_DELTAK_GADGET);

induce:
    if(CONF(prr, remove_cosmic_variance)) {
        fastpm_info("Remove Cosmic variance from initial condition.\n");
        fastpm_ic_remove_variance(fastpm->basepm, delta_k);
    }

    if(CONF(prr, set_mode)) {
        int method = 0;
        /* FIXME: use enums */
        if(0 == strcmp(CONF(prr, set_mode_method), "add")) {
            method = 1;
            fastpm_info("SetMode is add\n");
        } else {
            fastpm_info("SetMode is override\n");
        }
        int i;
        double * c = CONF(prr, set_mode);
        for(i = 0; i < CONF(prr, n_set_mode); i ++) {
            ptrdiff_t mode[4] = {
                c[i * 5 + 0],
                c[i * 5 + 1],
                c[i * 5 + 2],
                c[i * 5 + 3],
            };
            double value = c[i * 5 + 4];
            fastpm_apply_set_mode_transfer(fastpm->basepm, delta_k, delta_k, mode, value, method);
            double result = fastpm_apply_get_mode_transfer(fastpm->basepm, delta_k, mode);
            fastpm_info("SetMode %d : %td %td %td %td value = %g, to = %g\n", i, mode[0], mode[1], mode[2], mode[3], value, result);
        }
    }

    if(CONF(prr, inverted_ic)) {
        fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, -1);
    }

    double variance = pm_compute_variance(fastpm->basepm, delta_k);
    fastpm_info("Variance of input white noise is %0.8f, expectation is %0.8f\n", variance, 1.0 - 1.0 / pm_norm(fastpm->basepm));

    if(CONF(prr, write_whitenoisek)) {
        fastpm_info("Writing Fourier white noise to file '%s'.\n", CONF(prr, write_whitenoisek));
        write_complex(fastpm->basepm, delta_k, CONF(prr, write_whitenoisek), "WhiteNoiseK", prr->Nwriters);
    }

    /* introduce correlation */
    if(CONF(prr, f_nl_type) == FASTPM_FNL_NONE) {
        fastpm_info("Inducing correlation to the white noise.\n");

        fastpm_ic_induce_correlation(fastpm->basepm, delta_k,
            (fastpm_fkfunc) fastpm_powerspectrum_eval2, &linear_powerspectrum);
    } else {
        double kmax_primordial;
        kmax_primordial = CONF(prr, nc) / 2.0 * 2.0*M_PI/CONF(prr, boxsize) * CONF(prr, kmax_primordial_over_knyquist);
        fastpm_info("Will set Phi_Gaussian(k)=0 for k>=%f.\n", kmax_primordial);
        FastPMPNGaussian png = {
            .fNL = CONF(prr, f_nl),
            .type = CONF(prr, f_nl_type),
            .kmax_primordial = kmax_primordial,
            .pkfunc = (fastpm_fkfunc) fastpm_powerspectrum_eval2,
            .pkdata = &linear_powerspectrum,
            .h = CONF(prr, h),
            .scalar_amp = CONF(prr, scalar_amp),
            .scalar_spectral_index = CONF(prr, scalar_spectral_index),
            .scalar_pivot = CONF(prr, scalar_pivot)
        };
        fastpm_info("Inducing non gaussian correlation to the white noise.\n");
        fastpm_png_induce_correlation(&png, fastpm->basepm, delta_k);
    }

    /* The linear density field is not redshift zero, then evolve it with the model cosmology to 
     * redshift zero.
     * This matches the linear power at the given redshift, not necessarily redshift 0. */
    {
        double linear_density_redshift = CONF(prr, linear_density_redshift);
        double linear_evolve = fastpm_solver_growth_factor(fastpm, 1.0) /
                               fastpm_solver_growth_factor(fastpm, 1 / (linear_density_redshift + 1));

        fastpm_info("Reference linear density is calibrated at redshift %g; multiply by %g to extract to redshift 0.\n", linear_density_redshift, linear_evolve);

        fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, linear_evolve);
    }

    /* set the mean to 1.0 */
    ptrdiff_t mode[4] = { 0, 0, 0, 0, };

    fastpm_apply_modify_mode_transfer(fastpm->basepm, delta_k, delta_k, mode, 1.0);

    /* add constraints */
    if(CONF(prr, constraints)) {
        FastPM2PCF xi;

        fastpm_2pcf_from_powerspectrum(&xi, (fastpm_fkfunc) fastpm_powerspectrum_eval2, &linear_powerspectrum, CONF(prr, boxsize), CONF(prr, nc));

        FastPMConstrainedGaussian cg = {
            .constraints = malloc(sizeof(FastPMConstraint) * (CONF(prr, n_constraints) + 1)),
        };
        fastpm_info("Applying %d constraints.\n", CONF(prr, n_constraints));
        int i;
        for(i = 0; i < CONF(prr, n_constraints); i ++) {
            double * c = CONF(prr, constraints);
            cg.constraints[i].x[0] = c[4 * i + 0];
            cg.constraints[i].x[1] = c[4 * i + 1];
            cg.constraints[i].x[2] = c[4 * i + 2];
            cg.constraints[i].c = c[4 * i + 3];
            fastpm_info("Constraint %d : %g %g %g peak-sigma = %g\n", i, c[4 * i + 0], c[4 * i + 1], c[4 * i + 2], c[4 * i + 3]);
        }
        cg.constraints[i].x[0] = -1;
        cg.constraints[i].x[1] = -1;
        cg.constraints[i].x[2] = -1;
        cg.constraints[i].c = -1;

        if(CONF(prr, write_lineark)) {
            fastpm_info("Writing fourier space linear field before constraints to %s\n", CONF(prr, write_lineark));
            write_complex(fastpm->basepm, delta_k, CONF(prr, write_lineark), "UnconstrainedLinearDensityK", prr->Nwriters);
        }
        fastpm_cg_apply_constraints(&cg, fastpm->basepm, &xi, delta_k);

        free(cg.constraints);
    }

    fastpm_powerspectrum_destroy(&linear_powerspectrum);

    /* our write out and clean up stuff.*/
produce:

    if(CONF(prr, write_lineark)) {
        fastpm_info("Writing fourier space linear field to %s\n", CONF(prr, write_lineark));
        write_complex(fastpm->basepm, delta_k, CONF(prr, write_lineark), "LinearDensityK", prr->Nwriters);
    }

    if(CONF(prr, write_powerspectrum)) {
        FastPMPowerSpectrum ps;
        /* calculate the power spectrum */
        fastpm_powerspectrum_init_from_delta(&ps, fastpm->basepm, delta_k, delta_k);

        char buf[1024];
        sprintf(buf, "%s_linear.txt", CONF(prr, write_powerspectrum));
        fastpm_info("writing linear power spectrum to %s\n", buf);
        if(fastpm->ThisTask == 0) {
            fastpm_path_ensure_dirname(CONF(prr, write_powerspectrum));
            fastpm_powerspectrum_write(&ps, buf, pow(fastpm->config->nc, 3.0));
        }
        fastpm_powerspectrum_destroy(&ps);
    }

    fastpm_solver_setup_ic(fastpm, delta_k, CONF(prr, time_step)[0]);

    pm_free(fastpm->basepm, delta_k);
}

static void
_usmesh_ready_handler_free(void * userdata) {
    struct usmesh_ready_handler_data * data = userdata;
    fastpm_store_destroy(data->tail);
    free(data->aedges);
    free(data->hist);
    free(data);
}


static void
_smesh_ready_handler_free(void * userdata) {
    struct smesh_ready_handler_data * data = userdata;
    free(data->aedges);
    free(data->hist);
    free(data);
}

double * _get_aedges(FastPMSMeshSlice * slices, size_t Nslices)
{
    /* use the centers of the layers to ensure each bin has exactly one layer */
    double * aedges = malloc(sizeof(double) * (Nslices - 1));
    int i;
    for (i = 0; i < Nslices - 1; i ++) {
        aedges[i] = (slices[i].aemit + slices[i+1].aemit) * 0.5;
    }
    return aedges;
}

static void
prepare_lc(FastPMSolver * fastpm, Parameters * prr,
        FastPMLightCone * lc, FastPMUSMesh ** usmesh,
        FastPMSMesh ** smesh)
{
    {
        if(CONF(prr, ndim_lc_glmatrix) != 2 ||
           CONF(prr, shape_lc_glmatrix)[0] != 4 ||
           CONF(prr, shape_lc_glmatrix)[1] != 4
        ) {
            fastpm_raise(-1, "GL Matrix must be a 4x4 matrix\n");
        }

        int i, j;
        double * c = CONF(prr, lc_glmatrix);

        for (i = 0; i < 4; i ++) {
            for (j = 0; j < 4; j ++) {
                lc->glmatrix[i][j] = *c;
                c ++;
            }
            fastpm_info("GLTransformation [%d] : %g %g %g %g\n", i,
                lc->glmatrix[i][0], lc->glmatrix[i][1], lc->glmatrix[i][2], lc->glmatrix[i][3]);
        }
    }

    {
        int i;
        for(i = 0; i < CONF(prr, n_lc_octants); i ++) {
            int oct = CONF(prr, lc_octants)[i];
            lc->octants[oct % 8] = 1;
            fastpm_info("Using Octant %d\n", oct);
        }
    }

    fastpm_lc_init(lc);


    double lc_amin = HAS(prr, lc_amin)?CONF(prr, lc_amin):CONF(prr, time_step)[0];

    double lc_amax = HAS(prr, lc_amax)?CONF(prr, lc_amax):CONF(prr, time_step)[CONF(prr, n_time_step) - 1];

    fastpm_info("Unstructured Lightcone amin= %g amax=%g\n", lc_amin, lc_amax);

    *smesh = NULL;
    if(CONF(prr, lc_write_smesh)) {
        *smesh = malloc(sizeof(FastPMSMesh));

        double n = CONF(prr, nc) / CONF(prr, boxsize) * CONF(prr, lc_smesh_fraction);

        fastpm_smesh_init(*smesh, lc, fastpm->p->np_upper, 1 / n);

        if(lc->fov > 0) {
            fastpm_info("Creating healpix structured meshes for FOV=%g, with number density %g per (Mpc/h)**3. \n",
                lc->fov, n * n * n);
            fastpm_smesh_add_layers_healpix(*smesh,
                    n * n, n * n * n,
                    lc_amin, lc_amax,
                    CONF(prr, lc_smesh_max_nside),
                    fastpm->comm);
        } else {
            /* FIXME: use n, not nc */
            ptrdiff_t Nc1[3] = {CONF(prr, nc), CONF(prr, nc), CONF(prr, nc)};
            fastpm_smesh_add_layer_pm(*smesh, fastpm->basepm, NULL, Nc1, lc_amin, lc_amax);
        }

        fastpm_add_event_handler(&fastpm->event_handlers,
                FASTPM_EVENT_FORCE, FASTPM_EVENT_STAGE_AFTER,
                (FastPMEventHandlerFunction) smesh_force_handler,
                *smesh);

        size_t Nslices;
        FastPMSMeshSlice * slices = fastpm_smesh_get_aemit(*smesh, &Nslices);

        double * aedges = _get_aedges(slices, Nslices);

        /* use the centers of the layers to ensure each bin has exactly one layer */

        fastpm_info("Structured Mesh have %d layers, from a=%g to %g\n", Nslices, slices[0].aemit, slices[Nslices-1].aemit);

        struct smesh_ready_handler_data * data = malloc(sizeof(data[0]));
        data->fastpm = fastpm;
        data->prr = prr;
        data->slices = slices;
        data->Nslices = Nslices;

        data->aedges = aedges;
        data->Nedges = Nslices - 1;
        data->hist = calloc(data->Nedges + 1, sizeof(int64_t));

        fastpm_add_event_handler_free(&(*smesh)->event_handlers,
                FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
                (FastPMEventHandlerFunction) smesh_ready_handler,
                data, _smesh_ready_handler_free);

    }

    *usmesh = NULL;
    if(CONF(prr, lc_write_usmesh)) {
        *usmesh = malloc(sizeof(FastPMUSMesh));

        double (*tiles)[3];
        int ntiles;

        if(CONF(prr, ndim_lc_usmesh_tiles) != 2 ||
           CONF(prr, shape_lc_usmesh_tiles)[1] != 3
        ) {
            fastpm_raise(-1, "tiles must be a nx3 matrix, one row per tile.\n");
        }

        ntiles = CONF(prr, shape_lc_usmesh_tiles)[0];
        tiles = malloc(sizeof(tiles[0]) * ntiles);
        int i, j;
        double * c = CONF(prr, lc_usmesh_tiles);

        for (i = 0; i < ntiles; i ++) {
            for (j = 0; j < 3; j ++) {
                tiles[i][j] = (*c) * pm_boxsize(fastpm->basepm)[j];
                c ++;
            }
            fastpm_info("Lightcone tiles[%d] : %g %g %g\n", i,
                tiles[i][0], tiles[i][1], tiles[i][2]);
        }
        fastpm_usmesh_init(*usmesh, lc,
                CONF(prr, lc_usmesh_alloc_factor) * fastpm->p->np_upper,
                tiles, ntiles, lc_amin, lc_amax);

        fastpm_add_event_handler(&fastpm->event_handlers,
            FASTPM_EVENT_INTERPOLATION,
            FASTPM_EVENT_STAGE_BEFORE,
            (FastPMEventHandlerFunction) check_lightcone,
            *usmesh);

        free(tiles);

        struct usmesh_ready_handler_data * data = malloc(sizeof(data[0]));
        data->fastpm = fastpm;
        data->prr = prr;

        if (*smesh) {
            size_t Nslices;
            FastPMSMeshSlice * slices = fastpm_smesh_get_aemit(*smesh, &Nslices);
            double * aedges = _get_aedges(slices, Nslices);
            size_t Nedges = Nslices - 1;
            data->aedges = aedges;
            data->Nedges = Nedges;
            free(slices);
        } else {
            double amin = CONF(prr, lc_amin);
            double amax = CONF(prr, lc_amax);
            fastpm_info("No structured mesh requested, generating an AemitIndex with 128 layers for usmesh. \n");

            int nedges = 128;
            double * edges = malloc(sizeof(double) * (nedges));

            int i;

            for(i = 0; i < nedges - 1; i ++) {
                edges[i] = (amax - amin) * i / (nedges - 1) + amin;
            }

            edges[nedges - 1] = amax;
            data->aedges = edges;
            data->Nedges = nedges;
        }

        data->hist = calloc(data->Nedges + 1, sizeof(int64_t));

        fastpm_store_init(data->tail, 0, 0, FASTPM_MEMORY_FLOATING);

        fastpm_add_event_handler_free(&(*usmesh)->event_handlers,
                FASTPM_EVENT_LC_READY, FASTPM_EVENT_STAGE_AFTER,
                (FastPMEventHandlerFunction) usmesh_ready_handler,
                data, _usmesh_ready_handler_free);
    }

}

static void
smesh_ready_handler(FastPMSMesh * mesh, FastPMLCEvent * lcevent, struct smesh_ready_handler_data * data)
{
    CLOCK(io);
    CLOCK(sort);
    CLOCK(indexing);

    FastPMSolver * solver = data->fastpm;
    Parameters * prr = data->prr;

    int64_t np = lcevent->p->np;
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG, MPI_SUM, solver->comm);

    fastpm_info("Structured LightCone ready : a0 = %g a1 = %g, n = %td\n", lcevent->a0, lcevent->a1, np);

    char * fn = fastpm_strdup_printf(CONF(prr, lc_write_smesh));

    ENTER(sort);
    fastpm_sort_snapshot(lcevent->p, solver->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    if(lcevent->is_first) {
        fastpm_info("Creating smesh catalog in %s\n", fn);
        write_snapshot(solver, lcevent->p, fn, "1", "", prr->Nwriters);
        size_t na;
        double * a = fastpm_smesh_get_aemit(mesh, &na);
        write_snapshot_attr(fn, "Header", "LCLayersAemit", a, "f8", na, solver->comm);
        free(a);
    } else {
        fastpm_info("Appending smesh catalog to %s\n", fn);
        append_snapshot(solver, lcevent->p, fn, "1", "", prr->Nwriters);
    }

    LEAVE(io);
    ENTER(indexing);

    fastpm_store_histogram_aemit(lcevent->p, data->hist, data->aedges, data->Nedges, solver->comm);

    write_aemit_hist(fn, "1/.", data->hist, data->aedges, data->Nedges, solver->comm);

    LEAVE(indexing);

    free(fn);
}

static void
usmesh_ready_handler(FastPMUSMesh * mesh, FastPMLCEvent * lcevent, struct usmesh_ready_handler_data * data)
{
    CLOCK(io);
    CLOCK(sort);
    CLOCK(indexing);

    FastPMSolver * solver = data->fastpm;
    Parameters * prr = data->prr;
    FastPMStore * tail = data->tail;

    int64_t np = lcevent->p->np;
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG, MPI_SUM, solver->comm);

    fastpm_info("Unstructured LightCone ready : a0 = %g a1 = %g, n = %td\n", lcevent->a0, lcevent->a1, np);

    char * fn = fastpm_strdup_printf(CONF(prr, lc_write_usmesh));

    write_usmesh_fof(solver, lcevent->p, fn, prr, lcevent, tail, mesh->lc);

    /* subsample, this will remove the tail particles that were appended. */
    fastpm_store_subsample(lcevent->p, lcevent->p->mask, lcevent->p);

    ENTER(sort);
    fastpm_sort_snapshot(lcevent->p, solver->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    if(lcevent->is_first) {
        fastpm_info("Creating usmesh catalog in %s\n", fn);
        write_snapshot(solver, lcevent->p, fn, "1", "", prr->Nwriters);
        write_snapshot_attr(fn, "Header", "AemitGrid", data->aedges, "f8", data->Nedges, solver->comm);
    } else {
        fastpm_info("Appending usmesh catalog to %s\n", fn);
        append_snapshot(solver, lcevent->p, fn, "1", "", prr->Nwriters);
    }

    LEAVE(io);
    ENTER(indexing);

    fastpm_store_histogram_aemit(lcevent->p, data->hist, data->aedges, data->Nedges, solver->comm);

    write_aemit_hist(fn, "1/.", data->hist, data->aedges, data->Nedges, solver->comm);

    LEAVE(indexing);

    free(fn);
}

/* bridging force event to smesh interpolation */
static void
smesh_force_handler(FastPMSolver * solver, FastPMForceEvent * event, FastPMSMesh * smesh)
{
    CLOCK(potential);
    ENTER(potential);
    fastpm_smesh_compute_potential(smesh, event->pm, event->gravity, event->delta_k, event->a_f, event->a_n);
    LEAVE(potential);
}

static int
check_snapshots(FastPMSolver * fastpm, FastPMInterpolationEvent * event, Parameters * prr)
{
    CLOCK(io);

    fastpm_info("Checking Snapshots (%0.4f %0.4f) with K(%0.4f %0.4f %0.4f) D(%0.4f %0.4f %0.4f)\n",
        event->a1, event->a2,
        event->kick->af, event->kick->ai, event->kick->ac,
        event->drift->af, event->drift->ai, event->drift->ac
    );

    /* interpolate and write snapshots, assuming p 
     * is at time a_x and a_v. */
    FastPMStore * p = fastpm->p;
    int nout = CONF(prr, n_aout);
    double * aout= CONF(prr, aout);

    int iout;
    for(iout = 0; iout < nout; iout ++) {
        if(event->a1 == event->a2) {
            /* initial condition */
            if(event->a1 != aout[iout]) continue;
        } else {
            if(event->a1 >= aout[iout]) continue;
            if(event->a2 < aout[iout]) continue;
        }

        FastPMStore snapshot[1];

        fastpm_store_init(snapshot, p->np_upper,
                p->attributes & ~PACK_ACC,
                FASTPM_MEMORY_FLOATING
            );

        fastpm_info("Setting up snapshot at a = %6.4f (z=%6.4f)\n", aout[iout], 1.0f/aout[iout]-1);
        fastpm_info("Growth factor of snapshot %6.4f (a=%0.4f)\n", fastpm_solver_growth_factor(fastpm, aout[iout]), aout[iout]);

        fastpm_set_snapshot(fastpm, event->drift, event->kick, snapshot, aout[iout]);

        if(CONF(prr, write_fof)) {
            char filebase[1024];

            sprintf(filebase, "%s_%0.04f", CONF(prr, write_fof), aout[iout]);

            write_fof(fastpm, snapshot, filebase, prr);
        }

        /* in place subsampling to avoid creating another store, we trash it immediately anyways. */
        fastpm_store_subsample(snapshot, snapshot->mask, snapshot);

        take_a_snapshot(fastpm, snapshot, aout[iout], prr);

        fastpm_store_destroy(snapshot);

    }
    return 0;
}

static void
write_fof(FastPMSolver * fastpm, FastPMStore * snapshot, char * filebase, Parameters * prr)
{
    CLOCK(fof);
    CLOCK(io);
    CLOCK(sort);
    ENTER(fof);
    FastPMFOFFinder fof = {
        /* convert from fraction of mean separation to simulation distance units. */
        .linkinglength = CONF(prr, fof_linkinglength) * CONF(prr, boxsize) / CONF(prr, nc),
        .periodic = 1,
        .nmin = CONF(prr, fof_nmin),
        .kdtree_thresh = CONF(prr, fof_kdtree_thresh),
    };

    fastpm_fof_init(&fof, snapshot, fastpm->basepm);

    FastPMStore halos[1];

    ENTER(fof);

    fastpm_fof_execute(&fof, halos);

    LEAVE(fof);

    fastpm_info("Writing fof %s with %d writers\n", filebase, prr->Nwriters);

    ENTER(sort);
    fastpm_sort_snapshot(halos, fastpm->comm, FastPMSnapshotSortByLength, 0);
    LEAVE(sort);

    ENTER(io);

    char * dataset = fastpm_strdup_printf("LL-%05.3f", CONF(prr, fof_linkinglength));
    write_snapshot(fastpm, halos, filebase, dataset, prr->string, prr->Nwriters);

    free(dataset);
    LEAVE(io);

    fastpm_info("FOF Catalog %s written\n", filebase);

    fastpm_store_destroy(halos);
    fastpm_fof_destroy(&fof);
}

static void
_halos_ready (FastPMFOFFinder * finder,
    FastPMHaloEvent * event, void ** userdata)
{
    double rmin = *((double*) userdata[0]);
    double halosize = *((double*) userdata[1]);
    FastPMLightCone * lc = (FastPMLightCone*) userdata[2];
    char * keep_for_tail = (char *) userdata[3];
    FastPMStore * halos = event->halos;
    FastPMStore * p = event->p;

    ptrdiff_t i;

    uint8_t * halos_mask = malloc(halos->np);

    for(i = 0; i < halos->np; i ++) {
        double r = fastpm_lc_distance(lc, halos->x[i]);
        /* only keep reliable halos */
        if(r > rmin + halosize * 0.5) {
            halos_mask[i] = 1;
        } else {
            halos_mask[i] = 0;
        }
    }

    for(i = 0; i < p->np; i ++) {
        ptrdiff_t hid = event->ihalo[i];
        double r_p = fastpm_lc_distance(lc, p->x[i]);
        if (r_p > rmin + halosize) {
            /* not near the tail of the lightcone, do not keep for tail */
            keep_for_tail[i] = 0;
        } else {
            /* near the tail of the lightcone, do not keep those formed reliable halos */
            if(hid >= 0) {
                /* unreliable halos, keep them */
                keep_for_tail[i] = !halos_mask[hid];
            } else {
                /* not in halos, keep them too */
                keep_for_tail[i] = 1;
            }
        }
    }
    userdata[4] = halos_mask;
}


static void
write_usmesh_fof(FastPMSolver * fastpm,
        FastPMStore * snapshot,
        char * filebase, Parameters * prr,
        FastPMLCEvent * lcevent,
        FastPMStore * tail,
        FastPMLightCone * lc)
{
    CLOCK(fof);
    CLOCK(io);
    CLOCK(sort);

    int append = !lcevent->is_first;
    double maxhalosize = CONF(prr, lc_usmesh_fof_padding); /* MPC/h, used to cut along z direction. */
    FastPMStore * p = lcevent->p;
    ptrdiff_t i;

    /* not the first segment, add the left over particles to the FOF */

    /* avoid including the tail particles in the subsample;
     * by clearing their mask */
    for(i = 0; i < tail->np; i ++) {
        tail->mask[i] = 0;
    }
    fastpm_store_append(tail, lcevent->p);
    fastpm_store_destroy(tail);

    /* FIXME: register event to mask out particles*/
    uint8_t * keep_for_tail = fastpm_memory_alloc(p->mem, "keep", lcevent->p->np_upper, FASTPM_MEMORY_HEAP);

    ENTER(fof);

    FastPMFOFFinder fof = {
        /* convert from fraction of mean separation to simulation distance units. */
        .linkinglength = CONF(prr, fof_linkinglength) * CONF(prr, boxsize) / CONF(prr, nc),
        .periodic = 0,
        .nmin = CONF(prr, fof_nmin),
        .kdtree_thresh = CONF(prr, fof_kdtree_thresh),
    };

    fastpm_fof_init(&fof, snapshot, fastpm->basepm);

    double rmin = lc->speedfactor * HorizonDistance(lcevent->a1, lc->horizon);

    void * userdata[5];
    userdata[0] = & rmin;
    userdata[1] = & maxhalosize;
    userdata[2] = lc;
    userdata[3] = keep_for_tail;
    userdata[4] = NULL; /* halos_mask , return value of the handler */

    FastPMStore halos[1];
    fastpm_add_event_handler(&fof.event_handlers,
        FASTPM_EVENT_HALO,
        FASTPM_EVENT_STAGE_AFTER,
        (FastPMEventHandlerFunction) _halos_ready,
        userdata);

    ENTER(fof);

    fastpm_fof_execute(&fof, halos);

    LEAVE(fof);

    fastpm_info("Writing fof %s with %d writers\n", filebase, prr->Nwriters);

    uint8_t * halos_mask = userdata[4];
    fastpm_store_subsample(halos, halos_mask, halos);
    free(halos_mask);

    ENTER(sort);
    fastpm_sort_snapshot(halos, fastpm->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    char * dataset = fastpm_strdup_printf("LL-%05.3f", CONF(prr, fof_linkinglength));
    if(!append) {
        write_snapshot(fastpm, halos, filebase, dataset, prr->string, prr->Nwriters);
    } else {
        append_snapshot(fastpm, halos, filebase, dataset, prr->string, prr->Nwriters);
    }
    free(dataset);
    LEAVE(io);

    uint64_t nhalos = halos->np;
    MPI_Allreduce(MPI_IN_PLACE, &nhalos, 1, MPI_LONG, MPI_SUM, fastpm->comm);

    fastpm_info("FOF Catalog %s written with %td halos\n", filebase, nhalos);

    fastpm_store_destroy(halos);

    fastpm_fof_destroy(&fof);

    uint64_t ntail = 0;
    for(i = 0; i < p->np; i ++) {
        if(keep_for_tail[i]) ntail ++;
    }

    fastpm_store_init(tail, ntail, p->attributes, FASTPM_MEMORY_FLOATING);
    fastpm_store_subsample(p, keep_for_tail, tail);

    MPI_Allreduce(MPI_IN_PLACE, &ntail, 1, MPI_LONG, MPI_SUM, fastpm->comm);
    fastpm_info("%td particles will be reused in next batch\n", ntail);
    fastpm_memory_free(p->mem, keep_for_tail);
}

static int 
take_a_snapshot(FastPMSolver * fastpm, FastPMStore * snapshot, double aout, Parameters * prr) 
{
    CLOCK(io);
    CLOCK(sort);

    /* do this before write_snapshot, because white_snapshot messes up with the domain decomposition. */
    if(CONF(prr, write_nonlineark)) {
        char * filename = fastpm_strdup_printf("%s_%0.04f", CONF(prr, write_nonlineark), aout);
        FastPMPainter painter[1];

        FastPMFloat * rho_x = pm_alloc(fastpm->basepm);
        FastPMFloat * rho_k = pm_alloc(fastpm->basepm);

        fastpm_painter_init(painter, fastpm->basepm, fastpm->config->PAINTER_TYPE, fastpm->config->painter_support);

        fastpm_paint(painter, rho_x, snapshot, 0);
        pm_r2c(fastpm->basepm, rho_x, rho_k);

        write_complex(fastpm->basepm, rho_k, filename, "DensityK", prr->Nwriters);

        pm_free(fastpm->basepm, rho_k);
        pm_free(fastpm->basepm, rho_x);
        free(filename);
    }

    if(CONF(prr, write_snapshot)) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;
        sprintf(filebase, "%s_%0.04f", CONF(prr, write_snapshot), aout);

        fastpm_info("Writing snapshot %s at z = %6.4f a = %6.4f with %d writers\n", 
                filebase, z_out, aout, prr->Nwriters);

        ENTER(sort);
        if(CONF(prr, sort_snapshot)) {
            fastpm_info("Snapshot is sorted by ID.\n");
            fastpm_sort_snapshot(snapshot, fastpm->comm, FastPMSnapshotSortByID, 1);
        } else {
            fastpm_info("Snapshot is not sorted by ID.\n");
        }
        LEAVE(sort);

        ENTER(io);
        write_snapshot(fastpm, snapshot, filebase, "1", prr->string, prr->Nwriters);
        LEAVE(io);

        fastpm_info("snapshot %s written\n", filebase);
    }
    if(CONF(prr, write_runpb_snapshot)) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;

        sprintf(filebase, "%s_%0.04f.bin", CONF(prr, write_runpb_snapshot), aout);

        ENTER(io);
        write_runpb_snapshot(fastpm, snapshot, filebase);

        LEAVE(io);

        fastpm_info("snapshot %s written z = %6.4f a = %6.4f\n", 
                filebase, z_out, aout);

    }
    return 0;
}

static int 
check_lightcone(FastPMSolver * fastpm, FastPMInterpolationEvent * event, FastPMUSMesh * usmesh)
{
    fastpm_usmesh_intersect(usmesh, event->drift, event->kick, fastpm);

    int64_t np = usmesh->p->np;

    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_LONG, MPI_SUM, fastpm->comm);

    fastpm_info("Total number of particles in light cone: %ld\n", np);

    return 0;
}

static int 
print_transition(FastPMSolver * fastpm, FastPMTransitionEvent * event, Parameters * prr)
{
    FastPMTransition * trans = event->transition;
    char * action;
    switch (trans->action) {
        case FASTPM_ACTION_FORCE:
            action = "FORCE";
        break;
        case FASTPM_ACTION_KICK:
            action = "KICK";
        break;
        case FASTPM_ACTION_DRIFT:
            action = "DRIFT";
        break;
        default:
            action = "Unknown";
        break;
    }
    fastpm_info("==== -> %03d [%03d %03d %03d] a_i = %6.4f a_f = %6.4f a_r = %6.4f Action = %s(%d) ====\n",
            trans->iend,
            trans->end->x,
            trans->end->v,
            trans->end->force,
            trans->a.i, trans->a.f, trans->a.r, action, trans->action);
    return 0;
}

static void
_memory_peak_handler(FastPMMemory * mem, void * userdata)
{
    MPI_Comm comm = *(MPI_Comm*) userdata;

    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    if(ThisTask == 0) {
        fastpm_ilog(INFO, "Peak memory usage on rank %d: %g MB\n", ThisTask, mem->peak_bytes / 1024. / 1024);
        fastpm_ilog(INFO, "Allocation Table\n");
        fastpm_memory_dump_status(mem, 1);
    }
}

static int
write_powerspectrum(FastPMSolver * fastpm, FastPMForceEvent * event, Parameters * prr) 
{

    int K_LINEAR = CONF(prr, enforce_broadband_kmax);

    CLOCK(compute);
    CLOCK(io);

    fastpm_info("Force Calculation Nmesh = %d ====\n", pm_nmesh(event->pm)[0]);

    fastpm_info("Load imbalance: min = %g max = %g std = %g\n",
        fastpm->info.imbalance.min, fastpm->info.imbalance.max, fastpm->info.imbalance.std);

    MPI_Barrier(fastpm->comm);
    ENTER(compute);

    FastPMPowerSpectrum ps;
    /* calculate the power spectrum */
    fastpm_powerspectrum_init_from_delta(&ps, event->pm, event->delta_k, event->delta_k);

    double Plin = fastpm_powerspectrum_large_scale(&ps, K_LINEAR);

    double Sigma8 = fastpm_powerspectrum_sigma(&ps, 8);

    Plin /= pow(fastpm_solver_growth_factor(fastpm, event->a_f), 2.0);
    Sigma8 /= pow(fastpm_solver_growth_factor(fastpm, event->a_f), 2.0);

    fastpm_info("D^2(%g, 1.0) P(k<%g) = %g Sigma8 = %g\n", event->a_f, K_LINEAR * 6.28 / pm_boxsize(event->pm)[0], Plin, Sigma8);

    LEAVE(compute);

    MPI_Barrier(fastpm->comm);

    ENTER(io);
    if(CONF(prr, write_powerspectrum)) {
        char buf[1024];
        sprintf(buf, "%s_%0.04f.txt", CONF(prr, write_powerspectrum), event->a_f);
        fastpm_info("writing power spectrum to %s\n", buf);
        if(fastpm->ThisTask == 0) {
            fastpm_path_ensure_dirname(CONF(prr, write_powerspectrum));
            fastpm_powerspectrum_write(&ps, buf, event->N);
        }
    }
    LEAVE(io);

    fastpm_powerspectrum_destroy(&ps);

    return 0;
}

int
read_powerspectrum(FastPMPowerSpectrum * ps, const char filename[], const double sigma8, MPI_Comm comm)
{
    fastpm_info("Powerspecectrum file: %s\n", filename);

    int myrank;
    MPI_Comm_rank(comm, &myrank);
    char * content;
    if(myrank == 0) {
        content = fastpm_file_get_content(filename);
        if(content == NULL) {
            fastpm_raise(-1, "Failed to read powerspectrum from file %s\n", filename);
        }
        int size = strlen(content);
        MPI_Bcast(&size, 1, MPI_INT, 0, comm);
        MPI_Bcast(content, size + 1, MPI_BYTE, 0, comm);
    } else {
        int size = 0;
        MPI_Bcast(&size, 1, MPI_INT, 0, comm);
        content = malloc(size + 1);
        MPI_Bcast(content, size + 1, MPI_BYTE, 0, comm);
    }
    if (0 != fastpm_powerspectrum_init_from_string(ps, content)) {
        fastpm_raise(-1, "Failed to parse the powerspectrum\n");
    }
    free(content);

    fastpm_info("Found %d pairs of values in input spectrum table\n", ps->size);

    double sigma8_input= fastpm_powerspectrum_sigma(ps, 8);
    fastpm_info("Input power spectrum sigma8 %f\n", sigma8_input);

    if(sigma8 > 0) {
        fastpm_info("Expected power spectrum sigma8 %g; correction applied. \n", sigma8);
        fastpm_powerspectrum_scale(ps, pow(sigma8 / sigma8_input, 2));
    }
    return 0;
}


static void 
parse_args(int * argc, char *** argv, Parameters * prr) 
{
    int opt;
    extern int optind;
    extern char * optarg;
    prr->UseFFTW = 0;
    ParamFileName = NULL;
    prr->NprocY = 0;
    prr->Nwriters = 0;
    prr->MemoryPerRank = 0;
    while ((opt = getopt(*argc, *argv, "h?y:fW:m:")) != -1) {
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
            case 'm':
                prr->MemoryPerRank = atoi(optarg);
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
    printf("Usage: fastpm [-W Nwriters] [-f] [-y NprocY] [-m MemoryBoundInMB] paramfile\n"
    "-f Use FFTW \n"
    "-y Set the number of processes in the 2D mesh\n"
    "-n Throttle IO (bigfile only) \n"
);
    MPI_Finalize();
    exit(1);
}

int read_parameters(char * filename, Parameters * param, int argc, char ** argv, MPI_Comm comm)
{
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    /* read the configuration file */
    char * confstr;
    int confstr_len;

    /* run the parameter file on root rank.
     * other ranks use the serialized string to avoid duplicated
     * error reports */
    if(ThisTask == 0) {
        char * error;
        confstr = lua_config_parse("_parse", filename, argc, argv, &error);
        if(confstr == NULL) {
            fastpm_raise(-1, "%s\n", error);
        }
        confstr_len = strlen(confstr) + 1;
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    } else {
        MPI_Bcast(&confstr_len, 1, MPI_INT, 0, comm);
        confstr = malloc(confstr_len);
        MPI_Bcast(confstr, confstr_len, MPI_BYTE, 0, comm);
    }

    fastpm_info("Configuration %s\n", confstr);

    param->config = lua_config_new(confstr);
    param->string = confstr;
    if(lua_config_error(param->config)) {
        fastpm_raise(-1, "error: %s\n", lua_config_error(param->config));
    }
    return 0;
}

static void
free_parameters(Parameters * param)
{
    lua_config_free(param->config);
    free(param->string);
}


