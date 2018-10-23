#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>
#include <signal.h>
#include <limits.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>
#include <fastpm/lightcone.h>
#include <fastpm/constrainedgaussian.h>
#include <fastpm/io.h>
#include <fastpm/fof.h>
#include <bigfile-mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lua-config.h"
#include "param.h"

/* c99 has no pi. */
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

typedef struct {
    CLIParameters * cli;
    LUAParameters * lua;
} Parameters;


extern void
init_stacktrace();

static void
_memory_peak_handler(FastPMMemory * mem, void * userdata);

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

static void
_write_parameters(const char * filebase, const char * dataset, Parameters * prr, MPI_Comm comm)
{
    BigFile bf;
    BigBlock bb;
    if(0 != big_file_mpi_open(&bf, filebase, comm)) {
        fastpm_raise(-1, "Failed to open the file: %s\n", big_file_get_error_message());
    }

    if(0 != big_file_mpi_open_block(&bf, &bb, dataset, comm)) {
        fastpm_raise(-1, "Failed to open the dataset : %s\n", big_file_get_error_message());
    }

    big_block_set_attr(&bb, "ParamFile", prr->lua->string, "S1", strlen(prr->lua->string) + 1);
    double ParticleFraction = CONF(prr->lua, particle_fraction);
    big_block_set_attr(&bb, "ParticleFraction", &ParticleFraction, "f8", 1);

    big_block_mpi_close(&bb, comm);
    big_file_mpi_close(&bf, comm);

}

int run_fastpm(FastPMConfig * config, Parameters * prr, MPI_Comm comm);

int main(int argc, char ** argv) {

    init_stacktrace();

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD; 

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    fastpm_info("This is FastPM, with libfastpm version %s.\n", LIBFASTPM_VERSION);

    char * error;
    CLIParameters * cli = parse_cli_args_mpi(argc, argv, comm);
    LUAParameters * lua = parse_config_mpi(cli->argv[0], cli->argc, cli->argv, &error, comm);

    Parameters prr[1] = {{cli, lua}};

    if(prr->lua) {
        fastpm_info("Configuration %s\n", prr->lua->string);
    } else {
        fastpm_info("Parsing configuration failed with error: %s\n", error);
        MPI_Finalize();
        exit(1);
    }

    libfastpm_set_memory_bound(prr->cli->MemoryPerRank * 1024 * 1024);
    fastpm_memory_set_handlers(_libfastpm_get_gmem(), NULL, _memory_peak_handler, &comm);

    /* convert parameter files pm_nc_factor into VPMInit */
    VPMInit * vpminit = NULL;
    if(CONF(prr->lua, ndim_pm_nc_factor) == 0) {
        vpminit = (VPMInit[]) {
            {.a_start = 0, .pm_nc_factor = CONF(prr->lua, pm_nc_factor)[0] },
            {.a_start = 1, .pm_nc_factor = 0},
            };
    } else
    if(CONF(prr->lua, ndim_pm_nc_factor) == 2) {
        vpminit = alloca(sizeof(VPMInit) * (CONF(prr->lua, shape_pm_nc_factor)[0] + 1));
        int i;
        for(i = 0; i < CONF(prr->lua, n_pm_nc_factor); i ++) {
            vpminit[i].a_start = CONF(prr->lua, pm_nc_factor)[2 * i];
            vpminit[i].pm_nc_factor = CONF(prr->lua, pm_nc_factor)[2 * i + 1];
        }
        /* mark the end */
        vpminit[i].pm_nc_factor = 0;
    } else {
        fastpm_raise(-1, "Unknown format of pm_nc_factor, either a scalar or a 2d array. ");
    }

    fastpm_info("np_alloc_factor = %g\n", CONF(prr->lua, np_alloc_factor));

    FastPMConfig * config = & (FastPMConfig) {
        .nc = CONF(prr->lua, nc),
        .alloc_factor = CONF(prr->lua, np_alloc_factor),
        .vpminit = vpminit,
        .boxsize = CONF(prr->lua, boxsize),
        .omega_m = CONF(prr->lua, omega_m),
        .hubble_param = CONF(prr->lua, h),
        .USE_DX1_ONLY = CONF(prr->lua, za),
        .nLPT = -2.5f,
        .USE_SHIFT = CONF(prr->lua, shift),
        .FORCE_TYPE = CONF(prr->lua, force_mode),
        .KERNEL_TYPE = CONF(prr->lua, kernel_type),
        .DEALIASING_TYPE = CONF(prr->lua, dealiasing_type),
        .PAINTER_TYPE = CONF(prr->lua, painter_type),
        .painter_support = CONF(prr->lua, painter_support),
        .NprocY = prr->cli->NprocY,
        .UseFFTW = prr->cli->UseFFTW,
        .ExtraAttributes = 0,
    };

    if(CONF(prr->lua, compute_potential)) {
        config->ExtraAttributes |= COLUMN_POTENTIAL;
    }

    run_fastpm(config, prr, comm);

    free_lua_parameters(prr->lua);
    free_cli_parameters(prr->cli);

    libfastpm_cleanup();

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
report_memory(MPI_Comm);

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
    const double M0 = CONF(prr->lua, omega_m) * rho_crit
                    * pow(CONF(prr->lua, boxsize) / CONF(prr->lua, nc), 3.0);
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
        .speedfactor = CONF(prr->lua, dh_factor),
        .cosmology = fastpm->cosmology,
        .fov = CONF(prr->lua, lc_fov),
        .octants = {0, 0, 0, 0, 0, 0, 0, 0},
        .tol = 2. / CONF(prr->lua, nc) / CONF(prr->lua, lc_smesh_fraction),
    }};

    FastPMUSMesh * usmesh = NULL;
    FastPMSMesh * smesh = NULL;

    prepare_lc(fastpm, prr, lc, &usmesh, &smesh);

    MPI_Barrier(comm);

    ENTER(ic);
    prepare_ic(fastpm, prr, comm);

    fastpm_store_fill_subsample_mask(fastpm->p, CONF(prr->lua, particle_fraction), fastpm->p->mask, comm);

    fastpm_info("dx1  : %g %g %g %g\n", 
            fastpm->info.dx1[0], fastpm->info.dx1[1], fastpm->info.dx1[2],
            (fastpm->info.dx1[0] + fastpm->info.dx1[1] + fastpm->info.dx1[2]) / 3.0);
    fastpm_info("dx2  : %g %g %g %g\n", 
            fastpm->info.dx2[0], fastpm->info.dx2[1], fastpm->info.dx2[2],
            (fastpm->info.dx2[0] + fastpm->info.dx2[1] + fastpm->info.dx2[2]) / 3.0);

    LEAVE(ic);

    MPI_Barrier(comm);
    ENTER(evolve);
    fastpm_solver_evolve(fastpm, CONF(prr->lua, time_step), CONF(prr->lua, n_time_step));
    LEAVE(evolve);

    if(smesh)
        fastpm_smesh_destroy(smesh);

    if(usmesh)
        fastpm_usmesh_destroy(usmesh);

    free(smesh);
    free(usmesh);

    fastpm_lc_destroy(lc);

    fastpm_solver_destroy(fastpm);

    report_memory(comm);

    fastpm_clock_stat(comm);

    return 0;
}

static void 
prepare_ic(FastPMSolver * fastpm, Parameters * prr, MPI_Comm comm) 
{
    /* we may need a read gadget ic here too */
    if(CONF(prr->lua, read_runpbic)) {
        FastPMStore * p = fastpm->p;
        int temp_dx1 = 0;
        int temp_dx2 = 0;
        if(p->dx1 == NULL) {
            p->dx1 = fastpm_memory_alloc(p->mem, "DX1", sizeof(p->dx1[0]) * p->np_upper, FASTPM_MEMORY_STACK);
            temp_dx1 = 1;
        }
        if(p->dx2 == NULL) {
            p->dx2 = fastpm_memory_alloc(p->mem, "DX2", sizeof(p->dx2[0]) * p->np_upper, FASTPM_MEMORY_STACK);
            temp_dx2 = 1;
        }

        read_runpb_ic(fastpm, fastpm->p, CONF(prr->lua, read_runpbic));
        fastpm_solver_setup_ic(fastpm, NULL, CONF(prr->lua, time_step)[0]);
        if(temp_dx2) {
            fastpm_memory_free(p->mem, p->dx2);
            p->dx2 = NULL;
        }
        if(temp_dx1) {
            fastpm_memory_free(p->mem, p->dx1);
            p->dx1 = NULL;
        }
        return;
    } 

    /* at this point generating the ic involves delta_k */
    FastPMFloat * delta_k = pm_alloc(fastpm->basepm);

    if(CONF(prr->lua, read_lineark)) {
        fastpm_info("Reading Fourier space linear overdensity from %s\n", CONF(prr->lua, read_lineark));
        read_complex(fastpm->basepm, delta_k, CONF(prr->lua, read_lineark), "LinearDensityK", prr->cli->Nwriters);

        if(CONF(prr->lua, inverted_ic)) {
            fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, -1);
        }
        goto produce;
    }

    /* at this power we need a powerspectrum file to convolve the guassian */
    if(!CONF(prr->lua, read_powerspectrum)) {
        fastpm_raise(-1, "Need a power spectrum to start the simulation.\n");
    }

    FastPMPowerSpectrum linear_powerspectrum;

    read_powerspectrum(&linear_powerspectrum, CONF(prr->lua, read_powerspectrum), CONF(prr->lua, sigma8), comm);

    if(CONF(prr->lua, read_grafic)) {
        fastpm_info("Reading grafic white noise file from '%s'.\n", CONF(prr->lua, read_grafic));
        fastpm_info("GrafIC noise is Fortran ordering. FastPMSolver is in C ordering.\n");
        fastpm_info("The simulation will be transformed x->z y->y z->x.\n");

        FastPMFloat * g_x = pm_alloc(fastpm->basepm);

        read_grafic_gaussian(fastpm->basepm, g_x, CONF(prr->lua, read_grafic));

        /* r2c will reduce the variance. Compensate here.*/
        fastpm_apply_multiply_transfer(fastpm->basepm, g_x, g_x, sqrt(pm_norm(fastpm->basepm)));
        pm_r2c(fastpm->basepm, g_x, delta_k);

        pm_free(fastpm->basepm, g_x);

        goto induce;
    }

    if(CONF(prr->lua, read_whitenoisek)) {
        fastpm_info("Reading Fourier white noise file from '%s'.\n", CONF(prr->lua, read_whitenoisek));

        read_complex(fastpm->basepm, delta_k, CONF(prr->lua, read_whitenoisek), "WhiteNoiseK", prr->cli->Nwriters);
        goto induce;
    }

    /* Nothing to read from, just generate a gadget IC with the seed. */
    fastpm_ic_fill_gaussiank(fastpm->basepm, delta_k, CONF(prr->lua, random_seed), FASTPM_DELTAK_GADGET);

induce:
    if(CONF(prr->lua, remove_cosmic_variance)) {
        fastpm_info("Remove Cosmic variance from initial condition.\n");
        fastpm_ic_remove_variance(fastpm->basepm, delta_k);
    }

    if(CONF(prr->lua, set_mode)) {
        int method = 0;
        /* FIXME: use enums */
        if(0 == strcmp(CONF(prr->lua, set_mode_method), "add")) {
            method = 1;
            fastpm_info("SetMode is add\n");
        } else {
            fastpm_info("SetMode is override\n");
        }
        int i;
        double * c = CONF(prr->lua, set_mode);
        for(i = 0; i < CONF(prr->lua, n_set_mode); i ++) {
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

    if(CONF(prr->lua, inverted_ic)) {
        fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, -1);
    }

    double variance = pm_compute_variance(fastpm->basepm, delta_k);
    fastpm_info("Variance of input white noise is %0.8f, expectation is %0.8f\n", variance, 1.0 - 1.0 / pm_norm(fastpm->basepm));

    if(CONF(prr->lua, write_whitenoisek)) {
        fastpm_info("Writing Fourier white noise to file '%s'.\n", CONF(prr->lua, write_whitenoisek));
        write_complex(fastpm->basepm, delta_k, CONF(prr->lua, write_whitenoisek), "WhiteNoiseK", prr->cli->Nwriters);
    }

    /* introduce correlation */
    if(CONF(prr->lua, f_nl_type) == FASTPM_FNL_NONE) {
        fastpm_info("Inducing correlation to the white noise.\n");

        fastpm_ic_induce_correlation(fastpm->basepm, delta_k,
            (fastpm_fkfunc) fastpm_powerspectrum_eval2, &linear_powerspectrum);
    } else {
        double kmax_primordial;
        kmax_primordial = CONF(prr->lua, nc) / 2.0 * 2.0*M_PI/CONF(prr->lua, boxsize) * CONF(prr->lua, kmax_primordial_over_knyquist);
        fastpm_info("Will set Phi_Gaussian(k)=0 for k>=%f.\n", kmax_primordial);
        FastPMPNGaussian png = {
            .fNL = CONF(prr->lua, f_nl),
            .type = CONF(prr->lua, f_nl_type),
            .kmax_primordial = kmax_primordial,
            .pkfunc = (fastpm_fkfunc) fastpm_powerspectrum_eval2,
            .pkdata = &linear_powerspectrum,
            .h = CONF(prr->lua, h),
            .scalar_amp = CONF(prr->lua, scalar_amp),
            .scalar_spectral_index = CONF(prr->lua, scalar_spectral_index),
            .scalar_pivot = CONF(prr->lua, scalar_pivot)
        };
        fastpm_info("Inducing non gaussian correlation to the white noise.\n");
        fastpm_png_induce_correlation(&png, fastpm->basepm, delta_k);
    }

    /* The linear density field is not redshift zero, then evolve it with the model cosmology to 
     * redshift zero.
     * This matches the linear power at the given redshift, not necessarily redshift 0. */
    {
        double linear_density_redshift = CONF(prr->lua, linear_density_redshift);
        double linear_evolve = fastpm_solver_growth_factor(fastpm, 1.0) /
                               fastpm_solver_growth_factor(fastpm, 1 / (linear_density_redshift + 1));

        fastpm_info("Reference linear density is calibrated at redshift %g; multiply by %g to extract to redshift 0.\n", linear_density_redshift, linear_evolve);

        fastpm_apply_multiply_transfer(fastpm->basepm, delta_k, delta_k, linear_evolve);
    }

    /* set the mean to 1.0 */
    ptrdiff_t mode[4] = { 0, 0, 0, 0, };

    fastpm_apply_modify_mode_transfer(fastpm->basepm, delta_k, delta_k, mode, 1.0);

    /* add constraints */
    if(CONF(prr->lua, constraints)) {
        FastPM2PCF xi;

        fastpm_2pcf_from_powerspectrum(&xi, (fastpm_fkfunc) fastpm_powerspectrum_eval2, &linear_powerspectrum, CONF(prr->lua, boxsize), CONF(prr->lua, nc));

        FastPMConstrainedGaussian cg = {
            .constraints = malloc(sizeof(FastPMConstraint) * (CONF(prr->lua, n_constraints) + 1)),
        };
        fastpm_info("Applying %d constraints.\n", CONF(prr->lua, n_constraints));
        int i;
        for(i = 0; i < CONF(prr->lua, n_constraints); i ++) {
            double * c = CONF(prr->lua, constraints);
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

        if(CONF(prr->lua, write_lineark)) {
            fastpm_info("Writing fourier space linear field before constraints to %s\n", CONF(prr->lua, write_lineark));
            write_complex(fastpm->basepm, delta_k, CONF(prr->lua, write_lineark), "UnconstrainedLinearDensityK", prr->cli->Nwriters);
        }
        fastpm_cg_apply_constraints(&cg, fastpm->basepm, &xi, delta_k);

        free(cg.constraints);
    }

    fastpm_powerspectrum_destroy(&linear_powerspectrum);

    /* our write out and clean up stuff.*/
produce:

    if(CONF(prr->lua, write_lineark)) {
        fastpm_info("Writing fourier space linear field to %s\n", CONF(prr->lua, write_lineark));
        write_complex(fastpm->basepm, delta_k, CONF(prr->lua, write_lineark), "LinearDensityK", prr->cli->Nwriters);
    }

    if(CONF(prr->lua, write_powerspectrum)) {
        FastPMPowerSpectrum ps;
        /* calculate the power spectrum */
        fastpm_powerspectrum_init_from_delta(&ps, fastpm->basepm, delta_k, delta_k);

        char buf[1024];
        sprintf(buf, "%s_linear.txt", CONF(prr->lua, write_powerspectrum));
        fastpm_info("writing linear power spectrum to %s\n", buf);
        if(fastpm->ThisTask == 0) {
            fastpm_path_ensure_dirname(CONF(prr->lua, write_powerspectrum));
            fastpm_powerspectrum_write(&ps, buf, pow(fastpm->config->nc, 3.0));
        }
        fastpm_powerspectrum_destroy(&ps);
    }

    fastpm_solver_setup_ic(fastpm, delta_k, CONF(prr->lua, time_step)[0]);

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
    free(data->slices);
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
        if(CONF(prr->lua, ndim_lc_glmatrix) != 2 ||
           CONF(prr->lua, shape_lc_glmatrix)[0] != 4 ||
           CONF(prr->lua, shape_lc_glmatrix)[1] != 4
        ) {
            fastpm_raise(-1, "GL Matrix must be a 4x4 matrix\n");
        }

        int i, j;
        double * c = CONF(prr->lua, lc_glmatrix);

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
        for(i = 0; i < CONF(prr->lua, n_lc_octants); i ++) {
            int oct = CONF(prr->lua, lc_octants)[i];
            lc->octants[oct % 8] = 1;
            fastpm_info("Using Octant %d\n", oct);
        }
    }

    fastpm_lc_init(lc);


    double lc_amin = HAS(prr->lua, lc_amin)?CONF(prr->lua, lc_amin):CONF(prr->lua, time_step)[0];

    double lc_amax = HAS(prr->lua, lc_amax)?CONF(prr->lua, lc_amax):CONF(prr->lua, time_step)[CONF(prr->lua, n_time_step) - 1];

    fastpm_info("Unstructured Lightcone amin= %g amax=%g\n", lc_amin, lc_amax);

    *smesh = NULL;
    if(CONF(prr->lua, lc_write_smesh)) {
        *smesh = malloc(sizeof(FastPMSMesh));

        double n = CONF(prr->lua, nc) / CONF(prr->lua, boxsize) * CONF(prr->lua, lc_smesh_fraction);

        fastpm_smesh_init(*smesh, lc, fastpm->p->np_upper, 1 / n);

        if(lc->fov > 0) {
            fastpm_info("Creating healpix structured meshes for FOV=%g, with number density %g per (Mpc/h)**3. \n",
                lc->fov, n * n * n);
            fastpm_smesh_add_layers_healpix(*smesh,
                    n * n, n * n * n,
                    lc_amin, lc_amax,
                    CONF(prr->lua, lc_smesh_max_nside),
                    fastpm->comm);
        } else {
            /* FIXME: use n, not nc */
            ptrdiff_t Nc1[3] = {CONF(prr->lua, nc), CONF(prr->lua, nc), CONF(prr->lua, nc)};
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
    if(CONF(prr->lua, lc_write_usmesh)) {
        *usmesh = malloc(sizeof(FastPMUSMesh));

        double (*tiles)[3];
        int ntiles;

        if(CONF(prr->lua, ndim_lc_usmesh_tiles) != 2 ||
           CONF(prr->lua, shape_lc_usmesh_tiles)[1] != 3
        ) {
            fastpm_raise(-1, "tiles must be a nx3 matrix, one row per tile.\n");
        }

        ntiles = CONF(prr->lua, shape_lc_usmesh_tiles)[0];
        tiles = malloc(sizeof(tiles[0]) * ntiles);
        int i, j;
        double * c = CONF(prr->lua, lc_usmesh_tiles);

        for (i = 0; i < ntiles; i ++) {
            for (j = 0; j < 3; j ++) {
                tiles[i][j] = (*c) * pm_boxsize(fastpm->basepm)[j];
                c ++;
            }
            fastpm_info("Lightcone tiles[%d] : %g %g %g\n", i,
                tiles[i][0], tiles[i][1], tiles[i][2]);
        }
        fastpm_usmesh_init(*usmesh, lc,
                CONF(prr->lua, lc_usmesh_alloc_factor) * fastpm->p->np_upper,
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
            double amin = CONF(prr->lua, lc_amin);
            double amax = CONF(prr->lua, lc_amax);
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

    char * fn = fastpm_strdup_printf(CONF(prr->lua, lc_write_smesh));

    ENTER(sort);
    fastpm_sort_snapshot(lcevent->p, solver->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    if(lcevent->is_first) {
        fastpm_info("Creating smesh catalog in %s\n", fn);

        write_snapshot(solver, lcevent->p, fn, "1", prr->cli->Nwriters);
        _write_parameters(fn, "Header", prr, solver->comm);

        double * aemit = malloc(sizeof(double) * (data->Nslices));
        double * r = malloc(sizeof(double) * (data->Nslices));
        int * nside = malloc(sizeof(int) * (data->Nslices));

        int i;
        for(i = 0; i < data->Nslices; i ++) {
            aemit[i] = data->slices[i].aemit;
            nside[i] = data->slices[i].nside;
            r[i] = data->slices[i].distance;
        }
        write_snapshot_attr(fn, "Header", "SMeshLayers.aemit", aemit, "f8", data->Nslices, solver->comm);
        write_snapshot_attr(fn, "Header", "SMeshLayers.nside", nside, "i4", data->Nslices, solver->comm);
        write_snapshot_attr(fn, "Header", "SMeshLayers.distance", r, "f8", data->Nslices, solver->comm);

        free(nside);
        free(r);
        free(aemit);
    } else {
        fastpm_info("Appending smesh catalog to %s\n", fn);
        append_snapshot(solver, lcevent->p, fn, "1", prr->cli->Nwriters);
    }

    LEAVE(io);
    ENTER(indexing);

    fastpm_store_histogram_aemit_sorted(lcevent->p, data->hist, data->aedges, data->Nedges, solver->comm);

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

    char * fn = fastpm_strdup_printf(CONF(prr->lua, lc_write_usmesh));

    write_usmesh_fof(solver, lcevent->p, fn, prr, lcevent, tail, mesh->lc);

    /* subsample, this will remove the tail particles that were appended. */
    fastpm_store_subsample(lcevent->p, lcevent->p->mask, lcevent->p);

    ENTER(sort);
    fastpm_sort_snapshot(lcevent->p, solver->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    if(lcevent->is_first) {
        fastpm_info("Creating usmesh catalog in %s\n", fn);
        write_snapshot(solver, lcevent->p, fn, "1", prr->cli->Nwriters);
        _write_parameters(fn, "Header", prr, solver->comm);
    } else {
        fastpm_info("Appending usmesh catalog to %s\n", fn);
        append_snapshot(solver, lcevent->p, fn, "1", prr->cli->Nwriters);
    }

    LEAVE(io);
    ENTER(indexing);

    fastpm_store_histogram_aemit_sorted(lcevent->p, data->hist, data->aedges, data->Nedges, solver->comm);

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
    int nout = CONF(prr->lua, n_aout);
    double * aout= CONF(prr->lua, aout);

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
                p->attributes & ~COLUMN_ACC,
                FASTPM_MEMORY_FLOATING
            );

        fastpm_info("Setting up snapshot at a = %6.4f (z=%6.4f)\n", aout[iout], 1.0f/aout[iout]-1);
        fastpm_info("Growth factor of snapshot %6.4f (a=%0.4f)\n", fastpm_solver_growth_factor(fastpm, aout[iout]), aout[iout]);

        fastpm_set_snapshot(fastpm, event->drift, event->kick, snapshot, aout[iout]);

        if(CONF(prr->lua, write_fof)) {
            char filebase[1024];

            sprintf(filebase, "%s_%0.04f", CONF(prr->lua, write_fof), aout[iout]);

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
        .linkinglength = CONF(prr->lua, fof_linkinglength) * CONF(prr->lua, boxsize) / CONF(prr->lua, nc),
        .periodic = 1,
        .nmin = CONF(prr->lua, fof_nmin),
        .kdtree_thresh = CONF(prr->lua, fof_kdtree_thresh),
    };

    fastpm_fof_init(&fof, snapshot, fastpm->basepm);

    FastPMStore halos[1];

    ENTER(fof);

    fastpm_fof_execute(&fof, halos);

    LEAVE(fof);

    ENTER(sort);
    fastpm_sort_snapshot(halos, fastpm->comm, FastPMSnapshotSortByLength, 0);
    LEAVE(sort);

    ENTER(io);
    char * dataset = fastpm_strdup_printf("LL-%05.3f", CONF(prr->lua, fof_linkinglength));
    write_snapshot(fastpm, halos, filebase, dataset, prr->cli->Nwriters);
    _write_parameters(filebase, "Header", prr, fastpm->comm);

    free(dataset);
    LEAVE(io);

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

    uint8_t * halos_established = malloc(halos->np);
    fastpm_info("halos_ready sees %td halos\n", halos->np);
    for(i = 0; i < halos->np; i ++) {
        double r = fastpm_lc_distance(lc, halos->x[i]);

        /* this shall not happen */
        if (halos->length[i] < finder->nmin) abort();

        /* only keep reliable halos */
        if(r > rmin + halosize * 0.5) {
            halos_established[i] = 1;
        } else {
            halos_established[i] = 0;
        }
        halos->mask[i] &= halos_established[i];
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
                /* particles in unreliable halos, keep them */
                keep_for_tail[i] = !halos_established[hid];
            } else {
                /* not in halos, keep them too */
                keep_for_tail[i] = 1;
            }
        }
    }
    free(halos_established);
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
    double maxhalosize = CONF(prr->lua, lc_usmesh_fof_padding); /* MPC/h, used to cut along z direction. */
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
        .linkinglength = CONF(prr->lua, fof_linkinglength) * CONF(prr->lua, boxsize) / CONF(prr->lua, nc),
        .periodic = 0,
        .nmin = CONF(prr->lua, fof_nmin),
        .kdtree_thresh = CONF(prr->lua, fof_kdtree_thresh),
    };

    fastpm_fof_init(&fof, snapshot, fastpm->basepm);

    double rmin = lc->speedfactor * HorizonDistance(lcevent->a1, lc->horizon);

    void * userdata[5];
    userdata[0] = & rmin;
    userdata[1] = & maxhalosize;
    userdata[2] = lc;
    userdata[3] = keep_for_tail;

    FastPMStore halos[1];
    fastpm_add_event_handler(&fof.event_handlers,
        FASTPM_EVENT_HALO,
        FASTPM_EVENT_STAGE_AFTER,
        (FastPMEventHandlerFunction) _halos_ready,
        userdata);

    ENTER(fof);

    fastpm_fof_execute(&fof, halos);

    LEAVE(fof);

    ENTER(sort);
    fastpm_sort_snapshot(halos, fastpm->comm, FastPMSnapshotSortByAEmit, 1);
    LEAVE(sort);

    ENTER(io);
    char * dataset = fastpm_strdup_printf("LL-%05.3f", CONF(prr->lua, fof_linkinglength));
    if(!append) {
        write_snapshot(fastpm, halos, filebase, dataset, prr->cli->Nwriters);
        _write_parameters(filebase, "Header", prr, fastpm->comm);
    } else {
        append_snapshot(fastpm, halos, filebase, dataset, prr->cli->Nwriters);
    }
    free(dataset);
    LEAVE(io);

    fastpm_store_destroy(halos);

    fastpm_fof_destroy(&fof);

    uint64_t ntail = 0;
    for(i = 0; i < p->np; i ++) {
        if(keep_for_tail[i]) ntail ++;
    }

    fastpm_store_init(tail, ntail, p->attributes, FASTPM_MEMORY_FLOATING);
    fastpm_store_subsample(p, keep_for_tail, tail);

    MPI_Allreduce(MPI_IN_PLACE, &ntail, 1, MPI_LONG, MPI_SUM, fastpm->comm);

    fastpm_info("%td particles will be reused in next batch for usmesh FOF\n", ntail);

    fastpm_memory_free(p->mem, keep_for_tail);
}

static int 
take_a_snapshot(FastPMSolver * fastpm, FastPMStore * snapshot, double aout, Parameters * prr) 
{
    CLOCK(io);
    CLOCK(sort);

    /* do this before write_snapshot, because white_snapshot messes up with the domain decomposition. */
    if(CONF(prr->lua, write_nonlineark)) {
        char * filename = fastpm_strdup_printf("%s_%0.04f", CONF(prr->lua, write_nonlineark), aout);
        FastPMPainter painter[1];

        FastPMFloat * rho_x = pm_alloc(fastpm->basepm);
        FastPMFloat * rho_k = pm_alloc(fastpm->basepm);

        fastpm_painter_init(painter, fastpm->basepm, fastpm->config->PAINTER_TYPE, fastpm->config->painter_support);

        fastpm_paint(painter, rho_x, snapshot, FASTPM_FIELD_DESCR_NONE);
        pm_r2c(fastpm->basepm, rho_x, rho_k);

        write_complex(fastpm->basepm, rho_k, filename, "DensityK", prr->cli->Nwriters);

        pm_free(fastpm->basepm, rho_k);
        pm_free(fastpm->basepm, rho_x);
        free(filename);
    }

    if(CONF(prr->lua, write_snapshot)) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;
        sprintf(filebase, "%s_%0.04f", CONF(prr->lua, write_snapshot), aout);

        fastpm_info("Writing snapshot %s at z = %6.4f a = %6.4f with %d writers\n", 
                filebase, z_out, aout, prr->cli->Nwriters);

        ENTER(sort);
        if(CONF(prr->lua, sort_snapshot)) {
            fastpm_info("Snapshot is sorted by ID.\n");
            fastpm_sort_snapshot(snapshot, fastpm->comm, FastPMSnapshotSortByID, 1);
        } else {
            fastpm_info("Snapshot is not sorted by ID.\n");
        }
        LEAVE(sort);
        ENTER(io);
        write_snapshot(fastpm, snapshot, filebase, "1", prr->cli->Nwriters);
        _write_parameters(filebase, "Header", prr, fastpm->comm);
        LEAVE(io);

        fastpm_info("snapshot %s written\n", filebase);
    }
    if(CONF(prr->lua, write_runpb_snapshot)) {
        char filebase[1024];
        double z_out= 1.0/aout - 1.0;

        sprintf(filebase, "%s_%0.04f.bin", CONF(prr->lua, write_runpb_snapshot), aout);

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
    report_memory(fastpm->comm);
    return 0;
}

static char PEAK_ALLOC_STATUS[8192];
static size_t PEAK_ALLOC_BYTES;
static void
_memory_peak_handler(FastPMMemory * mem, void * userdata)
{
    PEAK_ALLOC_BYTES = mem->peak_bytes;

    fastpm_memory_dump_status_str(mem, PEAK_ALLOC_STATUS, 8192);
}

static void
report_memory(MPI_Comm comm)
{
    static double oldpeak = 0;
    static int oldrank = -1;

    double max_used_bytes;
    double min_used_bytes;
    int max_rank;
    int min_rank;

    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);


    MPIU_stats(comm, PEAK_ALLOC_BYTES, "<,>.",
            &min_used_bytes,
            &min_rank,
            &max_used_bytes,
            &max_rank);

    if(max_used_bytes != oldpeak || max_rank != oldrank) {
        if(max_rank == ThisTask) {
            fastpm_ilog(INFO, "Task %d Peak memory usage max: %g MB min : %g MB\n",
                ThisTask,
                max_used_bytes / 1024. / 1024,
                min_used_bytes / 1024. / 1024
            );
            fastpm_ilog(INFO, PEAK_ALLOC_STATUS);
        }
    }
    oldpeak = max_used_bytes;
    oldrank = max_rank;
}

static int
write_powerspectrum(FastPMSolver * fastpm, FastPMForceEvent * event, Parameters * prr) 
{

    int K_LINEAR = CONF(prr->lua, enforce_broadband_kmax);

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
    if(CONF(prr->lua, write_powerspectrum)) {
        char buf[1024];
        sprintf(buf, "%s_%0.04f.txt", CONF(prr->lua, write_powerspectrum), event->a_f);
        fastpm_info("writing power spectrum to %s\n", buf);
        if(fastpm->ThisTask == 0) {
            fastpm_path_ensure_dirname(CONF(prr->lua, write_powerspectrum));
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


