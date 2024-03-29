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
#include <fastpm/io.h>
#include <fastpm/fof.h>
#include <fastpm/rfof.h>
#include <bigfile.h>
#include <bigfile-mpi.h>

#include "lua-config.h"
#include "param.h"
#include "prepare.h"

extern void
init_stacktrace();

LUAParameters * read_lua_parameters_mpi(const char * filebase, char ** error, MPI_Comm comm)
{
    *error = NULL;
    LUAParameters * prr = malloc(sizeof(prr[0]));

    BigFile bf[1];
    BigBlock bb[1];
    BigAttr * attr;

    big_file_mpi_open(bf, filebase, comm);

    big_file_mpi_open_block(bf, bb, "Header", comm);

    attr = big_block_lookup_attr(bb, "ParamFile");

    char * confstr = malloc(attr->nmemb + 1);
    big_block_get_attr(bb, "ParamFile", confstr, "S1", attr->nmemb);

    prr->config = lua_config_new(confstr);

    if(lua_config_error(prr->config)) {
        *error = fastpm_strdup_printf("%s", lua_config_error(prr->config));
        free(prr);
        return NULL;
    }

    prr->string = confstr;

    big_block_mpi_close(bb, comm);
    big_file_mpi_close(bf, comm);
    return prr;
}

int
main(int argc, char * argv[])
{
    init_stacktrace();

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);
    fastpm_info("This is FastPM, with libfastpm version %s.\n", LIBFASTPM_VERSION);

    char * error;

    CLIParameters * cli = parse_cli_args_mpi(argc, argv, comm);
    if(cli->argc < 2) {
        fastpm_raise(-1, "Must supply a snapshot file name and a linking length scaling factor.\n");
    }
    char * filebase = cli->argv[0];
    double scale = atof(cli->argv[1]);
    char * dataset = fastpm_strdup_printf("RFOF-%05.3f", scale);

    fastpm_info("Running RFOF on %s; writing to %s\n", filebase, dataset);
    LUAParameters * lua = read_lua_parameters_mpi(filebase, &error, comm);

    if(lua) {
        fastpm_info("Configuration %s\n", lua->string);
    } else {
        fastpm_info("Parsing configuration failed with error: %s\n", error);
        MPI_Finalize();
        exit(1);
    }

    libfastpm_set_memory_bound(cli->MemoryPerRank * 1024 * 1024);

    FastPMStore source[1];

    /* load the CDM species from the snapshot */
    fastpm_store_init_evenly(source,
            fastpm_species_get_name(FASTPM_SPECIES_CDM),
            pow(1.0 * CONF(lua, nc), 3),
            COLUMN_POS | COLUMN_VEL | COLUMN_ID,
            CONF(lua, np_alloc_factor), comm);

    PM * basepm = fastpm_create_pm(CONF(lua, nc), cli->NprocY, 1, CONF(lua, boxsize), comm);

    fastpm_store_write(source, filebase, "r", cli->Nwriters, comm);

    /* convert from fraction of mean separation to simulation distance units. */
    double sep = CONF(lua, boxsize) / CONF(lua, nc);
    double z = 1. / source->meta.a_x - 1.;

    CLOCK(fof);
    CLOCK(io);
    CLOCK(sort);
    ENTER(fof);
    FastPMRFOFFinder rfof = {
        .periodic = 1,
        .nmin = CONF(lua, rfof_nmin),
        .kdtree_thresh = CONF(lua, rfof_kdtree_thresh),
        .linkinglength = CONF(lua, rfof_linkinglength) * sep * scale,
        .l1 = CONF(lua, rfof_l1) * sep * scale,
        .l6 = CONF(lua, rfof_l6) * sep * scale,
        .A1 = CONF(lua, rfof_a1) * sep * scale,
        .A2 = CONF(lua, rfof_a2) * sep * scale,
        .B1 = CONF(lua, rfof_b1) * sep * scale,
        .B2 = CONF(lua, rfof_b2) * sep * scale,
    };

    FastPMCosmology cosmology[1] = {{0}};

    prepare_cosmology(cosmology, lua);
    fastpm_cosmology_init(cosmology);
    fastpm_rfof_init(&rfof, cosmology, source, basepm);

    FastPMStore halos[1];

    ENTER(fof);

    fastpm_store_set_name(halos, dataset);
    ptrdiff_t * ihalo = fastpm_rfof_execute(&rfof, halos, z);
    fastpm_memory_free(halos->mem, ihalo);
    fastpm_store_subsample(halos, halos->mask, halos);
    LEAVE(fof);

    ENTER(sort);
    fastpm_sort_snapshot(halos, comm, FastPMSnapshotSortByLength, 0);
    LEAVE(sort);

    ENTER(io);
    fastpm_store_write(halos, filebase, "w", cli->Nwriters, comm);

    LEAVE(io);

    fastpm_store_destroy(halos);
    fastpm_rfof_destroy(&rfof);
    free_lua_parameters(lua);
    free_cli_parameters(cli);

    fastpm_store_destroy(source);
    fastpm_free_pm(basepm);

    free(dataset);
    fastpm_clock_stat(comm);

    libfastpm_cleanup();

    MPI_Finalize();

    return 0;

}
