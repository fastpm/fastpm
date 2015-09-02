#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "pmpfft.h"

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void rungdb(const char* fmt, ...);

static MPI_Datatype MPI_PTRDIFF = NULL;


static void module_init() {
    if(MPI_PTRDIFF != NULL) return;
        
    if(sizeof(ptrdiff_t) == 8) {
        MPI_PTRDIFF = MPI_LONG;
    } else {
        MPI_PTRDIFF = MPI_INT;
    }
}

void pm_pfft_init(PM * pm, PMInit * init, PMIFace * iface, MPI_Comm comm) {

    module_init();

    pm->init = *init;
    pm->iface = *iface;

    /* initialize the domain */    
    MPI_Comm_rank(comm, &pm->ThisTask);
    MPI_Comm_size(comm, &pm->NTask);

    int Nx = 1;
    int Ny = pm->NTask;
    for(; Nx * Nx < pm->NTask; Nx ++) continue;
    for(; Nx >= 1; Nx--) {
        if (pm->NTask % Nx == 0) break;
        continue;
    }

    Ny = pm->NTask / Nx;
    pm->Nproc[0] = Nx;
    pm->Nproc[1] = Ny;

    pm->Nmesh[0] = init->Nmesh;
    pm->Nmesh[1] = init->Nmesh;
    pm->Nmesh[2] = init->Nmesh;

    pm->BoxSize[0] = init->BoxSize;
    pm->BoxSize[1] = init->BoxSize;
    pm->BoxSize[2] = init->BoxSize;

    pm->Below[0] = 0;
    pm->Below[1] = 0;
    pm->Below[2] = 0;

    pm->Above[0] = 1;
    pm->Above[1] = 1;
    pm->Above[2] = 1;

    pfft_create_procmesh(2, comm, pm->Nproc, &pm->Comm2D);
    pm->allocsize = 2 * pfft_local_size_dft_r2c(
                3, pm->Nmesh, pm->Comm2D, PFFT_TRANSPOSED_OUT, 
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);

    int d;
    for(d = 0; d < 2; d ++) {
        MPI_Comm projected;
        int remain_dims[2] = {0, 0};
        remain_dims[d] = 1; 

        pm->Grid.edges_int[d] = 
            malloc(sizeof(pm->Grid.edges_int[0][0]) * (pm->Nproc[d] + 1));
        pm->Grid.edges_float[d] = 
            malloc(sizeof(pm->Grid.edges_float[0][0]) * (pm->Nproc[d] + 1));

        pm->Grid.MeshtoCart[d] = malloc(sizeof(int) * pm->Nmesh[d]);

        MPI_Cart_sub(pm->Comm2D, remain_dims, &projected);
        MPI_Allgather(&pm->IRegion.start[d], 1, MPI_PTRDIFF, 
            pm->Grid.edges_int[d], 1, MPI_PTRDIFF, projected);
        int ntask;
        MPI_Comm_size(projected, &ntask);

        MPI_Comm_free(&projected);
        int j;
        for(j = 0; j < pm->Nproc[d]; j ++) {
            pm->Grid.edges_float[d][j] = 1.0 * pm->Grid.edges_int[d][j] / pm->Nmesh[d] * pm->BoxSize[d];
        }
        /* Last edge is at the edge of the box */
        pm->Grid.edges_float[d][j] = pm->BoxSize[d];
        pm->Grid.edges_int[d][j] = pm->Nmesh[d];
        /* fill in the look up table */
        for(j = 0; j < pm->Nproc[d]; j ++) {
            int i;
            for(i = pm->Grid.edges_int[d][j]; i < pm->Grid.edges_int[d][j+1]; i ++) {
                pm->Grid.MeshtoCart[d][i] = j;
            }
        }
    }

    void * buffer = pm->iface.malloc(pm->allocsize * sizeof(double));

    pm->r2c = pfft_plan_dft_r2c(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_FORWARD, PFFT_TRANSPOSED_OUT | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);

    pm->c2r = pfft_plan_dft_c2r(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_BACKWARD, PFFT_TRANSPOSED_IN | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);
    pm->iface.free(buffer);

    for(d = 0; d < 3; d++) {
        pm->MeshtoK[d] = malloc(pm->Nmesh[d] * sizeof(double));
        int i;
        for(i = 0; i < pm->Nmesh[d]; i++) {
            int ii = i;
            if(ii > pm->Nmesh[d] / 2) {
                ii -= pm->Nmesh[d];
            }
            pm->MeshtoK[d][i] = i * 2 * M_PI / pm->BoxSize[d];
        }
    }
}

int pm_pos_to_rank(PM * pm, double pos[3]) {
    int d;
    int rank2d[2];
    for(d = 0; d < 2; d ++) {
        int ipos = floor(pos[d] / pm->BoxSize[d] * pm->Nmesh[d]);
        while(ipos < 0) ipos += pm->Nmesh[d];
        while(ipos >= pm->Nmesh[d]) ipos -= pm->Nmesh[d];
        rank2d[d] = pm->Grid.MeshtoCart[d][ipos];
    }
    return rank2d[0] * pm->Nproc[1] + rank2d[1];
}
void pm_start(PM * pm) {
    pm->canvas = pm->iface.malloc(sizeof(double) * pm->allocsize);
    pm->workspace = pm->iface.malloc(sizeof(double) * pm->allocsize);
}
void pm_stop(PM * pm) {
    pm->iface.free(pm->canvas);
    pm->iface.free(pm->workspace);
    pm->canvas = NULL;
    pm->workspace = NULL;
}

void pm_r2c(PM * pm) {
    pfft_execute_dft_r2c(pm->r2c, pm->canvas, (pfft_complex*)pm->canvas);
}

void pm_c2r(PM * pm) {
    pfft_execute_dft_c2r(pm->c2r, (pfft_complex*) pm->workspace, pm->workspace);
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    pfft_init();

    PMStore pdata;
    memset(&pdata, 0, sizeof(pdata));
    pm_store_init(&pdata, 100);

    int i;
    for(i = 0; i < 100; i ++) {
        pdata.x[i][0] = i + 0.2;
        pdata.x[i][1] = i;
        pdata.x[i][2] = i;
        pdata.acc[0][i] = 1;
    }
    
    PMInit pminit = {
        .Nmesh = 4,
        .BoxSize = 4.,
        .AllAttributes = 0xffff,
        .GhostAttributes = PACK_POS,
    };
    PM pm;
    PMGhostData pgd;

    pm_pfft_init(&pm, &pminit, &pdata.iface, MPI_COMM_WORLD);
    pm_ghost_data_init(&pm, &pdata, 4, &pgd);
    pm_append_ghosts(&pm, 100, &pgd);
    
    pm_reduce_ghosts(&pm, &pgd, PACK_ACC_X); 
    pm_reduce_ghosts(&pm, &pgd, PACK_ACC_Y); 
    pm_reduce_ghosts(&pm, &pgd, PACK_ACC_Z); 

    pm_start(&pm);
    pm_paint(&pm, &pdata, 4);
    rungdb("p ((double*)%p)[0]@32", &pm.canvas[0]);
    pm_r2c(&pm);
    memcpy(pm.workspace, pm.canvas, pm.allocsize*sizeof(double));
    pm_c2r(&pm);
    pm_readout_one(&pm, &pdata, 0);

    rungdb("p ((double*)%p)[0]@32", &pm.workspace[0]);
    pm_stop(&pm);
    //rungdb("p *((PMGhostData*)%p)", &pgd);
    //rungdb("p ((double*) %p)[0]@12", pdata.pos);
    //rungdb("p ((double*) %p)[0]@12", pdata.acc);
    //
    //rungdb("p ((PM*)%p)->Grid.edges_float[0][0]@%d", &pm, pm.Nproc[0]+1);
    //rungdb("p ((PM*)%p)->Grid.edges_float[1][0]@%d", &pm, pm.Nproc[1]+1);


    pfft_cleanup();
    MPI_Finalize();
}


static void rungdb(const char* fmt, ...){
    /* dumpstack(void) Got this routine from http://www.whitefang.com/unix/faq_toc.html
 *     ** Section 6.5. Modified to redirect to file to prevent clutter
 *         */
    /* This needs to be changed... */
    char dbx[160];
    char cmd[160];
    char * tmpfilename;
    extern const char *__progname;
    va_list va;
    va_start(va, fmt);
    
    vsprintf(cmd, fmt, va);
    va_end(va);

    tmpfilename = tempnam(NULL, NULL);

    sprintf(dbx, "echo '%s\n' > %s", cmd, tmpfilename);
    system(dbx);

    sprintf(dbx, "echo 'where\ndetach' | gdb -batch --command=%s %s %d", tmpfilename, __progname, getpid() );
    system(dbx);
    unlink(tmpfilename);
    free(tmpfilename);

    return;
}

