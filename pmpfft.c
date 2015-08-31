#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <pfft.h>

#define MAX(a, b) (a)>(b)?(a):(b)

static void rungdb(const char* fmt, ...);

typedef struct {
    void * (*malloc )(size_t);
    void   (*free   )(void *);
    MPI_Datatype PType;
    void   (*get_position)(void * pdata, ptrdiff_t * index, double pos[3]);
    ptrdiff_t Nmesh;
    double BoxSize;
} PMInit;

typedef struct {
    ptrdiff_t start[3];
    ptrdiff_t size[3];
} PMRegion;

static MPI_Datatype MPI_PTRDIFF = NULL;

typedef struct {
    ptrdiff_t * edges_int[2];
    double * edges_float[2];
} PMGrid;

typedef struct {
    PMInit init;
    int NTask;
    int ThisTask;
    pfft_plan r2c;
    pfft_plan c2r;

    int Nproc[2];
    MPI_Comm Comm2D;

    ptrdiff_t Nmesh[3];
    double    BoxSize[3];
    ptrdiff_t GhostSize[3]; 

    ptrdiff_t allocsize;

    PMRegion IRegion;
    PMRegion IRegionGhost;
    PMRegion ORegion;
    
    PMGrid Grid;
} PM;

static void module_init() {
    if(MPI_PTRDIFF != NULL) return;
        
    if(sizeof(ptrdiff_t) == 8) {
        MPI_PTRDIFF = MPI_LONG;
    } else {
        MPI_PTRDIFF = MPI_INT;
    }
}
void pm_pfft_init(PM * pm, PMInit * init, MPI_Comm comm) {

    module_init();

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

    pm->GhostSize[0] = 1;
    pm->GhostSize[1] = 1;
    pm->GhostSize[2] = 1;

    pfft_create_procmesh(2, comm, pm->Nproc, &pm->Comm2D);

    pm->allocsize = 2 * pfft_local_size_dft_r2c(
                3, pm->Nmesh, comm, PFFT_TRANSPOSED_OUT, 
                pm->IRegion.size, pm->IRegion.start,
                pm->ORegion.size, pm->ORegion.start);
    
    pm->allocsize = MAX(pm->allocsize, 
         pfft_local_size_gc(3, 
             pm->IRegion.size, pm->IRegion.start, 
             pm->GhostSize, pm->GhostSize,
             pm->IRegionGhost.size, pm->IRegionGhost.start)
         );

    int d;
    for(d = 0; d < 2; d ++) {
        MPI_Comm projected;
        int remain_dims[2] = {0};
        remain_dims[d] = 1; 

        pm->Grid.edges_int[d] = 
            malloc(sizeof(pm->Grid.edges_int[0][0]) * (pm->Nproc[d] + 1));
        pm->Grid.edges_float[d] = 
            malloc(sizeof(pm->Grid.edges_float[0][0]) * (pm->Nproc[d] + 1));

        MPI_Cart_sub(pm->Comm2D, remain_dims, &projected);
        MPI_Allgather(&pm->IRegion.start[d], 1, MPI_PTRDIFF, 
            pm->Grid.edges_int[d], 1, MPI_PTRDIFF, projected);
        MPI_Comm_free(&projected);
        int j;
        for(j = 0; j < pm->Nproc[d]; j ++) {
            pm->Grid.edges_float[d][j] = 1.0 * pm->Grid.edges_int[d][j] / pm->Nmesh[d] * pm->BoxSize[d];
        }
        /* Last edge is at the edge of the box */
        pm->Grid.edges_float[d][j] = pm->BoxSize[d];
        pm->Grid.edges_int[d][j] = pm->Nmesh[d];
    }
    pm->init = *init;
    void * buffer = pm->init.malloc(pm->allocsize * sizeof(double));

    pm->r2c = pfft_plan_dft_r2c(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_FORWARD, PFFT_TRANSPOSED_OUT | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);

    pm->c2r = pfft_plan_dft_c2r(
            3, pm->Nmesh, buffer, buffer, 
            pm->Comm2D,
            PFFT_BACKWARD, PFFT_TRANSPOSED_IN | PFFT_ESTIMATE | PFFT_DESTROY_INPUT);

    pm->init.free(buffer);
}

static int pm_pos_to_rank(PM * pm, double pos[3]) {
    int d;
    int rank2d[2];
    static int lastrank2d[2] = {0, 0};
    for(d = 0; d < 2; d ++) {
        int ipos = ceil(pos[d] / pm->BoxSize[d] * pm->Nmesh[d]);
        while(ipos < 0) ipos += pm->Nmesh[d];
        while(ipos >= pm->Nmesh[d]) ipos -= pm->Nmesh[d];
        int a = pm->Grid.edges_int[d][lastrank2d[d]];
        int b = pm->Grid.edges_int[d][lastrank2d[d] + 1];
        if(a <= ipos && b > ipos) {
            rank2d[d] = lastrank2d[d];
        } else {
            int j = 0;
            /* FIXME: use bsearch and add a hint ! */
            while(ipos > pm->Grid.edges_int[d][j])
                j++;
            rank2d[d] = j - 1;
            lastrank2d[d] = rank2d[d];
        }
    }
    return rank2d[0] * pm->Nproc[1] + rank2d[1];
};

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    pfft_init();
    PMInit pminit = {
        .malloc = malloc,
        .free = free,
        .Nmesh = 128,
        .BoxSize = 32.,
    };
    PM pm;
    pm_pfft_init(&pm, &pminit, MPI_COMM_WORLD);

    //rungdb("p ((PM*)%p)->Grid.edges_float[0][0]@%d", &pm, pm.Nproc[0]+1);
    //rungdb("p ((PM*)%p)->Grid.edges_float[1][0]@%d", &pm, pm.Nproc[1]+1);

    double pos[100][3];
    double vel[100][3];

    {
        double pos[3];
        int i;
        for (i = 0; i < 32; i ++) {
            pos[0] = i;
            pos[1] = i;
            pos[2] = i;
            printf("---- %d ----\n", pm_pos_to_rank(&pm, pos));
        }
    }
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

