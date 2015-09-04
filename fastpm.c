#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "pmpfft.h"

#define MAX(a, b) (a)>(b)?(a):(b)
#define BREAKPOINT raise(SIGTRAP);

static void rungdb(const char* fmt, ...);
static int to_rank(PMStore * pdata, ptrdiff_t i, void * data) {
    PM * pm = (PM*) data;
    double pos[3];
    pdata->iface.get_position(pdata, i, pos);
    return pm_pos_to_rank(pm, pos);
}
int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    pfft_init();

    PMStore pdata;
    memset(&pdata, 0, sizeof(pdata));
    pm_store_init(&pdata, 100);

    int i;
    for(i = 0; i < 100; i ++) {
        pdata.x[i][1] = fmod(i + 0.5, 4);
        pdata.x[i][0] = 0;
        pdata.x[i][2] = 0;
        pdata.acc[0][i] = 1;
    }
    pdata.np = 4;
 
    PMInit pminit = {
        .Nmesh = 4,
        .BoxSize = 4.,
        .NprocX = 0, /* 0 for auto, 1 for slabs */
    };
    PM pm;

    pm_pfft_init(&pm, &pminit, &pdata.iface, MPI_COMM_WORLD);
    pm_store_wrap(&pdata, pm.BoxSize);
    pm_store_decompose(&pdata, to_rank, &pm, MPI_COMM_WORLD);

    PMGhostData pgd = {
        .pm = &pm,
        .pdata = &pdata,
        .np = pdata.np,
        .np_upper = pdata.np_upper,
        .attributes = PACK_POS,
    };
    pm_append_ghosts(&pgd);

    pm_start(&pm);
    pm_paint(&pm, &pdata, pdata.np + pgd.nghosts);
    pm_r2c(&pm);
    memcpy(pm.workspace, pm.canvas, pm.allocsize*sizeof(double));
    pm_c2r(&pm);
    pm_readout_one(&pm, &pdata, 0);

    pm_reduce_ghosts(&pgd, PACK_ACC_X); 
    pm_reduce_ghosts(&pgd, PACK_ACC_Y); 
    pm_reduce_ghosts(&pgd, PACK_ACC_Z); 

    sleep(pm.ThisTask * 2);
    printf("----rank = %d\n", pm.ThisTask);
    rungdb("p ((double*)%p)[0]@32", &pm.workspace[0]);
    rungdb("p *((PMGhostData*)%p)", &pgd);
    rungdb("p ((double*) %p)[0]@12", pdata.x);
    //rungdb("p ((double*) %p)[0]@12", pdata.acc);
    //
    rungdb("p ((PM*)%p)->Grid.edges_float[0][0]@%d", &pm, pm.Nproc[0]+1);
    rungdb("p ((PM*)%p)->Grid.edges_float[1][0]@%d", &pm, pm.Nproc[1]+1);


    pm_stop(&pm);
    pm_destroy_ghosts(&pgd);
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

