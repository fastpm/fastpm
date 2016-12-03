#include <string.h>
#include <stdlib.h>

#include <mpi.h>
#include <pfft.h>
#include <math.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmesh/_window_imp.h"

#if FASTPM_FFT_PRECISION == 64
#define _pfft_plan pfft_plan
#define _pfft_create_procmesh pfft_create_procmesh
#define _pfft_local_size_dft_r2c pfft_local_size_dft_r2c
#define _pfft_plan_dft_r2c pfft_plan_dft_r2c
#define _pfft_plan_dft_c2r pfft_plan_dft_c2r
#define _pfft_execute_dft_r2c pfft_execute_dft_r2c
#define _pfft_execute_dft_c2r pfft_execute_dft_c2r
#define _pfft_init pfft_init
#define _pfft_cleanup pfft_cleanup
#define _pfft_destroy_plan pfft_destroy_plan
#else
#define _pfft_plan pfftf_plan
#define _pfft_create_procmesh pfftf_create_procmesh
#define _pfft_local_size_dft_r2c pfftf_local_size_dft_r2c
#define _pfft_plan_dft_r2c pfftf_plan_dft_r2c
#define _pfft_plan_dft_c2r pfftf_plan_dft_c2r
#define _pfft_execute_dft_r2c pfftf_execute_dft_r2c
#define _pfft_execute_dft_c2r pfftf_execute_dft_c2r
#define _pfft_init pfftf_init
#define _pfft_cleanup pfftf_cleanup
#define _pfft_destroy_plan pfftf_destroy_plan
#endif

struct FastPMMeshPrivate {
    MPI_Datatype ptrtype;
    MPI_Comm cart;
    int * Mesh2Rank[2];
    size_t allocsize;
    _pfft_plan plan_r2c;
    _pfft_plan plan_r2c_inplace;
    _pfft_plan plan_c2r;
    _pfft_plan plan_c2r_inplace;
    int proc_strides[3];
};

void
fastpm_mesh_init(FastPMMesh * self,
        int ndim,
        double BoxSize,
        ptrdiff_t Nmesh,
        ptrdiff_t Nproc[],
        MPI_Comm comm)
{
    self->mem = _libfastpm_get_gmem();
    self->priv = malloc(sizeof(self->priv[0]));
    FastPMMeshPrivate * priv = self->priv;
    if(sizeof(ptrdiff_t) == 8) {
        priv->ptrtype = MPI_LONG;
    } else {
        priv->ptrtype = MPI_INT;
    }

    self->ndim = ndim;
    self->ral.ndim = ndim;
    self->cal.ndim = ndim;
    int d = 0;
    self->Volume = 1;
    self->Norm = 1;
    for(d = 0; d < ndim; d ++) {
        self->BoxSize[d] = BoxSize;
        self->Nmesh[d] = Nmesh;
        self->CellSize[d] = BoxSize / Nmesh;
        self->InvCellSize[d] = Nmesh / BoxSize;
        self->Volume *= BoxSize;
        self->Norm *= Nmesh;
    }
    self->comm = comm;
    int NTask;
    MPI_Comm_size(comm, &NTask);

    if(Nproc == NULL) {
        if(self->ndim == 3) {
            /* 2d decomposition */
            int x = 1;
            for(x = 1; x * x < NTask; x ++) continue;
            while(x >= 1 && (NTask % x != 0)) {
                x --;
            }
            self->Nproc[1] = x;
            self->Nproc[0] = NTask / x;
        } else if (self->ndim == 2) {
            self->Nproc[0] = NTask;
        } else if (self->ndim == 1) {
            self->Nproc[0] = 0;
        } else {
            fastpm_raise(-1, "Only support at most 3 dimensions\n");
        }
    } else {
        for(d = 0; d < ndim - 1; d ++) {
            self->Nproc[d] = Nproc[d];
        }
    }
    self->NTask = 1;
    for(d = 0; d < ndim - 1; d ++) {
        self->NTask *= self->Nproc[d];
    }
    priv->proc_strides[ndim - 2] = 1;
    for(d = ndim - 3; d >=0; d--) {
        priv->proc_strides[d] = priv->proc_strides[d + 1] * self->Nproc[d + 1];
    }
    fastpm_info("proc_strides = %d %d\n", priv->proc_strides[0], priv->proc_strides[1]);
    if(NTask != self->NTask) {
        fastpm_raise(-1, "NTask and Nproc mismatch\n");
    }

    _pfft_create_procmesh(self->ndim - 1, comm, self->Nproc, &priv->cart);

    priv->allocsize = 2 * _pfft_local_size_dft_r2c(
            self->ndim, self->Nmesh, priv->cart, 
            PFFT_TRANSPOSED_OUT | PFFT_PADDED_R2C, 
            self->ral.shape, self->ral.start,
            self->cal.shape, self->cal.start);

    if(self->ndim == 3) {
        self->ral.strides[2] = 1;
        self->ral.strides[1] = self->ral.shape[2];
        self->ral.strides[0] = self->ral.shape[1] * self->ral.strides[1];
        self->ral.size = self->ral.shape[0] * self->ral.strides[0];

        /* remove padding from the view */
        self->ral.shape[2] = self->Nmesh[2];

        /* PFFT transposed, y, z, x */
        self->cal.strides[0] = 1;
        self->cal.strides[2] = self->cal.shape[0];
        self->cal.strides[1] = self->cal.shape[2] * self->cal.strides[2];
        self->cal.size = self->cal.shape[1] * self->cal.strides[1];
    } else {
        fastpm_raise(-1, "only three dimensions are supported for now\n");
    }

    for(d = 0; d < self->ndim - 1; d ++) {
        MPI_Comm projected;
        int remain_dims[100];
        int j;
        for(j = 0; j < self->ndim - 1; j ++) {
            remain_dims[j] = 0;
        }
        remain_dims[d] = 1; 
        priv->Mesh2Rank[d] = malloc(sizeof(priv->Mesh2Rank[0]) * self->Nmesh[d]);
        ptrdiff_t * edges_int = malloc(sizeof(edges_int[0]) * (self->Nproc[d] + 1));

        MPI_Cart_sub(priv->cart, remain_dims, &projected);
        MPI_Allgather(&self->ral.start[d], 1, priv->ptrtype, 
            edges_int, 1, priv->ptrtype, projected);
        int ntask;
        MPI_Comm_size(projected, &ntask);

        MPI_Comm_free(&projected);
        /* Last edge is at the edge of the box */
        edges_int[ntask] = self->Nmesh[d];
        /* fill in the look up table */
        for(j = 0; j < self->Nproc[d]; j ++) {
            int i;
            for(i = edges_int[j]; i < edges_int[j+1]; i ++) {
                priv->Mesh2Rank[d][i] = j;
            }
        }
        free(edges_int);
    }

    FastPMFloat * workspace = fastpm_mesh_alloc(self);
    FastPMFloat * canvas = fastpm_mesh_alloc(self);
    priv->plan_r2c = _pfft_plan_dft_r2c(
            self->ndim, self->Nmesh, (void*) workspace, (void*) canvas, 
            priv->cart,
            PFFT_FORWARD, 
            PFFT_TRANSPOSED_OUT
            | PFFT_PADDED_R2C 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_r2c_inplace = _pfft_plan_dft_r2c(
            self->ndim, self->Nmesh, (void*) workspace, (void*) workspace, 
            priv->cart,
            PFFT_FORWARD, 
            PFFT_TRANSPOSED_OUT
            | PFFT_PADDED_R2C 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_c2r = _pfft_plan_dft_c2r(
            self->ndim, self->Nmesh, (void*) workspace, (void*) canvas, 
            priv->cart,
            PFFT_BACKWARD, 
            PFFT_TRANSPOSED_IN
            | PFFT_PADDED_C2R 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    priv->plan_c2r_inplace = _pfft_plan_dft_c2r(
            self->ndim, self->Nmesh, (void*) workspace, (void*) workspace, 
            priv->cart,
            PFFT_BACKWARD, 
            PFFT_TRANSPOSED_IN
            | PFFT_PADDED_C2R 
            | PFFT_ESTIMATE 
            | PFFT_TUNE
            );

    fastpm_mesh_free(self, canvas);
    fastpm_mesh_free(self, workspace);
}

void
fastpm_mesh_destroy(FastPMMesh * self)
{
    int d;
    for(d = 0; d < self->ndim - 1; d ++) {
        free(self->priv->Mesh2Rank[d]);
    }
    MPI_Comm_free(&self->priv->cart);
    _pfft_destroy_plan(self->priv->plan_r2c);
    _pfft_destroy_plan(self->priv->plan_c2r);
    _pfft_destroy_plan(self->priv->plan_r2c_inplace);
    _pfft_destroy_plan(self->priv->plan_c2r_inplace);
    free(self->priv);
}

FastPMFloat *
fastpm_mesh_alloc(FastPMMesh * self)
{
    void * p = fastpm_memory_alloc(self->mem, sizeof(FastPMFloat) * self->priv->allocsize, FASTPM_MEMORY_HEAP);
    fastpm_memory_tag(self->mem, p, "FastPMFloat");
    return p;
}

void
fastpm_mesh_free(FastPMMesh * self, FastPMFloat * array)
{
    fastpm_memory_free(self->mem, array);
}

void
fastpm_mesh_r2c(FastPMMesh * self, FastPMFloat * r, FastPMFloat * c)
{
    if(c == r) {
        pfft_execute_dft_r2c(self->priv->plan_r2c_inplace, (void*)r, (void*)c);
    } else {
        pfft_execute_dft_r2c(self->priv->plan_r2c, (void*) r, (void*) c);
    }
}

void
fastpm_mesh_c2r(FastPMMesh * self, FastPMFloat * c, FastPMFloat * r)
{
    if(c == r) {
        pfft_execute_dft_c2r(self->priv->plan_c2r_inplace, (void*) r, (void*) c);
    } else {
        pfft_execute_dft_c2r(self->priv->plan_c2r, (void*) r, (void*) c);
    }
}

void
fastpm_mesh_copy(FastPMMesh * self, FastPMFloat * from, FastPMFloat * to)
{
    if(from == to) {
        return;
    } else {
        memcpy(to, from, sizeof(FastPMFloat) * self->priv->allocsize);
    }
}

int
fastpm_mesh_ipos_to_rank(FastPMMesh * self, int ipos[])
{
    int rank = 0;
    int d;
    for(d = 0; d < self->ndim - 1; d ++) {
        int i = ipos[d];
        while(i >= self->Nmesh[d]) i -= self->Nmesh[d];
        while(i < 0) i += self->Nmesh[d];
        rank += self->priv->Mesh2Rank[d][i] * self->priv->proc_strides[d];
    }
    return rank;
}

int
fastpm_mesh_pos_to_ipos(FastPMMesh * self, double pos, int dim)
{
    return (int) floor(self->InvCellSize[dim] * pos);
}
