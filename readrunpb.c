/**
 * IO of RunPB format
 * 
 * Authors: Yu Feng <rainwoodman@gmail.com>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>
#include <glob.h>
#include <math.h>
#include <alloca.h>
#include "parameters.h"
#include "power.h"
#include "particle.h"
#include "msg.h"
#include "heap.h"

#define FILENAME  "%s.%02d"

typedef struct {
  int   npart;          /* Total number of p. */
  int   nsph;           /* Number of gas p.   */
  int   nstar;          /* Number of star p.  */
  float aa;             /* Scale factor. */
  float eps;            /* Gravitational softening    */
} FileHeader;

int read_runpb_ic(Parameters * param, double a_init, Particles * p) {
    int ThisTask;
    int NTask;
    size_t scratch_bytes = ((size_t) 4) * param->nc * param->nc * param->nc;
    void * scratch = heap_allocate(scratch_bytes);
    float * fscratch = (float*) scratch;
    int64_t * lscratch = (int64_t *) scratch;

    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    size_t chunksize = scratch_bytes;

    size_t Ntot = 0;

    int * NperFile = NULL;
    size_t * NcumFile = NULL;
    int Nfile = 0;
    double aa = 0;

    if (ThisTask == 0) {
       char buf[1024];
        for(int i = 0; ; i ++) {
            sprintf(buf, FILENAME, param->readic_filename, i);
            FILE * fp = fopen(buf, "r");
            if(!fp) {
                Nfile = i;
                break;
            }
            fclose(fp);
        }
        if (Nfile == 0) {
            msg_abort(0030, "No snapshot files were found.\n");
        }
        printf("Total number of files is %d", Nfile);

        MPI_Bcast(&Nfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
        NperFile = malloc(sizeof(int) * Nfile);

        for(int i = 0; i < Nfile; i ++) {
            sprintf(buf, FILENAME, param->readic_filename, i);
            FILE * fp = fopen(buf, "r");
            FileHeader header;
            int eflag, hsize;
            fread(&eflag, sizeof(int), 1, fp);
            fread(&hsize, sizeof(int), 1, fp);
            if(hsize != sizeof(header)) {
                msg_abort(0030, "Unable to read from %s\n", buf);
            }
            fread(&header, sizeof(FileHeader), 1, fp);
            aa = header.aa;
            printf("reading from file %s", buf);
            printf(" npart=%d", header.npart);
            printf(" aa=%g \n", header.aa);
            NperFile[i] = header.npart;
            Ntot += header.npart;        
            fclose(fp);
        }
        if (Ntot != param->nc * param->nc * param->nc) {
            msg_abort(0030, "Number of p does not match nc\n");
        }
        MPI_Bcast(NperFile, Nfile, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Ntot, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&aa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        /* other ranks receive */
        MPI_Bcast(&Nfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
        NperFile = malloc(sizeof(int) * Nfile);
        MPI_Bcast(NperFile, Nfile, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Ntot, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&aa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    NcumFile = malloc(sizeof(size_t) * Nfile);
    NcumFile[0] = 0;
    for(int i = 1; i < Nfile; i ++) {
        NcumFile[i] = NcumFile[i - 1] + NperFile[i - 1];
    }

    size_t start = ThisTask * Ntot / NTask;
    size_t end   = (ThisTask + 1) * Ntot / NTask;
    p->np_local = end - start;

    
    int offset = 0;
    int chunknpart = chunksize / (sizeof(float) * 3);
    msg_printf(verbose, "chunknpart = %d\n", chunknpart);
    for(int i = 0; i < Nfile; i ++) {
        ptrdiff_t mystart = start - NcumFile[i];
        ptrdiff_t myend = end - NcumFile[i];
        size_t nread;
        /* not overlapping */
        if(myend <= 0) continue;
        if(mystart >= NperFile[i]) continue;

        printf("Task %d reading at %d \n", ThisTask, offset);

        /* cut to this file */
        if(myend > NperFile[i]) {
            myend = NperFile[i];
        }
        if(mystart < 0) {
            mystart =0;
        }
        /* at this point, read mystart:myend from current file */
        char buf[1024];
        int eflag, hsize;
        FileHeader header;
        sprintf(buf, FILENAME, param->readic_filename, i);

        FILE * fp = fopen(buf, "r");
        /* skip these */
        fread(&eflag, sizeof(int), 1, fp);
        fread(&hsize, sizeof(int), 1, fp);
        fread(&header, sizeof(FileHeader), 1, fp);
        /* pos */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            fread(scratch, sizeof(float) * 3, nbatch, fp);
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                for(int d = 0; d < 3; d ++)
                    p->x[offset + nread + ip][d] = fscratch[q++];
            }
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* vel */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            fread(scratch, sizeof(float) * 3, nbatch, fp);
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                for(int d = 0; d < 3; d ++)
                    p->v[offset + nread + ip][d] = fscratch[q++];
            }
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* ID */
        fseek(fp, mystart * sizeof(int64_t), SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            fread(scratch, sizeof(int64_t), nbatch, fp);
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                p->id[offset + nread + ip] = lscratch[q++];
            }
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(int64_t), SEEK_CUR);
        fclose(fp);
        offset += myend - mystart;
    }
    if(offset != p->np_local) {
        msg_abort(0030, "mismatch %d != %d\n", offset, p->np_local);
    }
    const double Omega_m = param->omega_m;
    const double omega=Omega_m/(Omega_m + (1.0 - Omega_m)*aa*aa*aa);
    const float DplusIC = 1.0/GrowthFactor(aa, 1.0);
    const float Dplus = 1.0/GrowthFactor(a_init, 1.0);
    const double D2= Dplus*Dplus*pow(omega/Omega_m, -1.0/143.0);
    const double D20= pow(Omega_m, -1.0/143.0);
    const double f1 = pow(omega, (4./7));
    const double f2 = pow(omega, (6./11));

    int64_t strides[] = {param->nc * param->nc, param->nc, 1};

    /* RUN PB ic global shifting */
    const double offset0 = 0.5 * 1.0 / param->nc;
    double dx1max = 0;
    double dx2max = 0;
    for(int ip = 0; ip < offset; ip ++) {
        float * x = p->x[ip];
        float * v = p->v[ip];
        float * dx1 = p->dx1[ip];
        float * dx2 = p->dx2[ip];
        
        int64_t id = p->id[ip];
        int64_t id0 = id;

        for(int d = 0; d < 3; d ++ ) {
            double opos = (id / strides[d]) * (1.0 / param->nc) + offset0;
            id %= strides[d];
            double disp = x[d] - opos;
            if(disp < -0.5) disp += 1.0;
            if(disp > 0.5) disp -= 1.0;
            dx1[d] = (v[d] - disp * (2 * f2)) / (f1 - 2 * f2) / DplusIC;
            dx2[d] = (v[d] - disp * f1) / (2 * f2 - f1) / (DplusIC * DplusIC);
            /* evolve to a_init with 2lpt */
            double tmp = opos + dx1[d] * Dplus + dx2[d] * (D20 * D2);
            x[d] = tmp * param->boxsize;
            while(x[d] < 0.0) x[d] += param->boxsize;
            while(x[d] >= param->boxsize) x[d] -= param->boxsize;
            dx1[d] *= param->boxsize;

            if(dx1[d] > 100) {
                printf("id = %ld dx1[d] = %g v = %g pos = %g disp = %g opos=%g f1=%g f1=%g Dplus=%g, D20=%g, D2=%g, DplusIC=%g\n", 
                    id0, dx1[d], v[d], x[d], disp, opos, f1, f2, Dplus, D20, D2, DplusIC);
            }
            dx2[d] *= param->boxsize;

            if(fabs(dx1[d]) > dx1max) dx1max = fabs(dx1[d]);
            if(fabs(dx2[d]) > dx2max) dx2max = fabs(dx2[d]);
            v[d] = 0.0;
        }
        
    }
    double dx1maxg, dx2maxg;
    MPI_Allreduce(&dx1max, &dx1maxg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&dx2max, &dx2maxg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    msg_printf(verbose, "dx1 max = %g, dx2 max = %g", dx1maxg, dx2maxg);
    free(NcumFile);
    free(NperFile);
    heap_return(scratch);
    return 0;
}

int write_runpb_snapshot(Particles * p,  
        char * filebase){
    int ThisTask;
    int NTask;
    size_t scratch_bytes = p->np_total * 8;
    void * scratch = heap_allocate(scratch_bytes);

    float * fscratch = (float*) scratch;
    int64_t * lscratch = (int64_t *) scratch;

    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    size_t chunksize = scratch_bytes;

    int np_local = p->np_local;

    int * NperTask = alloca(sizeof(int) * NTask);
    size_t * NcumTask = alloca(sizeof(size_t) * (NTask + 1));

    MPI_Allgather(&np_local, 1, MPI_INT, NperTask, 1, MPI_INT, MPI_COMM_WORLD);

    NcumTask[0] = 0;
    for(int i = 1; i <= NTask; i ++) {
        NcumTask[i] = NcumTask[i - 1] + NperTask[i - 1];
    }
    size_t Ntot = NcumTask[NTask];

    double aa = p->a;
    int Nfile = (Ntot + (1024 * 1024 * 128 - 1)) / (1024 * 1024 * 128);

    msg_printf(verbose, "Writing to %td paritlces to  %d files\n", Ntot, Nfile);

    double vfac = 100. / aa;
    double RSD = aa / p->qfactor / vfac;

    int * NperFile = NperFile = alloca(sizeof(int) * Nfile);

    size_t * NcumFile = alloca(sizeof(size_t) * Nfile);

    NcumFile[0] = 0;
    for(int i = 0; i < Nfile; i ++) {
        NperFile[i] = (i+1) * Ntot / Nfile -i * Ntot/ Nfile;
    }
    for(int i = 1; i < Nfile; i ++) {
        NcumFile[i] = NcumFile[i - 1] + NperFile[i - 1];
    }


    size_t start = NcumTask[ThisTask];
    size_t end   = NcumTask[ThisTask + 1];

    int offset = 0;
    int chunknpart = chunksize / (sizeof(float) * 3);

    for(int i = 0; i < Nfile; i ++) {
        char buf[1024];
        sprintf(buf, FILENAME, filebase, i);
        if(ThisTask == 0) {
            FILE * fp = fopen(buf, "w");
            fclose(fp);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < Nfile; i ++) {
        ptrdiff_t mystart = start - NcumFile[i];
        ptrdiff_t myend = end - NcumFile[i];
        size_t nread;
        /* not overlapping */
        if(myend <= 0) continue;
        if(mystart >= NperFile[i]) continue;

        /* cut to this file */
        if(myend > NperFile[i]) {
            myend = NperFile[i];
        }
        if(mystart < 0) {
            mystart =0;
        }
        /* at this point, write mystart:myend from current file */
        char buf[1024];
        int eflag = 1, hsize = sizeof(FileHeader);
        FileHeader header;
        header.npart = NperFile[i];
        header.nsph = 0;
        header.nstar = 0;
        header.aa    = aa;
        header.eps =  0.1 / pow(Ntot, 1./3);

        sprintf(buf, FILENAME, filebase, i);

        FILE * fp = fopen(buf, "r+");
        /* skip these */
        fwrite(&eflag, sizeof(int), 1, fp);
        fwrite(&hsize, sizeof(int), 1, fp);
        fwrite(&header, sizeof(FileHeader), 1, fp);
        /* pos */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                for(int d = 0; d < 3; d++) 
                    fscratch[q++] = p->x[offset + nread + ip][d] / p->boxsize;
            }
            fwrite(scratch, sizeof(float) * 3, nbatch, fp);
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* vel */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                for(int d = 0; d < 3; d++) 
                    fscratch[q++] = p->v[offset + nread + ip][d] * RSD / p->boxsize;
            }
            fwrite(scratch, sizeof(float) * 3, nbatch, fp);
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* ID */
        fseek(fp, mystart * sizeof(int64_t), SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            for(int ip = 0, q = 0; ip < nbatch; ip ++) {
                lscratch[q++] = p->id[offset + nread + ip];
            }
            fwrite(scratch, sizeof(int64_t), nbatch, fp);
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(int64_t), SEEK_CUR);
        fclose(fp);
        offset += myend - mystart;
    }
    if(offset != p->np_local) {
        msg_abort(0030, "mismatch %d != %d\n", offset, p->np_local);
    }
    heap_return(scratch);
    return 0;
}
