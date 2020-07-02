/**
 * IO of RunPB format
 * 
 * Authors: Yu Feng <rainwoodman@gmail.com>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <mpi.h>
#include <glob.h>
#include <math.h>
#include <alloca.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/prof.h>
#include <fastpm/string.h>

#define FILENAME  "%s.%02d"

typedef struct {
  int   npart;          /* Total number of p. */
  int   nsph;           /* Number of gas p.   */
  int   nstar;          /* Number of star p.  */
  float aa;             /* Scale factor. */
  float eps;            /* Gravitational softening    */
} FileHeader;

void
fwrite_or_die(void * ptr, size_t elsize, size_t size, FILE * fp,
    const char * message
    )
{
    if(size != fwrite(ptr, elsize, size, fp)) {
        fastpm_raise(0030, "failed to write %td bytes for %s. \n", size, message);
    }
}

int
read_runpb_ic(FastPMSolver * fastpm, FastPMStore * p, const char * filename)
{
    int ThisTask = fastpm->ThisTask;
    int NTask = fastpm->NTask;
    MPI_Comm comm = fastpm->comm;

    size_t scratch_bytes = 32 * 1024 * 1024;
    void * scratch = malloc(scratch_bytes);
    float * fscratch = (float*) scratch;
    int64_t * lscratch = (int64_t *) scratch;

    size_t chunksize = scratch_bytes;

    size_t Ntot = 0;

    int * NperFile = NULL;
    size_t * NcumFile = NULL;
    int Nfile = 0;
    double aa = 0;

    if (ThisTask == 0) {
        char buf[1024];
        int i;
        for(i = 0; ; i ++) {
            sprintf(buf, FILENAME, filename, i);
            FILE * fp = fopen(buf, "r");
            if(!fp) {
                Nfile = i;
                break;
            }
            fclose(fp);
        }
        if (Nfile == 0) {
            fastpm_raise(0030, "No snapshot files were found.\n");
        }
        fastpm_ilog(INFO, "Total number of files is %d\n", Nfile);

        MPI_Bcast(&Nfile, 1, MPI_INT, 0, comm);
        NperFile = malloc(sizeof(int) * Nfile);

        for(i = 0; i < Nfile; i ++) {
            sprintf(buf, FILENAME, filename, i);
            FILE * fp = fopen(buf, "r");
            FileHeader header;
            int eflag, hsize;
            if(1 != fread(&eflag, sizeof(int), 1, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            if(1 != fread(&hsize, sizeof(int), 1, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            if(hsize != sizeof(header)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            if(1 != fread(&header, sizeof(FileHeader), 1, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            aa = header.aa;
            fastpm_ilog(INFO, "reading from file %s npart=%d aa=%g \n", buf, header.npart, header.aa);
            NperFile[i] = header.npart;
            Ntot += header.npart;        
            fclose(fp);
        }
        if (Ntot != fastpm->config->nc * fastpm->config->nc * fastpm->config->nc) {
            fastpm_raise(0030, "Number of p does not match nc\n");
        }
        MPI_Bcast(NperFile, Nfile, MPI_INT, 0, comm);
        MPI_Bcast(&Ntot, 1, MPI_LONG_LONG, 0, comm);
        MPI_Bcast(&aa, 1, MPI_DOUBLE, 0, comm);
    } else {
        /* other ranks receive */
        MPI_Bcast(&Nfile, 1, MPI_INT, 0, comm);
        NperFile = malloc(sizeof(int) * Nfile);
        MPI_Bcast(NperFile, Nfile, MPI_INT, 0, comm);
        MPI_Bcast(&Ntot, 1, MPI_LONG_LONG, 0, comm);
        MPI_Bcast(&aa, 1, MPI_DOUBLE, 0, comm);
    }

    fastpm_info("Ntot = %td aa=%g\n", Ntot, aa);


    NcumFile = malloc(sizeof(size_t) * Nfile);
    NcumFile[0] = 0;
    int i;
    for(i = 1; i < Nfile; i ++) {
        NcumFile[i] = NcumFile[i - 1] + NperFile[i - 1];
    }

    size_t start = ThisTask * Ntot / NTask;
    size_t end   = (ThisTask + 1) * Ntot / NTask;
    p->np = end - start;

    int offset = 0;
    int chunknpart = chunksize / (sizeof(float) * 3);
    fastpm_info("chunknpart = %d\n", chunknpart);
    for(i = 0; i < Nfile; i ++) {
        ptrdiff_t mystart = start - NcumFile[i];
        ptrdiff_t myend = end - NcumFile[i];
        size_t nread;
        /* not overlapping */
        if(myend <= 0) continue;
        if(mystart >= NperFile[i]) continue;

        //printf("Task %d reading at %d \n", ThisTask, offset);

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
        sprintf(buf, FILENAME, filename, i);

        FILE * fp = fopen(buf, "r");
        /* skip these */
        if(1 != fread(&eflag, sizeof(int), 1, fp)) {
            fastpm_raise(0030, "Unable to read from %s\n", buf);
        }
        if(1 != fread(&hsize, sizeof(int), 1, fp)) {
            fastpm_raise(0030, "Unable to read from %s\n", buf);
        }

        if(1 != fread(&header, sizeof(FileHeader), 1, fp)) {
            fastpm_raise(0030, "Unable to read from %s\n", buf);
        }
        /* pos */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            if (nbatch != fread(scratch, sizeof(float) * 3, nbatch, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                int d;
                for(d = 0; d < 3; d ++)
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
            if (nbatch != fread(scratch, sizeof(float) * 3, nbatch, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                int d;
                for(d = 0; d < 3; d ++)
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
            if (nbatch != fread(scratch, sizeof(int64_t), nbatch, fp)) {
                fastpm_raise(0030, "Unable to read from %s\n", buf);
            }
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                p->id[offset + nread + ip] = lscratch[q++];
            }
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(int64_t), SEEK_CUR);
        fclose(fp);
        offset += myend - mystart;
    }
    if(offset != p->np) {
        fastpm_raise(0030, "mismatch %d != %d\n", offset, p->np);
    }

    const double omega = Omega_cdm_a(aa, fastpm->cosmology);

    FastPMGrowthInfo gi;
    fastpm_growth_info_init(&gi, aa, fastpm->cosmology);
    const float DplusIC = gi.D1;
    const double f1 = pow(omega, (4./7));
    const double f2 = pow(omega, (6./11));

    int64_t strides[] = {fastpm->config->nc * fastpm->config->nc, fastpm->config->nc, 1};

    int d;
    for(d = 0; d < 3; d++){
        p->meta._q_strides[d] = strides[d];
        p->meta._q_scale[d] = (1.0 * fastpm->config->boxsize) / fastpm->config->nc;
        p->meta._q_shift[d] = 0.5 * (1.0 * fastpm->config->boxsize) / fastpm->config->nc;
    }

    /* RUN PB ic global shifting */
    const double offset0 = 0.5 * 1.0 / fastpm->config->nc;
    double dx1disp[3] = {0};
    double dx2disp[3] = {0};
    int ip;
    for(ip = 0; ip < offset; ip ++) {
        double * x = p->x[ip];
        float * v = p->v[ip];
        float * dx1 = p->dx1[ip];
        float * dx2 = p->dx2[ip];

        int64_t id = p->id[ip];
        //int64_t id0 = id;
        int d;
        for(d = 0; d < 3; d ++ ) {
            double opos = (id / strides[d]) * (1.0 / fastpm->config->nc) + offset0;
            id %= strides[d];
            double disp = x[d] - opos;
            if(disp < -0.5) disp += 1.0;
            if(disp > 0.5) disp -= 1.0;
            dx1[d] = (v[d] - disp * (2 * f2)) / (f1 - 2 * f2) / DplusIC;
            /* no 7/3 here, this ensures  = x0 + dx1 + dx2; we shift the position in pm_2lpt_evolve  */
            dx2[d] = (v[d] - disp * f1) / (2 * f2 - f1) / (DplusIC * DplusIC);
            double boxsize = fastpm->config->boxsize;
            double tmp = opos; 
            x[d] = tmp * boxsize;
            while(x[d] < 0.0) x[d] += boxsize;
            while(x[d] >= boxsize) x[d] -= boxsize;
            dx1[d] *= boxsize;

            /*
            if(dx1[d] > 100) {
                printf("id = %ld dx1[d] = %g v = %g pos = %g disp = %g opos=%g f1=%g f1=%g DplusIC=%g\n", 
                    id0, dx1[d], v[d], x[d], disp, opos, f1, f2, DplusIC);
            } */
            dx2[d] *= boxsize;

            v[d] = 0.0;
            dx1disp[d] += dx1[d] * dx1[d];
            dx2disp[d] += dx2[d] * dx2[d];
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, dx1disp, 3, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, dx2disp, 3, MPI_DOUBLE, MPI_SUM, comm);
    for(d =0; d < 3; d++) {
        dx1disp[d] /= Ntot;
        dx1disp[d] = sqrt(dx1disp[d]);
        dx2disp[d] /= Ntot;
        dx2disp[d] = sqrt(dx2disp[d]);
    }
    fastpm_info("dx1 disp : %g %g %g %g\n", 
            dx1disp[0], dx1disp[1], dx1disp[2],
            (dx1disp[0] + dx1disp[1] + dx1disp[2]) / 3.0);
    fastpm_info("dx2 disp : %g %g %g %g\n", 
            dx2disp[0], dx2disp[1], dx2disp[2],
            (dx2disp[0] + dx2disp[1] + dx2disp[2]) / 3.0);

    free(NcumFile);
    free(NperFile);
    free(scratch);

    return 0;
}

static void write_mine(const char * filebase,
            FastPMStore * p, double aa, FastPMCosmology * c, double boxsize, size_t Ntot,
            size_t * NcumFile, int * NperFile, int Nfile, 
            ptrdiff_t start, ptrdiff_t end) {
    size_t scratch_bytes = 32 * 1024 * 1024;
    void * scratch = malloc(scratch_bytes);

    float * fscratch = (float*) scratch;
    int64_t * lscratch = (int64_t *) scratch;

    double H0 = 100.;
    double RSD = 1 / (aa * HubbleEa(aa, c) * H0);


    size_t chunksize = scratch_bytes;

    int offset = 0; 
    int chunknpart = chunksize / (sizeof(float) * 3);

    int i;
    for(i = 0; i < Nfile; i ++) {
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
        fwrite_or_die(&eflag, sizeof(int), 1, fp, "eflag");
        fwrite_or_die(&hsize, sizeof(int), 1, fp, "hsize");
        fwrite_or_die(&header, sizeof(FileHeader), 1, fp, "header");
        /* pos */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                int d;
                for(d = 0; d < 3; d++) 
                    fscratch[q++] = p->x[offset + nread + ip][d] / boxsize;
            }
            fwrite_or_die(scratch, sizeof(float) * 3, nbatch, fp, "pos");
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* vel */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                int d;
                for(d = 0; d < 3; d++) 
                    fscratch[q++] = p->v[offset + nread + ip][d] * RSD / boxsize;
            }
            fwrite_or_die(scratch, sizeof(float) * 3, nbatch, fp, "vel");
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(float) * 3, SEEK_CUR);
        /* ID */
        fseek(fp, mystart * sizeof(int64_t), SEEK_CUR);
        nread = 0;
        while(nread != myend - mystart) {
            size_t nbatch = chunknpart;
            if (nbatch + nread > myend - mystart) nbatch = myend - mystart - nread;
            int ip, q;
            for(ip = 0, q = 0; ip < nbatch; ip ++) {
                lscratch[q++] = p->id[offset + nread + ip];
            }
            fwrite_or_die(scratch, sizeof(int64_t), nbatch, fp, "id");
            nread += nbatch;
        }
        fseek(fp, (NperFile[i] - myend) * sizeof(int64_t), SEEK_CUR);
        fclose(fp);
        offset += myend - mystart;
    }
    if(offset != p->np) {
        fastpm_raise(0030, "mismatch %d != %d\n", offset, p->np);
    }
    free(scratch);
}

int 
write_runpb_snapshot(FastPMSolver * fastpm, FastPMStore * p, const char * filebase)
{
    CLOCK(meta);

    int ThisTask = fastpm->ThisTask;
    int NTask = fastpm->NTask;
    MPI_Comm comm = fastpm->comm;

    ENTER(meta);
    fastpm_path_ensure_dirname(filebase);
    MPI_Barrier(fastpm->comm);
    LEAVE(meta);

    double aa = p->meta.a_x;

    int np = p->np;
    int i;
    int * NperTask = alloca(sizeof(int) * NTask);
    size_t * NcumTask = alloca(sizeof(size_t) * (NTask + 1));

    MPI_Allgather(&np, 1, MPI_INT, NperTask, 1, MPI_INT, comm);

    NcumTask[0] = 0;
    for(i = 1; i <= NTask; i ++) {
        NcumTask[i] = NcumTask[i - 1] + NperTask[i - 1];
    }
    size_t Ntot = NcumTask[NTask];

    int Nfile = (Ntot + (1024 * 1024 * 128 - 1)) / (1024 * 1024 * 128);

    fastpm_info("Writing to %td paritlces to  %d files\n", Ntot, Nfile);

    int * NperFile = NperFile = alloca(sizeof(int) * Nfile);

    size_t * NcumFile = alloca(sizeof(size_t) * Nfile);

    NcumFile[0] = 0;
    for(i = 0; i < Nfile; i ++) {
        NperFile[i] = (i+1) * Ntot / Nfile -i * Ntot/ Nfile;
    }
    for(i = 1; i < Nfile; i ++) {
        NcumFile[i] = NcumFile[i - 1] + NperFile[i - 1];
    }


    size_t start = NcumTask[ThisTask];
    size_t end   = NcumTask[ThisTask + 1];

    for(i = 0; i < Nfile; i ++) {
        char buf[1024];
        sprintf(buf, FILENAME, filebase, i);
        if(ThisTask == 0) {
            FILE * fp = fopen(buf, "w");
            fclose(fp);
        }
    }
    MPI_Barrier(comm);

    write_mine(filebase, p, aa, fastpm->cosmology, fastpm->config->boxsize, Ntot, NcumFile, NperFile, Nfile, start, end);
    return 0;
}

