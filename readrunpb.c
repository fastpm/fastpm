#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>
#include <glob.h>
#include <math.h>
#include "parameters.h"
#include "power.h"
#include "particle.h"

#define FILENAME  "%s.%02d"

typedef struct {
  int   npart;          /* Total number of particles. */
  int   nsph;           /* Number of gas particles.   */
  int   nstar;          /* Number of star particles.  */
  float aa;             /* Scale factor. */
  float eps;            /* Gravitational softening    */
} FileHeader;

int read_runpb_ic(Parameters * param, double a_init, Particles * particles, 
        void * scratch) {
    int ThisTask;
    int NTask;
    float * fscratch = (float*) scratch;
    long long * lscratch = (long long *) scratch;

    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

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
            msg_abort(0030, "Number of particles does not match nc\n");
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
    particles->np_local = end - start;

    
    Particle * par = particles->p;
    int offset = 0;
    for(int i = 0; i < Nfile; i ++) {
        ptrdiff_t mystart = start - NcumFile[i];
        ptrdiff_t myend = end - NcumFile[i];
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
        fread(scratch, sizeof(float) * 3, (myend - mystart), fp);
        for(int p = 0, q = 0; p < myend - mystart; p ++) {
            par[offset + p].x[0] = fscratch[q++];
            par[offset + p].x[1] = fscratch[q++];
            par[offset + p].x[2] = fscratch[q++];
        }
        /* vel */
        fseek(fp, mystart * sizeof(float) * 3, SEEK_CUR);
        fread(scratch, sizeof(float) * 3, (myend - mystart), fp);
        for(int p = 0, q = 0; p < myend - mystart; p ++) {
            par[offset + p].v[0] = fscratch[q++];
            par[offset + p].v[1] = fscratch[q++];
            par[offset + p].v[2] = fscratch[q++];
        }

        /* ID */
        fseek(fp, mystart * sizeof(long long), SEEK_CUR);
        fread(scratch, sizeof(long long), (myend - mystart), fp);
        for(int p = 0, q = 0; p < myend - mystart; p ++) {
            par[offset + p].id = lscratch[q++];
        }
        fclose(fp);
        offset += myend - mystart;
    }
    if(offset != particles->np_local) {
        msg_abort(0030, "mismatch %d != %d\n", offset, particles->np_local);
    }
    const double Omega_m = param->omega_m;
    const double omega=Omega_m/(Omega_m + (1.0 - Omega_m)*aa*aa*aa);
    const float DplusIC = 1.0/GrowthFactor(aa, 1.0);
    const float Dplus = 1.0/GrowthFactor(a_init, 1.0);
    const double D2= Dplus*Dplus*pow(omega/Omega_m, -1.0/143.0);
    const double D20= pow(Omega_m, -1.0/143.0);
    const double f1 = pow(omega, (4./7));
    const double f2 = pow(omega, (6./11));

    long long strides[] = {param->nc * param->nc, param->nc, 1};

    /* RUN PB ic global shifting */
    const double offset0 = 0.5 * 1.0 / param->nc;
    for(int p = 0; p < offset; p ++) {
        float * x = par[p].x;
        float * v = par[p].v;
        float * dx1 = par[p].dx1;
        float * dx2 = par[p].dx2;
        
        long long id = par[p].id;

        for(int d = 0; d < 3; d ++ ) {
            double opos = (id / strides[d]) * (1.0 / param->nc) + offset0;
            id %= strides[d];
            double disp = x[d] - opos;
            if(disp < -0.5) disp += 1.0;
            if(disp > 0.5) disp -= 1.0;
            dx1[d] = (v[d] - disp * (2 * f2)) / (f1 - 2 * f2) / DplusIC;
            dx2[d] = (v[d] - disp * f1) / (2 * f2 - f1) / (DplusIC * DplusIC);
            /* evolve to a_init with 2lpt */
            x[d] = opos + dx1[d] * Dplus + dx2[d] * (D20 * D2);
            while(x[d] >= 1.0) x[d] -= 1.0;
            while(x[d] < 0.0) x[d] += 1.0;
            x[d] *= param->boxsize;
            dx1[d] *= param->boxsize;
            dx2[d] *= param->boxsize;
            v[d] = 0.0;
        }
        
        
    }
    free(NcumFile);
    free(NperFile);
    return 0;
}
