#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/string.h>
#include <fastpm/transfer.h>
#include "chealpix.h"
#include <mpsort.h>

#include "pmpfft.h"
#include "pmghosts.h"

static double RNDTABLE[8192];

gsl_rng * random_generator;
void fastpm_utils_init_randtable() {
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    gsl_rng_set(random_generator, 37 * (rank + 1) + 1);  /* start-up seed */
    int i;
    for(i = 0; i < 8192; i ++) {
        RNDTABLE[i] = gsl_rng_uniform(random_generator);
    }
    //gsl_rng_free(random_generator);
}

double 
fastpm_utils_get_random(uint64_t id) 
{
    return gsl_rng_uniform(random_generator);
    uint64_t ind = 0;
    ind = id;
    
    while(id != 0) {
        ind += id;
        id /= 8192;
    } 
    ind %= 8192;
    return RNDTABLE[ind];
}

void 
fastpm_utils_dump(PM * pm , const char * filename, FastPMFloat *data) 
{
    char fn1[1024];
    char fn2[1024];
    if(pm->NTask > 1) {
        sprintf(fn1, "%s.%03d.geometry", filename, pm->ThisTask);
        sprintf(fn2, "%s.%03d", filename, pm->ThisTask);
    } else {
        sprintf(fn1, "%s.geometry", filename);
        sprintf(fn2, "%s", filename);
    }

    fastpm_path_ensure_dirname(filename);

    FILE * fp;
    fp = fopen(fn2, "w");
    if(pm->allocsize != fwrite(data, sizeof(FastPMFloat), pm->allocsize, fp)) {
        fastpm_raise(-1, "IO failed.\n");
    }
    fclose(fp);
    fp = fopen(fn1, "w");
    fprintf(fp, "# real\n");
    fprintf(fp, "start: %td %td %td\n", 
                    pm->IRegion.start[0],
                    pm->IRegion.start[1],
                    pm->IRegion.start[2]);
    fprintf(fp, "size: %td %td %td\n", 
                    pm->IRegion.size[0],
                    pm->IRegion.size[1],
                    pm->IRegion.size[2]);
    fprintf(fp, "strides: %td %td %td\n", 
                    pm->IRegion.strides[0],
                    pm->IRegion.strides[1],
                    pm->IRegion.strides[2]);

    fprintf(fp, "# complex\n");
    fprintf(fp, "start: %td %td %td\n", 
                    pm->ORegion.start[0],
                    pm->ORegion.start[1],
                    pm->ORegion.start[2]);
    fprintf(fp, "size: %td %td %td\n", 
                    pm->ORegion.size[0],
                    pm->ORegion.size[1],
                    pm->ORegion.size[2]);
    fprintf(fp, "strides: %td %td %td\n", 
                    pm->ORegion.strides[0],
                    pm->ORegion.strides[1],
                    pm->ORegion.strides[2]);

    fclose(fp);
}

void 
fastpm_utils_load(PM * pm , const char * filename, FastPMFloat *data) 
{
    char fn1[1024];
    char fn2[1024];
    if(pm->NTask > 1) {
        sprintf(fn1, "%s.%03d.geometry", filename, pm->ThisTask);
        sprintf(fn2, "%s.%03d", filename, pm->ThisTask);
    } else {
        sprintf(fn1, "%s.geometry", filename);
        sprintf(fn2, "%s", filename);
    }
    FILE * fp;
    fp = fopen(fn2, "r");
    if(pm->allocsize != fread(data, sizeof(FastPMFloat), pm->allocsize, fp)) {
        fastpm_raise(-1, "File was bad\n");
    };
    fclose(fp);
}

static double 
tk_eh(double k, struct fastpm_powerspec_eh_params * params)		/* from Martin White */
{
    double q, theta, ommh2, a, s, gamma, L0, C0;
    double tmp;
    double omegam, ombh2, hubble;

    /* other input parameters */
    hubble = params->hubble_param;
    omegam = params->omegam;
    ombh2 = params->omegab * hubble * hubble;

    theta = 2.728 / 2.7;
    ommh2 = omegam * hubble * hubble;
    s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
    a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
        + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
    gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
    gamma *= omegam * hubble;
    q = k * theta * theta / gamma;
    L0 = log(2. * exp(1.) + 1.8 * q);
    C0 = 14.2 + 731. / (1. + 62.5 * q);
    tmp = L0 / (L0 + C0 * q * q);
    return (tmp);
}

double 
fastpm_utils_powerspec_eh(double k, struct fastpm_powerspec_eh_params * param)	/* Eisenstein & Hu */
{
    return param->Norm * k * pow(tk_eh(k, param), 2);
}

double 
fastpm_utils_powerspec_white(double k, double * amplitude) 	/* White Noise*/
{
    return *amplitude;
}

static void
_sort_pix(const void * ptr, void * radix, void * arg)
{
    memcpy(radix, ptr, 8);
}

/* generate evenly balanced healpix positions that are inside the lightcone. */
void
fastpm_utils_healpix_ra_dec (
                int nside,
                double **ra,
                double **dec,
                uint64_t **pix,
                size_t * n,
                FastPMLightCone * lc,
                MPI_Comm comm
            )
{
    const double rad_to_degree = 180./M_PI;

    size_t npix = nside2npix (nside);
    int ThisTask, NTask;

    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NTask);

    //fastpm_info("healpix npix %ld \n",*npix);
    size_t i = 0;

    size_t pix_start = ThisTask * npix / NTask;
    size_t pix_end = (ThisTask + 1) * npix / NTask;

    uint64_t local_npix = 0;

    uint64_t * pixels = NULL;

    /* two iterations; estimate and fill */
    while(1) {
        size_t j = 0;

        for (i = pix_start; i < pix_end; i++)
        {
            double vec[3];

            pix2vec_ring(nside, i, vec);

            if(!fastpm_lc_inside(lc, vec)) continue;

            if(pixels != NULL) {
                pixels[j] = i;
            }

            j ++;
        }
        if(pixels == NULL) {
            local_npix = j;
            pixels = malloc(sizeof(uint64_t) * local_npix);

        } else {
            break;
        }
    }

    /* count total */
    uint64_t valid_npix = local_npix;
    MPI_Allreduce(MPI_IN_PLACE, &valid_npix, 1, MPI_LONG, MPI_SUM, comm);

    /* redistribute / balance */

    size_t localsize =
        (ThisTask + 1) * valid_npix / NTask
    -
        (ThisTask    ) * valid_npix / NTask;

    uint64_t * recv_buffer = malloc(sizeof(uint64_t) * localsize);

    mpsort_mpi_newarray(pixels, local_npix, recv_buffer, localsize, sizeof(uint64_t), _sort_pix, 8, NULL, comm);

    free(pixels);

    *ra = malloc(sizeof(double) * localsize);
    *dec = malloc(sizeof(double) * localsize);

    for(i = 0; i < localsize; i ++) {
        double phi, theta;

        pix2ang_ring(nside, recv_buffer[i], &theta, &phi);

        phi *= rad_to_degree;
        theta*= rad_to_degree;
        (*ra)[i] = phi;
        (*dec)[i]= 90 - theta;
        recv_buffer[i] += (((uint64_t)nside) << 48);
    }
    *n = localsize;
    *pix = recv_buffer;
}

void
fastpm_gldot(double glmatrix[4][4], double xi[4], double xo[4])
{
    int i, j;
    for(i = 0; i < 4; i ++) {
        xo[i] = 0;
        for(j = 0; j < 4; j ++) {
            xo[i] += glmatrix[i][j] * xi[j];
        }
    }
}

void
fastpm_gldotf(double glmatrix[4][4], float vi[4], float vo[4])
{
    int i, j;
    for(i = 0; i < 4; i ++) {
        vo[i] = 0;
        for(j = 0; j < 4; j ++) {
            vo[i] += glmatrix[i][j] * vi[j];
        }
    }
}

/* returns MPI_LOR reduction */
int
MPIU_Any(MPI_Comm comm, int value)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_LOR, comm);
    return value;
}

/* returns MPI_LAND reduction */
int
MPIU_All(MPI_Comm comm, int value)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_LAND, comm);
    return value;
}

/*
 * Compute the standard statistics of a value over communicator.
 * fmt is a string, describing the list of stats to return:
 *  < : minimum  (double)
 *  , : rank of minimum (int)
 *  > : maximum  (double)
 *  . : rank of maximum (int)
 *  - : mean  (double)
 *  s : std (biased standard derivation)
 *  v : variance (biased variance)
 *  S : unbiased std
 *  V : unbiased variance
 *
 * e.g.
 *
 *    double min, max;
 *    MPIU_stats(comm, r, "<>", &min, &max);
 *
 * returns number of items converted;
 * if a fmt char is unknown, the input variable is unmodified.
 *
 * */
int
MPIU_stats(MPI_Comm comm,
    const double value,
    const char * fmt, ...)
{
    
    va_list va;
    va_start(va, fmt);

    int NTask;
    MPI_Comm_size(comm, &NTask);
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);

    double sum = value;
    struct {
        double value;
        int rank;
    } min = {value, ThisTask},
      max = {value, ThisTask};

    double sumsq = value * value;

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &sumsq, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    double avg = (sum) / NTask;
    double avg_sq = (sumsq) / NTask;
    double var = avg_sq - avg * avg;
    double std = sqrt(avg_sq - avg * avg);

    int i;
    int c = 0;
    for(i = 0; i < strlen(fmt); i ++) {
        void * r = va_arg(va, void *);
        int * ir = (int * ) r;
        double * dr = (double * ) r;

        c++;
        switch(fmt[i]) {
            case '-':
                *dr = avg;
                break;
            case '<':
                *dr = min.value;
                break;
            case '>':
                *dr = max.value;
                break;
            case ',':
                *ir = min.rank;
                break;
            case '.':
                *ir = max.rank;
                break;
            case 's':
                *dr = std;
                break;
            case 'S':
                *dr = std * sqrt(1.0 * NTask / (NTask - 1.0));
                break;
            case 'v':
                *dr = var;
                break;
            case 'V':
                *dr = var * 1.0 * NTask / (NTask - 1.0);
                break;
            default:
                c--;
        }
    }
    va_end(va);

    return c;
}

/*
 *
 * create a copy of string src and broadcast it to all ranks.
 *
 * Every rank receives a copy of src as the return value, including the root.
 *
 * if free is not NULL, the function is called as free(src)
 * */
char *
MPIU_Bcast_string(MPI_Comm comm, char * src, int root, void(*free)(void * ptr))
{
    int srcNULL = src == NULL;
    MPI_Bcast(&srcNULL, 1, MPI_INT, root, comm);

    if(srcNULL) return NULL;

    int rank;
    int N;
    MPI_Comm_rank(comm, &rank);

    if(rank == root) {
        N = strlen(src) + 1;
    }
    MPI_Bcast(&N, 1, MPI_INT, root, comm);

    char * ret = malloc(N);
    if(rank == root) {
        strcpy(ret, src);
    }

    if(free && rank == root) free(src);
    MPI_Bcast(ret, N, MPI_BYTE, root, comm);

    return ret;
}
