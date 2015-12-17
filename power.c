//
// Reads CAMB matter power spectrum camb_matterpower.dat, and
// provides power spectrum to 2LPT calculation in lpt.c
//
// Originally written by Jun Kuda, based on 
//   N-GenIC power.c by Volker Springel
//   http://www.mpa-garching.mpg.de/gadget/right.html#ICcode
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <mpi.h>

#include "msg.h"
#include "power.h"

#define WORKSIZE 100000

static double PowerNorm;
static int NPowerTable;
static double Omega, OmegaLambda;

static struct pow_table
{
    double logk, logD;
} *PowerTable;

static void read_power_table_camb(const char filename[]);
static double normalize_power(const double a_init, const double sigma8);
static double TopHatSigma2(double R, double (*func)(double, void*), void*);

void power_init(const char filename[], const double a_init, const double sigma8, const double omega_m, const double omega_lambda,
            MPI_Comm comm)
{
    // CAMB matter power spectrum filename
    Omega= omega_m;
    OmegaLambda= omega_lambda;

    if(filename == NULL) return;
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if(myrank == 0) {
        read_power_table_camb(filename);
        PowerNorm= normalize_power(a_init, sigma8);
    }

    msg_printf(normal, "Powerspecectrum file: %s\n", filename);

    MPI_Bcast(&PowerNorm, 1, MPI_DOUBLE, 0, comm);

    MPI_Bcast(&NPowerTable, 1, MPI_INT, 0, comm);

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(struct pow_table), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Bcast(&PowerTable[0], NPowerTable, type, 0, comm);
    MPI_Type_free(&type);
}

void read_power_table_camb(const char filename[])
{
    //char buf[500];
    double k, p;
    double fac= 1.0/(2.0*M_PI*M_PI);
    PowerNorm= 1.0;

    FILE* fp= fopen(filename, "r");
    if(!fp)
        msg_abort(3000, "Error: unable to read input power spectrum: %s",
                filename);

    NPowerTable = 0;
    do {
        if(fscanf(fp, "%lg %lg ", &k, &p) == 2)
            NPowerTable++;
        else
            break;
    } while(1);

    msg_printf(verbose, 
            "Found %d pairs of values in input spectrum table\n", NPowerTable);

    PowerTable = malloc(NPowerTable * sizeof(struct pow_table));

    rewind(fp);

    int n = 0;
    do {
        if(fscanf(fp, " %lg %lg ", &k, &p) == 2) {
            PowerTable[n].logk = log10(k);
            PowerTable[n].logD = log10(fac*k*k*k*p);
            //printf("%e %e\n", PowerTable[n].logk, PowerTable[n].logD);
            n++;
        }
        else
            break;
    } while(1);
    assert(NPowerTable == n);

    fclose(fp);
}

double PowerSpecWithData(double k, void * data) {
    return PowerSpec(k);
}

double normalize_power(const double a_init, const double sigma8)
{
    // Assume that input power spectrum already has a proper sigma8
    const double R8 = 8.0; // 8 Mpc

    double res = TopHatSigma2(R8, PowerSpecWithData, NULL); 
    double sigma8_input= sqrt(res);

    msg_printf(info, "Input power spectrum sigma8 %f\n", sigma8_input);
    if (sigma8 == 0)
        return 1;
    msg_printf(info, "Expected power spectrum sigma8 %f\n", sigma8);

    if(fabs(sigma8_input - sigma8)/sigma8 > 0.05)
        msg_printf(info, "Input sigma8 %f is far from target sigma8 %f; correction applied\n",
                sigma8_input, sigma8);

    return sigma8 * sigma8/ res;
}

double PowerSpec(const double k)
{
    // k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
    // convert to h/Mpc

    const double logk = log10(k);

    if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
        return 0;

    int binlow = 0;
    int binhigh = NPowerTable - 1;

    while(binhigh - binlow > 1) {
        int binmid = (binhigh + binlow) / 2;
        if(logk < PowerTable[binmid].logk)
            binhigh = binmid;
        else
            binlow = binmid;
    }

    const double dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;
    assert(dlogk > 0.0);

    const double u = (logk - PowerTable[binlow].logk) / dlogk;

    const double logD= (1 - u) * PowerTable[binlow].logD + u * 
        PowerTable[binhigh].logD;

    const double Delta2 = pow(10.0, logD);

    double P = PowerNorm * Delta2 / (4.0*M_PI*k*k*k);

    // printf("%le %le\n", k, P);

    return P;
}



double sigma2_int(double k, void *param)
{
    double kr, kr3, kr2, w, x;
    void ** vparam = param;
    double (*func )(double, void *) = vparam[0];
    double r_tophat = *(double*) vparam[1];
    void * data = vparam[2];
    kr = r_tophat * k;
    kr2 = kr * kr;
    kr3 = kr2 * kr;

    if(kr < 1e-8)
        return 0;

    w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
    x = 4 * M_PI * k * k * w * w * func(k, data);

    return x;
}

double TopHatSigma2(double R, double (*func)(double, void*), void * param) {
    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);
    void * data [] = {
        (void*) func,
        (void*) &R,
        (void*) param,
    };

    F.function = &sigma2_int;
    F.params = data;
    void * handler = gsl_set_error_handler_off ();
    gsl_integration_qag(&F, 0, 500.0 * 1 / R,
            0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
    // high precision caused error
    gsl_integration_workspace_free(workspace);
    gsl_set_error_handler (handler);

    return result;

    // note: 500/R is here chosen as (effectively) infinity integration boundary
}
