/**\file
 * Contains calculations for the Fourier-space semi-linear neutrino method
 * described in Ali-Haimoud and Bird 2012.
 * delta_tot_table stores the state of the integrator, which includes the matter power spectrum over all past time.
 * This file contains routines for manipulating this structure; updating it by computing a new neutrino power spectrum,
 * from the non-linear CDM power.
 */

#include <math.h>
#include <string.h>
#include <stdint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>
#include <mpi.h>
#include <bigfile-mpi.h>

#include <fastpm/logging.h>
#include <fastpm/neutrinos_lra.h>

// #include "powerspectrum.h"
// #include "physconst.h"

/** Floating point accuracy*/
#define FLOAT_ACC   1e-6
/** Number of bins in integrations*/
#define GSL_VAL 400

#define LIGHT 9.715614e-15     // in units: [ h * (Mpc/h) * s^-1 ]
#define  HUBBLE          3.2407789e-18	/* in h/sec */
/** Ratio between the massless neutrino temperature and the CMB temperature.
 * Note there is a slight correction from 4/11
 * due to the neutrinos being slightly coupled at e+- annihilation.
 * See Mangano et al 2005 (hep-ph/0506164)
 * We use the CLASS default value, chosen so that omega_nu = m_nu / 93.14 h^2
 * At time of writing this is T_nu / T_gamma = 0.71611.
 * See https://github.com/lesgourg/class_public/blob/master/explanatory.ini
 * !!! For FastPM we actually use the Gamma_nu function from cosmology.c !!!
 */
//#define TNUCMB     (pow(4/11.,1/3.)*1.00328)
/** The Boltzmann constant in units of eV/K*/
#define BOLEVK 8.617333262145e-5

/** Allocates memory for delta_tot_table.
 * @param nk_in Number of bins stored in each power spectrum.
 * @param TimeTransfer Scale factor of the transfer functions.
 * @param TimeMax Final scale factor up to which we will need memory.
 * @param CP Cosmology parameters.*/
void init_neutrinos_lra(_delta_tot_table * d_tot, const int nk_in, const double TimeMax);

/* Get total neutrino mass, wrapper*/
double get_omega_nu(double Time, FastPMCosmology * cosmo)
{
    double hubble = HubbleEa(Time, cosmo);
    return Omega_ncdmTimesHubbleEaSq(Time, cosmo)/pow(hubble,2);
}

/* Get neutrino mass for single species, wrapper*/
double omega_nu_single(double Time, int i, FastPMCosmology * cosmo)
{
    double hubble = HubbleEa(Time, cosmo);
    return Omega_ncdm_iTimesHubbleEaSq(Time, i, cosmo)/pow(hubble, 2);
}

/** Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table * const d_tot, const double a, const double delta_cdm_curr[], const double delta_nu_curr[], const int overwrite);

/** Main function: given tables of wavenumbers, total delta at Na earlier times (< = a),
 * and initial conditions for neutrinos, computes the current delta_nu.
 * @param d_tot Initialised structure for storing total matter density.
 * @param a Current scale factor.
 * @param delta_nu_curr Pointer to array to store square root of neutrino power spectrum. Main output.
 * @param mnu Neutrino mass in eV.*/
void get_delta_nu(FastPMCosmology * CP, const _delta_tot_table * const d_tot, const double a, double delta_nu_curr[], const double mnu);

/** Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(FastPMCosmology * CP, const _delta_tot_table * const d_tot, const double a, double delta_nu_curr[]);

/** Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute*/
static double specialJ(const double x);

/** Free-streaming length (times Mnu/k_BT_nu, which is dimensionless) for a non-relativistic
particle of momentum q = T0, from scale factor ai to af.
Arguments:
@param logai log of initial scale factor
@param logaf log of final scale factor
@param mnu Neutrino mass in eV
@param light speed of light in internal length units.
@returns free-streaming length in Unit_Length/Unit_Time (same units as light parameter).
*/
double fslength(FastPMCosmology * CP, const double logai, const double logaf, const double light);

/** Combine the CDM and neutrino power spectra together to get the total power.
 * OmegaNua3 = OmegaNu(a) * a^3
 * Omeganonu = Omega0 - OmegaNu(1)
 * Omeganu1 = OmegaNu(1) */
static inline double get_delta_tot(const double delta_nu_curr, const double delta_cdm_curr, const double OmegaNua3, const double Omeganonu)
{
    const double fcdm = 1 - OmegaNua3/Omeganonu;
    return fcdm * (delta_cdm_curr + delta_nu_curr * OmegaNua3/Omeganonu);
}


/*Structure which holds the neutrino state*/
_delta_tot_table delta_tot_table;
FastPMFuncK t_init[1];

void load_transfer_data(const double TimeTransfer, FastPMFuncK *t_init_in)
{
    delta_tot_table.TimeTransfer = TimeTransfer;
    /** Structure to store the initial transfer functions from CAMB.
    * We store transfer functions because we want to use the
    * CDM + Baryon total matter power spectrum from the
    * first timestep of fastpm, so that possible scattering
    * in the initial conditions is included in the neutrino and radiation components. */
    /*This is T_nu / (T_not-nu), where T_not-nu is a weighted average of T_cdm and T_baryon*/
    t_init->size = t_init_in->size;
    t_init->k = malloc(t_init_in->size * sizeof(t_init->k[0]));
    t_init->f = malloc(t_init_in->size * sizeof(t_init->k[0]));
    /* k should be in log10 for below*/
    int ik;
    for(ik=0; ik < t_init->size; ik++) {
        t_init->k[ik] = log10(t_init_in->k[ik]);
        t_init->f[ik] = t_init_in->f[ik];
    }
}

/* Constructor. transfer_init_tabulate must be called before this function.
 * Initialises delta_tot (including from a file) and delta_nu_init from the transfer functions.
 * read_all_nu_state must be called before this if you want reloading from a snapshot to work
 * Note delta_cdm_curr includes baryons, and is only used if not resuming.*/
static void delta_tot_first_init(_delta_tot_table * const d_tot, const int nk_in, const double wavenum[], const double delta_cdm_curr[], const double Time, FastPMCosmology * cosmo)
{
    int ik;
    d_tot->nk=nk_in;
    /* Allocate memory: hard-coded upper limit of z=0*/
    init_neutrinos_lra(d_tot, nk_in, 1.0);
    /*Set the prefactor for delta_nu, and the units system, which is Mpc/s.*/
    d_tot->light = LIGHT;
    d_tot->delta_nu_prefac = 1.5 * cosmo->Omega_m * HUBBLE * HUBBLE /d_tot->light;
    /*Matter fraction excluding neutrinos*/
    d_tot->cosmo = cosmo;
    d_tot->Omeganonu = cosmo->Omega_m - get_omega_nu(1, d_tot->cosmo);
    const double OmegaNua3=get_omega_nu(d_tot->TimeTransfer, d_tot->cosmo)*pow(d_tot->TimeTransfer,3);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp * spline = NULL;
    if(t_init->size > 0) {
        if(t_init->size > 2)
            spline = gsl_interp_alloc(gsl_interp_cspline,t_init->size);
        else
            spline = gsl_interp_alloc(gsl_interp_linear,t_init->size);
        gsl_interp_init(spline,t_init->k,t_init->f, t_init->size);
        /*Check we have a long enough power table: power tables are in log_10*/
        if(log10(wavenum[d_tot->nk-1]) > t_init->k[t_init->size-1])
            fastpm_raise(2,"Want k = %g but maximum in CLASS table is %g\n",wavenum[d_tot->nk-1], pow(10, t_init->k[t_init->size-1]));
    }
    for(ik=0;ik<d_tot->nk;ik++) {
            /* f contains T_nu / T_cdm.*/
            double T_nubyT_nonu = 1;
            if(t_init->size > 0 && wavenum[ik] > 0)
                T_nubyT_nonu = gsl_interp_eval(spline,t_init->k,t_init->f, log10(wavenum[ik]),acc);
            /*Initialise delta_nu_init to use the first timestep's delta_cdm_curr
             * so that it includes potential Rayleigh scattering. */
            d_tot->delta_nu_init[ik] = delta_cdm_curr[ik]*T_nubyT_nonu;
            /*Initialise the first delta_tot*/
            d_tot->delta_tot[ik][0] = get_delta_tot(d_tot->delta_nu_init[ik], delta_cdm_curr[ik], OmegaNua3, d_tot->Omeganonu);
            /*Set up the wavenumber array*/
            d_tot->wavenum[ik] = wavenum[ik];
    }
    gsl_interp_accel_free(acc);
    gsl_interp_free(spline);
    /* Free initial memory*/
    free(t_init->f);
    free(t_init->k);
    /*If we are not restarting, make sure we set the scale factor*/
    d_tot->scalefact[0]=log(Time);
    d_tot->ia=1;
    return;
}

void delta_nu_from_power(nu_lra_power * nupow, FastPMFuncK* ps, FastPMCosmology * CP, const double Time)
{
    int i;
    /*This is done on the first timestep: we need nk_nonzero for it to work.*/
    if(!delta_tot_table.delta_tot_init_done) {
        if(delta_tot_table.ia == 0) {
            /* Compute delta_nu from the transfer functions if first entry.*/
            delta_tot_first_init(&delta_tot_table, ps->size, ps->k, ps->f, Time, CP);
        }
        /*Initialise the first delta_nu*/
        get_delta_nu_combined(CP, &delta_tot_table, exp(delta_tot_table.scalefact[delta_tot_table.ia-1]), delta_tot_table.delta_nu_last);
        delta_tot_table.delta_tot_init_done = 1;
    }
    for(i = 0; i < ps->size; i++)
        nupow->logknu[i] = log(ps->k[i]);

    double * Power_in = ps->f;
    /* Rebin the input power if necessary*/
    if(delta_tot_table.nk != ps->size) {
        Power_in = (double *) malloc(delta_tot_table.nk * sizeof(double));
        double * logPower = (double *) malloc(ps->size * sizeof(double));
        for(i = 0; i < ps->size; i++)
            logPower[i] = log(ps->f[i]);
        gsl_interp * pkint = gsl_interp_alloc(gsl_interp_linear, ps->size);
        gsl_interp_init(pkint, nupow->logknu, logPower, ps->size);
        gsl_interp_accel * pkacc = gsl_interp_accel_alloc();
        for(i = 0; i < delta_tot_table.nk; i++) {
            double logk = log(delta_tot_table.wavenum[i]);
            if(pkint->xmax < logk || pkint->xmin > logk)
                Power_in[i] = delta_tot_table.delta_tot[i][delta_tot_table.ia-1];
            else
                Power_in[i] = exp(gsl_interp_eval(pkint, nupow->logknu, logPower, logk, pkacc));
        }
        free(logPower);
        gsl_interp_accel_free(pkacc);
        gsl_interp_free(pkint);
    }

    /* If we get called twice with the same scale factor, do nothing: delta_nu
     * already stores the neutrino power from the current timestep.*/
    if(log(Time)-delta_tot_table.scalefact[delta_tot_table.ia-1] > FLOAT_ACC) {
        /*We need some estimate for delta_tot(current time) to obtain delta_nu(current time).
            Even though delta_tot(current time) is not directly used (the integrand vanishes at a = a(current)),
            it is indeed needed for interpolation */
        /*It was checked that using delta_tot(current time) = delta_cdm(current time) leads to no more than 2%
          error on delta_nu (and moreover for large k). Using the last timestep's delta_nu decreases error even more.
          So we only need one step. */
        /*This increments the number of stored spectra, although the last one is not yet final.*/
        update_delta_tot(&delta_tot_table, Time, Power_in, delta_tot_table.delta_nu_last, 0);
        /*Get the new delta_nu_curr*/
        get_delta_nu_combined(CP, &delta_tot_table, Time, delta_tot_table.delta_nu_last);
        /* Decide whether we save the current time or not */
        if (Time > exp(delta_tot_table.scalefact[delta_tot_table.ia-2]) + 0.009) {
            /* If so update delta_tot(a) correctly, overwriting current power spectrum */
            update_delta_tot(&delta_tot_table, Time, Power_in, delta_tot_table.delta_nu_last, 1);
        }
        /*Otherwise discard the last powerspectrum*/
        else
            delta_tot_table.ia--;

        fastpm_info("Done getting neutrino power: nk = %d, k = %g, delta_nu = %g, delta_cdm = %g,\n", delta_tot_table.nk, delta_tot_table.wavenum[1], delta_tot_table.delta_nu_last[1], Power_in[1]);
        /*kspace_prefac = M_nu (analytic) / M_particles */
        const double OmegaNu = get_omega_nu(Time, delta_tot_table.cosmo);
        /* Omega0 - Omega in neutrinos = Omega in particles*/
        nupow->nu_prefac = OmegaNu/(delta_tot_table.Omeganonu/pow(Time,3));
    }
    double * delta_nu_ratio = (double *) malloc(delta_tot_table.nk * sizeof(double));
    double * logwavenum = (double *) malloc(delta_tot_table.nk * sizeof(double));
    gsl_interp * pkint = gsl_interp_alloc(gsl_interp_linear, delta_tot_table.nk);
    gsl_interp_accel * pkacc = gsl_interp_accel_alloc();
    /*We want to interpolate in log space*/
    for(i=0; i < delta_tot_table.nk; i++) {
        /*Enforce positivity for sanity reasons*/
        if(delta_tot_table.delta_nu_last[i] < 0)
            delta_tot_table.delta_nu_last[i] = 0;
        delta_nu_ratio[i] = delta_tot_table.delta_nu_last[i]/ Power_in[i];
        logwavenum[i] = log(delta_tot_table.wavenum[i]);
    }
    if(delta_tot_table.nk != ps->size)
        free(Power_in);
    gsl_interp_init(pkint, logwavenum, delta_nu_ratio, delta_tot_table.nk);

    /*We want to interpolate in log space*/
    for(i=0; i < ps->size; i++) {
        if(ps->size == delta_tot_table.nk)
            nupow->delta_nu_ratio[i] = delta_nu_ratio[i];
        else {
            double logk = nupow->logknu[i];
            if(logk > pkint->xmax)
                logk = pkint->xmax;
            nupow->delta_nu_ratio[i] = gsl_interp_eval(pkint, logwavenum, delta_nu_ratio, logk, pkacc);
        }
    }

    gsl_interp_accel_free(pkacc);
    gsl_interp_free(pkint);
    free(logwavenum);
    free(delta_nu_ratio);
    // fastpm_info("Neutrino ratio: nk = %d, k = %g, delta_nu = %g,\n", ps->size, nupow->logknu[1], nupow->delta_nu_ratio[1]);
}

/*Save the neutrino power spectrum to a file*/

void powerspectrum_nu_save(FastPMPowerSpectrum * ps, char powerspectrum_file[])
{
    /* substitute the last neutrino power spectrum */
    int i;
    /* The k bins are the same because this is constructed from the same power spectrum*/
    for(i = 0; i < ps->base.size; i++)
        ps->base.f[i] = pow(delta_tot_table.delta_nu_last[i], 2);

    fastpm_powerspectrum_write(ps, powerspectrum_file, 1.0);
}

/* save a block to disk */
static void _save_block(BigFile * bf, const char * blockname, BigArray * array)
{
    BigBlock bb;
    BigBlockPtr ptr;
    size_t size = array->dims[0];
    int NumFiles = 1;

    /*Do not write empty files*/
    if(size == 0) {
        NumFiles = 0;
    }
    /* create the block */
    /* dims[1] is the number of members per item */
    if(0 != big_file_mpi_create_block(bf, &bb, blockname, array->dtype, array->dims[1], NumFiles, size, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to create block at %s:%s\n", blockname,
                    big_file_get_error_message());
    }
    if(0 != big_block_seek(&bb, &ptr, 0)) {
        fastpm_raise(-1, "Failed to seek:%s\n", big_file_get_error_message());
    }
    if(0 != big_block_mpi_write(&bb, &ptr, array, NumFiles, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to write :%s\n", big_file_get_error_message());
    }
    if(0 != big_block_mpi_close(&bb, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to close block at %s:%s\n", blockname,
                big_file_get_error_message());
    }
}

void ncdm_lr_save_neutrinos(BigFile * bf, int ThisTask)
{
    if(!delta_tot_table.delta_tot_init_done)
        return;
    double * scalefact = delta_tot_table.scalefact;
    size_t nk = delta_tot_table.nk, ia = delta_tot_table.ia;
    size_t ik, i;
    double * delta_tot = (double *) malloc(nk * ia * sizeof(double));
    /*Save a flat memory block*/
    for(ik=0;ik< nk;ik++)
        for(i=0;i< ia;i++)
            delta_tot[ik*ia+i] = delta_tot_table.delta_tot[ik][i];

    BigBlock bn;
    if(0 != big_file_mpi_create_block(bf, &bn, "Neutrino", NULL, 0, 0, 0, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to create block at %s:%s\n", "Neutrino",
                big_file_get_error_message());
    }
    if ( (0 != big_block_set_attr(&bn, "Nscale", &ia, "u8", 1)) ||
       (0 != big_block_set_attr(&bn, "scalefact", scalefact, "f8", ia)) ||
        (0 != big_block_set_attr(&bn, "Nkval", &nk, "u8", 1)) ) {
        fastpm_raise(-1, "Failed to write neutrino attributes %s\n",
                    big_file_get_error_message());
    }
    if(0 != big_block_mpi_close(&bn, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to close block %s\n",
                    big_file_get_error_message());
    }
    BigArray deltas = {0};
    size_t dims[2] = {0 , ia};
    /*The neutrino state is shared between all processors,
     *so only write on master task*/
    if(ThisTask == 0) {
        dims[0] = nk;
    }
    ptrdiff_t strides[2] = {(ptrdiff_t) (sizeof(double) * ia), (ptrdiff_t) sizeof(double)};
    big_array_init(&deltas, delta_tot, "=f8", 2, dims, strides);

    _save_block(bf, "Neutrino/Deltas", &deltas);
    free(delta_tot);
    /*Now write the initial neutrino power*/
    BigArray delta_nu = {0};
    dims[1] = 1;
    strides[0] = sizeof(double);
    big_array_init(&delta_nu, delta_tot_table.delta_nu_init, "=f8", 2, dims, strides);
    _save_block(bf, "Neutrino/DeltaNuInit", &delta_nu);
    /*Now write the k values*/
    BigArray kvalue = {0};
    big_array_init(&kvalue, delta_tot_table.wavenum, "=f8", 2, dims, strides);
    _save_block(bf, "Neutrino/kvalue", &kvalue);
}

/* read a block from disk, spread the values to memory with setters  */
static int _read_block(BigFile * bf, const char * blockname, BigArray * array) {
    BigBlock bb;
    BigBlockPtr ptr;

    /* open the block */
    if(0 != big_file_mpi_open_block(bf, &bb, blockname, MPI_COMM_WORLD)) {
            fastpm_raise(-1, "Failed to open block at %s:%s\n", blockname, big_file_get_error_message());
    }
    if(0 != big_block_seek(&bb, &ptr, 0)) {
            fastpm_raise(-1, "Failed to seek block %s: %s\n", blockname, big_file_get_error_message());
    }
    if(0 != big_block_mpi_read(&bb, &ptr, array, 1, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to read from block %s: %s\n", blockname, big_file_get_error_message());
    }
    if(0 != big_block_mpi_close(&bb, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to close block at %s:%s\n", blockname,
                    big_file_get_error_message());
    }
    return 0;
}


/*Read the neutrino data from the snapshot*/
int ncdm_lr_read_neutrinos(BigFile * bf, int ThisTask)
{
    size_t nk, ia, ik, i;
    BigBlock bn;
    if(0 != big_file_mpi_open_block(bf, &bn, "Neutrino", MPI_COMM_WORLD)) {
        return 1;
    }
    if(
    (0 != big_block_get_attr(&bn, "Nscale", &ia, "u8", 1)) ||
    (0 != big_block_get_attr(&bn, "Nkval", &nk, "u8", 1))) {
        fastpm_raise(-1, "Failed to read attr: %s\n",
                    big_file_get_error_message());
    }
    double *delta_tot = (double *) malloc(ia*nk*sizeof(double));
    /*Allocate list of scale factors, and space for delta_tot, in one operation.*/
    if(0 != big_block_get_attr(&bn, "scalefact", delta_tot_table.scalefact, "f8", ia))
        fastpm_raise(-1, "Failed to read attr: %s\n", big_file_get_error_message());
    if(0 != big_block_mpi_close(&bn, MPI_COMM_WORLD)) {
        fastpm_raise(-1, "Failed to close block %s\n",
                    big_file_get_error_message());
    }
    BigArray deltas = {0};
    size_t dims[2] = {0, ia};
    ptrdiff_t strides[2] = {(ptrdiff_t) (sizeof(double)*ia), (ptrdiff_t)sizeof(double)};
    /*The neutrino state is shared between all processors,
     *so only read on master task and broadcast*/
    if(ThisTask == 0) {
        dims[0] = nk;
    }
    big_array_init(&deltas, delta_tot, "=f8", 2, dims, strides);
    _read_block(bf, "Neutrino/Deltas", &deltas);
    if(nk > 1Lu*delta_tot_table.nk_allocated || ia > 1Lu*delta_tot_table.namax)
        fastpm_raise(-1, "Allocated nk %d na %d for neutrino power but need nk %d na %d\n", delta_tot_table.nk_allocated, delta_tot_table.namax, nk, ia);
    /*Save a flat memory block*/
    for(ik=0;ik<nk;ik++)
        for(i=0;i<ia;i++)
            delta_tot_table.delta_tot[ik][i] = delta_tot[ik*ia+i];
    delta_tot_table.nk = nk;
    delta_tot_table.ia = ia;
    free(delta_tot);
    /* Read the initial delta_nu. This is basically zero anyway,
     * so for backwards compatibility do not require it*/
    BigArray delta_nu = {0};
    dims[1] = 1;
    strides[0] = sizeof(double);
    memset(delta_tot_table.delta_nu_init, 0, delta_tot_table.nk);
    big_array_init(&delta_nu, delta_tot_table.delta_nu_init, "=f8", 2, dims, strides);
    _read_block(bf, "Neutrino/DeltaNuInit", &delta_nu);
    /* Read the k values*/
    BigArray kvalue = {0};
    memset(delta_tot_table.wavenum, 0, delta_tot_table.nk);
    big_array_init(&kvalue, delta_tot_table.wavenum, "=f8", 2, dims, strides);
    _read_block(bf, "Neutrino/kvalue", &kvalue);

    /*Broadcast the arrays.*/
    MPI_Bcast(&(delta_tot_table.ia), 1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&(delta_tot_table.nk), 1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(delta_tot_table.delta_nu_init,delta_tot_table.nk,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(delta_tot_table.wavenum,delta_tot_table.nk,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(delta_tot_table.ia > 0) {
        /*Broadcast data for scalefact and delta_tot, Delta_tot is allocated as the same block of memory as scalefact.
          Not all this memory will actually have been used, but it is easiest to bcast all of it.*/
        MPI_Bcast(delta_tot_table.scalefact,delta_tot_table.namax*(delta_tot_table.nk+1),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    return 0;
}

/*Allocate memory for delta_tot_table. This is separate from delta_tot_init because we need to allocate memory
 * before we have the information needed to initialise it*/
void init_neutrinos_lra(_delta_tot_table *d_tot, const int nk_in, const double TimeMax)
{
   int count;
   /*Memory allocations need to be done on all processors*/
   d_tot->nk_allocated=nk_in;
   d_tot->nk=nk_in;
   /* Allocate memory for delta_tot here, so that we can have further memory allocated and freed
    * before delta_tot_init is called. The number nk here should be larger than the actual value needed.*/
   /*Allocate pointers to each k-vector*/
   d_tot->namax=ceil((TimeMax-d_tot->TimeTransfer)/0.008)+4;
   d_tot->ia=0;
   d_tot->delta_tot =(double **) malloc(nk_in*sizeof(double *));
   /*Allocate list of scale factors, and space for delta_tot, in one operation.*/
   d_tot->scalefact = (double *) malloc(d_tot->namax*(nk_in+1)*sizeof(double));
   /*Allocate space for delta_nu and wavenumbers*/
   d_tot->delta_nu_last = (double *) malloc(sizeof(double) * 3 * nk_in);
   d_tot->delta_nu_init = d_tot->delta_nu_last + nk_in;
   d_tot->wavenum = d_tot->delta_nu_init + nk_in;
   /*Allocate actual data. Note that this means data can be accessed either as:
    * delta_tot[k][a] OR as
    * delta_tot[0][a+k*namax] */
   d_tot->delta_tot[0] = d_tot->scalefact+d_tot->namax;
   for(count=1; count< nk_in; count++)
        d_tot->delta_tot[count] = d_tot->delta_tot[0] + count*d_tot->namax;
   /*Allocate space for the initial neutrino power spectrum*/
}

/*Begin functions that do the actual computation of the neutrino power spectra.
 * The algorithms executed are explained in Ali-Haimoud & Bird 2012 and Bird, Ali-Haimoud, Feng & Liu 2018
 * arXiv:1209.0461 and arXiv:1803.09854.
 * This is a Fourier-space linear response method for computing neutrino overdensities from CDM overdensities.*/

/*Function which wraps three get_delta_nu calls to get delta_nu three times,
 * so that the final value is for all neutrino species*/
void get_delta_nu_combined(FastPMCosmology * CP, const _delta_tot_table * const d_tot, const double a, double delta_nu_curr[])
{
    const double Omega_nu_tot=get_omega_nu(a, d_tot->cosmo);
    int mi;
    /*Initialise delta_nu_curr*/
    memset(delta_nu_curr, 0, d_tot->nk*sizeof(double));
    /*Get each neutrinos species and density separately and add them to the total.
     * Neglect perturbations in massless neutrinos.*/
    for(mi=0; mi<d_tot->cosmo->N_ncdm; mi++) {
        int ik;
        double * delta_nu_single = (double *) malloc(sizeof(double) * d_tot->nk);
        const double omeganu = omega_nu_single(a, mi, d_tot->cosmo);
        get_delta_nu(CP, d_tot, a, delta_nu_single,d_tot->cosmo->m_ncdm[mi]);
        for(ik=0; ik<d_tot->nk; ik++)
            delta_nu_curr[ik]+=delta_nu_single[ik]*omeganu/Omega_nu_tot;
        free(delta_nu_single);
    }
    return;
}

/*Update the last value of delta_tot in the table with a new value computed
 from the given delta_cdm_curr and delta_nu_curr.
 If overwrite is true, overwrite the existing final entry.*/
void update_delta_tot(_delta_tot_table * const d_tot, const double a, const double delta_cdm_curr[], const double delta_nu_curr[], const int overwrite)
{
  const double OmegaNua3 = get_omega_nu(a, d_tot->cosmo)*pow(a,3);
  int ik;
  if(!overwrite)
    d_tot->ia++;
  /*Update the scale factor*/
  d_tot->scalefact[d_tot->ia-1] = log(a);
  /* Update delta_tot(a)*/
  for (ik = 0; ik < d_tot->nk; ik++){
    d_tot->delta_tot[ik][d_tot->ia-1] = get_delta_tot(delta_nu_curr[ik], delta_cdm_curr[ik], OmegaNua3, d_tot->Omeganonu);
  }
}

/*Kernel function for the fslength integration*/
double fslength_int(const double loga, void *params)
{
    FastPMCosmology * CP = (FastPMCosmology *) params;
    /*This should be M_nu / k_B T_nu (which is dimensionless)*/
    const double a = exp(loga);
    return 1./a/(a*HubbleEa(a, CP));
}

/******************************************************************************************************
Free-streaming length (times Mnu/k_BT_nu, which is dimensionless) for a non-relativistic
particle of momentum q = T0, from scale factor ai to af.
Arguments:
logai - log of initial scale factor
logaf - log of final scale factor
light - speed of light in internal units.
Result is in Unit_Length/Unit_Time.
******************************************************************************************************/
double fslength(FastPMCosmology * CP, const double logai, const double logaf, const double light)
{
  double abserr;
  double fslength_val;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
  F.function = &fslength_int;
  F.params = CP;
  if(logai >= logaf)
      return 0;
  gsl_integration_qag (&F, logai, logaf, 0, 1e-6,GSL_VAL,6,w,&(fslength_val), &abserr);
  gsl_integration_workspace_free (w);
  return light*fslength_val;
}

/**************************************************************************************************
Fit to the special function J(x) that is accurate to better than 3% relative and 0.07% absolute
    J(x) = Integrate[(Sin[q*x]/(q*x))*(q^2/(Exp[q] + 1)), {q, 0, Infinity}]
    and J(0) = 1.
    Mathematica gives this in terms of the PolyGamma function:
   (PolyGamma[1, 1/2 - i x/2] - PolyGamma[1, 1 - i x/2] -    PolyGamma[1, 1/2 + i x/2] +
   PolyGamma[1, 1 + i x/2])/(12 x Zeta[3]), which we could evaluate exactly if we wanted to.
***************************************************************************************************/
static inline double specialJ(const double x)
{

  double x2, x4, x8;
  if (x <= 0.)
      return 1.;
  x2 = x*x;
  x4 = x2*x2;
  x8 = x4*x4;

  return (1.+ 0.0168 * x2 + 0.0407* x4)/(1. + 2.1734 * x2 + 1.6787 * exp(4.1811*log(x)) +  0.1467 * x8);
}

/**A structure for the parameters for the below integration kernel*/
struct _delta_nu_int_params
{
    /**Current wavenumber*/
    double k;
    /**Neutrino mass divided by k_B T_nu*/
    double mnubykT;
    gsl_interp_accel *acc;
    gsl_interp *spline;
    FastPMCosmology * CP;
    /**Precomputed free-streaming lengths*/
    gsl_interp_accel *fs_acc;
    gsl_interp *fs_spline;
    double * fslengths;
    double * fsscales;
    /**Make sure this is at the same k as above*/
    double * delta_tot;
    double * scale;
};
typedef struct _delta_nu_int_params delta_nu_int_params;

/**GSL integration kernel for get_delta_nu*/
double get_delta_nu_int(double logai, void * params)
{
    delta_nu_int_params * p = (delta_nu_int_params *) params;
    double fsl_aia = gsl_interp_eval(p->fs_spline,p->fsscales,p->fslengths,logai,p->fs_acc);
    double delta_tot_at_a = gsl_interp_eval(p->spline,p->scale,p->delta_tot,logai,p->acc);
    double specJ = specialJ(p->k*fsl_aia/p->mnubykT);
    double ai = exp(logai);
    return fsl_aia/(ai*HubbleEa(ai, p->CP)) * specJ * delta_tot_at_a;
}

/*
Main function: given tables of wavenumbers, total delta at Na earlier times (<= a),
and initial conditions for neutrinos, computes the current delta_nu.
Na is the number of currently stored time steps.
*/
void get_delta_nu(FastPMCosmology * CP, const _delta_tot_table * const d_tot, const double a, double delta_nu_curr[],const double mnu)
{
    double fsl_A0a,deriv_prefac;
    int ik;
    /*Number of stored power spectra. This includes the initial guess for the next step*/
    const int Na = d_tot->ia;
    double TNUCMB = Gamma_nu(CP);
    double kBtnu = BOLEVK * TNUCMB * d_tot->cosmo->T_cmb;
    const double mnubykT = mnu / kBtnu;
    /*Tolerated integration error*/
    double relerr = 1e-6;
    //       message(0,"Start get_delta_nu: a=%g Na =%d wavenum[0]=%g delta_tot[0]=%g m_nu=%g\n",a,Na,wavenum[0],d_tot->delta_tot[0][Na-1],mnu);

    fsl_A0a = fslength(CP, log(d_tot->TimeTransfer), log(a),d_tot->light);
    /*Precompute factor used to get delta_nu_init. This assumes that delta ~ a, so delta-dot is roughly 1.*/
    deriv_prefac = d_tot->TimeTransfer*(HubbleEa(d_tot->TimeTransfer, CP)/d_tot->light)* d_tot->TimeTransfer;
    for (ik = 0; ik < d_tot->nk; ik++) {
      /* Initial condition piece, assuming linear evolution of delta with a up to startup redshift */
      /* This assumes that delta ~ a, so delta-dot is roughly 1. */
      /* Also ignores any difference in the transfer functions between species.
       * This will be good if all species have similar masses, or
       * if two species are massless.
       * Also, since at early times the clustering is tiny, it is very unlikely to matter.*/
      /*For zero mass neutrinos just use the initial conditions piece, modulating to zero inside the horizon*/
      const double specJ = specialJ(d_tot->wavenum[ik]*fsl_A0a/(mnubykT > 0 ? mnubykT : 1));
      delta_nu_curr[ik] = specJ*d_tot->delta_nu_init[ik] *(1.+ deriv_prefac*fsl_A0a);
    }
    /*If only one time given, we are still at the initial time*/
    /*If neutrino mass is zero, we are not accurate, just use the initial conditions piece*/
    if(Na > 1 && mnubykT > 0){
        delta_nu_int_params params;
        params.acc = gsl_interp_accel_alloc();
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (GSL_VAL);
        gsl_function F;
        F.function = &get_delta_nu_int;
        F.params=&params;
        /*Use cubic interpolation*/
        if(Na > 2) {
                params.spline=gsl_interp_alloc(gsl_interp_cspline,Na);
        }
        /*Unless we have only two points*/
        else {
                params.spline=gsl_interp_alloc(gsl_interp_linear,Na);
        }
        params.scale=d_tot->scalefact;
        params.mnubykT=mnubykT;
        /* Massively over-sample the free-streaming lengths.
         * Interpolation is least accurate where the free-streaming length -> 0,
         * which is exactly where it doesn't matter, but
         * we still want to be safe. */
        int Nfs = Na*16;
        params.fs_acc = gsl_interp_accel_alloc();
        params.fs_spline=gsl_interp_alloc(gsl_interp_cspline,Nfs);
        params.CP = CP;
        /*Pre-compute the free-streaming lengths, which are scale-independent*/
        double * fslengths = (double *) malloc(Nfs* sizeof(double));
        double * fsscales = (double *) malloc(Nfs* sizeof(double));
        for(ik=0; ik < Nfs; ik++) {
            fsscales[ik] = log(d_tot->TimeTransfer) + ik*(log(a) - log(d_tot->TimeTransfer))/(Nfs-1.);
            fslengths[ik] = fslength(CP, fsscales[ik], log(a),d_tot->light);
        }
        params.fslengths = fslengths;
        params.fsscales = fsscales;

        if(!params.spline || !params.acc || !w || !params.fs_spline || !params.fs_acc || !fslengths || !fsscales)
              fastpm_raise(2016,"Error initialising and allocating memory for gsl interpolator and integrator.\n");

        gsl_interp_init(params.fs_spline,params.fsscales,params.fslengths,Nfs);
        for (ik = 0; ik < d_tot->nk; ik++) {
            double abserr,d_nu_tmp;
            params.k=d_tot->wavenum[ik];
            params.delta_tot=d_tot->delta_tot[ik];
            gsl_interp_init(params.spline,params.scale,params.delta_tot,Na);
            gsl_integration_qag (&F, log(d_tot->TimeTransfer), log(a), 0, relerr,GSL_VAL,6,w,&d_nu_tmp, &abserr);
            delta_nu_curr[ik] += d_tot->delta_nu_prefac * d_nu_tmp;
         }
         gsl_integration_workspace_free (w);
         gsl_interp_free(params.spline);
         gsl_interp_accel_free(params.acc);
         free(fsscales);
         free(fslengths);
    }
//     for(ik=0; ik< 3; ik++)
//         message(0,"k %g d_nu %g\n",wavenum[d_tot->nk/8*ik], delta_nu_curr[d_tot->nk/8*ik]);
    return;
}
