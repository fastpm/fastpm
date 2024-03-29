#ifndef NEUTRINOS_LRA_H
#define NEUTRINOS_LRA_H

#include <fastpm/libfastpm.h>
#include <bigfile.h>

/** Now we want to define a static object to store all previous delta_tot.
 * This object needs a constructor, a few private data members, and a way to be read and written from disk.
 * nk is fixed, delta_tot, scalefact and ia are updated in get_delta_nu_update*/
typedef struct _delta_tot_table {
    /** Number of actually non-zero k values stored in each power spectrum*/
    int nk;
    /** Size of arrays allocated to store power spectra*/
    int nk_allocated;
    /** Maximum number of redshifts to store. Redshifts are stored every delta a = 0.01 */
    int namax;
    /** Number of already "recorded" time steps, i.e. scalefact[0...ia-1] is recorded.
    * Current time corresponds to index ia (but is only recorded if sufficiently far from previous time).
    * Caution: ia here is different from Na in get_delta_nu (Na = ia+1).*/
    int ia;
    /** Prefactor for use in get_delta_nu. Should be 3/2 Omega_m H^2 /c. Units are h^2 / sec / Mpc  */
    double delta_nu_prefac;
    /** Set to unity once the init routine has run.*/
    int delta_tot_init_done;
    /** Pointer to nk arrays of length namax containing the total power spectrum.*/
    double **delta_tot;
    /** Array of length namax containing scale factors at which the power spectrum is stored*/
    double * scalefact;
    /** Pointer to array of length nk storing initial neutrino power spectrum*/
    double * delta_nu_init;
    /** Last-seen neutrino power spectrum*/
    double * delta_nu_last;
    /**Pointer to array storing the effective wavenumbers for the above power spectra*/
    double * wavenum;
    /** Matter density excluding neutrinos*/
    double Omeganonu;
    /** Light speed in internal units, Mpc/s.*/
    double light;
    /** The time at which the simulation starts*/
    double TimeTransfer;
    /* Cosmology factors*/
    FastPMCosmology * cosmo;
} _delta_tot_table;

/* Structure for the computed neutrino data.*/
typedef struct nu_lra_power
{
    double * logknu;
    double * delta_nu_ratio;
    int size;
    double nu_prefac;
    gsl_interp *nu_spline;
} nu_lra_power;

/*Computes delta_nu from a CDM power spectrum.*/
void delta_nu_from_power(nu_lra_power * nupow, FastPMFuncK* ps, FastPMCosmology * CP, const double Time);

/* Load the neutrino transfer function*/
void load_transfer_data(const double TimeTransfer, FastPMFuncK t_init_in[]);

/*These functions save and load neutrino related data from the snapshots*/
void ncdm_lr_save_neutrinos(BigFile * bf, int ThisTask);
int ncdm_lr_read_neutrinos(BigFile * bf, int ThisTask);

/*Save the neutrino power spectrum to a file*/
void powerspectrum_nu_save(FastPMPowerSpectrum * ps, char powerspectrum_file[], double MtotbyMcdm);
#endif
