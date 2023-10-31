"""This module creates a neutrino transfer function file using Classylss,
a python interface to the CLASS Boltzmann code.

See:
http://classylss.readthedocs.io/en/stable/
http://class-code.net/

Cite CLASS paper:
 D. Blas, J. Lesgourgues, T. Tram, arXiv:1104.2933 [astro-ph.CO], JCAP 1107 (2011) 034
"""

import math
import os.path
import argparse
import numpy as np
#import classylss
import classylss.binding as CLASS
import configobj

def make_class_power(paramfile, tfile):
    """Main routine: parses a parameter file and makes a transfer function ratio.
    Will not over-write power spectra if already present.
    Options are loaded from the fastpm parameter file.
    Supported:
        - Omega_fld and DE parameters.
        - Massive neutrinos.
        - Using Sigma8 to set the power spectrum scale.
        - Different transfer functions.

    We use class velocity transfer functions to have accurate initial conditions
    even on superhorizon scales, and to properly support multiple species.
    The alternative is to use rescaling.

    Not supported:
        - Warm dark matter power spectra.
        - Rescaling with different transfer functions."""
    config = configobj.ConfigObj(infile=paramfile, file_error=True)
    #Input sanitisation

    #Precision
    pre_params = {'tol_background_integration': 1e-9, 'tol_perturb_integration' : 1.e-7, 'tol_thermo_integration':1.e-5, 'k_per_decade_for_pk': 50,
                  'k_bao_width': 8, 'k_per_decade_for_bao':  200, 'neglect_CMB_sources_below_visibility' : 1.e-30,
                  'transfer_neglect_late_source': 3000., 'l_max_g' : 50, 'l_max_ur':150}

    #Important! Densities are in synchronous gauge!
    pre_params['gauge'] = 'synchronous'

    h0 = config['h']
    omeganu = sum(config['m_ncdm'])/93.14/h0**2
    ocdm = config['Omega_m']- omeganu
    gparams = {'h':config['h'], 'Omega_cdm':ocdm, 'Omega_b':0.0486, 'Omega_k':0, 'n_s': 0.97}
    #Set up massive neutrinos
    gparams['m_ncdm'] = '%.8f,%.8f,%.8f' % (config['MNue'], config['MNum'], config['MNut'])
    gparams['N_ncdm'] = config['N_nu']
    gparams['N_ur'] = config['N_eff'] - config['N_nu']
    #Neutrino accuracy: Default pk_ref.pre has tol_ncdm_* = 1e-10,
    #which takes 45 minutes (!) on my laptop.
    #tol_ncdm_* = 1e-8 takes 20 minutes and is machine-accurate.
    #Default parameters are fast but off by 2%.
    #I chose 1e-4, which takes 20 minutes and is accurate to 1e-4
    # gparams['tol_ncdm_newtonian'] = 1e-4
    # gparams['tol_ncdm_synchronous'] = 1e-4
    # gparams['tol_ncdm_bg'] = 1e-10
    # gparams['l_max_ncdm'] = 50
    #For accurate, but very slow, P_nu, set ncdm_fluid_approximation = 3
    #CAMB does this better.
    #gparams['ncdm_fluid_approximation'] = 2
    #gparams['ncdm_fluid_trigger_tau_over_tau_k'] = 10000.
    #Power spectrum amplitude: sigma8 is ignored by classylss.
    gparams['A_s'] = 2.130624e-9
    pre_params.update(gparams)
    redshift = config['z_i']
    outputs = redshift
    #Pass options for the power spectrum
    boxmpc = config['boxsize']
    maxk = max(10, 2*math.pi/boxmpc*config['nc']*4)
    #CLASS needs the first redshift to be relatively high for some internal interpolation reasons
    maxz = max(1 + np.max(outputs), 99)
    powerparams = {'output': 'dTk vTk mPk', 'P_k_max_h/Mpc' : maxk, "z_max_pk" : maxz,'z_pk': outputs, 'extra metric transfer functions': 'y'}
    pre_params.update(powerparams)
    #Print the class parameters to terminal in a format
    #readable by the command line class.
    #Make the power spectra module
    engine = CLASS.ClassEngine(pre_params)
    powspec = CLASS.Spectra(engine)
    print("sigma_8(z=0) = ", powspec.sigma8, "A_s = ",powspec.A_s)
    #Get and save the transfer functions if needed
    trans = powspec.get_transfer(z=redshift)
    if os.path.exists(tfile):
        raise IOError("Refusing to write to existing file: ",tfile)
    #Get and save the transfer function ratio
    tcb = trans[2] + trans[3]
    tnu = trans[5] + trans[6] + trans[7]
    save_transfer(np.hstack([trans['k'], tnu/tcb]), tfile)

def save_transfer(transfer, transferfile):
    """Save a transfer function ratio between neutrinos and CDM."""
    header=r"""Transfer function T_i(k) for the ratio between neutrinos and CDM.
    This is T_i(k) / T_{cb}(k) where T_i is the total neutrinos and cb is the cdm and baryon combined transfer function.
    """
    #This format matches the default output by CLASS command line.
    np.savetxt(transferfile, transfer, header=header)

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('paramfile', type=str, help='fastpm paramfile')
    parser.add_argument('transfer', type=str, help='Output transfer function file')
    args = parser.parse_args()
    make_class_power(args.paramfile, tfile=args.transfer)
