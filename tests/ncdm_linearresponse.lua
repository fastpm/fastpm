---Example parameter file for a run with ncdm (massive neutrinos)
------ Size of the simulation --------
nc = 64
boxsize = 1000

-------- Time Sequence --------
-- take n_steps_log steps in log(a) from a_i to a_m,
-- then n_steps_lin steps in a from a_m to a_f.
n_steps_log = 5
n_steps_lin = 20

z_i = 99
z_m = 19
a_i = 1. / (1. + z_i)
a_m = 1. / (1. + z_m)
a_f = 1

time_step = loglinspace(a_i, a_m, a_f, n_steps_log, n_steps_lin)

-- Redshifts of output snapshots
output_redshifts = {2, 1, 0}

-------- Cosmology --------
Omega_m = 0.3175
h       = 0.6711
T_cmb   = 2.7255
-- neutrino (ncdm) parameters
N_eff   = 3.046
N_nu    = 3                 -- number of neutrinos species (including massless species)
m_ncdm  = {0.12, 0.06, 0.02}
n_shell = 0
every_ncdm = 0              -- this defines the ratio of the cdm grid number to the ncdm grid number
ncdm_freestreaming = true   -- choose whether to treat ncdm as free-streaming for the growth ODE and Poisson source terms
ncdm_matterlike = false     -- choose whether to approximate ncdm as matter-like in the background
--lra params
ncdm_linearresponse = true
ncdm_transfer_redshift = z_i
ncdm_transfer_nu_file = "lra_trans.txt"

-------- Perturbations --------
-- Input powerspectrum and growth rate
-- We advise using REPS: https://github.com/matteozennaro/reps
-- Power spectrum units: k P(k) in Mpc/h units
-- Must be compatible with the cosmology parameters above
read_powerspectrum = "Pcb.txt"
--read_powerspectrum_ncdm = "Pncdm.txt"
--read_linear_growth_rate = "fcb.txt"
--read_linear_growth_rate_ncdm = "fncdm.txt"

linear_density_redshift = z_i      -- the redshift of the input cdm files
--linear_density_redshift_ncdm = z_i -- the redshift of the input ncdm files

random_seed = 100
particle_fraction = 1.0
--sort_snapshot = false

-------- Approximation Method ---------
force_mode = "fastpm"
pm_nc_factor = {{0.0, 1}, {0.0001, 2}}   -- mesh size
remove_cosmic_variance = true
growth_mode = "ODE"
za = true

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------
-- destination of particle outputs (all particles)
write_snapshot= "ncdm_lra/fastpm"
-- destination of 1d power spectrum (raw)
write_powerspectrum = "ncdm_lra/powerspec"
