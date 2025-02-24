-- parameter file
------ Size of the simulation---------- 

-- For Testing
--nc = 1024
nc = 512
--nc = 256
--nc = 128
--boxsize = 1000.0
boxsize = 250.0

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
time_step = linspace(0.01, 1.0, 2)
--time_step = {0.01}

output_redshifts= {99}  -- redshifts of output

-- Cosmology --
Omega_m = 0.3089
h       = 0.6774

------- Initial power spectrum --------
-- See libfastpm/pngaussian.c for documentation
-- Amplitude of primordial power spectrum at pivot scale

scalar_amp = 2.142e-9   -- same as scalar_amp in CAMB

-- Pivot scale k_pivot in 1/Mpc
scalar_pivot = 0.05  -- same as pivot_scalar in CAMB
-- Primordial spectral index n_s
scalar_spectral_index = 0.9667  -- same as scalar_spectral_index in CAMB

sigma8 = 0.8147
--
-- Start with a power spectrum file
-- Linear power spectrum at z=0: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum = "parameter_files/Pk_Planck15_Table4.txt"
--read_powerspectrum = "parameter_files/Planck2015Table4LastColumn_matterpower.dat"

f_nl_type = 'local'

--f_nl = 0.
f_nl = 100.

random_seed = 1
--inverted_ic = true
remove_cosmic_variance = true


-------- Approximation Method ---------------
force_mode = "fastpm"

pm_nc_factor = 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm = 0.2            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "/home/adrian/UNIT_PNG/FastPM_new/output_tests/snaps/fnew_N512_fnl100_Pkcut_L250"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "/home/adrian/UNIT_PNG/FastPM_new/output_tests/pks/pk_fnew_N512_fnl100_Pkcut_L250"
--write_powerspectrum = "/home/adrian/UNIT_PNG/FastPM_new/output_tests/pks/pk_N256_gauss_as_2em9_kcut_100_kovnyq_0.25_new"
write_noisek = "/home/adrian/UNIT_PNG/FastPM_new/output_tests/snaps/noisek"
write_noise = "/home/adrian/UNIT_PNG/FastPM_new/output_tests/snaps/noise"
