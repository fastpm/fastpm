-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 1024

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
--time_step = linspace(0.1, 1, 10)
--time_step = linspace(0.01, 0.02, 10)
time_step = linspace(0.01, 0.02, 10)

output_redshifts= {99.0, 49.}  -- redshifts of output

-- Cosmology --
omega_m = 0.307494
omega_ncdm = 0.001404
h       = 0.6774

--ncdm split
m_ncdm = {0.15}
n_shell = 10
n_side = 1
lvk = false
every_ncdm = 1

-- Start with a linear density field
-- Power spectrum of the linear density field: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
linear_density_redshift = 0.0 -- the redshift of the linear density field.
random_seed= 100
particle_fraction = 1.0
--
-------- Approximation Method ---------------
force_mode = "fastpm"

pm_nc_factor = 2

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "ncdm/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "ncdm/powerspec"
-- write_fof = "ncdm/fastpm"

