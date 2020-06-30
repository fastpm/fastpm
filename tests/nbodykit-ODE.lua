-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 384.0

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 3)

output_redshifts= {9., 0.5, 0.0}  -- redshifts of output

-- Cosmology --
Omega_m = 0.307494
h       = 0.6774

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
kernel_type = "1_4"

growth_mode = "ODE"

pm_nc_factor = {{0.0, 1}, {0.01, 2}}

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_runpb_snapshot= "nbodykit/tpm"
write_snapshot= "nbodykit/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "nbodykit/powerspec"
write_fof = "nbodykit/fastpm"

