-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 1024
boxsize = 1890.0

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)

time_step = linspace(0.01, 1.0, 40)
aout = {0.020000, 0.250000, 0.344828, 0.370370, 0.400000, 0.434783, 0.476190, 0.526316, 0.588235, 0.666667, 0.769230, 0.909091, 1.0}  -- redshifts of output

-- Cosmology --
omega_m = 0.31320
h       = 0.67310

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

pm_nc_factor = {{0.0, 1}, {0.01, 2}}

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "/lustre/home/acct-phyzpj/phyzpj-minjioh/FASTPM/fastpm_massive_in_and_out/restart-large/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "/lustre/home/acct-phyzpj/phyzpj-minjioh/FASTPM/fastpm_massive_in_and_out/restart-large/powerspec"

