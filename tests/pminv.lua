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
time_step = linspace(0.1, 1, 10)

output_redshifts= {0.0}  -- redshifts of output

-- Cosmology --
Omega_m = 0.307494
h       = 0.6774

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
random_seed= 100
inverted_ic = true

-------- Approximation Method ---------------
force_mode = "pm"

pm_nc_factor = {2,   2,   2, }            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm =    {0, 0.2, 0.5, }            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "pminv/fastpm"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "pminv/powerspec"

