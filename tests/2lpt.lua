-- parameter file
------ Size of the simulation -------- 

-- RunPB size
--nc = 2048
--boxsize = 1380.0

-- For Testing
nc = 128
boxsize = 1280.0

-- Broken, always use 1
nrealization= 1

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = {1.0}

output_redshifts= {0.0}  -- redshifts of output

-- Cosmology --
omega_m = 1-0.708
h       = 0.69

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
random_seed= 100
sigma8  = 0.82 -- 0 if the power spectrum is already normalized

-- Alternatively, give a RunPB Initial conditions file
-- readic="ic/snp00100_0.1000.bin" 

-------- Approximation Method ---------------
force_mode = "pm"

pm_nc_factor = {1, }            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm =    {0,}            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "2lpt/snp"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "2lpt/powerspec"

