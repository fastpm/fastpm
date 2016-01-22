-- parameter file
------ Size of the simulation -------- 

-- RunPB size
--nc = 2048
--boxsize = 1380.0

-- For Testing
nc = 128
boxsize = 384.0

-- Broken, always use 1
nrealization= 1

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 10)

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
--readic="/scratch1/scratchdirs/yfeng1/PB01/tpmsph_ic.bin" 

-------- Approximation Method ---------------
force_mode = "cola"

pm_nc_factor = {2,   2,   2, }            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm =    {0, 0.2, 0.5, }            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= "cola/snp"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = "cola/powerspec"

