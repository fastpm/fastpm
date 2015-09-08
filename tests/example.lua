-- parameter file
------ Size of the simulation -------- 

-- RunPB size
--nc = 2048
--boxsize = 1380.0

-- For Testing
nc = 64
boxsize = 1280.0

-- Broken, always use 1
nrealization= 1

-------- Time Sequence ----
-- Uniform time steps in a
ntimestep= 10
a_init = 0.1
a_final= 1.0
-- Alternatively, specify the steps manually
--time_step = {0.1, 0.5, 1.0}

output_redshifts= {1.0, 0.55, 0.09, 0.0}  -- redshifts of output

-- Cosmology --
omega_m = 1-0.708
h       = 0.69

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
powerspectrum= "powerspec.txt"
random_seed= 100
sigma8  = 0.82 -- 0 if the power spectrum is already normalized

-- Alternatively, give a RunPB Initial conditions file
--readic="/scratch1/scratchdirs/yfeng1/PB01/tpmsph_ic.bin" 

-------- Approximation Method ---------------
force_mode = "pm"
cola_nonstdda = false
smoothing = 2.0 -- Unused
diff_order = 1
enforce_broadband = true

pm_nc_factor1= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm = 0.5            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1
pm_nc_factor2= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension after a > change_pm
np_alloc_factor= 1.5      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

-------- Output ---------------

-- Dark matter particle outputs (all particles)
snapshot= "example-output/snp"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
measure_power = "example-output/powerspec"

