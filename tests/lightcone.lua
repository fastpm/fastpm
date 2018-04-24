-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 256

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 10)

output_redshifts= {9.0, 0.0}  -- redshifts of output
compute_potential = true
compute_tidal = true

-- Cosmology --
omega_m = 0.307494
h       = 0.6774

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
random_seed= 100

-------- Approximation Method ---------------
force_mode = "fastpm"
-- force_mode = "cola"

pm_nc_factor = 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning

np_alloc_factor= 8.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
-- write_runpb_snapshot= "nbodykit/tpm"
write_snapshot = "lightcone/fastpm" 
write_fof = "lightcone/fof"
-- 1d power spectrum (raw), without shotnoise correction

dh_factor = 0.01
lc_fov = 0
lc_fov = 361

--s =[[glmatrix = { 
--        {1, 0, 0, -128},
--        {0, 1, 0, -128},
--        {0, 0, 1, -128},
--        {0, 0, 0, 1},
--        }
--]]

lc_glmatrix = fastpm.translation(-128, -128, -128)
lc_smesh_max_nside=512
lc_amin = 0.1
lc_amax = 1.0
--
-- lc_usmesh_tiles = fastpm.outerproduct({0}, {0}, {0, 1, 2, 3})
lc_write_usmesh = "lightcone/usmesh"
lc_write_smesh = "lightcone/smesh"
