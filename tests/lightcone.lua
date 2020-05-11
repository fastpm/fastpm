-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 512

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 8)

output_redshifts= {0.0}  -- redshifts of output
compute_potential = true
compute_tidal = true

-- Cosmology --
Omega_m = 0.307494
h       = 0.6774

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
random_seed= 100
remove_cosmic_variance=true
-------- Approximation Method ---------------
force_mode = "fastpm"
-- force_mode = "cola"
growth_mode = "LCDM"
pm_nc_factor = 2
lpt_nc_factor = 1
np_alloc_factor = 2.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
-- write_runpb_snapshot= "nbodykit/tpm"
write_snapshot = "lightcone/fastpm"
write_fof = "lightcone/fof"
-- 1d power spectrum (raw), without shotnoise correction

particle_fraction = 1.0
fof_linkinglength = 0.2
fof_nmin = 4
dh_factor = 0.1

lc_fov = 0

--lc_octants = {0}
lc_fov = 360

--s =[[glmatrix = { 
--        {1, 0, 0, -128},
--        {0, 1, 0, -128},
--        {0, 0, 1, -128},
--        {0, 0, 0, 1},
--        }
--]]

-- lc_glmatrix = fastpm.translation(-128, -128, -128)
lc_amin = 0.1
lc_amax = 1.0

lc_write_usmesh = "lightcone/usmesh"
lc_usmesh_tiles = fastpm.outerproduct({-2, -1, 0, 1}, {-2, -1, 0, 1}, {-2, -1, 0, 1})
lc_usmesh_fof_padding = 20.0
lc_usmesh_alloc_factor = 2.0
