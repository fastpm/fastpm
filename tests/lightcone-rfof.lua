-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 64
boxsize = 512

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 8)

output_redshifts= {0.0}  -- redshifts of output

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
growth_mode = "LCDM"
pm_nc_factor = 1
lpt_nc_factor = 1
np_alloc_factor = 2.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
-- write_runpb_snapshot= "nbodykit/tpm"
write_snapshot = "lightcone/fastpm"
-- FIXME: currently fof and rfof cannot be both enabled during lightcone,
-- because we append the 'tail' to lightcone event particles, and if both
-- are enabled then one of them is wrong!
-- write_fof = "lightcone/fof"
write_rfof = "lightcone/rfof"
-- 1d power spectrum (raw), without shotnoise correction

particle_fraction = 1.0
dh_factor = 0.1

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
lc_usmesh_ell_limit = 200  -- downsamples the lightcone particles.
