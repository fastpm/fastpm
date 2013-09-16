-- cola code parameter file
nc = 64
boxsize = 100.0

random_seed= 100
nrealization= 1

ntimestep= 10
a_final= 1.0
output_redshifts= {0.73, 0.6, 0.44, 0.0}  -- redshifts of output

omega_m = 0.273
h       = 0.705
sigma8  = 0.812

pm_nc_factor= 3            -- Particle Mesh grid pm_nc_factor*nc per dimension
np_alloc_factor= 1.25      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

powerspectrum= "camb0_matterpower.dat" -- Initial power spectrum: k P(k)

-- Options
--   Following outputs can be turned off by commenting out
--   fof, snapshot, subsample, coarse_grid

-- FoF halo catalogue
fof= "fof"                 -- base filename
linking_factor= 0.2        -- FoF linking length= linking_factor*mean_separation

-- Dark matter particle outputs (all particles)
snapshot= "snp"       -- comment out to suppress snapshot output

-- Dark matter particle subsample
subsample= "sub"           -- base filename
subsample_factor= 1/100    -- fraction of particles to output

-- Dark matter density grid
coarse_grid= "grid"        -- base filename
coarse_grid_nc= 16         -- number of grid per dimension

-- Initial condition for other N-body simulation (at a_init= a_final/ntimestep)
-- initial= "init"


-- Use 8-byte long id for GADGET snapshot
write_longid= false

