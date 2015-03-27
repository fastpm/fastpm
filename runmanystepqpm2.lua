-- cola code parameter file
nc = 256
boxsize = 256.0

random_seed= 100
nrealization= 1

ntimestep= 40
a_final= 1.0
output_redshifts= {0.55, 0.0}  -- redshifts of output

omega_m = 1-0.7168
h       = 0.697
sigma8  = 0.817

-- somewhat canonical qpm parameters
qpm = true
time_step = 'a'
stdda = true
smoothing = 2.0
diff_order = 2
nopm = false


pm_nc_factor1= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm = 0.7            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1
pm_nc_factor2= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension after a > change_pm
np_alloc_factor= 1.25      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

powerspectrum= "powerspec.txt" -- Initial power spectrum: k P(k)

-- Options
--   Following outputs can be turned off by commenting out
--   fof, snapshot, subsample, coarse_grid

-- FoF halo catalogue
fof= "runmanystepqpm2/fof"                 -- base filename
linking_factor= 0.2        -- FoF linking length= linking_factor*mean_separation

-- Dark matter particle outputs (all particles)
snapshot= "runmanystepqpm2/snp"       -- comment out to suppress snapshot output

-- Dark matter particle subsample
-- subsample= "runmanystepqpm2/sub"           -- base filename
-- subsample_factor= 1/100    -- fraction of particles to output

-- Dark matter density grid
-- coarse_grid= "runmanystepqpm2/grid"        -- base filename
-- coarse_grid_nc= 16         -- number of grid per dimension

-- Initial condition for other N-body simulation (at a_init= a_final/ntimestep)
-- initial= "init"


-- Use 8-byte long id for GADGET snapshot
write_longid= false

measure_power = "runmanystepqpm2/powerspec"
