-- cola code parameter file
nc = 256
boxsize = 400.0

random_seed= 100
nrealization= 1

ntimestep= 10
a_final= 1.0
output_redshifts= {0.73, 0.6, 0.44, 0.0}  -- redshifts of output

omega_m = 1-0.7168
h       = 0.697
sigma8  = 0.817

qpm = true
time_step = "loga"
stdda = true
smoothing = 2.0
diff_order = 1
nopm = false


pm_nc_factor1= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm = 0.5            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1
pm_nc_factor2= 3            -- Particle Mesh grid pm_nc_factor*nc per dimension after a > change_pm
np_alloc_factor= 1.25      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

-- use a RunPB IC file. nc and omega shall match!
-- readic="PMIC/tpmsph_ic.bin" 
powerspectrum= "/home/energy/cola_optim/cola_halo/camb0_matterpower.dat" -- Initial power spectrum: k P(k)

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

measure_power = "powerspec"
