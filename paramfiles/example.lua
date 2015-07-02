-- cola code parameter file
nc = 2048
boxsize = 1380.0

random_seed= 100
nrealization= 1

ntimestep= 10
a_init = 0.1
a_final= 1.0
--time_step = {0.1, 0.5, 1.0}
output_redshifts= {1.0, 0.55, 0.0}  -- redshifts of output

omega_m = 1-0.708
h       = 0.69
sigma8  = 0.82

force_mode = "cola"
cola_stdda = false
smoothing = 2.0
diff_order = 1
enforce_broadband = false

pm_nc_factor1= 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm = 0.5            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1
pm_nc_factor2= 3            -- Particle Mesh grid pm_nc_factor*nc per dimension after a > change_pm
np_alloc_factor= 1.5      -- Amount of memory allocated for particle
loglevel=0                 -- 0=verbose increase value to reduce output msgs

-- use a RunPB IC file. nc and omega shall match!
readic="/scratch1/scratchdirs/yfeng1/PB01/tpmsph_ic.bin" 
-- powerspectrum= "/home/energy/cola_optim/cola_halo/camb0_matterpower.dat" -- Initial power spectrum: k P(k)

-- Options
--   Following outputs can be turned off by commenting out
--   snapshot

-- Dark matter particle outputs (all particles)
snapshot= "/scratch3/scratchdirs/energy/qrpmc10/snp_PB01"       -- comment out to suppress snapshot output

-- Use 8-byte long id for GADGET snapshot
write_longid= true

measure_power = "/scratch3/scratchdirs/energy/qrpmc10/powerspec_PB01"
