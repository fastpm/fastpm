-- parameter file
------ Size of the simulation -------- 
-- ncdm.lua n_shell n_side n_p
-- For Testing
nc = 128
boxsize = 1000

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- logspace: Uniform time steps in loga
time_step = linspace(0.05,1.,41)

output_redshifts= {19.0, 9., 2., 1., 0.}  -- redshifts of output

-- Cosmology --
omega_m = 0.307494
omega_ncdm = 0.001404
h       = 0.6774

--ncdm split
m_ncdm = {0.05,0.05,0.05}
n_shell = 10
n_side = 2
lvk = true
every_ncdm = 1

-- Start with a linear density field
-- Power spectrum of the linear density field: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum = "Pcb0.txt"
linear_density_redshift = 0. -- the redshift of the linear density field.

read_powerspectrum_ncdm = "Pncdm99.txt" --comment out for bias run! Could defo add bias para or something
linear_density_redshift_ncdm = 99. -- the redshift of the linear density field.

random_seed= 42
particle_fraction = 1.0

--sort_snapshot = false
--
-------- Approximation Method ---------------
force_mode = "fastpm"

pm_nc_factor = 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

remove_cosmic_variance = true

dealiasing_type = "none"   --aggressive means GADGET now

-------- Output ---------------

filename = string.format("NEWsrunE_SH%d_NS%d_every%d_proc%d_nc%d_size%d_lvk%s_rcv%s_dt%s_pnf%d_z99", n_shell, n_side, every_ncdm, os.get_nprocs(), nc, boxsize, lvk, remove_cosmic_variance, dealiasing_type, pm_nc_factor)  --add time_step to fname?
--loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/Pncdm_init_test/"
loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/lowz/"

-- Dark matter particle outputs (all particles)
write_snapshot = loc .. filename .. "/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = loc .. filename .. "/powerspec"
--write_fof = loc .. filename .. "/fastpm"


