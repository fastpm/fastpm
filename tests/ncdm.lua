-- parameter file
------ Size of the simulation -------- 
-- ncdm.lua n_shell n_side n_p
-- For Testing
nc = 128
boxsize = 1000

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- logspace: Uniform time steps in loga
time_step = linspace(0.05,1.,11)

output_redshifts = {19.,0.}  --{19.0, 14., 9., 2., 1., 0.}  -- redshifts of output

-- Cosmology --
omega_m = 0.307494
omega_ncdm = 0.001404  --THIS IS HARDCODED ATM. SO THIS PARA DOESNT DO ANYTHING!
h       = 0.6774

--ncdm split
m_ncdm = {0.05}
n_shell = 10
n_side = 2
lvk = false
every_ncdm = 2

-- Start with a linear density field
-- Power spectrum of the linear density field: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
power_loc = "/global/cscratch1/sd/abayer/nbodykit_cosmology/"

read_powerspectrum = power_loc .. "Pcb0_Nncdm1.txt"
linear_density_redshift = 0. -- the redshift of the linear density field.

read_powerspectrum_ncdm = power_loc .. "Pncdm19_Nncdm1.txt"
linear_density_redshift_ncdm = 19. -- the redshift of the linear density field.

random_seed= 42
particle_fraction = 1.0

--sort_snapshot = false
--
-------- Approximation Method ---------------
force_mode = "fastpm"

pm_nc_factor = 2

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

remove_cosmic_variance = true

force_softening_type = "none"

-------- Output ---------------

filename = string.format("coriE_SH%d_NS%d_every%d_proc%d_nc%d_size%d_lvk%s_rcv%s_fst%s_pnf%d_z19_Nncdm1", n_shell, n_side, every_ncdm, os.get_nprocs(), nc, boxsize, lvk, remove_cosmic_variance, force_softening_type, pm_nc_factor)  --add time_step to fname?
--filename = string.format("LU2_SH10_NS2_Ndm1")
--loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/Pncdm_init_test/"
loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/lowz/"
--loc = "/global/cscratch1/sd/abayer/fastpm/trash/"

-- Dark matter particle outputs (all particles)
write_snapshot = loc .. filename .. "/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = loc .. filename .. "/powerspec"
--write_fof = loc .. filename .. "/fastpm"


