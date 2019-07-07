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
m_ncdm = {0.03}
n_shell = 10
n_side = 2
every_ncdm = 2
lvk = true

mass_string = ""
for i,m in pairs(m_ncdm) do
    mass_string = mass_string .. m .. "-"
end
mass_string = mass_string:sub(1, -2)    --not lua uses 1-indexing

-- Start with a linear density field
-- Power spectrum of the linear density field: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
power_loc = "/global/cscratch1/sd/abayer/nbodykit_cosmology/"

linear_density_redshift = 0. -- the redshift of the linear density field.
power_fname_cb = string.format("Pcb%d_mncdm%s.txt", linear_density_redshift, mass_string)
read_powerspectrum = power_loc .. power_fname_cb

linear_density_redshift_ncdm = 19. -- the redshift of the linear density field.
power_fname_ncdm = string.format("Pncdm%d_mncdm%s.txt", linear_density_redshift_ncdm, mass_string)
read_powerspectrum_ncdm = power_loc .. power_fname_ncdm

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

filename = string.format("coriE_SH%d_NS%d_every%d_proc%d_nc%d_size%d_lvk%s_rcv%s_fst%s_pnf%d_z%d_m%s", n_shell, n_side, every_ncdm, os.get_nprocs(), nc, boxsize, lvk, remove_cosmic_variance, force_softening_type, pm_nc_factor, linear_density_redshift_ncdm, mass_string)  --add time_step to fname?
--filename = string.format("LU2_SH10_NS2_Ndm1")
--loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/Pncdm_init_test/"
loc = "/global/cscratch1/sd/abayer/fastpm/ncdm/lowz/"
--loc = "/global/cscratch1/sd/abayer/fastpm/trash/"

-- Dark matter particle outputs (all particles)
write_snapshot = loc .. filename .. "/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = loc .. filename .. "/powerspec"
--write_fof = loc .. filename .. "/fastpm"


