-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 384.0
if args[2] == 'za' then
    za = true
    force_mode = "pm"
    time_step = {1.0}
elseif args[2] == '2lpt' then
    za = false
    force_mode = "pm"
    time_step = {1.0}
elseif args[2] == 'cola' then
    za = false
    force_mode = "cola"
    time_step = linspace(0.1, 1, 10)
elseif args[2] == 'pm' then
    za = false
    force_mode = "pm"
    time_step = linspace(0.1, 1, 10)
elseif args[2] == 'zola' then
    za = false
    force_mode = "zola"
    time_step = linspace(0.1, 1, 10)
elseif args[2] == 'ic' then
    za = false
    force_mode = "zola"
    time_step = {0.1}
end

local function has(keyword)
    for i,k in pairs(args) do
        if k == keyword then
            return true
        end
    end
    return false
end

if has('inverted') then
    inverted_ic = true
else
    inverted_ic = false
end

if has('remove_variance') then
    remove_cosmic_variance = true
else
    remove_cosmic_variance = false
end

prefix = 'results'
for i,k in pairs(args) do
    if i > 1 then
    prefix = prefix .. '-' .. k
    end
end

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
output_redshifts= {9.0, 0.0}  -- redshifts of output

-- Cosmology --
omega_m = 0.307494
h       = 0.6774

if has('lineark') then
    read_lineark = 'results-ic/lineark'
else
    --Start with a power spectrum file
    -- Initial power spectrum: k P(k) in Mpc/h units
    -- Must be compatible with the Cosmology parameter
    read_powerspectrum = "powerspec.txt"
    random_seed= 100
    if has('whitenoisek') then
        read_whitenoisek = 'results-ic/whitenoisek'
    end
end

-------- Approximation Method ---------------
force_mode = "pm"

pm_nc_factor = {2, }            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
change_pm =    {0,}            -- time(scaling factor) when the pm_nc_factor is changed, range from 0 to 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= prefix .. "/fastpm"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = prefix .. "/powerspec"
write_whitenoisek = prefix .. "/whitenoisek"

write_nonlineark = prefix .. "/nonlineark"
write_lineark = prefix .. "/lineark"
