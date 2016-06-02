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
    time_step = linspace(0.1, 1, 5)
elseif args[2] == 'pm' then
    za = false
    force_mode = "pm"
    time_step = linspace(0.1, 1, 5)
elseif args[2] == 'zola' then
    za = false
    force_mode = "zola"
    time_step = linspace(0.1, 1, 5)
elseif args[2] == 'fastpm' then
    za = false
    force_mode = "fastpm"
    time_step = linspace(0.1, 1, 5)
elseif args[2] == 'ic' then
    za = false
    force_mode = "zola"
    time_step = {0.1}
else
    error("wrong arg!")
end

local function has(keyword)
    for i,k in pairs(args) do
        if k == keyword then
            return true
        end
    end
    return false
end
if has('lanczos3') then
    painter_type = "lanczos"
    painter_support = 3
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

if has('shift') then
    shift = true
else
    shift = false
end
if has('gaussian36') then
    dealiasing_type = 'gaussian36'
end
if has('aggressive') then
    dealiasing_type = 'aggressive'
end
if has('gaussian') then
    dealiasing_type = 'gaussian'
end
if has('twothird') then
    dealiasing_type = 'twothird'
end
if has('eastwood') then
    kernel_type = 'eastwood'
end
if has('gadget') then
    kernel_type = 'gadget'
end
if has('naive') then
    kernel_type = 'naive'
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

if has('fnl') then
    -- FIXME: what?
    f_nl_type = "local"
    scalar_amp = 2.441e-9
    scalar_pivot = 0.002
    scalar_spectral_index = 0.9667
    f_nl = 0.1
else
    f_nl_type = "none"
end
-------- Approximation Method ---------------

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
