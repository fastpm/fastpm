-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 128
boxsize = 384.0
if args[1] == 'za' then
    za = true
    force_mode = "pm"
    time_step = {1.0}
elseif args[1] == '2lpt' then
    za = false
    force_mode = "pm"
    time_step = {1.0}
elseif args[1] == 'cola' then
    za = false
    force_mode = "cola"
    time_step = linspace(0.1, 1, 5)
elseif args[1] == 'pm' then
    za = false
    force_mode = "pm"
    time_step = linspace(0.1, 1, 5)
elseif args[1] == 'zola' then
    za = false
    force_mode = "zola"
    time_step = linspace(0.1, 1, 5)
elseif args[1] == 'fastpm' then
    za = false
    force_mode = "fastpm"
    time_step = linspace(0.1, 1, 5)
elseif args[1] == 'ic' then
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
    painter_support = 6
end
if has('lanczos2') then
    painter_type = "lanczos"
    painter_support = 4
end
if has('linear1') then
    painter_type = "linear"
    painter_support = 2
end
if has('linear2') then
    painter_type = "linear"
    painter_support = 4
end
if has('inverted') then
    inverted_ic = true
else
    inverted_ic = false
end

if has('fixed_mode') then
    set_mode_method = "add"
    set_mode = {
            {0, 1, 0, 0, 0.1},
            {0, 1, 0, 1, 0.0},
            }
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
    if i > 0 then
    prefix = prefix .. '-' .. k
    end
end

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
output_redshifts= {9.0, 2.0, 1.0, 0.0}  -- redshifts of output

-- Cosmology --
Omega_m = 0.307494
h       = 0.6774

if has('lineark') then
    read_lineark = 'results-ic/IC'
else
    --Start with a power spectrum file
    -- Initial power spectrum: k P(k) in Mpc/h units
    -- Must be compatible with the Cosmology parameter
    read_powerspectrum = "powerspec.txt"
    random_seed= 100
    if has('whitenoisek') then
        read_whitenoisek = 'results-ic/IC'
    end
end

if has('fnl') then
    -- FIXME: what?
    f_nl_type = "local"
    scalar_amp = 2.130624e-9
    scalar_pivot = 0.05
    scalar_spectral_index = 0.9667
    f_nl = 10.0
    kmax_primordial_over_knyquist = 0.25
else
    f_nl_type = "none"
end
-------- Approximation Method ---------------

pm_nc_factor = 3            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning
lpt_nc_factor = 1

np_alloc_factor= 4.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
write_snapshot= prefix .. "/fastpm"       -- comment out to suppress snapshot output
-- 1d power spectrum (raw), without shotnoise correction
write_powerspectrum = prefix .. "/powerspec"
write_whitenoisek = prefix .. "/IC"

write_nonlineark = prefix .. "/fastpm"
write_lineark = prefix .. "/IC"

if has('lightcone') then
    write_lightcone = prefix .. "/lightcone"
    dh_factor = 0.05
end

if has('constrain') then
    constraints = {
        {boxsize * 0.5, boxsize * 0.5, boxsize * 0.5, 100.},
        {boxsize * 0.25, boxsize * 0.25, boxsize * 0.5, 100.},
        {boxsize * 0.75, boxsize * 0.75, boxsize * 0.5, 100.},
    }
end
