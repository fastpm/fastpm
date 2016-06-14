----------------------------------------------------
-- This is the main LUA runtime library of FastPM.
--
-- Author: Yu Feng <rainwoodman@gmail.com> 2016
----------------------------------------------------

local config = require("lua-runtime-config")

local _NAME = ... or 'main'

local fastpm = {}

local schema = config.Schema()
schema.declare{name='nc',                type='int', required=true, help="Number of Particles Per side"}
schema.declare{name='boxsize',           type='number', required=true, help="Size of box in Mpc/h"}
schema.declare{name='time_step',         type='array:number',  required=true, help="Scaling factor of steps, can be linspace(start, end, Nsteps)." }
schema.declare{name='output_redshifts',  type='array:number',  required=true, help="Redshifts for outputs" }
schema.declare{name='aout',  type='array:number',  required=false}
-- set aout from output_redshifts
function schema.output_redshifts.action(output_redshifts)
    local aout = {}
    for i, z in pairs(output_redshifts) do
        aout[i] = 1.0 / (z + 1.)
    end
    schema.aout.default = aout
end

schema.declare{name='omega_m',           type='number', required=true, default=0.3 }
schema.declare{name='h',                 type='number', required=true, default=0.7, help="Dimensionless Hubble parameter"}
schema.declare{name='pm_nc_factor',      type='array:number',  required=true, help="A list of PM resolutions, must be the same length of change_pm"}
schema.declare{name='change_pm',         type='array:number',  required=true, help="A list of scaling factor that the PM resolution changes" }
schema.declare{name='np_alloc_factor',   type='number', required=true, help="Over allocation factor for load imbalance" }

-- Force calculation --
schema.declare{name='painter_type',        type='enum', default='cic', help="Type of painter."}
schema.declare{name='painter_support',     type='int', default=1, help="Support (size) of the painting kernel"}
schema.painter_type.choices = {
    cic = 'FASTPM_PAINTER_CIC',
    linear = 'FASTPM_PAINTER_LINEAR',
    lanczos = 'FASTPM_PAINTER_LANCZOS',
}
function schema.painter_type.action(painter_type)
    if painter_type ~= 'cic' then
        schema.painter_support.required = true
    end
end
schema.declare{name='force_mode',        type='enum', default='fastpm'}
schema.force_mode.choices = {
    cola = 'FASTPM_FORCE_COLA',
    zola = 'FASTPM_FORCE_FASTPM',
    fastpm = 'FASTPM_FORCE_FASTPM',
    pm   = 'FASTPM_FORCE_PM',
}
schema.declare{name='enforce_broadband_mode',  type='enum',
            default='none', choices={'pm', 'linear', '2lpt', 'za', 'none'}}
schema.enforce_broadband_mode.choices = {
    pm = 'FASTPM_MODEL_PM',
    linear = 'FASTPM_MODEL_LINEAR',
    ['2lpt'] = 'FASTPM_MODEL_2LPT',
    za = 'FASTPM_MODEL_ZA',
    none = 'FASTPM_MODEL_NONE',
}

schema.declare{name='enforce_broadband_kmax',  type='int', default=4}
schema.declare{name='cola_stdda',           type='boolean'}

function schema.force_mode.action (force_mode)
    if force_mode == "pm" then
        schema.cola_stdda.default = true
        schema.enforce_broadband_mode.default = "pm"
    elseif force_mode == "cola" then
        schema.cola_stdda.default = false
        schema.enforce_broadband_mode.default = "none"
    elseif force_mode == "zola" then
        schema.cola_stdda.default = true
        schema.enforce_broadband_mode.default = "none"
    elseif force_mode == "fastpm" then
        schema.cola_stdda.default = true
        schema.enforce_broadband_mode.default = "none"
    end
end

-- Primordial Non-Gaussianity --
schema.declare{name='f_nl_type', type='enum', default='none'}
schema.f_nl_type.choices = {
    ['local'] = 'FASTPM_FNL_LOCAL',
    ['none']  = 'FASTPM_FNL_NONE',
}
schema.declare{name='f_nl', type='number'}
schema.declare{name='kmax_primordial_over_knyquist', type='number', default=0.25}
schema.declare{name='scalar_amp', type='number'}
schema.declare{name='scalar_pivot', type='number'}
schema.declare{name='scalar_spectral_index', type='number'}
function schema.f_nl_type.action (f_nl_type)
    if f_nl_type ~= 'none' then
        schema.f_nl.required = true
        schema.scalar_amp.required = true
        schema.scalar_pivot.required = true
        schema.scalar_spectral_index.required = true
    end
end

-- Initial condition --
schema.declare{name='read_grafic',        type='string'}
schema.declare{name='read_lineark',        type='string'}
schema.declare{name='read_runpbic',       type='string'}
schema.declare{name='read_whitenoisek',         type='string'}

schema.declare{name='fix_ic_mode', type='array:int'}
schema.declare{name='fix_ic_value', type='number'}

function schema.fix_ic_mode.action (fix_ic_mode)
    if fix_ic_mode ~= nil then
        schema.fix_ic_value.required = true
        if #fix_ic_mode ~= 4 then
            error("fix_ic_mode must be four integers (x, y, z, real_or_imag)")
        end
        if fix_ic_mode[4] ~= 0 and
           fix_ic_mode[4] ~= 1 then
            error("the fourth component of the mode must be 0 (real) or 1 (imaginary)")
        end
    end
end

schema.declare{name='read_powerspectrum', type='file'}
schema.declare{name='sigma8',             type='number', default=0}
schema.declare{name='random_seed',         type='int'}
schema.declare{name='shift',             type='boolean', default=false}
schema.declare{name='inverted_ic',             type='boolean', default=false}
schema.declare{name='remove_cosmic_variance',  type='boolean', default=false}

function schema.read_grafic.action (read_grafic)
    if read_grafic ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = false
    end
end

function schema.read_lineark.action (read_lineark)
    if read_lineark ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = false
    end
end
function schema.read_runpbic.action (read_runpbic)
    if read_runpbic ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = false
    end
end
function schema.read_whitenoisek.action (read_whitenoisek)
    if read_whitenoisek ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = true
    end
end

schema.declare{name='write_lineark',         type='string'}
schema.declare{name='write_whitenoisek',         type='string'}
schema.declare{name='write_runpbic',       type='string'}
schema.declare{name='write_powerspectrum', type='string'}
schema.declare{name='write_snapshot',      type='string'}
schema.declare{name='write_nonlineark',      type='string'}
schema.declare{name='write_runpb_snapshot', type='string'}

schema.declare{name='za',                      type='boolean', default=false}
schema.declare{name='kernel_type',             type='enum', default="3_4"}
schema.kernel_type.choices = {
    ['3_4'] = 'FASTPM_KERNEL_3_4',
    ['5_4'] = 'FASTPM_KERNEL_5_4',
    ['eastwood'] = 'FASTPM_KERNEL_EASTWOOD',
    ['gadget'] = 'FASTPM_KERNEL_GADGET',
    ['naive'] = 'FASTPM_KERNEL_NAIVE',
    ['3_2'] = 'FASTPM_KERNEL_3_2',
}
schema.declare{name='dealiasing_type',             type='enum', default="none"}
schema.dealiasing_type.choices = {
    none = 'FASTPM_DEALIASING_NONE',
    gaussian = 'FASTPM_DEALIASING_GAUSSIAN',
    aggressive = 'FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN',
    gaussian36 = 'FASTPM_DEALIASING_GAUSSIAN36',
    twothird = 'FASTPM_DEALIASING_TWO_THIRD',
}

function fastpm.linspace(a, e, N)
-- Similar to numpy.linspace, but always append the end
-- point to the result, returning N + 1 elements.
--
-- https://mail.scipy.org/pipermail/numpy-discussion/2016-February/075065.html
    local r = {}
    N1 = N + 1
    for i=1,N1 do
        r[i] = 1.0 * (e - a) * (i - 1) / N + a
    end
    r[N1] = e
    return r
end

function fastpm.logspace(a, e, N)
-- a and end are in log10.
-- Returns N+1 elements, including e.
    local r
    r = linspace(a, e, N)
    for i, j in pairs(r) do
        r[i] = math.pow(10, j)
    end
    return r
end
function fastpm.blendspace(a, e, a1, a2)
    local r = {}
    local a = a
    local i = 1
    while a < e do
        r[i] = a
        dlna = math.pow(math.pow(1/a1, 2) + math.pow(a/a2, 2), -0.5)
        a = math.exp(math.log(a) + dlna)
        i = i + 1
    end
    r[i] = e
    return r
end

function fastpm.test()
    ns = {
        __file__ = "standard.lua",
        boxsize = 384.0,
        cola_stdda = true,
        dealiasing_type = "none",
        enforce_broadband_kmax = 4,
        enforce_broadband_mode = "pm",
        f_nl = 0.100000000000000006,
        f_nl_type = "local",
        force_mode = "pm",
        h = 0.677400000000000002,
        inverted_ic = false,
        kernel_type = "3_4",
        nc = 128,
        np_alloc_factor = 4.0,
        omega_m = 0.30749399999999999,
        prefix = "results-za-nongaussian",
        random_seed = 100,
        read_powerspectrum = "powerspec.txt",
        remove_cosmic_variance = false,
        scalar_amp = 0.000000002441,
        scalar_pivot = 0.002,
        scalar_spectral_index = 0.966700000000000004,
        sigma8 = 0,
        write_lineark = "results-za-nongaussian/lineark",
        write_nonlineark = "results-za-nongaussian/nonlineark",
        write_powerspectrum = "results-za-nongaussian/powerspec",
        write_snapshot = "results-za-nongaussian/fastpm",
        write_whitenoisek = "results-za-nongaussian/whitenoisek",
        za = true,
        args = {
            "standard.lua",
            "za",
            "nongaussian",
        },
        change_pm = {
            0,
        },
        output_redshifts = {
            9.0,
            0.0,
        },
        pm_nc_factor = {
            2,
        },
        time_step = {
            1.0,
        },
    }
    schema.print()
    for i,k in pairs(schema.bind(ns)) do
        print(i, k)
    end
end

fastpm.schema = schema

--
-- The main functions, must be a global symbol
--
function _parse_runmain(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    logspace = fastpm.logspace
    linspace = fastpm.linspace

    return config.parse(fastpm.schema, filename, true, {...})
end

function _parse(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    logspace = fastpm.logspace
    linspace = fastpm.linspace

    return config.parse(fastpm.schema, filename, false, {...})
end

function _help(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    return fastpm.schema.format_help()
end

return fastpm
