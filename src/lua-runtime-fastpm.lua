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
schema.declare{name='pm_nc_factor',      type='array:number',  required=true, help="A list of {a, PM resolution}, "}
schema.declare{name='np_alloc_factor',   type='number', required=true, help="Over allocation factor for load imbalance" }
schema.declare{name='compute_potential',   type='boolean', required=false, default=false, help="Calculate the gravitional potential."}

-- Force calculation --
schema.declare{name='painter_type',        type='enum', default='cic', help="Type of painter."}
schema.declare{name='painter_support',     type='int', default=2, help="Support (size) of the painting kernel"}
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
schema.declare{name='enforce_broadband_kmax',  type='int', default=4}

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

schema.declare{name='read_powerspectrum', type='file'}
schema.declare{name='sigma8',             type='number', default=0}
schema.declare{name='random_seed',         type='int'}
schema.declare{name='shift',             type='boolean', default=false}
schema.declare{name='inverted_ic',             type='boolean', default=false}
schema.declare{name='remove_cosmic_variance',  type='boolean', default=false}

function schema.read_grafic.action (read_grafic)
    if read_grafic ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = true
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
schema.declare{name='write_lightcone',         type='string'}
schema.declare{name='write_lightcone_potential',         type='string'}

schema.declare{name='dh_factor',    type='number', default=1.0, help='Scale Hubble distance to amplify the lightcone effect'}
schema.declare{name='fov',     type='number', default=0.0, help=' field of view of the sky. 0 for flat sky. the beam is along the z-direction after glmatrix.'}
schema.declare{name='glmatrix',     type='array:number',
        default={
            {1, 0, 0, 0,},
            {0, 1, 0, 0,},
            {0, 0, 1, 0,},
            {0, 0, 0, 0,},
        },
        help=[[transformation matrix to move (x, y, z, 1) to the observer coordinate with a left dot product.
               The observer is sitting at z=0. The last column of the matrix is the translation in Mpc/h.]]}

schema.declare{name='tiles',     type='array:number',
        default={
            {0, 0, 0},
        },
        help=[[tiling of the simulation box, in units of box edges.
              all tiles will be considered during lightcone construction.
              tiling occurs before the glmatrix.]]
        }

schema.declare{name='za',                      type='boolean', default=false, help='use ZA initial condition not 2LPT'}

schema.declare{name='kernel_type',             type='enum', default="3_4", help='Force kernel; very little effect.'}
schema.kernel_type.choices = {
    ['3_4'] = 'FASTPM_KERNEL_3_4',
    ['5_4'] = 'FASTPM_KERNEL_5_4',
    ['eastwood'] = 'FASTPM_KERNEL_EASTWOOD',
    ['gadget'] = 'FASTPM_KERNEL_GADGET',
    ['naive'] = 'FASTPM_KERNEL_NAIVE',
    ['3_2'] = 'FASTPM_KERNEL_3_2',
}
schema.declare{name='dealiasing_type',             type='enum', default="none", help='Dealiasing kernel (wipes out small scale force), very litle effect)'}
schema.dealiasing_type.choices = {
    none = 'FASTPM_DEALIASING_NONE',
    gaussian = 'FASTPM_DEALIASING_GAUSSIAN',
    aggressive = 'FASTPM_DEALIASING_AGGRESSIVE_GAUSSIAN',
    gaussian36 = 'FASTPM_DEALIASING_GAUSSIAN36',
    twothird = 'FASTPM_DEALIASING_TWO_THIRD',
}

schema.declare{name='constraints',      type='array:number',  help="A list of {x, y, z, overdensity}, giving the constraints in MPC/h units. "}
function schema.constraints.action (constraints)
    if constraints == nil then
        return
    end
    for i,v in pairs(constraints) do
        if #v ~= 4 then
            error("contraints must be a list of 4-vectors (x, y, z, real_or_imag, value)")
        end
    end
end
schema.declare{name='set_mode_method',  type='string', default="override", help="override or relative"}
schema.declare{name='set_mode',         type='array:number', help="A list of {kix, kiy, kiz, ri, value}, set the IC mode at integer k (ri for real and imag) to value"}

function schema.set_mode.action (set_mode)
    if set_mode == nil then
        return
    end
    for i,v in pairs(set_mode) do
        if #v ~= 5 then
            error("set_mode must be a list of 5-vectors (x, y, z, real_or_imag, value)")
        end
        if v[4] ~= 1 and v[4] ~=0 then
            error("the fourth component specifies real or imag part of the mode. must be 0 or 1")
        end
    end
end

function fastpm.translation(dx, dy, dz)
-- generate a translation gl matrix that shifts the coordinates
    return {
        {1, 0, 0, dx},
        {0, 1, 0, dy},
        {0, 0, 1, dz},
        {0, 0, 0, 0 },
        }
end

function fastpm.outerproduct(a, b, c)
-- generate a list that is outer product of elements of a, b, c
    local r = {}
    for i=1, #a do
        for j=1, #b do
            for k=1, #c do
                r[#r + 1] = {a[i], b[j], c[k]}
            end
        end
    end
    return r
end

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

    local globals = setmetatable({}, {__index=_G})

    globals.fastpm = fastpm
    globals.logspace = fastpm.logspace
    globals.linspace = fastpm.linspace

    return config.parse(fastpm.schema, filename, true, globals, {...})
end

function _parse(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    local globals = setmetatable({}, {__index=_G})
    globals.fastpm = fastpm
    globals.logspace = fastpm.logspace
    globals.linspace = fastpm.linspace

    return config.parse(fastpm.schema, filename, false, globals, {...})
end

function _help(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    return fastpm.schema.format_help()
end

return fastpm
