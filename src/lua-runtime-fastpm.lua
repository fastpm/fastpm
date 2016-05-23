local function Schema()
    local self = {}

    local names = {}
    local trace_actions = {}

    function self.declare(options)
        local name = options.name
        local type = options.type
        local required = options.required or false
        local default = options.default
        local choices = options.choices
        local action = options.action

        self[name] = {type=type, required=required, default=default, choices=choice, action=action}
        names[#names + 1] = name
    end

    function self.print()
        for i,k in pairs(names) do
            v = self[k]
            local header = string.format('name: %s ', k)
            local body = ''
            for kk,vv in pairs(v) do
                body = string.format('%s %s:%s', body, kk,vv)
            end
            print(string.format('%s - %s', header, body))
        end
    end

    local function exists(filename)
        local file = io.open(filename, 'r')
        if file == nil then
            return true
            --return false
        end
        file:close()
        return true
    end

    local function check_array(name, entry, value)
        eletype = entry.type:match('array:(%a+)')
        for i, v in pairs(value) do
            if type(v) ~= eletype then
                error(string.format("entry %d of key `%s` is of type `%s`, but `%s` is expected",
                    i, name, type(v), eletype))
            end
        end
    end

    local function check_file(name, entry, value)
        if not exists(value) then
            error(string.format("file `%s' is not readable.", value))
        end
    end

    local function check_enum(name, entry, value)
        if entry.choices[value] == nil then
            local s = ''
            for v,_ in pairs(entry.choices) do
                s = s .. string.format('`%s`', v) .. ' '
            end
            print(s)
            error(string.format("value `%s` of key `%s` is not one of %s",
                value, name, s))
        end
    end

    local function check_type(name, entry, value)
        if entry.type == 'file' then
            -- File?
            check_file(name, entry, value)
        elseif entry.type:match('array:(%a+)') ~= nil then
            -- Array?
            check_array(name, entry, value)
        elseif entry.type == 'enum' then
            -- Enum?
            check_enum(name, entry, value)
        elseif type(value) ~= entry.type then
            -- Anything else
            error(string.format("key `%s` is of type `%s`, but `%s` is expected",
                name, type(value), entry.type))
        end
    end

    local function bind_one(name, entry, value)
        if value == nil then
            if entry.required then
                error(string.format("`%s` is required but undefined.", name))
            else
                return entry.default
            end
        else

            check_type(name, entry, value)

            if entry.action ~= nil then
                entry.action(value)
            end
            return value
        end
    end

    function self.bind(namespace)
        result = {}
        for i,k in pairs(names) do
            result[k] = bind_one(k, self[k], namespace[k])
        end
        return result
    end

    function self.set_action(name, action)
        self[name].action = action
    end

    return self
end

schema = Schema()
schema.declare{name='nc',                type='number', required=true}
schema.declare{name='boxsize',           type='number', required=true}
schema.declare{name='time_step',         type='array:number',  required=true }
schema.declare{name='output_redshifts',  type='array:number',  required=true }
schema.declare{name='omega_m',           type='number', required= true }
schema.declare{name='h',                 type='number', required=true }
schema.declare{name='pm_nc_factor',      type='array:number',  required=true }
schema.declare{name='change_pm',         type='array:number',  required=true }
schema.declare{name='np_alloc_factor',   type='number', required=true }

-- Force calculation --
schema.declare{name='force_mode',        type='enum', required=true}
schema.force_mode.choices = {
    cola = 'FASTPM_FORCE_MODE_COLA',
    zola = 'FASTPM_FORCE_MODE_ZOLA',
    pm   = 'FASTPM_FORCE_MODE_PM',
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

schema.declare{name='enforce_broadband_kmax',  type='number', default=4}
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
    end
end

-- Primordial Non-Gaussianity --
schema.declare{name='f_nl_type', type='enum', default='none'}
schema.f_nl_type.choices = {
    ['local'] = 'FASTPM_FNL_TYPE_LOCAL',
    ['none']  = 'FASTPM_FNL_TYPE_NONE',
}
schema.declare{name='f_nl', type='number'}
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
schema.declare{name='read_whitenoise',         type='string'}

schema.declare{name='read_powerspectrum', type='file'}
schema.declare{name='sigma8',             type='number', default=0}
schema.declare{name='random_seed',         type='number'}
schema.declare{name='inverted_ic',             type='boolean', default=false}
schema.declare{name='remove_cosmic_variance',  type='boolean', default=false}

function schema.read_grafic.action (read_grafic)
    schema.random_seed.required = false
    schema.read_powerspectrum.required = false
end

function schema.read_lineark.action (read_lineark)
    schema.random_seed.required = false
    schema.read_powerspectrum.required = false
end
function schema.read_runpbic.action (read_lineark)
    schema.random_seed.required = false
    schema.read_powerspectrum.required = false
end
function schema.read_whitenoise.action (read_whitenoise)
    schema.random_seed.required = false
    schema.read_powerspectrum.required = true
end

schema.declare{name='write_lineark',         type='string'}
schema.declare{name='write_whitenoise',         type='string'}
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
    ['3_2'] = 'FASTPM_KERNEL_3_2',
}
schema.declare{name='dealiasing_type',             type='enum', default="none"}
schema.dealiasing_type.choices = {
    none = 'FASTPM_DEALIASING_NONE',
    gaussian = 'FASTPM_DEALIASING_GAUSSIAN',
    twothird = 'FASTPM_DEALIASING_TWO_THIRD',
}

local fastpm = {}

fastpm.schema = schema

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
    r = fastpm.linspace(a, e, N)
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

return "fastpm", fastpm
