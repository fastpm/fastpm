local fastpm = {}

local Schema = {}
Schema.__index = Schema

function io.exists(filename)
    local file = io.open(filename, 'r')
    if file == nil then
        return false
    end
    file:close()
    return true
end

function Schema.new()
    local self = setmetatable({}, Schema)
    return self
end

function Schema.print(self)
    for k,v in pairs(self) do
        local header = string.format('name: %s ', k)
        local body = ''
        for kk,vv in pairs(v) do
            body = string.format('%s %s:%s', body, kk,vv)
        end
        print(string.format('%s - %s', header, body))
    end
end

function Schema.add(self, options)

    local name = options.name
    local type = options.type
    local required = options.required or false
    local default = options.default
    local choices = options.choices

    local choiceset
    if choices ~= nil then
        choiceset = {}
        for _,v in pairs(choices) do
            choiceset[v] = true
        end
    else
        choiceset = nil
    end
    self[name] = {type=type, required=required, default=default, choices=choiceset}
end

function Schema.validate(self, namespace, required_only)
    for name, entry in pairs(self) do
        -- Request all fields or just the required fields.
        -- After bind() is called, the default values would have
        -- been all set.
        value = namespace[name]
        if value == nil then
            if entry.required then
                error(string.format("key `%s` is requested but undefined", name))
            end
            if not entry.required and not required_only and not entry.default == nil then
                error(string.format("optional key `%s` is requested but undefined", name))
            end
        else
            local truetype = entry.type
            if entry.type == 'file' then
                truetype = 'string'
            end
            if type(value) ~= truetype then
                error(string.format("key `%s` is of type `%s`, but `%s` is expected",
                    name, type(value), entry.type))
            end
            if entry.choices and entry.choices[value] == nil then
                local s = ''
                for v,_ in pairs(entry.choices) do
                    s = s .. string.format('`%s`', v) .. ' '
                end
                error(string.format("value `%s` of key `%s` is not one of %s",
                    value, name, s))
            end
            if entry.type == 'file' then
                if not io.exists(value) then
                    error(string.format("file `%s' is not readable.", value))
                end
            end
        end
    end
end

function Schema.bind(self, namespace)
    local ret = {}
    for k, v in pairs(namespace) do
        ret[k] = v
    end
    local functionals = {}
    -- First resolve all plain value defaults
    for name, entry in pairs(self) do
        if (not entry.required) and ret[name] == nil then
            if type(entry.default) == 'function' then
                table.insert(functionals, name)
            else
                ret[name] = entry.default
            end
        end
    end

    -- Now iterate and try to resolve the functional default values
    -- If there is no circular dependency this will converge

    local errors = {}
    local N
    repeat
        N = #functionals
        for i, name in pairs(functionals) do
            f = self[name].default
            type = self[name].type
            status, v = pcall(f, ret)
            -- must be a good function call and a right return type
            -- (e.g. not nil)
            if status and type(v) == type then
                ret[name] = v
                table.remove(functionals, i)
                errors[name] = nil
            else
                if v == nil then
                    v = 'Type error'
                end
                errors[name] = v
            end
        end
    -- if this iteration failed to deduce any functional default values
    -- we are done
    until N == #functionals

    -- remaining items is an error
    if N ~= 0 then
        local names = ''
        for k,v in pairs(errors) do
            names = string.format('%s %s', names, k)
        end
        error("Cannot resolve functional default values : " .. names)
    end
    return ret
end

function Schema.dependency(self, ns)
    if #ns.change_pm ~= #ns.pm_nc_factor then
        error("change_pm and pm_nc_factor must be of the same length")
    end
    if not ns.read_lineark and not ns.read_linear and not ns.read_runpbic then
        -- need a power spectrum and gaussian random field
        self.read_powerspectrum.required = true
        if not ns.read_grafic then
            self.random_seed.required = true
        else
            self.random_seed.required = false
        end
    end

    if ns.force_mode == "pm" then
        self.cola_stdda.default = true
        self.enforce_broadband_mode.default = "pm"
    elseif ns.force_mode == "cola" then
        self.cola_stdda.default = false
        self.enforce_broadband_mode.default = "none"
    elseif ns.force_mode == "zola" then
        self.cola_stdda.default = true
        self.enforce_broadband_mode.default = "none"
    end

    --if self.time_step[0] <= 0.0 then
     --   error("Initial time must be greater than a = 0.0")
    --end
end

schema = Schema.new()
schema:add{name='nc',                type='number', required=true}
schema:add{name='boxsize',           type='number', required=true}
schema:add{name='time_step',         type='table',  required=true }
schema:add{name='output_redshifts',  type='table',  required=true }
schema:add{name='omega_m',           type='number', required= true }
schema:add{name='h',                 type='number', required=true }
schema:add{name='pm_nc_factor',      type='table',  required=true }
schema:add{name='change_pm',         type='table',  required=true }
schema:add{name='np_alloc_factor',   type='number', required=true } schema:add{name='force_mode',        type='string', required=true,
          choices={'cola', 'zola', 'pm'}}

schema:add{name='read_grafic',        type='string'}
schema:add{name='read_lineark',        type='string'}
schema:add{name='write_lineark',         type='string'}
schema:add{name='read_runpbic',       type='string'}
schema:add{name='read_powerspectrum', type='file'}
schema:add{name='scalar_amp, type='number', required=true}
schema:add{name='pivot_scalar, type='number', required=true}
schema:add{name='scalar_spectral_index', type='number', required=true}
schema:add{name='f_nl, type='number', default=0.0}
schema:add{name='f_nl_type, type='string, default='local',
           choices={'local'}}
schema:add{name='sigma8',             type='number', default=0}
schema:add{name='random_seed',         type='number'}
schema:add{name='read_whitenoise',         type='string'}
schema:add{name='write_whitenoise',         type='string'}
schema:add{name='write_runpbic',       type='string'}
schema:add{name='write_powerspectrum', type='string'}
schema:add{name='write_snapshot',      type='string'}
schema:add{name='write_nonlineark',      type='string'}
schema:add{name='write_runpb_snapshot', type='string'}
schema:add{name='cola_stdda',           type='boolean'}
schema:add{name='enforce_broadband_mode',  type='string',
            default='none', choices={'pm', 'linear', '2lpt', 'za', 'none'}}
schema:add{name='enforce_broadband_kmax',  type='number', default=4}
schema:add{name='za',                      type='boolean', default=false}
schema:add{name='kernel_type',             type='string',
        default="3_4", choices={'3_4', '5_4', 'eastwood', '3_2'}}
schema:add{name='dealiasing_type',             type='string',
        default="none", choices={'none', 'gaussian', 'twothird'}}
schema:add{name='inverted_ic',             type='boolean', default=false}
schema:add{name='remove_cosmic_variance',  type='boolean', default=false}

fastpm.FORCE_MODE_PM = "pm"
fastpm.FORCE_MODE_ZOLA = "zola"
fastpm.FORCE_MODE_COLA = "cola"

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

return "fastpm", fastpm
