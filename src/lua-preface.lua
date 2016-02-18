-- Define a few time stepping schemes

function linspace(start, e, N)
    local r = {} 
    N1 = N + 1 
    for i=1,N1 do
        r[i] = 1.0 * (e - start) * (i - 1) / N + start
    end 
    r[N1] = e 
    return r
end 
function logspace(start, e, N)
    local r 
    r = linspace(start, e, N)
    for i, j in pairs(r) do
        r[i] = math.pow(10, j)
    end 
    return r
end
function blendspace(start, e, a1, a2)
    local r = {}
    a = start
    i = 1 
    while a < e do 
        r[i] = a 
        dlna = math.pow(math.pow(1/a1, 2) + math.pow(a/a2, 2), -0.5)
        a = math.exp(math.log(a) + dlna) 
        i = i + 1 
    end 
    r[i] = e 
    return r
end 

function parse_file(filename, mode)
    local namespace = setmetatable({}, {__index=_G})
    namespace['__file__'] = filename
    namespace['__mode__'] = mode
    local param = loadfile(filename, 'bt', namespace)
    if param == nil then
        error(string.format("Could not open file %s", filename))
    end
    param()
    local required = {
        nc = 'number',
        boxsize = 'number',
        time_step = 'table',
        output_redshifts = 'table',
        omega_m = 'number',
        h = 'number',
        pm_nc_factor = 'table',
        change_pm = 'table',
        np_alloc_factor = 'number',
        force_mode = 'string',
    }
    local optional = {
        read_grafic = 'string',
        read_noisek = 'string',
        read_noise = 'string',
        read_runpbic = 'string',
        read_powerspectrum = 'file',
        read_sigma8 = 'number',
        random_seed = 'number',
        read_write_powerspectrum = 'string',
        read_write_runpb_snapshot = 'string',
        read_write_snapshot = 'string',
        read_write_noisek = 'string',
        read_write_noise = 'string',
        cola_stdda = 'boolean',
        enforce_broadband = 'boolean',
        enforce_broadband_kmax = 'number',
        za = 'boolean',
    }
    check_schema(namespace, required, true)
    check_schema(namespace, optional, false)
    return dump.tostring(namespace)
end

function check_schema(namespace, schema, required)
    for key, t in pairs(schema) do
        if namespace[key] == nil then
            if required then
                error(string.format("key `%s` is required but not defined", key))
            end 
        else
            if t == 'file' then
                t1 = 'string'
            else
                t1 = t
            end
            if type(namespace[key]) ~= t1 then
                error(string.format("key `%s' is not of type `%s', but `%s' ", key, t1, type(namespace[t])))
            end
            if t == 'file' then
                file = io.open(namespace[key], 'r')
                if file == nil then
                    error(string.format("file `%s' is not readable", namespace[key]))
                end
                file:close()
            end
        end
    end 
end
