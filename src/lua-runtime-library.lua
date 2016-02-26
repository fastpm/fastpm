----------------------------------------------------
-- This is the main LUA runtime library of FastPM.
--
-- Author: Yu Feng <rainwoodman@gmail.com> 2016
----------------------------------------------------

function linspace(a, e, N)
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
function logspace(a, e, N)
-- a and end are in log10.
-- Returns N+1 elements, including e.
    local r 
    r = linspace(a, e, N)
    for i, j in pairs(r) do
        r[i] = math.pow(10, j)
    end 
    return r
end
function blendspace(a, e, a1, a2)
    local r = {}
    a = a
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

function parse_file(filename, runmain, ...)
-- Parse(run) a lua file
--
-- This is the first function we land off from the
-- C part of fastpm/fastpm-lua.
--
-- filename: the filename of the lua file
-- runmain: if true, the main function in the file is executed
-- ... : args

    local namespace = setmetatable({}, {__index=_G})
    namespace['__file__'] = filename
    namespace['args'] = {...}
    local param, err = loadfile(filename, 'bt', namespace)
    if param == nil then
        error(err)
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
    local main = namespace['main']
    if main ~= nil then
        if runmain then
            main()
        end
        namespace['main'] = nil
    end    
    for k,t in pairs(namespace) do
        if type(t) == 'function' then
            error(string.format("in `%s`, function `%s` shall be local", filename, k))
        end
    end
    local ret, err = dump.tostring(namespace)
    if ret == nil then
        error(err)
    end
    return ret
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
