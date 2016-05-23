----------------------------------------------------
-- This is the main LUA runtime library of FastPM.
--
-- Author: Yu Feng <rainwoodman@gmail.com> 2016
----------------------------------------------------

local function _runmain(filename, runmain, ...)

-- Parse(run) a lua file
--
-- This is the first function we land off from the
-- C part of fastpm/fastpm-lua.
--
-- filename: the filename of the lua file
-- runmain: if true, the main function in the file is executed
-- ... : args

    local namespace = setmetatable({}, {__index=_G})

    -- alias for easy access to logspace and linspace
    -- note they are added as global variables in _G

    logspace = fastpm.logspace
    linspace = fastpm.linspace

    namespace['__file__'] = filename
    namespace['args'] = {...}

    local param, err = loadfile(filename, 'bt', namespace)
    if param == nil then
        error(err)
    end
    param()

    -- prune main function from the parameter file
    local main = namespace['main']
    namespace['main'] = nil

    if main ~= nil then
        if runmain then
            main()
        end
    end

    local namespace2 = fastpm.schema.bind(namespace)

    for k,t in pairs(namespace2) do
        if type(t) == 'function' then
            error(string.format("in `%s`, function `%s` shall be local", filename, k))
        end
    end

    local ret, err = dump.tostring(namespace2)
    if ret == nil then
        error(err)
    end
    return ret
end

return "_runmain", _runmain
