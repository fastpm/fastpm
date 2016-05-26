
----------------------------------------------------
-- This is the main LUA runtime library of FastPM.
--
-- Author: Yu Feng <rainwoodman@gmail.com> 2016
----------------------------------------------------

function _main(filename, runmain, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    logspace = fastpm.logspace
    linspace = fastpm.linspace

    return config.parse(fastpm.schema, filename, runmain, {...})
end
