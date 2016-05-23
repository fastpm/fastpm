a, config = dofile('lua-runtime-config.lua')
a, fastpm = dofile('lua-runtime-fastpm.lua')

module = arg[1]

filename_c = module .. '.c'
filename_h = module .. '.h'

stream_h, stream_c = config.compile(fastpm.schema,
{
    prefix ="lua_config",
    global_headers = {
        'fastpm/libfastpm.h',
    },
    local_headers = {
    }
}
)

file_c = io.open(filename_c, 'w')
file_c:write(stream_c)
file_c:close()

file_h = io.open(filename_h, 'w')
file_h:write(stream_h)
file_h:close()
