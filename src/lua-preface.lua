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

function parse_file(filename)
    env = setmetatable({}, {__index=_G})
    param = loadfile(filename, 'bt', env)
    param()
    return dump.toscript(env)
end
