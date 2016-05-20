import numpy
import os

def read_real(filename, size):
    data = numpy.zeros((size, size, size), dtype='f4')
    i = 0
    while True:
        fn = '%s.%03d' % (filename, i)
        geofn = '%s.%03d.geometry' % (filename, i)
        if not os.path.exists(fn):
            if i == 0:
                fn = filename
                geofn = '%s.geometry' % filename
                if not os.path.exists(fn):
                    raise OSError("File not found")
            else:
                break
        d = numpy.fromfile(fn, dtype='f4')
        strides = numpy.loadtxt(open(geofn).readlines()[3].split()[1:], dtype=int)
        offset = numpy.loadtxt(open(geofn).readlines()[1].split()[1:], dtype=int)
        shape = numpy.loadtxt(open(geofn).readlines()[2].split()[1:], dtype=int)    
        ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])        
        d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * 4)
        #print ind, d.shape, d.max()
        data[ind] = d
        i = i + 1
    return data

def read_complex(filename, size):
    data = numpy.zeros((size, size, size // 2 + 1), dtype='complex64')
    i = 0
    while True:
        fn = '%s.%03d' % (filename, i)
        geofn = '%s.%03d.geometry' % (filename, i)
        #print fn
        if not os.path.exists(fn):
            if i == 0:
                fn = filename
                geofn = '%s.geometry' % filename
                if not os.path.exists(fn):
                    raise OSError("File not found")
            else:
                break
        d = numpy.fromfile(fn, dtype='complex64')
        strides = numpy.loadtxt(open(geofn).readlines()[7].split()[1:], dtype=int)
        offset = numpy.loadtxt(open(geofn).readlines()[5].split()[1:], dtype=int)
        shape = numpy.loadtxt(open(geofn).readlines()[6].split()[1:], dtype=int)    
        ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])        
        d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * 8)
        #d = d.view(dtype='complex64')
        #print ind, d.shape, d.max()
        data[ind] = d
        i = i + 1
    return data
