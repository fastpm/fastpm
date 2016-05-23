import numpy
import os

class DumpFile(object):
    def __init__(self, path):
        self.path = path
        self.filenames = []
        i = 0
        while True:
            fn = '%s.%03d' % (path, i)
            if not os.path.exists(fn):
                if i == 0:
                    fn = path
                    if not os.path.exists(fn):
                        raise OSError("File not found")
                else:
                    break
            self.filenames.append(fn)
            i = i + 1

    def as_real(self):
        shape = self._guess_size('real')
        data = numpy.zeros(shape, dtype='f4')
        for fn in self.filenames:
            geo = fn + '.geometry'
            strides, offset, shape = self._parse_geo(geo, 'real')
            d = numpy.fromfile(fn, dtype='f4')
            ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])
            d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * 4)
            data[ind] = d
        return data

    def as_complex(self):
        shape = self._guess_size('complex')
        data = numpy.zeros(shape, dtype='complex64')
        for fn in self.filenames:
            geo = fn + '.geometry'
            strides, offset, shape = self._parse_geo(geo, 'complex')
            d = numpy.fromfile(fn, dtype='complex64')
            ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])
            d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * 8)
            data[ind] = d
        return data

    def _parse_geo(self, geofn, mode):
        if mode == 'real':
            strides = numpy.loadtxt(open(geofn).readlines()[3].split()[1:], dtype=int)
            offset = numpy.loadtxt(open(geofn).readlines()[1].split()[1:], dtype=int)
            shape = numpy.loadtxt(open(geofn).readlines()[2].split()[1:], dtype=int)
        elif mode == 'complex':
            strides = numpy.loadtxt(open(geofn).readlines()[7].split()[1:], dtype=int)
            offset = numpy.loadtxt(open(geofn).readlines()[5].split()[1:], dtype=int)
            shape = numpy.loadtxt(open(geofn).readlines()[6].split()[1:], dtype=int)    

        return strides, offset, shape

    def _guess_size(self, mode):
        size = None
        for fn in self.filenames:
            geo = fn + '.geometry'
            strides, offset, shape = self._parse_geo(geo, mode)
            last = shape + offset
            if size is None:
                size = last
            else:
                size = numpy.maximum(size, last)
        return size
