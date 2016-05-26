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

def power(f1, f2=None, boxsize=1.0, return_modes=False):
    """ stupid power spectrum calculator.

        f1 f2 must be density fields in configuration or fourier space.

        For convenience if f1 is strictly overdensity in fourier space,
        (zero mode is zero) the code still works.

        Returns k, p or k, p, N if return_modes is True

    """
    def tocomplex(f1):
        if f1.dtype.kind == 'c':
            return f1
        else:
            return numpy.rfftn(f1)
    f1 = tocomplex(f1)
    if f1[0, 0, 0] != 0.0:
        f1 /= abs(f1[0, 0, 0])

    if f2 is None:
        f2 = f1

    if f2 is not f1:
        f2 = tocomplex(f2)
        if f2[0, 0, 0] != 0.0:
            f2 /= abs(f2[0, 0, 0])

    def fftk(shape, boxsize, symmetric=True):
        k = []
        for d in range(len(shape)):
            kd = numpy.fft.fftfreq(shape[d])
            kd *= 2 * numpy.pi / boxsize * shape[d]
            kdshape = numpy.ones(len(shape), dtype='int')
            if symmetric and d == len(shape) -1:
                kd = kd[:shape[d]//2 + 1]
            kdshape[d] = len(kd)
            kd = kd.reshape(kdshape)
            k.append(kd)
        kk = sum([i ** 2 for i in k])
        return kk ** 0.5

    x = (f1 * f2.conjugate()).real

    k = fftk(f1.shape, boxsize, symmetric=False)
    w = numpy.ones(k.shape, dtype='f4') * 2
    w[:, :, 0] = 1
    x *= w
    edges = numpy.arange(f1.shape[-1]) * 2 * numpy.pi / boxsize
    H, edges = numpy.histogram(k.flat, weights=x.flat, bins=edges)
    N, edges = numpy.histogram(k.flat, weights=w.flat, bins=edges)
    center = 0.5 * (edges[1:] + edges[:-1])
    if return_modes:
        return center, H / N *boxsize ** 3, N
    else:
        return center, H / N *boxsize ** 3

