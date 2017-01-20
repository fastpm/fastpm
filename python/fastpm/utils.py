import numpy
import os

class DumpFile(object):
    def __init__(self, path, dtype):
        self.path = path
        self.filenames = []
        dtype = numpy.dtype(dtype)
        if dtype == numpy.dtype('f8'):
            self.rdtype = numpy.dtype('f8')
            self.cdtype = numpy.dtype('complex128')
        else:
            self.rdtype = numpy.dtype('f4')
            self.cdtype = numpy.dtype('complex64')

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
        data = numpy.zeros(shape, dtype=self.rdtype)
        for fn in self.filenames:
            geo = fn + '.geometry'
            strides, offset, shape = self._parse_geo(geo, 'real')
            d = numpy.fromfile(fn, dtype=self.rdtype)
            ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])
            d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * self.rdtype.itemsize)
            data[ind] = d
        return data

    def as_complex(self):
        shape = self._guess_size('complex')
        data = numpy.zeros(shape, dtype='complex64')
        for fn in self.filenames:
            geo = fn + '.geometry'
            strides, offset, shape = self._parse_geo(geo, 'complex')
            d = numpy.fromfile(fn, dtype=self.cdtype)
            ind = tuple([slice(x, x+o) for x, o in zip(offset, shape)])
            d = numpy.lib.stride_tricks.as_strided(d, shape=shape, strides=strides * self.cdtype.itemsize)
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

def power(f1, f2=None, boxsize=1.0, average=True):
    """ stupid power spectrum calculator.

        f1 f2 must be density fields in configuration or fourier space.

        For convenience if f1 is strictly overdensity in fourier space,
        (zero mode is zero) the code still works.

        Returns k, p or k, p * n, N if average is False

    """
    def tocomplex(f1):
        if f1.dtype.kind == 'c':
            return f1
        else:
            return numpy.fft.rfftn(f1)

    f1 = tocomplex(f1)
    if f1[0, 0, 0] != 0.0:
        f1 /= abs(f1[0, 0, 0])

    if f2 is None:
        f2 = f1

    if f2 is not f1:
        f2 = tocomplex(f2)
        if f2[0, 0, 0] != 0.0:
            f2 /= abs(f2[0, 0, 0])

    def fftk(shape, boxsize):
        k = []
        for d in range(len(shape)):
            kd = numpy.arange(shape[d])

            if d != len(shape) - 1:
                kd[kd > shape[d] // 2] -= shape[d] 
            else:
                kd = kd[:shape[d]]

            kdshape = numpy.ones(len(shape), dtype='int')
            kdshape[d] = len(kd)
            kd = kd.reshape(kdshape)

            k.append(kd)
        return k

    k = fftk(f1.shape, boxsize)

    def find_root(kk):
        solution = numpy.int64(numpy.sqrt(kk) - 2).clip(0)
        solution[solution < 0] = 0
        mask = (solution + 1) ** 2 < kk
        while(mask.any()):
            solution[mask] += 1
            mask = (solution + 1) ** 2 <= kk

        return solution

    ksum = numpy.zeros(f1.shape[0] //2, 'f8')
    wsum = numpy.zeros(f1.shape[0] //2, 'f8')
    xsum = numpy.zeros(f1.shape[0] //2, 'f8')

    for i in range(f1.shape[0]):
        kk = k[0][i] ** 2 + k[1] ** 2 + k[2] ** 2

        # remove unused dimension
        kk = kk[0]

        d = find_root(kk)

        w = numpy.ones(d.shape, dtype='f4') * 2
        w[..., 0] = 1
        w[..., -1] = 1

        xw = abs(f1[i] * f2[i].conjugate()) * w

        kw = kk ** 0.5 * 2 * numpy.pi / boxsize * w

        ksum += numpy.bincount(d.flat, weights=kw.flat, minlength=f1.shape[0])[:f1.shape[0] // 2]
        wsum += numpy.bincount(d.flat, weights=w.flat, minlength=f1.shape[0])[:f1.shape[0] // 2]
        xsum += numpy.bincount(d.flat, weights=xw.flat, minlength=f1.shape[0])[:f1.shape[0] // 2]

    center = ksum / wsum

    if not average:
        return center, xsum * boxsize**3, wsum
    else:
        return center, xsum / wsum * boxsize **3

def fftdown(field, size):
    """ Down samples a fourier space field. Size can be scalar or vector.

        Hermitian should be handled correctly. But I have only tested this
        on hermitian compressed fields.
    """

    ashape = numpy.array(field.shape)
    asize = ashape.copy()
    fullaxes = ashape == ashape[0]
    asize[:] = size
    if numpy.isscalar(size):
        asize[~fullaxes] = size // 2 + 1

    t = numpy.fft.fftshift(field, axes=fullaxes.nonzero()[0])
    left = (ashape - asize) // 2
    left[~fullaxes] = 0
    right = left + asize

    sl = [ slice(left[i], right[i]) for i in range(len(ashape))]
    t = t[sl]

    return numpy.fft.ifftshift(t, axes=fullaxes.nonzero()[0])

def complex_to_fastpm(fn, ds, complex, BoxSize):
    """
        >>> c = numpy.fft.rfftn(numpy.random.normal(size=(128, 128, 128)))

        >>> complex_to_fastpm("IC", "Copy", c, 128.)

        >>> a = bigfile.BigFile("IC")['Copy'][:]
        >>> print(a[1], c.ravel()[1])
        >>> print(a[300], c.ravel()[300])
    """
    import bigfile
    import numpy

    bf = bigfile.BigFile(fn, create=True)
    bb = bf.create(ds, complex.dtype, complex.size, Nfile=1)
    Nmesh = complex.shape[0]

    bb.attrs['Nmesh'] = Nmesh
    bb.attrs['BoxSize'] = BoxSize
    # This ensures we write the array as C-Contiguous
    bb.write(0, complex)
    # Thus we set the strides and shape accordingly for FastPM to pick up.
    bb.attrs['ndarray.ndim'] = 3
    bb.attrs['ndarray.shape'] = (Nmesh, Nmesh, Nmesh // 2 + 1)
    bb.attrs['ndarray.strides'] = (Nmesh * (Nmesh // 2 + 1), Nmesh // 2 + 1, 1)
