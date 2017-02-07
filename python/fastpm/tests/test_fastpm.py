from __future__ import print_function
from mpi4py_test import MPITest
from numpy.testing import assert_allclose
import numpy

import fastpm
from astropy.cosmology import FlatLambdaCDM
from pmesh.pm import ParticleMesh

import fastpm.operators as operators

def _addampl(white):
    """ returns dlin_k and amplitude """
    ampl = white.copy()

    def kernel(k, v):
        kk = sum(ki ** 2 for ki in k)
        ka = kk ** 0.5
        p = (ka / 0.1) ** -2 * .4e4 * (1.0 / v.BoxSize).prod()
        p[ka == 0] = 0
        return p ** 0.5 + p ** 0.5 * 1j

    ampl.apply(kernel, out=Ellipsis)
    white[...] *= ampl.real

    return white, ampl

@MPITest([1, 4])
def test_lpt(comm):

    from pmesh.pm import ParticleMesh

    pm = ParticleMesh(BoxSize=128.0, Nmesh=(4, 4), comm=comm)
    vm = fastpm.Evolution(pm, shift=0.5)

    dlink = pm.generate_whitenoise(1234, mode='complex')
    code = vm.code()
    code.Displace(D1=1.0, v1=0, D2=0.0, v2=0.0)
    code.Chi2(variable='s')

    dlink, ampl = _addampl(dlink)

    _test_model(code, dlink, ampl)

@MPITest([1, 4])
def test_prior(comm):

    from pmesh.pm import ParticleMesh

    pm = ParticleMesh(BoxSize=128.0, Nmesh=(4, 4), comm=comm)
    vm = fastpm.Evolution(pm, shift=0.5)

    dlink = pm.generate_whitenoise(1234, mode='complex')
    code = vm.code()
    code.Prior(powerspectrum=lambda k : 1.0)
    code.CopyVariable(a='prior', b='chi2')
    dlink, ampl = _addampl(dlink)

    _test_model(code, dlink, ampl)

@MPITest([1, 4])
def test_gravity(comm):
    from pmesh.pm import ParticleMesh
    import fastpm.operators as operators

    pm = ParticleMesh(BoxSize=4.0, Nmesh=(4, 4), comm=comm, method='cic', dtype='f8')

    vm = fastpm.Evolution(pm, shift=0.5)

    dlink = pm.generate_whitenoise(12345, mode='complex')

    code = vm.code()
    code.Displace(D1=1.0, v1=0, D2=0.0, v2=0.0)
    code.Force(factor=0.1)
    code.Chi2(variable='f')

    # FIXME: without the shift some particles have near zero dx1.
    # or near 1 dx1.
    # the gradient is not well approximated by the numerical if
    # any of the left or right value shifts beyond the support of
    # the window.
    #

    dlink, ampl = _addampl(dlink)
    _test_model(code, dlink, ampl)

def _test_model(code, dlin_k, ampl):

    def objective(dlin_k):
        init = {'dlin_k': dlin_k}
        return code.compute('chi2', init, monitor=None)

    def gradient(dlin_k):
        tape = code.vm.tape()
        init = {'dlin_k': dlin_k}
        code.compute(['chi2'], init, tape=tape, monitor=None)
        gcode = code.vm.gradient(tape)
        init = {'_chi2' : 1}
        return gcode.compute('_dlin_k', init, monitor=None)

    y0 = objective(dlin_k)
    yprime = gradient(dlin_k)

    num = []
    ana = []
    print('------')
    for ind1 in numpy.ndindex(*(list(dlin_k.cshape) + [2])):
        dlinkl = dlin_k.copy()
        dlinkr = dlin_k.copy()
        old = dlin_k.cgetitem(ind1)
        pert = ampl.cgetitem(ind1) * 1e-5
        left = dlinkl.csetitem(ind1, old - pert)
        right = dlinkr.csetitem(ind1, old + pert)
        diff = right - left
        yl = objective(dlinkl)
        yr = objective(dlinkr)
        grad = yprime.cgetitem(ind1)
        print(ind1, old, pert, yl, yr, grad * diff, yr - yl)
        ana.append(grad * diff)
        num.append(yr - yl)
    print('------')

    assert_allclose(num, ana, rtol=1e-3, atol=1e-10)

@MPITest([1, 4])
def test_kdk(comm):
    cosmo = FlatLambdaCDM(Om0=0.3, H0=70, Tcmb0=0)

    pm = ParticleMesh(BoxSize=128.0, Nmesh=(4,4,4), comm=comm, dtype='f8')
    vm = fastpm.KickDriftKick(pm, shift=0.5)
    dlink = pm.generate_whitenoise(12345, mode='complex', unitary=True)

    dlink, ampl = _addampl(dlink)

    data = dlink.c2r()
    data[...] = 0
    sigma = data.copy()
    sigma[...] = 1.0

    code = vm.simulation(cosmo, 0.1, 1.0, 5)
    code.Paint()
    code.Subtract(data_x=data, sigma_x=sigma)
    code.Chi2(variable='mesh')

    _test_model(code, dlink, ampl)

@MPITest([1, 4])
def test_resample(comm):
    if comm.size == 4:
        comm.barrier()
    from pmesh.pm import ParticleMesh

    # testing the gradient of down sampling, which we use the most.
    pm = ParticleMesh(BoxSize=128.0, Nmesh=(8, 8), comm=comm)
    pm1 = ParticleMesh(BoxSize=128.0, Nmesh=(4, 4), comm=comm)
    vm = fastpm.Evolution(pm, shift=0.5)

    dlink = pm.generate_whitenoise(1234, mode='complex')

    dlink, ampl = _addampl(dlink)
    print(dlink.value.round(2))
    dlink1 = pm1.create(mode='complex')
    dlink.resample(dlink1)
    print(dlink1.value.round(2))
    code = vm.code()
    code.C2R()
    code.Resample(pm=pm1)
    code.Chi2(variable='mesh')

    _test_model(code, dlink, ampl)
