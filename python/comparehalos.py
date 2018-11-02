from nbodykit.lab import FFTPower, BigFileCatalog
from nbodykit import setup_logging
import numpy
import argparse
import warnings
from mpi4py import MPI
import os

# usage:
#
# python halobias.py output --nmin x --nmax x --nn x dmcatalog [...] -- halocatalog [...]
#
# measure the bias and growth rate from kaiser model by cross correlating halos of different size
# to the dark matter.
#
# example
#
# python halobias.py test.json --with-plot ../tests/nbodykit/fastpm_1.0000/ --dataset=1 -- ../tests/nbodykit/fof_1.0000/ --dataset=LL-0.200

# for matter
# python halobias.py test.json --with-plot ../tests/nbodykit/fastpm_1.0000/ --dataset=1 -- ../tests/nbodykit/fof_1.0000/ --dataset=1

ap = argparse.ArgumentParser()
ap.add_argument("output", help='e.g. power.json (FFTPower.load) or power.txt (numpy.loadtxt)')
ap.add_argument("--nmin", default=8, type=int)
ap.add_argument("--kmax", default=None, type=float, help="cut to stop using kmax, scale where kaiser is bad")
ap.add_argument("--nmax", default=1000, type=int)
ap.add_argument("--nn", default=10, type=int)
ap.add_argument("--unique-k", action='store_true', default=False, help='compute for all unique k values.')
ap.add_argument("--nmesh", type=int, default=256, help='mesh resolution')
ap.add_argument("--verbose", action='store_true', default=False, help='print progress')

cat_ap = argparse.ArgumentParser()

cat_ap.add_argument("catalog", help='e.g. fastpm_1.0000 or fof_1.0000')
cat_ap.add_argument("--dataset", default='LL-0.200', help='data set to select; for a dm catalog, use 1 for a halo catalog, usually LL-0.200')

ns, args = ap.parse_known_args()

if '--' in args:
    split = args.index('--')
    ns1 = cat_ap.parse_args(args[:split])
    ns2 = cat_ap.parse_args(args[split+1:])
else:
    ns1 = cat_ap.parse_args(args)
    ns2 = ns1

def read_cat(ns, nmin=None):
    cat = BigFileCatalog(ns.catalog, header='Header', dataset=ns.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3

    if nmin is not None and nmin != 0:
        sel = True
        sel = sel & (cat['Length'] >= nmin)

        cat['Selection'] = sel
#        cat = cat[sel]

    cat['RSDPosition'] = cat['Position'] + cat.attrs['RSDFactor'] * cat['Velocity'] * [0, 0, 1]
    return cat

# use bisection to find the nmin to match the number of nsel

def read_cat_nsel(ns, nsel, nmin0, nmin1):
    cat = BigFileCatalog(ns.catalog, header='Header', dataset=ns.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3

    if 'Length' in cat.columns:
        while nmin1 - nmin0 > 1:
            nminc = (nmin1 + nmin0) / 2

            sel = True
            sel = sel & (cat['Length'] >= nminc)

            nsel1 = cat.comm.allreduce(sel.sum().compute())
            if nsel1 < nsel: # too few
                nmin1 = nminc
            else:
                nmin0 = nminc

        if cat.comm.rank == 0:
            print('found nmin', nmin1, nmin0, 'nsel is', nsel1, 'target is', nsel)

        cat['Selection'] = sel
    cat['RSDPosition'] = cat['Position'] + cat.attrs['RSDFactor'] * cat['Velocity'] * [0, 0, 1]
    return cat

def main(ns, ns1, ns2):
    if ns.verbose:
        setup_logging('info')

    if ns.unique_k:
        dk = 0
    else:
        dk = None

    cat1 = read_cat(ns1)
    cat2 = read_cat(ns2)

    nmin = numpy.unique(numpy.int32(numpy.logspace(numpy.log10(ns.nmin), numpy.log10(ns.nmax), ns.nn, endpoint=True)))
    if 'Length' in cat1.columns:
        nmin0 = cat1.comm.allreduce(cat1['Length'].min().compute() if cat1.size > 0 else 10000000, MPI.MIN)
        nmax0 = cat1.comm.allreduce(cat1['Length'].max().compute() if cat1.size > 0 else 0, MPI.MAX)
        nmin = nmin[nmin >= nmin0]
        nmin = nmin[nmin < nmax0]
    else:
        nmin = [0]

    if 'Length' in cat2.columns:
        nmin2 = cat2.comm.allreduce(cat2['Length'].min().compute() if cat2.size > 0 else 10000000, MPI.MIN)
        nmax2 = cat2.comm.allreduce(cat2['Length'].max().compute() if cat2.size > 0 else 0, MPI.MAX)
    else:
        nmin2 = 0
        nmax2 = 1

    if cat1.comm.rank == 0:
        os.makedirs(os.path.dirname(ns.output), exist_ok=True)

    for nmin1 in nmin:
        cat1 = read_cat(ns1, nmin1)
        nsel = cat1.comm.allreduce(cat1['Selection'].sum().compute())
        cat2 = read_cat_nsel(ns2, nsel, nmin2, nmax2)

        mesh1 = cat1.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh, position='RSDPosition')
        mesh2 = cat2.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh, position='RSDPosition')

        r1 = FFTPower(mesh1, second=mesh1, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax)
        r2 = FFTPower(mesh2, second=mesh2, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax)
        rx = FFTPower(mesh1, second=mesh2, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax)

        save_bs(ns.output, 'nmin-%05d-r1' % nmin1, r1)
        save_bs(ns.output, 'nmin-%05d-r2' % nmin1, r2)
        save_bs(ns.output, 'nmin-%05d-rx' % nmin1, rx)
        if cat1.comm.rank == 0:
            print("nmin = ", nmin1, "finished")

def save_bs(filename, dataset, r):

    basename = filename.rsplit('.', 1)[0]
    if filename.endswith('.json'):
        r.save(basename + '-%s.json' % dataset)
    elif filename.endswith('.txt'):
        if r.comm.rank == 0:
            for var in r.power.data.dtype.names:
                numpy.savetxt(basename + '-%s-%s.txt' % (dataset, var),
                    r.power[var].real
                )

main(ns, ns1, ns2)
