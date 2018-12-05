from nbodykit.lab import FFTPower, BigFileCatalog, BigFileMesh
from nbodykit import setup_logging
import numpy
import argparse
import warnings
from mpi4py import MPI
import os

# usage:
#
# python cutslice.py output halocatalog [...] -- catalog [...]
#
# Cut a slice around a halo in halocatalog from catalog. Save all particles in the slice to output

# example
#
# python cutslice.py slices ../tests/nbodykit/fof_1.0000/ --dataset=LL-0.200 -- ../tests/nbodykit/fastpm_1.0000/ --dataset=1

ap = argparse.ArgumentParser()
ap.add_argument("output", help='')
ap.add_argument("--output-dataset", default=None, help='')
ap.add_argument("--haloid", default=5, type=int)
ap.add_argument("--los", default='z', help='los direction')
ap.add_argument("--verbose", action='store_true', default=False, help='print progress')
ap.add_argument("--thickness", default=10, type=float, help='thickness in Mpc/h')

cat_ap = argparse.ArgumentParser()

cat_ap.add_argument("catalog", help='e.g. fastpm_1.0000 or fof_1.0000')
cat_ap.add_argument("--dataset", default='LL-0.200', help='data set to select; for a dm catalog, use 1 for a halo catalog, usually LL-0.200')

ns, args = ap.parse_known_args()

split = args.index('--')
ns1 = cat_ap.parse_args(args[:split])
ns2 = cat_ap.parse_args(args[split+1:])

def read_cat(ns):
    cat = BigFileCatalog(ns.catalog, header='Header', dataset=ns.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3

    return cat

def main(ns, ns1, ns2):

    if ns.output_dataset is None:
        ns.output_dataset = '%sS-HID-%04d' % (ns2.dataset, ns.haloid)

    los = dict(x=[1, 0, 0], y=[0, 1, 0], z=[0, 0, 1])[ns.los]

    if ns.verbose:
        setup_logging('info')

    cat1 = read_cat(ns1)
    cat2 = read_cat(ns2)

    mask = (cat1.Index == ns.haloid).nonzero()[0]
    pos = (cat1['Position'][mask]).compute()
    if len(pos) > 0:
        pos = pos[0]
    else:
        pos = numpy.zeros((3))


    BoxSize = cat2.attrs['BoxSize']

    pos = cat1.comm.allreduce(pos) % BoxSize

    if cat1.comm.rank == 0:
        cat1.logger.info('Center position for halo %d is at %s' % (ns.haloid, str(pos)))

    r = (cat2['Position'] - pos)

    r = r + ((r > BoxSize * 0.5) * -BoxSize
          + (r < -BoxSize * 0.5) * BoxSize)

    r = (r * los).sum(axis=-1) ** 2
    sel = r < (ns.thickness * 0.5) ** 2

    catsel = cat2[sel]

    if cat2.comm.rank == 0:
        cat2.logger.info('Selected %d objects' % catsel.csize)

    catsel.attrs['BoxCenter'] = pos

    columns = sorted(set(catsel.columns) - set(['Weight', 'Selection', 'Value']))

    catsel.save(ns.output, columns=columns, dataset=ns.output_dataset)

    if catsel.comm.rank == 0:
        catsel.logger.info('saved to %s : %s' % (ns.output, ns.output_dataset))
main(ns, ns1, ns2)
