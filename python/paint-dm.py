from nbodykit.lab import FFTPower, BigFileCatalog
from nbodykit import setup_logging
import numpy
import argparse
import warnings
import os

# usage:
#
# python paint-dm.py output --nmesh x --nn x dmcatalog

# for matter
# python paint-dm.py ../fastpm_1.0000/1-mesh ../tests/nbodykit/fastpm_1.0000/ --dataset=1

ap = argparse.ArgumentParser()
ap.add_argument("output", help='e.g. power.json (FFTPower.load) or power.txt (numpy.loadtxt)')
ap.add_argument("--output-dataset", default=None, help='output dataset, default to N%04d where %04d is nmesh ')
ap.add_argument("--nmesh", type=int, default=256, help='mesh resolution')
ap.add_argument("--verbose", action='store_true', default=False, help='print progress')

cat_ap = argparse.ArgumentParser()

cat_ap.add_argument("catalog", help='e.g. fastpm_1.0000')
cat_ap.add_argument("--dataset", default='1', help='data set to select; for a dm catalog, use 1')

ns, args = ap.parse_known_args()

ns1 = cat_ap.parse_args(args)

def read_cat(ns):
    cat = BigFileCatalog(ns.catalog, header='Header', dataset=ns.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3
    return cat

def main(ns, ns1):
    if ns.verbose:
        setup_logging('info')

    cat1 = read_cat(ns1)
    mesh1 = cat1.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)

    if ns.output_dataset is None:
        ns.output_dataset = 'N%04d' % ns.nmesh

    mesh1.save(ns.output, dataset=ns.output_dataset)

main(ns, ns1)
