from nbodykit.lab import BigFileCatalog
from nbodykit.lab import FOF
from nbodykit.lab import HaloCatalog
from nbodykit.lab import KDDensity
from argparse import ArgumentParser 
from nbodykit.cosmology import Planck15
import numpy

def main():
    ap = ArgumentParser()
    ap.add_argument('fpm', help='e.g. /scratch/fpm_0.1000/')
    ap.add_argument('ll', type=float, help='e.g. 0.2 or 0.168') 
    ap.add_argument('--with-peak', help='Find Peaks KDDensity estimation (slow)', default=True)
    ap.add_argument('fof', help='e.g. /scratch/fpm_0.1000/fof . Will write to {fof}/{ll}')
    ap.add_argument('--nmin', type=int, default=20, help='min number of particles to be in the catalogue')

    ns = ap.parse_args()

    cat = BigFileCatalog(ns.fpm, header='Header', dataset='1/')

    cat.attrs['BoxSize']  = numpy.ones(3) * cat.attrs['BoxSize'][0]
    cat.attrs['Nmesh']  = numpy.ones(3) * cat.attrs['NC'][0]

    cosmo = Planck15.match(Omega0_m=cat.attrs['OmegaM'])

    M0 = cat.attrs['OmegaM'][0] * 27.75 * 1e10 * cat.attrs['BoxSize'].prod() / cat.csize

    if cat.comm.rank == 0:
        print('BoxSize', cat.attrs['BoxSize'])
        print('Nmesh', cat.attrs['Nmesh'])
        print('Mass of a particle', M0)
        print('OmegaM', cosmo.Om0)


    if ns.with_peak:
        cat['Density'] = KDDensity(cat).density

    fof = FOF(cat, linking_length=ns.ll, nmin=ns.nmin)

    if ns.with_peak:
        features = fof.find_features(peakcolumn='Density')
    else:
        features = fof.find_features(peakcolumn=None)

    features['Mass'] = M0 * features['Length']
    if fof.comm.rank == 0:
        print('Total number of features found', features.csize)
        print('Saving columns', features.columns)

    features.save(ns.fof + '/%0.3f' % ns.ll, features.columns)


main()
