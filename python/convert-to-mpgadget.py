""" This script converts a FastPM snapshot to a MP-Gadget snapshot.

    Three blocks are converted, position, velocity and ID. A new Mass column is added.

"""

from argparse import ArgumentParser
import numpy
import bigfile

ap = ArgumentParser()
ap.add_argument('source', help='FastPM snapshot (bigfile)')
ap.add_argument('dest', help='Gadget filename base; dir will be created on the fly. ')

def main(ns):
    def copy(bbi, bbo, chunksize=1024*1024):
        for i in range(0, bbi.size, chunksize):
            bbo.write(i, bbi[i:i+chunksize])

    with bigfile.File(ns.source) as bfi, \
         bigfile.File(ns.dest, create=True) as bfo:
        header = bfi['Header']
        with bfi['1/Position'] as bbi:
            with bfo.create('1/Position', dtype=bbi.dtype, size=bbi.size, Nfile=bbi.Nfile) as bbo:
                copy(bbi, bbo)
            npart = bbi.size
            Nfile = bbi.Nfile
        with bfi['1/Velocity'] as bbi:
            with bfo.create('1/Velocity', dtype=bbi.dtype, size=bbi.size, Nfile=bbi.Nfile) as bbo:
                copy(bbi, bbo)
        with bfi['1/ID'] as bbi:
            with bfo.create('1/ID', dtype=bbi.dtype, size=bbi.size, Nfile=bbi.Nfile) as bbo:
                copy(bbi, bbo)

        with bfo.create('1/Mass', dtype='f4', size=npart, Nfile=Nfile) as bbo:
            mass = numpy.broadcast_to(header.attrs['MassTable'][1], npart)
            copy(mass, bbo)

        with bfo.create('Header') as bbo:
            for key in header.attrs:
                bbo.attrs[key] = header.attrs[key]

if __name__ == '__main__':
    ns = ap.parse_args()
    main(ns)
