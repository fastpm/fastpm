"""Read and save linear density field at lagrangian position of particles.

Examples

Using fastpm_1.0000's settings:

- read_powerspectrum = "powerspectrum.txt"
- nc = 128
- random seed = 100
- remove_cosmic variance = false

Save to 'linear' column of linear density of Lagrangian position computed from ID of each
particle, where powerspectrum.txt is at /home/yfeng1/source/fastpm/tests
(aka the current working directory of the original simulation). The dataset "1" is implied.

  python read-linear-q.py ../tests/nbodykit/fastpm_1.0000/ linear --cwd /home/yfeng1/source/fastpm/tests/

If we do not want to write to ../tests/nbodykit/fastpm_1.0000, then add --ocatalog=here.

For a read-only CrowCanyon-II dataset, something like this would work:

  python read-linear-q.py /.../readonly/5001-8192/lightcone/usmesh/ linear \
       --ocatalog=/global/cscratch1/sd/yfeng1/m3404/desi/5001-8192/lightcone/usmesh/

There was no need for cwd because we use absolute paths in CrowCanyon2.


This script will likely not work with ncdm, because the way ncdm particle ID relates to Q. 
A fix should be trivial.

"""

from nbodykit.lab import BigFileCatalog
import json
import numpy
import argparse
import warnings
from mpi4py import MPI
import os
from pmesh.pm import ParticleMesh
from nbodykit import setup_logging
import dask.array as da

ap = argparse.ArgumentParser()
ap.add_argument("catalog", help="e.g. fastpm_1.0000 or fof_1.0000")
ap.add_argument("ocolumn", help="e.g. linear_rho")
ap.add_argument("--nmesh", default=None, type=int, help="if not given, use nc of the params")
ap.add_argument("--cwd", default="", help="working directory to resolve relative paths in the fastpm parameters.")
ap.add_argument("--ocatalog", default=None, help="if given, do not write to the catalog")
ap.add_argument("--dataset", default="1", help="data set to use.")
ap.add_argument("--verbose", default=True, help="show logging info")

# FIXME: move this somewhere -- nbodykit or fastpm-python?
class FastPMParams(object):
    def __init__(self, paramstr, cwd=""):
        """
        Create fastpm parameters from a serialized parameter string.

        Parameters
        ----------
        paramstr: a python string, produced by fastpm using lua's dump module.
        cwd: the working directory for relative paths in the file.
        """

        parser = self.get_parser()
        self.params = parser.parseString(paramstr, parseAll=True).asDict()
        self.cwd = cwd
        
    def __getitem__(self, key):
        return self.params[key]

    def read_powerspectrum(self):
        filename = os.path.join(self.cwd, self.params['read_powerspectrum'])
        k, p = numpy.loadtxt(filename, unpack=True)
        if (k < 0).any():  # logscaled power spectrum
            warnings.warn("File %s is likely in log10, converting to usual k, p pairs." % filename,
                RuntimeWarning)
            return 10**k, 10**p
        return k, p
 
    @staticmethod
    def get_parser():
        """A Parser that parses the dumped ParamFile attribute by FastPM.

        This must be a result produced by the lua dump module. Must be
        a lua table. e.g.
        { a = 3, b = {0, 1, 2,} }

        (modified from jsonParser.py example under pyparsing)

        When using the parser, convert the attribute from an array of U1 to a string first.
        """
        import pyparsing as pp
        from pyparsing import pyparsing_common as ppc

        def make_keyword(kwd_str, kwd_value):
            return pp.Keyword(kwd_str).setParseAction(pp.replaceWith(kwd_value))

        TRUE = make_keyword("true", True)
        FALSE = make_keyword("false", False)
        NULL = make_keyword("nil", None)

        LBRACE, RBRACE, ASSIGN, COMMA = map(pp.Suppress, "{}=,")

        luaName = pp.Word(pp.alphas + "_", pp.alphanums + "_")
        luaString = pp.dblQuotedString().setParseAction(pp.removeQuotes)
        luaNumber = ppc.number()

        luaObject = pp.Forward()
        luaValue = pp.Forward()
        luaElements = pp.delimitedList(luaValue) + pp.Optional(COMMA)
        luaArray = pp.Group(LBRACE + pp.Optional(luaElements, []) + RBRACE)
        luaValue << (
            luaString | luaNumber | pp.Group(luaObject) | luaArray | TRUE | FALSE | NULL
        )
        memberDef = pp.Group(luaName + ASSIGN + luaValue)
        luaMembers = pp.delimitedList(memberDef) + pp.Optional(COMMA)
        luaObject << pp.Dict(LBRACE + pp.Optional(luaMembers) + RBRACE)
        return luaObject

def id2q(pid, strides, scale, shift):
    ipos = pid[:, None] // strides[None, :]
    # the size of the highest dimension is unknown.
    sizes = (strides[:-1] // strides[1:])
    ipos[:, 1:] = ipos[:, 1:] % sizes
    pos = (ipos + shift) * scale
    return pos.astype("f8")

def main():
    ns = ap.parse_args()
    if ns.verbose:
        setup_logging('info')
    cat = BigFileCatalog(ns.catalog, dataset="1")
    params = FastPMParams(''.join(cat.attrs['ParamFile']), cwd=ns.cwd)

    if params['sigma8'] != 0:
        raise ValueError("overriding sigma8 is not supported")

    nmesh = ns.nmesh or params['nc']
    if cat.comm.rank == 0:
        cat.logger.info("Nmesh = %d", nmesh)

    strides = cat.attrs["q.strides"]
    scale = cat.attrs["q.scale"]
    shift = cat.attrs["q.shift"]

    pm = ParticleMesh([nmesh] * 3, cat.attrs['BoxSize'], comm=cat.comm)

    k, Pk = params.read_powerspectrum()
    def pklin(k_):
        return  numpy.interp(k_, k, Pk)
    def tf(k_):
        return (pklin(k_.normp(2, zeromode=1.0)**0.5) / pm.BoxSize.prod()**0.5)

    dlin = pm.generate_whitenoise(params['random_seed'],
                                unitary=params['remove_cosmic_variance'],
                                type='untransposedcomplex')
    dlin = dlin.apply(lambda k, v: tf(k) * v).c2r()

    if cat.comm.rank == 0:
        cat.logger.info("linear field generated.")

    Nchunks = max(cat.comm.allgather(cat.size)) // (4 * 1024 * 1024)
    Nchunks = max(1, Nchunks)

    if cat.comm.rank == 0:
        cat.logger.info("Nchunks = %d", Nchunks)

    ID = cat['ID']

    delta = numpy.empty(len(ID), dtype='f4')
    
    if cat.comm.rank == 0:
        cat.logger.info("delta allocated for %d particles.", cat.csize)

    def work(ID, i):
        if cat.comm.rank == 0:
            cat.logger.info("%d / %d", i, Nchunks)
        Q = id2q(ID, strides=strides, scale=scale, shift=shift)
        csize = cat.comm.allreduce(ID.size)
        if csize == 0:
            return numpy.zeros(len(Q), dtype='f4')
        cmin = numpy.min(cat.comm.allgather(Q.min(axis=0)), axis=0)
        cmax = numpy.max(cat.comm.allgather(Q.max(axis=0)), axis=0)
        if cat.comm.rank == 0:
            cat.logger.info("Q = [ %s ] - [ %s ], len(Q) = %d", cmin, cmax, csize)

        layout = pm.decompose(Q)
        if cat.comm.rank == 0:
            cat.logger.info("decompose finished.")
        delta = dlin.readout(Q, layout=layout)
        if cat.comm.rank == 0:
            cat.logger.info("readout done.")
        csum1 = cat.comm.allreduce(delta.sum())
        csum2 = cat.comm.allreduce((delta ** 2).sum())
        cmean = csum1 / csize
        cstd = (csum2 / csize - (csum1 / csize)**2) ** 0.5
        if cat.comm.rank == 0:
            cat.logger.info("On rank0, <delta> = %g, std(delta) = %g", cmean, cstd)
        return delta

    for i in range(Nchunks):
        chunk = slice(i * len(ID) // Nchunks, (i + 1) * len(ID) // Nchunks)
        delta[chunk] = work(ID[chunk].compute(), i)

    cat[ns.ocolumn] = delta

    if ns.ocatalog is None:
        ns.ocatalog = ns.catalog

    cat.save(ns.ocatalog, columns=[ns.ocolumn], header=None, dataset=ns.dataset)
    cat.comm.barrier()
    if cat.comm.rank == 0:
        cat.logger.info("done")

if __name__ == "__main__":
    main()
