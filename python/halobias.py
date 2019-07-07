from nbodykit.lab import FFTPower, BigFileCatalog, BigFileMesh
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
ap.add_argument("--kmax", default=0.04, type=float, help="cut to stop using kmax, scale where kaiser is bad")
ap.add_argument("--nmax", default=1000, type=int)
ap.add_argument("--nn", default=10, type=int)
ap.add_argument("--with-plot", action='store_true', default=False, help='make a plot (will be same basename as output')
ap.add_argument("--unique-k", action='store_true', default=False, help='compute for all unique k values.')
ap.add_argument("--nmesh", type=int, default=256, help='mesh resolution')
ap.add_argument("--verbose", action='store_true', default=False, help='print progress')

cat_ap = argparse.ArgumentParser()

cat_ap.add_argument("catalog", help='e.g. fastpm_1.0000 or fof_1.0000')
cat_ap.add_argument("--dataset", default='LL-0.200', help='data set to select; for a dm catalog, use 1 for a halo catalog, usually LL-0.200')
cat_ap.add_argument("--mesh", action='store_true', default=False, help='data set is a BigFileMesh')

ns, args = ap.parse_known_args()

if '--' in args:
    split = args.index('--')
    ns1 = cat_ap.parse_args(args[:split])
    ns2 = cat_ap.parse_args(args[split+1:])
else:
    ns1 = cat_ap.parse_args(args)
    ns2 = ns1

def read_cat1(ns, ns1, nmin=None):
    if ns1.mesh:
        mesh = BigFileMesh(ns1.catalog, dataset=ns1.dataset)
    else:
        cat = BigFileCatalog(ns1.catalog, header='Header', dataset=ns1.dataset)
        volume = cat.attrs['BoxSize'][0] ** 3

        if nmin is not None and nmin != 0:
            sel = True
            sel = sel & (cat['Length'] >= nmin)

            cat['Selection'] = sel
        cat['RSDPosition'] = cat['Position'] + cat.attrs['RSDFactor'] * cat['Velocity'] * [0, 0, 1]

        mesh = cat.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)

    return mesh

def read_cat2(ns, ns1, nmin=None):
    cat = BigFileCatalog(ns1.catalog, header='Header', dataset=ns1.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3

    if nmin is not None and nmin != 0:
        sel = True
        sel = sel & (cat['Length'] >= nmin)

        cat['Selection'] = sel
    cat['RSDPosition'] = cat['Position'] + cat.attrs['RSDFactor'] * cat['Velocity'] * [0, 0, 1]
    return cat

def main(ns, ns1, ns2):
    if ns.verbose:
        setup_logging('info')

    cat1 = read_cat1(ns, ns1)

    mesh1 = cat1.paint(mode='complex')

    cat2 = read_cat2(ns, ns2)

    if ns.unique_k:
        dk = 0
    else:
        dk = None

    rm = FFTPower(mesh1, second=mesh1, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax * 10)
    nmin = numpy.unique(numpy.int32(numpy.logspace(numpy.log10(ns.nmin), numpy.log10(ns.nmax), ns.nn, endpoint=True)))
    if 'Length' in cat2.columns:
        nmin0 = cat1.comm.allreduce(cat2['Length'].min().compute() if cat2.size > 0 else 10000000, MPI.MIN)
        nmax0 = cat1.comm.allreduce(cat2['Length'].max().compute() if cat2.size > 0 else 0, MPI.MAX)
        nmin = nmin[nmin >= nmin0]
        nmin = nmin[nmin < nmax0]
    else:
        nmin = [0]

    Nmodes = (rm.power['modes'] * (rm.power['k'] < ns.kmax)).sum()

    if cat1.comm.rank == 0:
        print('Using %d modes to estimate bias and growth rate' % Nmodes)

        dirname = os.path.dirname(ns.output)
        if len(dirname) > 0:
            os.makedirs(dirname, exist_ok=True)

    save_bs(ns.output, 'a-matter', rm)

    r = []
    b = []
    a = []
    f = []
    if cat1.comm.rank == 0:
        print('# Nmin bias growthrate abundance')
    for nmin1 in nmin:
        cat2 = read_cat2(ns, ns2, nmin1)
        mesh2 = cat2.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh, position='RSDPosition')
        mesh3 = cat2.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh, position='Position')

        r_rsd = FFTPower(mesh1, second=mesh2, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax * 10)
        r_real = FFTPower(mesh1, second=mesh3, mode='2d', dk=dk, Nmu=10, kmax=ns.kmax * 10)

        r.append(r_rsd)

        save_bs(ns.output, 'x-nmin-%05d' % nmin1, r[-1])
        bias, gr = fit_bias(r_rsd, r_real, rm, ns.kmax)
        abundance = r[-1].attrs['N2'] / cat2.attrs['BoxSize'][0] ** 3
        b.append(bias)
        a.append(abundance)
        f.append(gr)
        if cat1.comm.rank == 0:
            print(nmin1, bias, gr, abundance)

    basename = ns.output.rsplit('.', 1)[0]

    if cat1.comm.rank == 0:
        numpy.savetxt(basename + '-bias.txt', numpy.array([nmin, b, f, a]).T)

    if ns.with_plot:
        if cat1.comm.rank == 0:
            figure = make_plot(rm, r, nmin, b, f, ns.kmax)
            figure.savefig(basename + '.png')

def fit_bias(r_rsd, r_real, rm, kmax):
    # fit bias with real cross power
    #
    # then fit RSD with
    #
    # d_h^rsd = d_h^real + f mu^2 d_m
    #
    # this enhances variance cancallation, at the non-linear order.

    mu = rm.power['mu']
    mask = rm.power['k'] < kmax

    def loss_b(b):
        with numpy.errstate(all='ignore'):
            model = b * rm.power['power']
            res = (r_real.power['power'].real - model.real)

        res[numpy.isnan(res)] = 0
        res *= rm.power['modes']
        res *= mask

        return numpy.sum(res ** 2)

    def loss_f(f):
        with numpy.errstate(all='ignore'):
            model = r_real.power['power'] + (f * mu ** 2) * rm.power['power']
            res = (r_rsd.power['power'].real - model.real)

        res[numpy.isnan(res)] = 0
        res *= rm.power['modes']
        res *= mask

        return numpy.sum(res ** 2)

    from scipy.optimize import minimize
    res = minimize(lambda x: loss_b(x[0]), x0=(1,), method='Nelder-Mead')

    b = res.x[0]

    res = minimize(lambda x: loss_f(x[0]), x0=(0.), method='Nelder-Mead')
    f = res.x[0]
    return b, f

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

def make_plot(rm, r, nmin, b, f, kmax):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    figure = Figure()
    canvas = FigureCanvasAgg(figure)
    ax = figure.add_subplot(111)

    for r1, nmin1, bias, gr in zip(r, nmin, b, f):
        lines = None
        for i in range(rm.power['k'].shape[1]):
            mask = rm.power['modes'][:, i] > 0
            if i == 0:
                label='nmin = %d' % nmin1
            else:
                label= None
            if lines is None:
                color = None
            else:
                color = lines.get_color()
            lines, = ax.plot(rm.power['k'][:, i][mask],
                    (r1.power['power'].real / rm.power['power'].real)[:, i][mask], '.',
                        color=color, alpha= 0.3 - 0.7 * i / rm.power['k'].shape[1])
            lines, = ax.plot(rm.power['k'][:, i][mask],
                    (bias + gr * r1.power['mu'] ** 2)[:, i][mask],
                    label=label, color=lines.get_color())

    ax.axvline(kmax, ls='--', color='k')
    ax.legend()
    ax.set_ylim(0, numpy.max(b) * 2)
    ax.set_xlim(1e-3, kmax * 8)
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel('k h/Mpc')
    ax.set_ylabel('Px / Pa')

    return figure

main(ns, ns1, ns2)
