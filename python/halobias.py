from nbodykit.lab import FFTPower, BigFileCatalog
from nbodykit import setup_logging
import numpy
import argparse
import warnings

# usage:
#
# python halobias.py output --nmin x --nmax x --nn x dmcatalog [...] -- halocatalog [...]
#
# measure the bias by cross correlating halos of different size to the dark matter.
#
# example
#
# python halobias.py test.json --with-plot ../tests/nbodykit/fastpm_1.0000/ --dataset=1 -- ../tests/nbodykit/fof_1.0000/ --dataset=LL-0.200

ap = argparse.ArgumentParser()
ap.add_argument("output", help='e.g. power.json (FFTPower.load) or power.txt (numpy.loadtxt)')
ap.add_argument("--nmin", default=8, type=int)
ap.add_argument("--nmax", default=1000, type=int)
ap.add_argument("--nn", default=10, type=int)
ap.add_argument("--with-plot", action='store_true', default=False, help='make a plot (will be same basename as output')
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

    if nmin is not None:
        sel = True
        sel = sel & (cat['Length'] >= nmin)

        cat['Selection'] = sel

    return cat

def main(ns, ns1, ns2):
    if ns.verbose:
        setup_logging('info')

    cat1 = read_cat(ns1)
    mesh1 = cat1.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)
    cat2 = read_cat(ns2)

    if ns.unique_k:
        dk = 0
    else:
        dk = None

    rm = FFTPower(mesh1, second=mesh1, mode='1d', dk=dk)
    nmin = numpy.unique(numpy.int32(numpy.logspace(numpy.log10(ns.nmin), numpy.log10(ns.nmax), ns.nn, endpoint=True)))
    nmin0 = cat2['Length'].min().compute()
    nmax0 = cat2['Length'].max().compute()
    nmin = nmin[nmin >= nmin0]
    nmin = nmin[nmin < nmax0]

    save_bs(ns.output, 'a-matter', rm)

    r = []
    b = []
    for nmin1 in nmin:
        cat2 = read_cat(ns2, nmin1)
        mesh2 = cat2.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)

        r.append(FFTPower(mesh1, second=mesh2, mode='1d', dk=dk))

        save_bs(ns.output, 'x-nmin-%05d' % nmin1, r[-1])
        bias = fit_bias(r[-1], rm)
        b.append(bias)
        print(nmin1, bias)

    basename = ns.output.rsplit('.', 1)[0]
    numpy.savetxt(basename + '-bias.txt', numpy.array([nmin, b]).T)

    if ns.with_plot:
        figure = make_plot(rm, r, nmin)
        figure.savefig(basename + '.png')

def fit_bias(r, rm):
    rat = r.power['power'][:20].real / rm.power['power'][:20].real
    good = ~numpy.isnan(rat)
    return rat[good].mean()

def save_bs(filename, dataset, r):

    basename = filename.rsplit('.', 1)[0]
    if filename.endswith('.json'):
        r.save(basename + '-%s.json' % dataset)
    elif filename.endswith('.txt'):
        for var in r.power.data.dtype.names:
            numpy.savetxt(basename + '-%s-%s.txt' % (dataset, var),
                r.power[var]
            )

def make_plot(rm, r, nmin):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    figure = Figure()
    canvas = FigureCanvasAgg(figure)
    ax = figure.add_subplot(111)

    for r1, nmin1 in zip(r, nmin):
        lines = ax.plot(rm.power['k'], r1.power['power'].real / rm.power['power'].real, label='nmin = %d' % nmin1)

    ax.legend()

    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel('k h/Mpc')
    ax.set_ylabel('Px / Pa')

    return figure

main(ns, ns1, ns2)
