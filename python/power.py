from nbodykit.lab import FFTPower, BigFileCatalog
import numpy
import argparse
import warnings

# usage:
#
# python power.py output [...] catalog1 [...] -- catalog2 [...]
#
# for cross correlation the second catalog is given after '--'

ap = argparse.ArgumentParser()
ap.add_argument("output", help='e.g. power.json (FFTPower.load) or power.txt (numpy.loadtxt)')
ap.add_argument("--mode", choices=['1d', '2d'], default=None)
ap.add_argument("--with-plot", action='store_true', default=False, help='make a plot (will be same basename as output')
ap.add_argument("--unique-k", action='store_true', default=False, help='compute for all unique k values.')
ap.add_argument("--nmesh", type=int, default=256, help='mesh resolution')

cat_ap = argparse.ArgumentParser()

cat_ap.add_argument("catalog", help='e.g. fastpm_1.0000 or fof_1.0000')
cat_ap.add_argument("--dataset", default='LL-0.200', help='data set to select; for a dm catalog, use 1 for a halo catalog, usually LL-0.200')
cat_ap.add_argument("--with-rsd", action='store_true', default=False, help='Apply RSD')
cat_ap.add_argument("--nmax", type=int, default=None, help='max number of particles, inclusive')
cat_ap.add_argument("--nmin", type=int, default=None, help='min number of particles, inclusive')
cat_ap.add_argument("--abundance", type=float, default=None, help='cut with abundance per (Mpc/h)^3), overrides nmax or nmin')

ns, args = ap.parse_known_args()

if '--' in args:
    split = args.index('--')
    ns1 = cat_ap.parse_args(args[:split])
    ns2 = cat_ap.parse_args(args[split+1:])
else:
    ns1 = cat_ap.parse_args(args)
    ns2 = ns1

def read_cat(ns):
    cat = BigFileCatalog(ns.catalog, header='Header', dataset=ns.dataset)
    volume = cat.attrs['BoxSize'][0] ** 3
    sel = True
    if ns.abundance is not None:
        sel = sel | (cat['Index'] < volume * ns.abundance)
    if ns.nmin is not None:
        sel = sel | (cat['Length'] >= ns.nmin)
    if ns.nmax is not None:
        sel = sel | (cat['Length'] <= ns.nmax)

    cat['VelocityOffset'] = cat['Velocity'] * cat.attrs['RSDFactor']

    cat['Selection'] = sel
    if ns.with_rsd:
        cat['Position'] += cat['VelocityOffset'] * [0, 0, 1.]
    return cat

def main(ns, ns1, ns2):
    cat1 = read_cat(ns1)
    cat2 = read_cat(ns2)

    mesh1 = cat1.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)
    mesh2 = cat2.to_mesh(interlaced=True, compensated=True, window='tsc', Nmesh=ns.nmesh)

    if ns1.with_rsd != ns2.with_rsd:
        warnings.warn("Two catalogs have different with-rsd settings, this may not be intended.")

    if ns.mode is None:
        if ns1.with_rsd or ns2.with_rsd:
            ns.mode = '2d'
        else:
            ns.mode = '1d'

    if ns.unique_k:
        dk = 0
    else:
        dk = None

    r = FFTPower(mesh1, second=mesh2, mode=ns.mode, dk=dk)

    if ns.output.endswith('.json'):
        basename = ns.output[:-5]
        r.save(ns.output)
    elif ns.output.endswith('.txt'):
        basename = ns.output[:-4]

        for var in r.power.data.dtype.names:
            numpy.savetxt(basename + '-%s.txt' % var,
                r.power[var]
            )

    if ns.with_plot:
        figure = make_plot(r)
        figure.savefig(basename + '.png')

def make_plot(r):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    figure = Figure()
    canvas = FigureCanvasAgg(figure)
    ax = figure.add_subplot(111)

    if 'mu' in r.power:
        lines = ax.plot(r.power['k'], r.power['power'].real)
        labels = ['power-%.2g' % mu for mu in r.power.coords['mu']]
    else:
        lines = ax.plot(r.power['k'], r.power['power'].real)
        labels = ['power']

    ax.axhline(r.attrs['shotnoise'], label='shotnoise', ls='--')
    ax.legend(lines, labels)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('k h/Mpc')
    ax.set_ylabel('P')

    return figure

main(ns, ns1, ns2)
