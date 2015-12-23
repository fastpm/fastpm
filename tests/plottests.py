from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy
import glob

def plotps(filename, ax):
    klin, plin = numpy.loadtxt('powerspec.txt', unpack=True)
    k, p, j = numpy.loadtxt(filename, unpack=True)
    plin = numpy.interp(k, klin, plin)
    ax.plot(k, p / plin, label=filename)

fig = Figure()
ax = fig.add_subplot(121)
for fn in sorted(glob.glob('cola/powerspec*.txt')):
    plotps(fn, ax)
ax.set_ylim(1, 1e5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('cola')

ax = fig.add_subplot(122)
for fn in sorted(glob.glob('pm/powerspec*.txt')):
    plotps(fn, ax)

ax.set_ylim(1, 1e5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('pm')

canvas = FigureCanvasAgg(fig)
fig.savefig('tests-result.svgz')
