from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy
import glob
import os
TRAVIS_BUILD_NUMBER = os.environ.get('TRAVIS_BUILD_NUMBER', 'local')
TRAVIS_COMMIT = os.environ.get('TRAVIS_COMMIT', 'local')

def plotps(filename, ax):
    klin, plin = numpy.loadtxt('powerspec.txt', unpack=True)
    k, p, j = numpy.loadtxt(filename, unpack=True)
    p[0] = numpy.nan
    plin = numpy.interp(k, klin, plin)
    ax.plot(k, p / plin, label=filename)

fig = Figure(figsize=(12, 6))
ax = fig.add_subplot(121)
for fn in sorted(glob.glob('cola/powerspec_*.txt')):
    plotps(fn, ax)
ax.set_ylim(0, 2)
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylabel('P/P_lin')
ax.set_title('cola: ' + TRAVIS_BUILD_NUMBER + '\n' + TRAVIS_COMMIT)

ax = fig.add_subplot(122)
for fn in sorted(glob.glob('pm/powerspec_*.txt')):
    plotps(fn, ax)

ax.set_ylim(0, 2)
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylabel('P/P_lin')
ax.set_title('pm: ' + TRAVIS_BUILD_NUMBER + '\n' + TRAVIS_COMMIT)

canvas = FigureCanvasAgg(fig)
fig.savefig('tests-result.png', dpi=72)
