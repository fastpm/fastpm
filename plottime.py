from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy

def getboth(aa):
    kgroundtruth, pgroundtruth = numpy.loadtxt('rungroundtruth/powerspec-%0.4f.txt' % aa, unpack=True) 
    k0 = kgroundtruth
    pmanystepqpm = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm1/powerspec-%0.4f.txt' % aa, unpack=True))
    return k0, pmanystepqpm / pgroundtruth
figure = Figure()
ax = figure.add_subplot(111)
ax.plot(*getboth(1.0), label='hack QPM 40 1.0')
ax.plot(*getboth(0.3250), label='hack QPM 40 0.3250')
ax.plot(*getboth(0.5050), label='hack QPM 40 0.5050')
ax.plot(*getboth(0.8200), label='hack QPM 40 0.8200')
ax.plot(*getboth(0.7300), label='hack QPM 40 0.7300')

ax.set_yscale('linear')
ax.set_ylim(0.85, 1.15)
ax.set_xlim(0.01, 5)
ax.set_xscale('log')
ax.legend(loc='upper left', fontsize='small')
ax.grid()
canvas = FigureCanvasAgg(figure)
figure.savefig('plotpowerbytime.png')

