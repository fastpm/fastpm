from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy

kgroundtruth, pgroundtruth = numpy.loadtxt('rungroundtruth/powerspec-1.0000.txt', unpack=True) 
k0 = kgroundtruth
pqpm = numpy.interp(k0, *numpy.loadtxt('runqpm/powerspec-1.0000.txt', unpack=True))
pcola = numpy.interp(k0, *numpy.loadtxt('runcola/powerspec-1.0000.txt', unpack=True))
pmanystepqpm = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm/powerspec-1.0000.txt', unpack=True))
pmanystepqpm1 = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm1/powerspec-1.0000.txt', unpack=True))
pmanystepqpm2 = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm2/powerspec-1.0000.txt', unpack=True))
pmanystepqpm3 = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm3/powerspec-1.0000.txt', unpack=True))

figure = Figure()
ax = figure.add_subplot(111)
ax.plot(k0, pcola / pgroundtruth, label='COLA')
ax.plot(k0, pqpm / pgroundtruth, label='QPM')
ax.plot(k0, pmanystepqpm / pgroundtruth, label='hack QPM 40 steps 0')
ax.plot(k0, pmanystepqpm1 / pgroundtruth, label='hack QPM 40 steps1')
ax.plot(k0, pmanystepqpm2 / pgroundtruth, label='hack QPM 40 steps2')
ax.plot(k0, pmanystepqpm3 / pgroundtruth, label='hack QPM 40 steps3')

ax.set_yscale('linear')
ax.set_ylim(0.85, 1.15)
ax.set_xlim(0.01, 5)
ax.set_xscale('log')
ax.legend(loc='upper left', fontsize='small')
ax.grid()
canvas = FigureCanvasAgg(figure)
figure.savefig('plotpowerbyorder.png')

