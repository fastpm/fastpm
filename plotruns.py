from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy

kgroundtruth, pgroundtruth = numpy.loadtxt('rungroundtruth/powerspec-1.0000.txt', unpack=True) 
k0 = kgroundtruth
pqpm = numpy.interp(k0, *numpy.loadtxt('runqpm/powerspec-1.0000.txt', unpack=True))
pcola = numpy.interp(k0, *numpy.loadtxt('runcola/powerspec-1.0000.txt', unpack=True))
pmanystepqpm = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm/powerspec-1.0000.txt', unpack=True))
pmanystepqpm2 = numpy.interp(k0, *numpy.loadtxt('runmanystepqpm2/powerspec-1.0000.txt', unpack=True))
pnewqpm = numpy.interp(k0, *numpy.loadtxt('runnewqpm/powerspec-1.0000.txt', unpack=True))
pnewbland= numpy.interp(k0, *numpy.loadtxt('runnewbland/powerspec-1.0000.txt', unpack=True))
pnewexhaust= numpy.interp(k0, *numpy.loadtxt('runnewexhaust/powerspec-1.0000.txt', unpack=True))
pnewluxury= numpy.interp(k0, *numpy.loadtxt('runnewluxury/powerspec-1.0000.txt', unpack=True))
pgrowthexhaust= numpy.interp(k0, *numpy.loadtxt('rungrowthexhaust/powerspec-1.0000.txt', unpack=True))
p2lpt = numpy.interp(k0, *numpy.loadtxt('run2lpt/powerspec-1.0000.txt', unpack=True))

figure = Figure()
ax = figure.add_subplot(111)
ax.plot(k0, pcola / pgroundtruth, label='COLA')
ax.plot(k0, pqpm / pgroundtruth, label='QPM')
ax.plot(k0, pmanystepqpm / pgroundtruth, label='hack QPM 40 steps')
ax.plot(k0, pmanystepqpm2 / pgroundtruth, label='hack QPM 40 steps highorder')
ax.plot(k0, pnewqpm / pgroundtruth, label='NEW QPM')
ax.plot(k0, pnewbland / pgroundtruth, label='NEW BLAND')
ax.plot(k0, pnewexhaust / pgroundtruth, label='NEW EXHAUST')
#ax.plot(k0, pnewluxury/ pgroundtruth, label='NEW LUXURY')
ax.plot(k0, pgrowthexhaust / pgroundtruth, label='QPM GROWTH')
ax.plot(k0, p2lpt / pgroundtruth, label='2LPT')

ax.set_yscale('linear')
ax.set_ylim(0.85, 1.15)
ax.set_xlim(0.01, 5)
ax.set_xscale('log')
ax.legend(loc='upper left', fontsize='small')
ax.grid()
canvas = FigureCanvasAgg(figure)
figure.savefig('plotpower.png')

