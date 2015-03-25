from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy

kcola, pcola = numpy.loadtxt('runcola/powerspec-1.0000.txt', unpack=True) 
k0 = kcola
pqpm = numpy.interp(k0, *numpy.loadtxt('runqpm/powerspec-1.0000.txt', unpack=True))
pnewqpm = numpy.interp(k0, *numpy.loadtxt('runnewqpm/powerspec-1.0000.txt', unpack=True))
pnewbland= numpy.interp(k0, *numpy.loadtxt('runnewbland/powerspec-1.0000.txt', unpack=True))
pnewexhaust= numpy.interp(k0, *numpy.loadtxt('runnewexhaust/powerspec-1.0000.txt', unpack=True))
pgrowthexhaust= numpy.interp(k0, *numpy.loadtxt('rungrowthexhaust/powerspec-1.0000.txt', unpack=True))
p2lpt = numpy.interp(k0, *numpy.loadtxt('run2lpt/powerspec-1.0000.txt', unpack=True))

figure = Figure()
ax = figure.add_subplot(111)
ax.plot(k0, pcola / pcola, label='COLA')
ax.plot(k0, pqpm / pcola, label='QPM')
ax.plot(k0, pnewqpm / pcola, label='NEW QPM')
ax.plot(k0, pnewbland / pcola, label='NEW BLAND')
ax.plot(k0, pnewexhaust / pcola, label='NEW EXHAUST')
ax.plot(k0, pgrowthexhaust / pcola, label='QPM GROWTH')
ax.plot(k0, p2lpt / pcola, label='2LPT')

ax.set_yscale('linear')
ax.set_ylim(0.8, 1.2)
ax.set_xscale('log')
ax.legend(loc='upper left')
canvas = FigureCanvasAgg(figure)
figure.savefig('plotpower.png')

