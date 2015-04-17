QrPM
====

Quick and bRutal Particle Mesh (QrPM) 
is a naive and fast gravity solver, for 
formation of large scale structures.

QrPM solves the gravity Possion equation with
an augumented particle mesh that is 2x to 3x of the 
particle grid.

Halos are identified with a Friend-of-Friend algorithm[3].

QrPM is extremely efficient, with 80% of wall time spent 
in Fast Fourier Transform.

To solver a 8 billion particle system over the 13.7 billion years
of the Universe (from z=9 to z=0), it takes

- with 20 steps,
- on 1024 cores,
- in 30 minutes.

The infrastrcture code in cola_halo by Dr. Jun Koda have greatly helped
us to speed up the development. QrPM still maintains a legacy 
COLA-compatible mode.

This code is based on publicly avaiable codes ([1], [2], [3], [4]),
and previously privated codes (QPM by Martin White), modified and redistributed 
under GPLv3.


[1]: http://cosmo.nyu.edu/roman/2LPT/
[2]: https://bitbucket.org/tassev/colacode/
[3]: https://github.com/junkuda/cola_halo/
[4]: http://www-hpcc.astro.washington.edu/tools/fof.html

