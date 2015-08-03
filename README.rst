fastPM
======

.. image:: https://api.travis-ci.org/rainwoodman/fastPM.svg
    :alt: Build Status
    :target: https://travis-ci.org/rainwoodman/fastPM/

A simple Particle Mesh (fastPM) solver for N-body cosmological simulations,
for the formation of large scale structures.

fastPM solves the gravity Possion equation with an augumented particle mesh 
that is 2x to 3x of the particle grid. 

The post analysis chain is nbodykit [3]_, consisting tools for calculating the
power spectrum, and halo identification.

fastPM is extremely efficient, with 80% of wall time spent 
in Fast Fourier Transform.

To solver a 8 billion particle system over the 13.7 billion years
of the Universe (from z=9 to z=0), it takes

- with 20 steps,
- on 1024 cores,
- in 30 minutes.

The infrastrcture code in cola_halo by Dr. Jun Koda have greatly helped
us to speed up the development. fastPM still maintains a legacy 
COLA-compatible mode.

This code is based on publicly avaiable codes ([1]_, [2]_)
and previously privated codes (QPM by Martin White), modified and redistributed 
under GPLv3.


.. [1] http://cosmo.nyu.edu/roman/2LPT/
.. [2] https://bitbucket.org/tassev/colacode/
.. [3] http://github.com/bccp/nbodykit/

