fastPM
======

.. image:: https://api.travis-ci.org/rainwoodman/fastPM.svg
    :alt: Build Status
    :target: https://travis-ci.org/rainwoodman/fastPM/

A simple Particle Mesh (fastPM) solver for N-body cosmological simulations.
The code is indented to study the formation of large scale structure.

fastPM solves the gravity Possion equation with an augumented particle mesh 
that is 2x to 3x of the particle grid. 

2 parallel fourier transform backends are
supported, slab-based FFTW and pencil-based PFFT. 
The PFFT backend enables fastPM to scale to hundred thousand MPI ranks, allows
runs with significantly more than 2048**3 particles. 
The FFTW backend allows one to run smaller simulations (<= 2048**3 particles on 1024 ranks
more efficiently.

fastPM supports plain PM (correctly implemented) and Comoving-Lagranian (COLA)
solvers. Several optimizations in fastPM allows one to use arbitrary differential kernel 
in fastPM at little additional cost.

The post analysis chain is nbodykit [3]_, consisting tools for calculating 
power spectrum and identifying haloes and subhalos.

fastPM is extremely efficient, with 80% of wall time spent 
in Fast Fourier Transform.

To solver a 8 billion particle system over the 13.7 billion years
of the Universe (from z=9 to z=0), it takes

- 20 steps,
- on 1024 cores,
- 30 minutes.

The infrastrcture code in cola_halo [4]_ by Dr. Jun Koda have greatly helped
us to speed up the development. 

This code is based on publicly avaiable codes ([1]_, [2]_)
and previously privated codes (QPM by Martin White), modified and redistributed 
under GPLv3.


.. [1] http://cosmo.nyu.edu/roman/2LPT/
.. [2] https://bitbucket.org/tassev/colacode/
.. [3] http://github.com/bccp/nbodykit/
.. [4] http://github.com/junkoda/cola_halo

Installation
------------

First gsl. Most super-computing have these already installed.

FFTW and PFFT are bundled and installed 

.. code::

    depends/install_pfft.sh $PWD/depends/install

Some minor tweaks to install_pfft.sh on the configure scripts may be needed.
Especially the --enable-avx and --enable-sse / --enable-sse2 flags if compliation
fails with strange errors.

The automatical dependency requires a working version of gcc.

Point GSL_DIR to the installation dir and use

.. code::

    make GSL_DIR= .... MPICC=....


Examples
--------

Refer to tests/example.lua and tests/runtest.sh



