fastPM
======

.. image:: https://api.travis-ci.org/rainwoodman/fastPM.svg
    :alt: Build Status
    :target: https://travis-ci.org/rainwoodman/fastPM/

A simple Particle Mesh (fastPM) solver for N-body cosmological simulations.
The code is indented to study the formation of large scale structure.

fastPM solves the gravity Possion equation with an variable augumented particle mesh.

2 parallel fourier transform backends are supported, slab-based FFTW and pencil-based PFFT. 

- The PFFT backend enables fastPM to scale to hundred thousand MPI ranks, allows
runs with significantly more than 2048**3 particles. 

- The FFTW backend allows one to run smaller simulations (<= 2048**3 particles on 1024 ranks
more efficiently.

fastPM supports plain PM and Comoving-Lagranian (COLA) solvers. The finite differentiation kernel
in fastPM is the 4 point low-noise super-lanzcos kernel. A discrete laplacian operator is used to solve
Poisson's equation.

The IO format of fastPM is identical to TPM by Martin White.  
The post analysis chain nbodykit [3]_ natively support these formats. nbodykit [3]_ provides
 tools for calculating 2 point functions, making QPM mocks, and identifying Friend-of-Friend 
haloes and calculating spherical overdensity properties of subhalos.

The Particle Mesh solver and 2LPT initial condition generator are written from scratch.

The following files in fastPM are originally from cola_halo [4] by Jun Koda:

- msg.c :  Provides logging infrastructure

- pmtimer.c : Provides time profiling

- pmsteps.c : COLA stepping

- power.c : Cosmology functions

fastPM is based on publicly avaiable codes ([1]_, [2]_)
and previously privated codes (QPM and ic_2lpt by Martin White), modified and redistributed 
under GPLv3.

.. [1] http://cosmo.nyu.edu/roman/2LPT/
.. [2] https://bitbucket.org/tassev/colacode/
.. [3] http://github.com/bccp/nbodykit/
.. [4] http://github.com/junkoda/cola_halo

Installation
------------

First gsl. Most super-computing have these already installed.

FFTW and PFFT are bundled and installed .

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

