fastPM
======

.. image:: https://api.travis-ci.org/rainwoodman/fastPM.svg
    :alt: Build Status
    :target: https://travis-ci.org/rainwoodman/fastPM/

Introduction
------------

fastPM is a simple Particle Mesh solver for N-body cosmological simulations.
The code is indented to study the formation of large scale structure.

fastPM solves the gravity Possion equation with an variable augumented particle mesh.

2 parallel fourier transform backends are supported, slab-based FFTW and pencil-based PFFT. 

- The PFFT backend enables fastPM to scale to hundred thousand MPI ranks, allows
  runs with significantly more than 2048**3 particles. 

- The FFTW backend allows one to run smaller simulations 
   (<= 2048**3 particles on 1024 ranks) more efficiently.


The size of mesh in fastPM can vary with time, allowing one to use coarse force mesh at high redshift
with increase temporal resolution for accurate large scale modes.

fastPM supports plain PM and Comoving-Lagranian (COLA) solvers. The finite differentiation kernel
in fastPM is the 4 point low-noise super-lanzcos kernel. A discrete laplacian operator is used to solve
Poisson's equation. These kernels, which suppresses artifical noise at small scales, should have been 
used in all PM codes.

IO and Compatibility
--------------------

The IO format of fastPM is identical that of TPM by Martin White.
The post analysis chain nbodykit [NBODYKIT]_ natively support this format. 
nbodykit [NBODYKIT]_ provides tools for calculating two point functions, generating QPM mocks, 
and identifying Friends-of-Friends (DBSCAN)
haloes and calculating spherical overdensity properties of subhalos.

In addition to the snapshots, fastPM calculates and writes the power-spectrum at each time step. These
power spectrum files are compatible with numpy plain text files, yet contains the same meta data present
in power spectrum files written by nbodykit.

Acknowledgement
---------------

fastPM is based on publicly avaiable codes (
[2LPT]_, [COLA]_, [LUA]_, [NBODYKIT]_, [MP-GADGET]_)
and previously privated codes (QPM and ic_2lpt by Martin White). 

The Particle Mesh solver and 2LPT initial condition generator in fastPM are written from scratch
to properly support pencil / stencil domain decomposition schemes.

The following files in fastPM are originally from other projects:

- fastpm-steps.c : COLA stepping [COLAHALO]_ by Jun Koda.

- power.c : Reading CAMB power spectrum [COLAHALO]_

- cosmology.c : Cosmology functions [COLAHALO]_ 

The source code of fastPM is distributed under GPLv3.

.. [2LPT] http://cosmo.nyu.edu/roman/2LPT/
.. [COLA] https://bitbucket.org/tassev/colacode/
.. [NBODYKIT] http://github.com/bccp/nbodykit/
.. [COLAHALO] http://github.com/junkoda/cola_halo
.. [LUA] http://lua.org/
.. [MP-GADGET] http://bluetides-project.org/code

Installation
------------

Set up the compilers and location of files in Makefile.local. An example
is provided in Makefile.local.example which shall work on a recent version of
Fedora .

gsl : Most super-computing facility have these already installed. Locate the
path.  Point GSL_DIR to the installation dir.

pfft : bundled and installed automatically via 

.. code::

    make -f Makefile.pfft

Some minor tweaks to Makefile.pfft on the configure scripts may be needed.
Especially the --enable-avx and --enable-sse / --enable-sse2 flags 
if compliation fails with strange errors.

The automatical dependency requires a working version of gcc.

.. code::

    make GSL_DIR= .... CC=....


For example, on Edison

.. code::

    # this will set GSL_DIR automatically
    module load gsl
    make CC=cc

Examples
--------

Refer to tests/example.lua and tests/runtest.sh

