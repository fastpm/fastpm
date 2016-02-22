FastPM
======

CI-Status

.. image:: https://api.travis-ci.org/rainwoodman/fastpm.svg
    :alt: Build Status
    :target: https://travis-ci.org/rainwoodman/fastpm/

Introduction
------------

FastPM solves the gravity Possion equation with a boosted particle mesh. Arbitrary
time steps can be used.  
The code is indented to study the formation of large scale structure.

Serious attempts are made to ensure the source code of FastPM is relatively 
clear to read (mirroring the complicity of the algorithms).
FastPM provides a library and a C-API (still unstable) to be reused at binary level
in another application.

FastPM supports plain PM and Comoving-Lagranian (COLA) solvers. 
A broadband correction enforces the linear theory model growth
factor at large scale. See the code paper [FIXME].

Thanks to the PFFT Fourier Transform library, FastPM scales extremely well, 
to hundred thousand MPI ranks. 

The size of mesh in FastPM can vary with time, allowing one to use coarse force mesh at high redshift
with increase temporal resolution for accurate large scale modes.

The finite differentiation kernel in FastPM is the 4 point low-noise super-lanzcos kernel. 
A discrete laplacian operator is used to solve Poisson's equation. 

A parameter file interpreter is provided to validate and execute the configuration 
files without running the simulation, allowing creative usages of the configuration files.


IO and Compatibility
--------------------

The snapshot format of FastPM is [bigfile]_. The format can be easily accessed by python, C, or Fortran.
Alternatively, the snapshot can be written as TPM by Martin White, which can be subtly accessed by 
Python, C, or Fortran.

Due to unjustified complicity, snapshots in legacy GADGET format is not supported by FastPM. 

The nbody post-analysis package [NBODYKIT]_ natively supports bigfile and TPM formats.

[NBODYKIT]_ provides tools for calculating two point functions, 
generating QPM mocks, and identifying Friends-of-Friends (DBSCAN)
haloes and calculating spherical overdensity properties of subhalos.

In addition to the snapshots, FastPM calculates and writes the power-spectrum at each time step. 
These power spectrum files are compatible with numpy plain text files. 

Acknowledgement
---------------

FastPM uses or references publicly avaiable codes ([PFFT]_
[2LPT]_, [COLAHALO]_, [LUA]_, [NBODYKIT]_, [MP-GADGET]_)
and private codes (mpipm and ic_2lpt by Martin White). 

The Particle Mesh solver and 2LPT initial condition generator in FastPM are written from scratch
to properly support pencil / stencil domain decomposition schemes.

The following files distributed in FastPM are originally from other projects:

- fastpm-steps.c : COLA stepping [COLAHALO]_ [COLA]_ 

- power.c : Reading CAMB power spectrum [COLAHALO]_ [2LPT]_ [MP-Gadget]_ N-GenIC

- cosmology.c : Cosmology functions [COLAHALO]_ [2LPT]_ [MP-Gadget]_ N-GenIC

- lua : Lua library [LUA]_

The source code of FastPM is distributed under GPLv3, with the exception files in
lua directory distributed under any appropriate license of lua. 

.. [PFFT] http://github.com/mpip/pfft
.. [2LPT] http://cosmo.nyu.edu/roman/2LPT/
.. [COLA] https://bitbucket.org/tassev/colacode/
.. [NBODYKIT] http://github.com/bccp/nbodykit/
.. [COLAHALO] http://github.com/junkoda/cola_halo
.. [LUA] http://lua.org/
.. [MP-GADGET] http://bluetides-project.org/code
.. [bigfile] https://github.com/rainwoodman/bigfile

Installation
------------

Set up the compilers and location of files in Makefile.local. An example
is provided in Makefile.local.example which shall work on a recent version of
Fedora .

gsl : Most super-computing facility have these already installed. Locate the
path.  Point GSL_DIR to the installation dir. (parent directory of lib and include)

pfft : bundled and built statically via 

.. code::

    make -f Makefile.pfft

Some minor tweaks to Makefile.pfft on the configure scripts may be needed.
Especially the --enable-avx and --enable-sse / --enable-sse2 flags 
if compliation fails with strange errors.

The automatical dependency requires a working version of gcc.

For example, on Edison

.. code::

    # this will set GSL_DIR automatically
    module load gsl
    cp Makefile.edison Makefile.local
    make

Examples
--------

Refer to tests/example.lua and tests/runtest.sh

CI
--

Lastest power spectrum from TravisCI: 

.. image:: https://rainwoodman.github.io/fastpm/tests-result.png

