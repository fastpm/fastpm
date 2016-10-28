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
factor at large scale. See the code paper [FastPMPaper]_

Thanks to the PFFT Fourier Transform library, FastPM scales extremely well,
to hundred thousand MPI ranks.

The size of mesh in FastPM can vary with time, allowing one to use coarse force mesh at high redshift
with increase temporal resolution for accurate large scale modes.

FastPM supports a huge varieties of Greens function and differentiation kernels, though for most practical
simulations the choice of kernels does not make a difference.

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

IO Units
----------

*FastPM format*

- Position is in units of Mpc/h.
- Velocity is in peculiar Km/s, :math:`v_p = \frac{a}\dot{x}`

To convert velocity to RSD, use
:math:`\delta x_{rsd} = \frac{1}{aH} v_p`
This number is saved in the header.

*RunPB format*

- Position is between 0 and 1.
- Velocity is in RSD units, between 0 and 1.

*Interanal Units*

- Position is in units of Mpc/h.
- Velocity is :math:`v_i = \frac{a^2}{H_0}\dot{x}`

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
.. [FastPMPaper] http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1603.00476].


Installation
------------

FastPM works out-of-the-box on many GNU/Linux and compatible systems.
Notably the development is performed on a Cray system, a Fedora Workstation,
and the integrated continuous integration tests are performed on a Ubuntu based Linux distribution.
The recommended compiler is `gcc`. FastPM is built with the GNU `make` tool.

Set up the compilers and location of files in Makefile.local. An example
is provided in Makefile.local.example which shall work on a recent version of
Fedora .

- gsl : Most super-computing facility have these already installed. Locate the
  path.  Point GSL_DIR to the installation dir. (parent directory of lib and include)

- pfft : bundled and built statically in depends directory  `Makefile.pfft`.
  Some minor tweaks to Makefile.pfft on the configure scripts may be needed.
  Especially the `--enable-avx` and `--enable-sse` / `--enable-sse2` flags 
  if compliation fails with strange errors about invalid instructions.

The automatical dependency requires a working version of gcc, so its the best
to compile with the gnu compilers.

The make process requires a `Makefile.local` file, which sets the variables
like compiler (`MPICC`). A few examples are provided, but you shall customize
it based on the example for your site.

.. code::

    # the following example works at NERSC
    # this will set GSL_DIR automatically

    module load gsl

    # copy the edison example file to Makefile.local

    cp Makefile.edison Makefile.local

    # the rest is just make. It may take a while.
    make

Examples
--------

Refer to tests/standard.lua and tests/runtest.sh

Commandline Interface
---------------------
The CLI consists of two main executable files:

 - `fastpm` is the main executable file of FastPM.
 - `fastpm-lua` is an interpreter that executes the `main` function defined in a parameter file.

A parameter file instructs the run of FastPM. The parameter file is written in the LUA programming language.
We refer the readers to the Lua Reference manual for syntax and run-time libraries of the LUA programming language.
In a parameter file, the command-line arguments to fastpm can be accessed by the `args` variable, allowing dynamic generation of parameters during run-time. 
The interpreter `fastpm-lua` can be used to process the parameter file and generate job script files.
The example parameter file `standard.lua` is distributed with the software in the code repository.

FastPM use the initial condition from a 3-dimensional white-noise, a linear density field `read_lineark`, 
or initial position and velocity of particles `read_runpbic`.

- The white noise field requires a linear theory power spectrum input. The white noise can be retrieved from
a Fourier space dump from FastPM (`read_whitenoisek`), or a configuration space dump from GRAFIC.
The GRAFIC file contains a set of FORTRAN 77 unformatted data blocks, one per each slab in z-y plane. 
The size of the GRAFIC mesh must match with the number of particles in FastPM. 
It is important to be aware that the coordinates in FastPM is transposed from GRAFIC, 
with the transformation :math:`x \to z, y \to y, z \to x`.
(`read_grafic`),
or generated from a random seed (`random_seed`) based on the scale invariant Gadget N-GenIC sequence.

- A linear density field in Fourier space (`read_lineark`). The field shall have the correct linear theory power at z=0.

- Particle position and velocity evolved with 2LPT initial condition generator. (`read_runpbic`).
  The Lagrangian position of the particles are assumed to be on a regular grid,
  and the :math:`s_1`, :math:`s_2` terms are recovered from velocity and
  displacement according to the cosmology specified in the parameter file. This
  type of input is used for the comparison with RunPB TreePM simulations.

An arbitrary list of time steps can be specified in the parameter
file(`time_steps`). We provide functions the create three commonly
used time stepping: 

- `linspace(a_0, a_1, N)`: N + 1 steps linear in scaling factor :math:`a \in [a_0, a_1]`.
- `logspace(log a_0, log a_1, N)`: N + 1 steps linear in :math:`\log a \in [\lg a_0, \lg a_1]`.

The names are inspired from similar functions to
generate sequences in numpy, but be aware of the subtle differences.
Functions here always includes an additional `end` point, while those in numpy do not.

FastPM measures and stores the dark matter power spectrum at each Kick step to
a path specified in the parameter file(`write_powerspectrum`). The
measurement is performed on the density field that produces the gravitational
force; no correction for aliasing or shot noise is applied.

At selected redshifts (`output_redshift`), FastPM writes snapshot in [bigfile]_ format to a path (`write_snapshot`). 
The bigfile format stores data in a sequence of plain binary files and meta data in plain text files. 

The snapshots can be read by nbodykit via the `FastPM` data source plugin,
or directly accessed as binary files. Position is stored in the unit of Mpc/h, 
and velocity is stored as peculiar velocity in km/s. The attributes of the header block stores the meta data about the simulation.

C Application Programming Interface
-----------------------------------

The FastPM CLI is built on top of `libfastpm`. The core functionality of
`libfastpm` is to evolve a linear theory over-density field to a non-linear
density field and a list of particle displacement and velocities. There are
also tools for measurement of power spectrum and generating Gaussian
realizations of initial linear density field.

The library is built as `libfastpm/libfastpm.a`. To use the library,
include `fastpm/libfastpm.h` from the `api` directory. 
Two solver classes are provided,

- `FastPM` : for multi-step particle mesh simulations)
- `FastPM2LPT` : for 1/2LPT particle mesh simulations).

We refer interested users to `src/test2lpt.c` and `src/testpm.c` for example uses of the C-API.
We make the best effort to ensure the API is compatible with C++. If not, please report an issue.

