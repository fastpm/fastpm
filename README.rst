FastPM
======

CI-Status

.. image:: https://github.com/fastpm/fastpm/workflows/main/badge.svg
    :alt: Build Status
    :target: https://github.com/fastpm/fastpm/actions?query=workflow%3Amain

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
factor at large scale. See the code paper [FastPMPaper]_. 
For details on the neutrino implementation see [FastPMNeutrinoPaper]_.

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

The snapshot container format of FastPM is [bigfile]_.
The format can be easily accessed by python, C, or Fortran.

The snapshots can be read by nbodykit via the `BigFileCatalog`, or via
the `DataSet` interface of `bigfile` python package, or directly accessed
as binary files from Fortran, C, or with the unix command `od`.

Position is stored in the unit of Mpc/h as double precision 3 vectors.
Velocity is stored in peculiar velocity km/s as single precision 3 vectors.
ID is stored as 8 bit unsigned integers.

The attributes of the `Header` column/dataset stores the meta data
about the simulation, including the fully evaluated parameter file.

A rarely used feature is to store the snapshot in the TPM container format by
Martin White, which can be subtly accessed by Python, C, or Fortran.

The FastPM container can be converted to a legacy Gadget-1/2 container format with
the supplied python script in python directory.

The nbody post-analysis package [NBODYKIT]_ natively supports bigfile
and TPM formats via the `BigFileCatalog` class and `TPMBinaryCatalog` interfaces.
Refer to nbodykit documentation (search for the class name) for examples.

[NBODYKIT]_ provides tools for calculating two point functions,
generating QPM mocks, and identifying Friends-of-Friends (also know as DBSCAN)
halos and calculating spherical overdensity properties of subhalos.

In addition to the snapshots, FastPM calculates and writes
the power-spectrum at each time step.
These power spectrum files are compatible with
`k, p = numpy.loadtxt(filename, unpack=True)` Note that these are the power spectrum of the density field that goes
into the force calculation, thus contain the mass assignment window effects.

FastPM can also write the Linear Density field, white noise, or non-linear density.
Use `write_lineark`, `write_whitenoise`, `white_nonlineark` in the parameter file.
The output is also stored as columns in a `bigfile` container. Read these via
the `BigFileMesh` class from [NBBODYKIT]_, or directly access the bigfile
container via the `bigfile` module.

Units
-----

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

FastPM shared a majority of the particle mesh code with the `pmesh` python package.

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
.. [FastPMPaper] http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1603.00476]
.. [FastPMNeutrinoPaper] https://ui.adsabs.harvard.edu/abs/arXiv:2007.13394].


Source Code Installation
------------------------

FastPM works out-of-the-box on many GNU/Linux and compatible systems.
The integrated continuous integration tests are performed on a Ubuntu
based Linux distribution via Travis-ci.org.

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

    cp Makefile.local.example Makefile.local

    # the rest is just make. It may take a while.
    make

Binary installation via Anaconda
--------------------------------

Anaconda is a popular Python distribution that provides portable
binary distributions of software on most x86-64 and Linux platforms.
FastPM compiles cleanly under the MPI provided by Anaconda.

Binaries for Linux-64 and OSX-64 are provided. Sorry we do not have
enough expertise on Windows builds.

The following command will install FastPM and nbodykit to the cfastpm
environment. 

.. code::

    conda create -n cfastpm
    conda activate cfastpm

    conda install -c bccp cfastpm nbodykit

Notice that there is a package called `fastpm` from Python,
which is a Python rewrite of FastPM that provides a playground for
different ParticleMesh based Poisson solvers.

For now, openmp does not seem to work with Anaconda, unless the
anaconda compiler is used (installed via gcc_linux-64), but this
currently interferes with the MPI compiler provided by the
mpich2 package. Most problems we solve with FastPM are small enough
that hybrid with threads is not necessary; for real large problems
we likely will need to recompiler from source code on the super-computer
anyways.

Anaconda Development Environment
--------------------------------

The current development is mainly performed on a Anaconda Linux-64 environment.

The following command creates the conda environment for development.

Notice that we install Anaconda's generic linux-64 gcc compiler and use the mpich
provided by the BCCP channel, which is a special version of mpich-3.2 that produces
correct binaries with the anaconda compiler.

.. code::

    conda create -n cfastpm
    conda activate cfastpm

    conda install -c bccp mpich gcc_linux-64 gsl

Please also refer to the file Makefile.dev.example.

Docker
------

There is a basic docker configuration file to set up a container for FastPM. 

To build it, run:

.. code::

    # first remove all prebuilt binary files

    make deep-clean

    sudo docker build -t fastpm .

To start the docker container in interactive mode, 
with port 8888 exposed and linking ``/my/file/directory`` to ``/worksapce``, run

.. code::

    sudo docker run -it -v /my/file/directory:/workspace -p 8888:8888 fastpm

We install a jupyter notebook service in the docker image, which listens on the
forwarded port of 8888.

.. code::

    jupyter notebook --ip=* --allow-root

As of now, proper set up of docker needs root access.
It may be necesssary to prepend `su -c` or `sudo` in docker command line, see [docker-root]_.

.. [docker-root] http://www.projectatomic.io/blog/2015/08/why-we-dont-let-non-root-users-run-docker-in-centos-fedora-or-rhel/

Examples
--------

- refer to tests/nbodykit.lua for a basic parameter file.
- refer to python/make-pklin.py for generation a linear power spectrum to start the simulation.
- refer to python/fof.py for halo finding. It is a MPI capable script that we
  frequently use on a few thousand cores!
- refer to python/convert-to-gadget-1.py for conversion from FastPM's bigfile to
  Gadget container format.
  The result can be used as an 2LPT or non-linear
  intial condition for Gadget.
  The script is currently sequential and takes about 6 hours
  to convert a `4096**3` simulation. 

Feel free to copy and modify these files to fit your own need, especially if you
have strong opinions on the choice data containers.

*Massive Neutrino Simulations*

- massive neutrinos are referred to as ncdm (not-cold dark matter) in the code.
- see [FastPMNeutrinoPaper]_ for details on the implementation.
- refer to tests/ncdm.lua for an example parameter file.

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

