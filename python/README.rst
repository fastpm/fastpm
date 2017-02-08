fastpm
======

This is the Python implementation of the FastPM numerical scheme for quasi-nbody simulations.

Install
-------

The best result is obtained by using anaconda.

Python 3 is the current development environment.

First set up the basic requirements,

.. code::

    conda install cython numpy scipy mpi4py nose

    env LD_LIBRARY_PATH=$CONDA_PREFIX/lib pip install pmesh bigfile

    conda install dask h5py

FastPM is build as a forward model with `abopt`, so we need a recent version
of that:

.. code::

    pip install https://github.com/bccp/abopt/archive/master.zip


For running the nonlinear reconstruction code,
we need a recent version of nbodykit and matplotlib

.. code::

    pip install https://github.com/bccp/nbodykit/archive/master.zip

    conda install matplotlib


Finally, we can install fastpm, either from PYPI (the latest release, we have none yet)

.. code::

    pip install fastpm

or from the git clone :

.. code::

    python setup.py install


`LD_LIBRARY_PATH` hack
----------------------

When installing pmesh, prefix `LD_LIBRARY_PATH` helps
the compilation of a package named `pfft-python`, the parallel
FFT software we use.


Development
-----------

To run the tests

.. code::

    python runtests.py

.. code::

    python runtests.py --mpirun

To run a single test (e.g. `test_fastpm.py:test_name`) :

.. code::

    python runtests.py --mpirun -t fastpm/tests/test_fastpm:test_name



