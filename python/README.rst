fastpm
======

This is the Python implementation of the FastPM numerical scheme for quasi-nbody simulations.

Install
-------

The best result is obtained by using anaconda.

Python 3 is the current development environment.

First the requirements:

.. code::

    conda install cython numpy scipy mpi4py nose

    pip install pmesh bigfile kdcount
    conda install dask h5py

    pip install https://github.com/bccp/abopt/archive/master.zip
    pip install https://github.com/bccp/nbodykit/archive/master.zip

Then install fastpm, either from PYPI (the latest release)

.. code::

    pip install fastpm

or from the git clone :

.. code::

    python setup.py install


MAC special notes
-----------------

When installing pmesh, prefix `LD_LIBRARY_PATH` helps
the compilation of a package named `pfft-python`, the parallel
FFT software we use.

.. code::

    env LD_LIBRARY_PATH=$CONDA_PREFIX/lib pip install pmesh bigfile kdcount


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



