from distutils.core import setup
import numpy
import os

def find_version(path):
    import re
    # path shall be a plain ascii text file.
    s = open(path, 'rt').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              s, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Version not found")


setup(
    name="fastpm", version=find_version("fastpm/version.py"),
    author="Yu Feng",
    description="FastPM in Python",
    package_dir = {'fastpm': 'fastpm'},
    packages= ['fastpm', 'fastpm.tests'],
    install_requires=['cython', 'numpy', 'scipy', 'pmesh', 'abopt', 'mpi4py_test'],
)

