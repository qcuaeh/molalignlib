import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from os import path

source_files = []
with open(path.join('molalignlib', 'source_files'), 'r') as file:
    for line in file.readlines():
        if not line.startswith('#'):
            source_files.append(path.join('molalignlib', line.rstrip()))

f2py_files = []
with open(path.join('molalignlib', 'f2py_files'), 'r') as file:
    for line in file.readlines():
        if not line.startswith('#'):
            f2py_files.append(path.join('molalignlib', line.rstrip()))

molalignlib = ('molalignlib', dict(sources=source_files, extra_f90_compile_args=['-O3', '-ffast-math']))
molalignlibext = Extension('molalignlibext', sources=f2py_files, libraries=['molalignlib', 'lapack'])

setup(
    libraries = [molalignlib],
    ext_modules = [molalignlibext]
)
