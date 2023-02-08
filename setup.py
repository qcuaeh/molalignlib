import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from os import path

f90_files = []
with open(path.join('molalignlib', 'f90_files'), 'r') as file:
    for line in file.readlines():
        if not line.startswith('#'):
            f90_files.append(path.join('molalignlib', line.rstrip()))

f2py_files = []
with open(path.join('molalignlib', 'f2py_files'), 'r') as file:
    for line in file.readlines():
        if not line.startswith('#'):
            f2py_files.append(path.join('molalignlib', line.rstrip()))

molalignlib = ('molalignlib', dict(sources=f90_files, extra_f90_compile_args=['-O3', '-ffast-math']))
molalignlibext = Extension('molalignlibext', sources=f2py_files, libraries=['molalignlib', 'lapack'])

setup(
    libraries = [molalignlib],
    ext_modules = [molalignlibext]
)
