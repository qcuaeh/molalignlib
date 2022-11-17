from os import path, symlink, remove
from setuptools import setup

setupdir = path.dirname(path.abspath(__file__))

symlink(path.join(setupdir, 'scripts', 'gnu.env'), path.join(setupdir, 'molalignlib', 'build.env'))
symlink(path.join(setupdir, 'scripts', 'compile.sh'), path.join(setupdir, 'molalignlib', 'compile.sh'))
symlink(path.join(setupdir, 'scripts', 'link_pythonlib.sh'), path.join(setupdir, 'molalignlib', 'link_pythonlib.sh'))

setup()

remove(path.join(setupdir, 'molalignlib', 'build.env'))
remove(path.join(setupdir, 'molalignlib', 'compile.sh'))
remove(path.join(setupdir, 'molalignlib', 'link_pythonlib.sh'))
