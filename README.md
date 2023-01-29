MolAlignLib
===========

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalignlib.git/HEAD?labpath=examples)

**MolAlignLib** is a fortran library based on random rotations and quasi-local RMSD minimizations to align rigid molecules and atom
clusters. It includes a command line and a Python interface.

![graphical abstract](abstract.png)

Before installing
-----------------

You will need a modern Fortran compiler and LAPACK to build the program from source, and additionally Python 3 to use
the python interface. It is recommended to install gfortran and LAPACK with your package manager:

in RHEL or Fedora use *yum*

    yum install gcc-gfortran lapack-devel

and in Debian, Ubuntu, etc. use *apt*

    apt install gfortran liblapack-dev

Install molalignlib with pip
----------------------------

The Python library only supports Python 3 so make sure that you are using the right version of *pip*:

    pip3 install molalignlib

It will install the *molalign* executable in your path and the *molalignlib* python module which provides the *assign_atoms* function
and the *align_to* method, which extends ASE's *Atoms* class (see example1.py):

    >>> from ase import Atoms
    >>> from molalignlib import assign_atoms

Build molalign from source 
--------------------------

Clone the repository:

    git clone https://github.com/qcuaeh/molalign.git

then enter its directory, rename *gnu.env* (or *intel.env* if you have the Intel fortran compiler and MKL installed) as *build.env*
and run:

    ./build.sh

It will create the *molalign* executable inside the *build* directory.

Program options
---------------

These options are supported by both, the native executable and the python script:

<code>-sort</code> Sort atoms by RMSD-minimizing order.  
<code>-trials *N*</code> Set maximum number of trials to *N*.  
<code>-fast</code> Enable biasing and iterative minimization steps.  
<code>-count *N*</code> Set the count threshold to *N* (defaults to 10).  
<code>-tol *TOL*</code> Set the biasing tolerance to *TOL* (defaults to 0.35 Å).  
<code>-out *NAME*</code> Set the output file name to *NAME* (defaults to *aligned.xyz*).  
<code>-rec *N*</code> Record up to *N* assignments (defaults to 1).  
<code>-enan</code> Align with enantiomer (mirrored coordinates).  
<code>-stats</code> Print detailed stats of the calculation.  
<code>-mass</code> Use mass weighted coordinates.  

These options are only supported by the native executable:

<code>-live</code> Show progress in real time.  
<code>-test</code> Produce repeatable results for testing.  
<code>-stdin *EXT*</code> Read coordinates from standard input in *EXT* format.  
<code>-stdout *EXT*</code> Write coordinates to standard output in *EXT* format.  
 
Notice that the input format is determined from the extensions of the input files. Also bear in mind that the format specifier
*EXT* must be an standard three or four letter extension, and that the native executable only supports *xyz* as input format and
*xyz* and *mol2* as output formats, while the python script supports many more input and oputput formats.

Basic usage
-----------

The syntax of the command is

    molalign [option[s]] file [file]

If only a file is specified then two sets of coordinates will be read from the file, otherwise a single set of coordinates will
be read from each file.

To align the atoms without reordering run the command without options and to align with optimal assigment add the `-sort` option.

Advanced usage
--------------

The performance of the computation when sorting is requested can be greatly improved by enabling biasing and iteration with the
`-fast` option, however if the tolerance is set too small the assignment will fail. In such cases the tolerance must be increased
with the `-tol` option to a value larger than the maximum expected displacement of the atoms.

By default a threshold of 10 counts is used to stop the computation, reducing this threshold with the `-count` option will
proportionally reduce the computation time but the probability of obtaining suboptimal assignmnets will increase. To avoid too
long computations the `-trials` option can be used to stop the computation before reaching the count threshold.

The algorithm always explores multiple possible assignments, but only the best one is stored by default. For symmetric clusters
multiple equivalent assignments are possible so it can be useful to print more than one result using the `-rec` option.

Examples
--------

For small atom displacements the default tolerance is enough:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast
    Optimized RMSD = 0.0506
    
    ./build/molalign examples/Co138_0.xyz examples/Co138_2.xyz -sort -fast
    Optimized RMSD = 0.0977

but for atom displacements larger than the tolerance the assignment will fail:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast
    Error: Assignment failed

Increasing the tolerance will fix the problem but the calculation will slow down:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast -tol 0.7
    Optimized RMSD = 0.1973

Sometimes there is more than one optimal assignment due to symmetry:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast -stats -rec 5
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10     10.0      74.4      54.5     0.0506
       2       10      9.0      69.8      48.5     0.0506
       3       16     11.3      83.2      63.6     0.0506
       4        1      2.0       8.4       7.1     0.6652
       5        1      7.0      30.5       9.8     0.6716
    -----------------------------------------------------
    Random trials = 81
    Minimization steps = 706
    Visited local minima > 5
    Optimized RMSD = 0.0506

The ouput shows that there are 3 equivalent optimal assignments.

Note: The stats can be different than those shown in the examples, due to the use of different random seeds on each run,
but the optimized RMSD should be the same.
