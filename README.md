MolAlign
========

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalign.git/HEAD?labpath=examples)

MolAlign is a fortran library based on random rotations and pseudo local minimizations to reliably find solutions to the
near-congruence problem for atom clusters.

![graphical abstract](abstract.png)

Before installing
-----------------

You will need a modern Fortran compiler and LAPACK to build the program from source or to install the python package with pip. It is
recommended to use your package manager to install them:

in RHEL or Fedora use *yum*

    yum install gcc-gfortran lapack-devel

and in Debian, Ubuntu, etc. use *apt*

    apt install gfortran liblapack-dev

Install molalign with pip
----------------------------

The python library only supports Python 3 so make sure that you are using the right version of *pip*:

    pip3 install molalign

It will install the *molalign* executable in your path and the *molalignlib* python module which provides the classes *Align* and
*Assign*:

    >>> from molalignlib import Atoms, assign_atoms

Build molalign from source 
--------------------------

Clone this repository

    git clone https://github.com/qcuaeh/molalign.git

then enter its directory, rename *gnu.env* (or *intel.env* if you have the Intel fortran compiler and MKL installed) as *build.env*
and run:

    ./build.sh

It will create the executable inside the *build* directory.

Program options
---------------

These options are supported by both, the native executable and the python script:

<code>-sort</code> Sort atoms by optimal assignment.  
<code>-fast</code> Enable biasing and iterative convergence.  
<code>-tol *TOL*</code> Set biasing tolerance to *TOL* (defaults to 0.35 Å).  
<code>-count *N*</code> Set count threshold to *N* (defaults to 10).  
<code>-trials *N*</code> Set maximum number of trials to *N*.  
<code>-rec *N*</code> Set number of recorded solutions to *N* (defaults to 1).  
<code>-out *EXT*</code> Set output file format to *EXT* (defaults to *xyz*).  
<code>-test</code> Generate the same sequence of random numbers on every run.  
<code>-mass</code> Use mass weighted coordinates.  
<code>-enan</code> Use mirrored coordinates.  

These options are only supported by the native executable:

<code>-live</code> Show progress in real time.  
<code>-stdin *EXT*</code> Read coordinates from standard input as an *EXT* file.  
 
Notice that the input format is determined from the extensions of the input files. Also bear in mind that the format specifier
*EXT* must be an standard three or four letter extension and that the native executable only supports *xyz* as input format and
*xyz* and *mol2* as output formats, while the python script supports many more input and oputput formats.

Basic usage
-----------

The syntax of the command is

    molalign [option[s]] file [file]

If only a file is specified then two sets of coordinates will be read from the file, otherwise a single set of coordinates will be
read from each file.

To align the atoms without reordering run the command without options and to align with optimal assigment add the `-sort` option.

Advanced usage
--------------

The performance of the optimal assignment algorithm can be improved by around three orders of magnitude by enabling biasing and
iteration with the `-fast` option, however, this option can potentially lead to suboptimal results when the difference between the
coordinates is too large.

If the RMSD increases when the `-fast` option is used then the tolerance must be increased from the default of 0.35 Å to a larger
value with the `-tol` option.

By default a threshold of 10 counts is used to stop the computation, but it can be changed with the `-count` option. Reducing
this threshold will proportionally reduce the computation time but the probability of suboptimal assignmnets will increase.

The algorithm always explores multiple possible assignments, but only the best one is stored by default. For symmetric clusters it
can be useful to report more than one assignment using the `-rec` option.

To avoid too long computations the `-trials` option can be used to stop the computation if the number of trials exceeds *N*.

To have reproducible results use the option `-test`, which will force the generation of the same stream of random numbers on every run.

Examples
--------

To reproduce exactly the same output as in these examples include the option `-test`.

For small distortions the default tolerance is enough:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10     10.0      74.4      54.5     0.0506
    -----------------------------------------------------
    
    ./build/molalign examples/Co138_0.xyz examples/Co138_2.xyz -sort -fast
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10     14.4      84.9      65.7     0.0977
    -----------------------------------------------------

but if the maximum distortion is larger than the tolerance then a wrong alignmnet can be obtained:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10      1.0     133.3     133.3     2.9315
    -----------------------------------------------------

Increasing the tolerance will fix the problem but will significatively slow the calculation:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast -tol 0.69
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10      2.4      15.1       8.2     0.1973
    -----------------------------------------------------

Sometimes it is necessary to record more than one alignment due to the cluster symmetry, for example:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast -rec 5
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10     10.0      74.4      54.5     0.0506
       2       10      9.0      69.8      48.5     0.0506
       3       16     11.3      83.2      63.6     0.0506
       4        1      2.0       8.4       7.1     0.6652
       5        1      7.0      30.5       9.8     0.6716
    -----------------------------------------------------

The ouput shows that there are 3 degenerated solutions due to the symmetry of the cluster.
