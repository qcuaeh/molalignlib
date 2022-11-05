molalign
========
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalign.git/HEAD?labpath=tests)

**molalign** is a fortran code which applies random rotations and pseudo
local minimizations of RMSD to find reliable solutions to the near-congruence
problem for atom clusters.

![graphic1](assets/graphic1.png)

Building
--------

To build the program you will need a machine with Bash, GNU or Intel
Fortran compiler, LAPACK, and optionally, with Python 3 and F2PY.

Clone the repository and enter in the directory, then create a symbolic
link named *build.env* pointing to either *gnu.env* or *intel.env*,
depending on your platform, and run `./build.sh` to  build the executable.
You may use the following options to customize the building:

-d : Compile the debug version of the executable.  
-q : Quick compile (compile only modified source files).  
-r[4|8] : Use 4-byte or 8-byte (default) floating point numbers.  
-l[s|d|py] : Build a static/dynamic/python library instead of an executable.  

After runnig the script the executable and/or libraries will be created
in the same directory.

Usage
-----

### Options

-sort : Sort atoms by optimal assignment.  
-bias : Enable biasing and iterative convergence.  
-tol *ϵ* : Set biasing tolerance to *ϵ* (defaults to 0.2 Å).  
-test : Use the same pseudo random numbers on every run.  
-count *N* : Set the count convergence threshold to *N* (defaults to 10).  
-out xyz|mol2 : Set the output format to XYZ or Mol2 (defaults to XYZ).  
-rec *N* : Set the number of recorded solutions to *N* (defaults to 1).  
-trials *MAX* : Set the maximum number of trials to *MAX*.  
-scale *α* : Set biasing scale to *α* (defaults to 1000).  
-stdin : Read coordinates from standard input.  
-weight : Use mass weighted RMSD.  
-live : Show progress in real time.  
 
### Basic usage

To align two molecules without reordering just run:

    ./molalign tests/r005/Co100.xyz

To align two molecules with reordering you can run:

    ./molalign tests/r005/Co100.xyz -sort

by default a convergence threshold of 10 counts is used but a lower threshold
can be used to reduce the computation time at the expense of less reliable
results:

    ./molalign tests/r005/Co100.xyz -sort -count 3

When reordering is requested performance can be improved by two to four orders
of magnitude by enabling biasing. For example to reorder and align two molecules
applying biasing run:

    ./molalign tests/r005/Co100.xyz -sort -bias

A biasing tolerance of 0.2 Å is used by default, which is enough to account for
small distortions on the geometry, however if the distortion is too large then the
tolerance must be increased, for example to increase the tolerance to 0.35 Å run:

    ./molalign tests/r01/Co100.xyz -sort -bias -tol 0.35

Tolerance must be chosen carefully because a value too small can lead to wrong
results while a value too large will not improve performance.

### Advanced usage

Most of the options not covered in the basic usage are intended for testing purposes.

The option `-test` forces the generation of the same stream of random numbers on
every run in order to have reproducible results:

    ./molalign tests/r005/Co100.xyz 10 -sort -bias -test

The algorithm by default generates several possible reorderings, but only one is
printed by default. To print more than one solution use the `-rec` option:

    ./molalign tests/r005/Co100.xyz 10 -sort -bias -rec 10

By default the biasing weight is much larger than the euclidean cost but
in can be reduced with the `-scale` option:

    ./molalign tests/r005/Co100.xyz 10 -sort -bias -scale 10

however it is rarely useful for clusters.

To avoid too long computations the option `-trials` can be used:

    ./molalign tests/r005/Co100.xyz -sort -bias -trials 1000

which will halt the computation if the number of trials exceeds 1000.

To use mass weighted coordinates instead of unweighted coordinates use the
`-weight` option:

    ./molalign tests/r005/Co100.xyz -sort -bias -weight

Finally the options `-stdin`, `-out` and `-live` just control the source and the
destination of the input and results, and the way that progress is displayed.

### Examples

Notice that because the random number generator is initialized with a different
seed each time the output will not be exactly the same as in the examples (to
reproduce exactly the same output use the option `-test`):

To run the program with biased costs run:

    ./molalign tests/r005/Co100.xyz -sort -bias
 
The ouput should look as follows:

     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10      2.0     124.7     119.4     0.0482
    -----------------------------------------------------
    Found more than 1 mapping(s) in 12 random trial(s)

To run the program with biased costs (recording up to 10 mappings) run:

    ./molalign tests/r005/Au161Pd40.xyz -rec 10 -sort -bias

The output should look as follows:

     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10      3.4      62.5      35.0     0.0497
       2       12      3.7      70.9      40.7     0.0497
       3        8      3.9      69.2      43.3     0.0497
       4        8      4.0      75.6      41.9     0.0497
       5       13      2.9      56.6      31.6     0.0497
       6       15      3.0      63.1      30.3     0.0497
       7        7      2.3      51.3      30.8     0.0497
       8       10      3.9      71.6      41.4     0.0497
       9        1      1.0       4.3       4.3     5.0285
      10        1      2.0      10.9       7.6     5.0285
    -----------------------------------------------------
    Found more than 10 mapping(s) in 97 random trial(s)

To run the program with biased costs and mass weighted coordinates (recording up to 10 mappings) run:

    ./molalign tests/r005/Au161Pd40.xyz -rec 10 -sort -bias -weight

The output should look as follows:

     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10      2.9      64.8      33.6     0.1010
       2       16      3.1      61.1      32.5     0.1010
       3       15      3.5      69.0      38.2     0.1010
       4       14      4.3      75.3      44.4     0.1010
       5       11      2.8      53.2      31.6     0.1010
       6       17      3.8      68.9      38.0     0.1010
       7       14      3.4      62.6      38.2     0.1010
       8       14      2.6      52.9      33.8     0.1010
       9        1      3.0      29.5      14.9     5.1384
      10        1      1.0       6.5       6.5     5.1384
    -----------------------------------------------------
    Found more than 10 mapping(s) in 125 random trial(s)

