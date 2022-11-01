molalign
========
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalign.git/HEAD?labpath=tests)

**molalign** is a fortran program based on random rotations and the Hungarian algorithm to solve the approximate congruence problem for atomic systems.

![graphic1](assets/graphic1.png)

Building
--------

Run `./build.sh` without arguments to build the program in double
precision with all optimizations enabled. The script accepts the
following options to change the default behaviour:

-q: Recompile only modified source files.  
-d: Compile non optimized code with debug info.  
-r single|double: Set the real precision to single or double (default).  
-l static|shared|python: Build the selected library instead of the program.  

After runnig the script the program will be created in the *bin* directory and the libraries in the *lib* directory.

Usage
-----

### Options

-sort: Sort atoms.  
-bias: Enable biasing and iterative convergence.  
-tol *ϵ*: Set biasing tolerance to *ϵ* (defaults to 0.2 Å).  
-count *N*: Set the count convergence threshold to *N* (defaults to 10).  
-trials *MAX*: Set the maximum number of trials to *MAX*.  
-out xyz|mol2: Set the output format to XYZ or Mol2 (defaults to XYZ).  
-rec *N*: Set the number of recorded solutions to *N* (defaults to 1).  
-test: Use the same pseudo random numbers on every run.  
-scale *α*: Set biasing scale to *α* (defaults to 1000).  
-stdin: Read coordinates from standard input.  
-weight: Use mass weighted distances.  
-live: Show progress in real time.  
 
### Basic usage

To align two molecules without reordering just run:

    ./bin/molalign tests/r005/Co100.xyz

To align two molecules with reordering you can run:

    ./bin/molalign tests/r005/Co100.xyz -sort

by default a convergence threshold of 10 counts is used but a lower threshold
can be used to reduce the computation time at the expense of less reliable
results:

    ./bin/molalign tests/r005/Co100.xyz -sort -count 3

When reordering is requested performance can be improved by two to four orders
of magnitude by enabling biasing. For example to reorder and align two molecules
applying biasing run:

    ./bin/molalign tests/r005/Co100.xyz -sort -bias

A tolerance of 0.2 Å is used by default, which is enough to account for numerical
errors, however if the clusters are only approximately congruent then a larger
tolerance is required, for example to apply biasing with a tolerance of 0.1 Å run:

    ./bin/molalign tests/r01/Co100.xyz -sort -bias -tol 0.35

Tolerance must be carefully determined because a too small value will lead to
wrong results while a too large value will not improve performance.

### Advanced usage

Most of the options not covered in the basic usage are intended for testing purposes.

The option `-test` forces the generation of the same stream of random numbers on
every run in order to have reproducible results:

    ./bin/molalign tests/r005/Co100.xyz 10 -sort -bias -test

The algorithm by default generates several possible reorderings, but only one is
printed by default. To print more than one solution use the `-rec` option:

    ./bin/molalign tests/r005/Co100.xyz 10 -sort -bias -rec 10

By default the biasing weight is much larger than the euclidean cost but
in can be reduced with the `-scale` option:

    ./bin/molalign tests/r005/Co100.xyz 10 -sort -bias -scale 10

however it is rarely useful for clusters.

To avoid too long computations the option `-trials` can be used:

    ./bin/molalign tests/r005/Co100.xyz -sort -bias -trials 1000

which will halt the computation if the number of trials exceeds 1000.

To use mass weighted coordinates instead of unweighted coordinates use the
`-weight` option:

    ./bin/molalign tests/r005/Co100.xyz -sort -bias -weight

Finally the options `-stdin`, `-out` and `-live` just control the source and the
destination of the input and results, and the way that progress is displayed.

### Examples

To run the program using mass weighted coordinates (recording up to 10 mappings) run:

    ./bin/molalign tests/r005/Co100.xyz -test -sort -bias
 
The ouput should look as follows:

     Map   Count   Cycles   Meanrot   Totalrot      RMSD
    ----------------------------------------------------
       1      10      2.0      58.5     104.2     0.0482
    ----------------------------------------------------
    Found 1 mapping(s) in 10 random trial(s)

To run the program limited up to 1000 trials (recording up to 10 mappings) run:

    ./bin/molalign tests/r005/Au161Pd40.xyz -test -rec 10 -sort -bias

The output should look as follows:

     Map   Count   Cycles   Meanrot   Totalrot      RMSD
    ----------------------------------------------------
       1      10      3.4      19.4      35.0     0.0497
       2      12      3.7      19.2      40.7     0.0497
       3       8      3.9      19.1      43.3     0.0497
       4       8      4.0      19.2      41.9     0.0497
       5      13      2.9      21.5      31.6     0.0497
       6      15      3.0      21.9      30.3     0.0497
       7       7      2.3      23.9      30.8     0.0497
       8      10      3.9      19.4      41.4     0.0497
       9       1      1.0       4.3       4.3     5.0285
      10       1      2.0       5.5       7.6     5.0285
    ----------------------------------------------------
    Found more than 10 mapping(s) in 97 random trial(s)

