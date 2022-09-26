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

### Program options

-live: Print live stats.  
-iter: Perform iterated trials.  
-tol *ϵ*: Set biasing tolerance to *ϵ*.  
-count *N*: Set the count convergence threshold to *N*.  
-trials *MAX*: Set the maximum number of trials to *MAX*.  
-bias: Use biased distances (must be used together with -tol).  
-rec *N*: Set the number of recorded solutions to *N* (defaullt is *N* = 1).  
-out xyz|mol2: Set the output format to XYZ or Mol2 (default is XYZ).  
-scale *α*: Set biasing scale to *α* (default is *α* = 1000).  
-test: Use the same pseudo random numbers on every run.  
-stdin: Read coordinates from standard input.  
-weight: Use mass weighted distances.  
 
### Examples

#### Example 1

The following command will run the program with
up to 10 recorded mappings,
a convergence threshold of 10 counts,
biasing with a tolerance of 0.17 Å,
iteration,
mass weighted distances
and repeatable pseudo random numbers:

    ./bin/molalign tests/r005/Co100.xyz -rec 10 -count 10 -bias -tol 0.17 -iter -test
 
The ouput should look as follows:

     Map   Trial   Count   Cycles   Meanrot   Totalrot      RMSD
    ------------------------------------------------------------
       1       2      10      2.0      65.2     116.0     0.0482
       2       1       4      1.5      78.0     110.0     2.1060
    ------------------------------------------------------------
    Found 2 mapping(s) in 14 random trial(s)


#### Example 2

The following command will run the program with
a maximum of 2000 trials,
biasing with a tolerance of 0.17 Å,
iteration,
mass weighted distances
and repeatable pseudo random numbers:

    ./bin/molalign tests/r005/Au161Pd40.xyz -trials 2000 -bias -tol 0.17 -iter -weight -test

The output should look as follows:

     Map   Trial   Count   Cycles   Meanrot   Totalrot      RMSD
    ------------------------------------------------------------
       1      14      17      1.0      21.9      21.9     0.1010
    ------------------------------------------------------------
    Found more than 1 mapping(s) in 2000 random trial(s)

