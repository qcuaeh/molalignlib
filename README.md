molalign
========
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalign.git/HEAD?labpath=tests)

**molalign** is a fortran program to solve the approximate congruence problem for atomic systems.

Building
--------

Run `./build.sh` without arguments to build the program in double precision with all
optimizations enabled. The script accepts the following options to build the library
instead of the program or to change the optimization level and the numeric precision:

&ensp; -o fast|debug : Set the optimization level.  
&ensp; -r single|double : Set the precision for real numbers.  
&ensp; -l static|shared|python : Build the library instead of the program.  
&ensp; -q : Quick compile (recompile only modified source files).  

The program and libraries will be created in the *bin* directory.

Usage
-----

### Program options

&ensp; -live : Print live stats.  
&ensp; -iter : Use iterative steps.  
&ensp; -mass : Use mass weighted distances.  
&ensp; -stdin : Read coordinates from standard input.  
&ensp; -out xyz|mol2 : Set output format to XYZ or Mol2.  
&ensp; -test : Use the same pseudo random numbers on every run.  
&ensp; -rec *NUM* : Set the number of recorded solutions to *NUM*.  
&ensp; -count *NUM* : Set count convergence threshold to *NUM*.  
&ensp; -trial *MAX* : Set maximum number of trials to *MAX*.  
&ensp; -bias *TOL* : Use biasing with tolerance *TOL*.  
&ensp; -scale *SCALE* : Set length scale to *SCALE*.  
 
### Examples

Te following command will run the program with
up to 10 recorded mappings,
a convergence threshold of 10 counts,
a maximum of 1000 trials,
biasing with a tolerance of 0.17 Å,
iteration,
mass weighted distances
and repeatable pseudo random numbers:

    ./bin/molalign tests/r005/Co100.xyz -rec 10 -count 10 -max 1000 -bias 0.17 -iter -mass -test
    
The ouput should look as follows:

     Map   Trial   Count   Cycles   Meanrot   Totalrot      RMSD
    ------------------------------------------------------------
       1       2      10      2.0      65.2     116.0     0.0482
       2       1       4      1.5      78.0     110.0     2.1060
    ------------------------------------------------------------
    Found 2 mapping(s) in 14 random trial(s)
