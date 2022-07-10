Compile the program
===================

Run `./build` without arguments to build the program with the highest
optimization level. The script accepts the following options to change the
defaults:

&ensp; -slow : Build without optimizations.  
&ensp; -debug : Build with complete checks and debug info.  
&ensp; -recomp : Recompile all source files from scratch.  

The executables will be created in the *bin* directory.

Compile the shared library
==========================

Run `./build-lib` without arguments to compile the shared library with the
highest optimization level. The script accepts the same options as the *build*
script and the libraries will be also created in the *bin* directory.

Run the program
===============

### Program options

&ensp; -maps MAX : Set number of recorded maps to MAX.  
&ensp; -out xyz|mol2 : Set output format to *XYZ* or *MOL2*.  
&ensp; -live : Print live stats.  
&ensp; -iter : Use iterative steps.  
&ensp; -bias TOL : Use biasing with tolerance TOL.  
&ensp; -weight none|mass : Set weighting property to *None* or *Atomic Mass*.  
&ensp; -count MAX : Set maximum map counting to MAX.  
&ensp; -trials MAX : Set maximum number of trials to MAX.  
&ensp; -scale FACTOR : Set length scale to FACTOR.  
&ensp; -test : Use a repeateble random number sequence for testing.  
 
### Examples
 
./bin/ralign tests/r005/100cobalt_j5.xyz -iter -weight mass -bias 0.17 -test -trials 1000 -count 10
