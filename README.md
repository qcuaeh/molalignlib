Compile the program
===================

Run `./build` without arguments to build the program with the highest
optimization level. The script accepts the following options:

&ensp; -slow : Build without optimizations.  
&ensp; -debug : Build with complete checks and debug info.  
&ensp; -recompile : Recompile all source files from scratch.  

The program files will be created in the *bin* directory.

Compile the shared library
==========================

Run `./build-lib` without arguments to compile the shared library with the
highest optimization level. This script accepts the same options as the *build*
script and the library files will be also created in the *bin* directory.

Run the program
===============

### Program options

&ensp; -live : Print live stats.  
&ensp; -iter : Use iterative steps.  
&ensp; -stdin : Read coordinates from standard input.  
&ensp; -test : Use the same pseudo random number sequence on every run.  
&ensp; -maps MAX : Set number of recorded maps to MAX.  
&ensp; -out xyz|mol2 : Set output format to *XYZ* or *MOL2*.  
&ensp; -weight none|mass : Set weighting property to *None* or *Atomic Mass*.  
&ensp; -count MAX : Set maximum map counting to MAX.  
&ensp; -trials MAX : Set maximum number of trials to MAX.  
&ensp; -bias TOL : Use biasing with tolerance TOL.  
&ensp; -scale SCALE : Set length scale to SCALE.  
 
### Examples
 
./bin/ralign tests/r005/100cobalt_j5.xyz -iter -weight mass -bias 0.17 -test -trials 1000 -count 10
