Building
========

Run `./build` without arguments to build the program in double precision
with optimizations enabled. The program will be created in the *bin*
directory. Moreover the script accepts the following options:

#### Build type:

&ensp; -program : Build the program (default).  
&ensp; -library : Build the shared library instead.  

#### Optimization level:

&ensp; -fast : Enable optimizations with fast math (default).  
&ensp; -debug : Disable optimizations and include debug info.  

#### Numeric precision:

&ensp; -single : Use single precision.  
&ensp; -double : Use double precision (default).  

#### Miscellaneous:

&ensp; -recompile : Recompile all sources from scratch.  

Running the program
===================

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
