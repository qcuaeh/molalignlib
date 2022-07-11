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

Running ralign
==============

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

Te following line will run the program using biases with a tolerance of 0.17 Å, iteration, mass weigthed coordinates, a stopping threshold of 10 counts, a maximum of 1000 trials and repeatable pseudo random numbers:

    ./bin/ralign tests/r005/100cobalt_j5.xyz -bias 0.17 -iter -weight mass -count 10 -trials 1000 -test
    
The ouput should look as follows:

     Map   Trial   Count   Cycles   Meanrot   Totalrot      RMSD
    ------------------------------------------------------------
       1       2      10      2.0      65.2     116.0     0.0482
       2       1       4      1.5      78.0     110.0     2.1060
    ------------------------------------------------------------
    Found 2 mapping(s) in 14 random trial(s)
