Building
========

Run `./build` without arguments to build the program in double precision with
all optimizations enabled. The script accepts the following options to build
the library instead of the program or to change the optimization level and the
numeric precision:

### Library type (-l)

&ensp; static : Build the static library
&ensp; shared : Build the shared library
&ensp; python : Build the python library.  

### Optimization level (-o)

&ensp; fast : Enable all optimizations (default).  
&ensp; debug : Enable debugging and disable optimizations.  

### Real number precision (-r)

&ensp; single : Use single precision.  
&ensp; double : Use double precision (default).  

### Quick recompile (-q)

The program and libraries will be created in the *bin* directory.

Running the program
===================

### Program options

&ensp; -live : Print live stats.  
&ensp; -iter : Use iterative steps.  
&ensp; -stdin : Read coordinates from standard input.  
&ensp; -test : Use the same pseudo random numbers on every run.  
&ensp; -rec *NUM* : Set number of printed solutions to *NUM*.  
&ensp; -out xyz|mol2 : Set output format to XYZ or Mol2.  
&ensp; -weight none|mass : Set weights to unity or atomic masses.  
&ensp; -count *MAX* : Set map counting threshold to *MAX*.  
&ensp; -trial *MAX* : Set maximum number of trials to *MAX*.  
&ensp; -bias *TOL* : Use biasing with tolerance *TOL*.  
&ensp; -scale *NUM* : Set length scale to *NUM*.  
 
### Examples

Te following line will run the program using biases with a tolerance of 0.17 Å,
iteration, mass weigthed coordinates, a stopping threshold of 10 counts, a maximum
of 1000 trials and repeatable pseudo random numbers:

    ./bin/ralign tests/r005/100cobalt_j5.xyz -bias 0.17 -iter -weight mass -count 10 -trial 1000 -test -rec 10
    
The ouput should look as follows:

     Map   Trial   Count   Cycles   Meanrot   Totalrot      RMSD
    ------------------------------------------------------------
       1       2      10      2.0      65.2     116.0     0.0482
       2       1       4      1.5      78.0     110.0     2.1060
    ------------------------------------------------------------
    Found 2 mapping(s) in 14 random trial(s)
