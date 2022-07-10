Compile the program
===================

Run ./build without arguments to compile the executable with the highest
optimization level. It accepts the following options to change the default
behaviour:

### Optimization level options

-fast : Build the optimized version (default)  
-debug : Build the debug version  
  
### Build type options

-exe : Build the executable program (default)  
-lib : Build the shared libraries  

Run ./rebuild instead of ./build to recompile all source files from scratch.
Both scripts accept the same options.

### Notes

* Executables are created in the "bin" directory
* Object files are created in the "_build_dir" directory
* The compilation list is located in ./source/compilelist
* The optimized version is compiled with fast math

Run the program
===============

### Program options

-maps <recorded maps>
-out <output format (xyz|mol2)>
-live
-iter
-bias <bias tolerance>
-weight <weighting property>
-count <maximum count>
-trials <maimum trials>
-scale <bias scale>
-test

### Examples

./bin/ralign tests/r005/100cobalt_j5.xyz -iter -weight mass -bias 0.17 -test -trials 1000 -count 10
