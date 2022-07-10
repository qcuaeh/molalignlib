Compile the program
===================

Run ./build without arguments to compile the executable with the highest
optimization level. Run ./rebuild instead of ./build to recompile all source files from scratch.
Both scripts accept the same options:

### Optimization level options

-fast : Build the optimized version (default)  
-debug : Build the debug version  
  
### Build type options

-exe : Build the executable program (default)  
-lib : Build the shared libraries  

### Notes

* Executables are created in the "bin" directory
* Object files are created in the "_build_dir" directory
* The compilation list is located in ./source/compilelist
* The optimized version is compiled with fast math

Run the program
===============

### Program options

-maps INTEGER : Number of recorded maps
-out xyz|mol2 : Output format (XYZ or MOL2)

-live : Print live stats
-iter : Iterative trial step
-bias REAL : Tolerance for biasing
-weight STRING : Weighting property
-count INTEGER : Maximum map counting
-trials INTEGER : Maximum number of trials
-scale REAL : Length scale
-test : Repeateble run for testing

### Examples

./bin/ralign tests/r005/100cobalt_j5.xyz -iter -weight mass -bias 0.17 -test -trials 1000 -count 10
