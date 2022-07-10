Compile the program
===================

Run ./build without arguments to build the program with the highest
optimization level.

Compile the library
===================

Run ./build-lib without arguments to compile the library with the highest
optimization level.

Both scripts accept the same options:

-debug : Build the debug version
-rebuild : Recompile all source files from scratch

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
