**MolAlignLib** is a Fortran and Python library based on random rotations and quasi-local RMSD minimizations to align rigid
molecules and clusters. The details of the method can be found in publication [[1]](#1).

![graphical abstract](abstract.png)

Before installing
-----------------

You will need gfortran and LAPACK to build the program from source.
It is recommended to install them with your package manager:

in RHEL, Fedora, etc. use *yum*

    yum install gcc-gfortran lapack-devel

or in Debian, Ubuntu, etc. use *apt*

    apt install gfortran liblapack-dev

Build from source 
-----------------

Clone the repository:

    git clone https://github.com/qcuaeh/molalignlib.git

then enter the cloned directory, edit the *build.env* file to suit your system, and run:

    ./build.sh

It will create the *molalign* executable inside the *build* directory.

The library also has Python interface, to install it run:

    python3 setup.py install --user

It will compile and install the *molalign* Python script and the *molalignlib* Python module in your path.

or install with pip
-------------------

You can also install the library with *pip*:

    pip3 install molalignlib

It will install the *molalign* Python script and the *molalignlib* Python module in your system path.

or try on Binder
----------------

You can try the Python interface without installing anything:

* [Open example1.ipynb on Binder](https://mybinder.org/v2/gh/qcuaeh/molalignlib.git/HEAD?labpath=examples/example1.ipynb)

Program options
---------------

These options are supported by both, the native executable and the python script:

<code>-sort</code>&nbsp; Reorder the atoms to minimize the RMSD.  
<code>-trials <em>N</em></code>&nbsp;  Set maximum number of trials to *N*.  
<code>-fast</code>&nbsp; Enable biasing and iterative minimization steps.  
<code>-count <em>N</em></code>&nbsp; Set the count threshold to *N* (defaults to 10).  
<code>-tol <em>TOL</em></code>&nbsp; Set the biasing tolerance to *TOL* (defaults to 0.35 Å).  
<code>-out <em>NAME</em></code>&nbsp; Set the output file name to *NAME* (defaults to *aligned.xyz*).  
<code>-rec <em>N</em></code>&nbsp; Record up to *N* assignments (defaults to 1).  
<code>-stats</code>&nbsp; Print detailed stats of the calculation.  
<code>-test</code>&nbsp; Use a fixed random seed for testing.  
<code>-mirror</code>&nbsp; Reflect aligned coordinates.  
<code>-mass</code>&nbsp; Use mass weighted RMSD.  

These options are only supported by the native executable:

<code>-live</code>&nbsp; Show progress in real time.  
<code>-stdin <em>EXT</em></code>&nbsp; Read coordinates from standard input in *EXT* format.  
<code>-stdout <em>EXT</em></code>&nbsp; Write coordinates to standard output in *EXT* format.  
 
The native executable only supports the *xyz* and *mol2* formats, while the python script supports many more.
Note that the format is determined from the file extension when reading from or writing to a file.

Basic usage
-----------

The syntax of the command is:

    molalign [option[s]] file [file]

If only a file is specified, two sets of coordinates will be read from the file. If two files are specified, a single
set of coordinates will be read from each file.

To align the atoms without reordering, run the command without options. To reorder and align the atoms run the
command with the `-sort` option.

Advanced usage
--------------

When reodering is performed, the computation time can be greatly reduced with the `-fast` option. However, the assignment
can fail if the clusters are too different. In such cases use the `-tol` option to increase the tolerance to a value larger
than the maximum expected displacement of the atoms.

By default a threshold of 10 counts is used to stop the computation, but it can be adjusted with the `-count` option.
Lower thresholds reduce the computation time but increase the probability of obtaining suboptimal results.

To avoid too long computations use the `-trials` option to stop the calculation when the number of trials reaches the
specified limit, regardless of the count threshold.

The algorithm always explores multiple possible assignments, but only the best one is recorded by default. Use the `-rec`
option to record and print more than one assignment.

Examples
--------

For small atom displacements the default tolerance is enough:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast
    Optimized RMSD = 0.0506
    
    ./build/molalign examples/Co138_0.xyz examples/Co138_2.xyz -sort -fast
    Optimized RMSD = 0.0977

but for atom displacements larger than the tolerance the assignment will fail:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast
    Error: Assignment failed

Increasing the tolerance will fix the problem but the calculation will slow down:

    ./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -sort -fast -tol 0.7
    Optimized RMSD = 0.1973

Sometimes there is more than one optimal assignment due to symmetry:

    ./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -sort -fast -stats -rec 5
     Map    Count    Steps     Total      Real       RMSD
    -----------------------------------------------------
       1       10     10.0      74.4      54.5     0.0506
       2       10      9.0      69.8      48.5     0.0506
       3       16     11.3      83.2      63.6     0.0506
       4        1      2.0       8.4       7.1     0.6652
       5        1      7.0      30.5       9.8     0.6716
    -----------------------------------------------------
    Random trials = 81
    Minimization steps = 706
    Visited local minima > 5
    Optimized RMSD = 0.0506

Note:
Due to the use of different random seeds, the stats will be different on each run, but the optimized RMSD should be always
the same.

References
----------

<a id="1">[1]</a>
J. M. Vasquez-Perez, L. A. Zarate-Hernandez, C. Z. Gomez-Castro, U. A. Nolasco-Hernandez.
A Practical Algorithm to Solve the Near-Congruence Problem for Rigid Molecules and Clusters,
Journal of Chemical Information and Modeling (2023), DOI: <https://doi.org/10.1021/acs.jcim.2c01187>
