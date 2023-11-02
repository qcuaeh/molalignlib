![graphical abstract](abstract.png) 

**MolAlignLib** is a Fortran and Python library based on random rotations and quasi-local
RMSD minimizations to align rigid molecules and clusters. The details of the method can be
found in publication [[1]](#1).

Before installing
-----------------

If you can't install stuff on your computer then
[try&nbsp;MolAlignLib&nbsp;on&nbsp;Binder](https://notebooks.gesis.org/binder/v2/gh/qcuaeh/molalignlib/HEAD?filepath=examples%2Fjupyter).

* To build the native executable you will need a Fortran 2008 compiler and LAPACK.

* To build the python module you will need GFortran 4.9 or higher, LAPACK, Python 3.6 or
higher, *NumPy* and *ASE*.

The easiest way to install the required packages is with your distro's package manager:

#### RHEL, Centos, Fedora, etc.

```
yum install git python3 python3-pip gcc-gfortran lapack-devel
```

#### Debian, Ubuntu, Mint, etc.

```
apt install git python3 python3-pip gfortran liblapack-dev
```

Build the native executable
---------------------------

Clone the repository:

```
git clone https://github.com/qcuaeh/molalignlib.git
```

then enter the cloned directory, edit the *build.env* file to suit your system and run:

```
./build.sh
```

It will create the *molalign* executable inside the *build* directory.

or build and install the python module
--------------------------------------

Simply run:

```
pip3 install git+https://github.com/qcuaeh/molalignlib.git
```

It will install *NumPy*, *ASE* and *MolAlignLib* in your site packages and the *molalign* script in your path.

Program options
---------------

#### Options supported by the native executable and the python script

<code>-reorder</code>&nbsp; Reorder atoms to minimize the RMSD.  
<code>-fast</code>&nbsp; Prune assignments that surpass the displacement tolerance.  
<code>-bond</code>&nbsp; Prioritize assignments that minimize bonding differences.  
<code>-tol <em>TOL</em></code>&nbsp; Set the displacement tolerance to *TOL* (defaults to 0.35 Å).  
<code>-out <em>NAME</em></code>&nbsp; Set the output file name to *NAME* (defaults to aligned.xyz).  
<code>-count <em>N</em></code>&nbsp; Set the count threshold to *N* (defaults to 10).  
<code>-trials <em>N</em></code>&nbsp; Set the maximum number of trials to *N*.  
<code>-rec <em>N</em></code>&nbsp; Record up to *N* assignments (defaults to 1).  
<code>-stats</code>&nbsp; Print detailed stats of the calculation.  
<code>-test</code>&nbsp; Use a fixed random seed for testing.  
<code>-mirror</code>&nbsp; Reflect aligned coordinates.  
<code>-mass</code>&nbsp; Use mass weighted RMSD.  

#### Options only supported by the native executable

<code>-live</code>&nbsp; Print stats in real time (if stats are enabled).  
<code>-stdin <em>FMT</em></code>&nbsp; Read coordinates from standard input with format *FMT*.  
<code>-stdout <em>FMT</em></code>&nbsp; Write coordinates to standard output with format *FMT*.  

Basic usage
-----------

The syntax of the command is:

```
molalign [option[s]] file1 file2
```

The coordinates in file2 will be aligned to the coordinates in file1. If there is
more than one set of coordinates in a file, only the first one will be read. The native
executable only reads *xyz* and *mol2* files while the python script reads all the file
formats supported by *ASE*.

* To align the atoms without reordering, run the command without options.

* To reorder the atoms to minimize the RMSD run the command with the `-reorder` option.

Advanced usage
--------------

When reordering is performed the computation can take a lot of time to complete but
can be considerably speeded up with the `-fast` option which enables pruning of any
assignment that surpass the displacement tolerance. However, if the atom displacements
are larger than this tolerance, the assignment will fail, or, if they are very close,
the assignment can be suboptimal. In such cases the displacement tolerance should be
increased with the `-tol` option.

The count threshold is used to decide if the procedure is converged. A threshold of 10 
counts is used by default, which works fine for almost all cases, but you can change it
with the `-count` option. To avoid too long computations you can set a maximum number of
random trials with the `-trials` option, the computation will be aborted if it is reach
before the count threshold.

The algorithm always explores multiple possible assignments, but only the best one is
recorded by default. To record additional assignments use the `-rec` option. To print
the stats of the computation use the `-trials` option, but notice that they will be
different on each repeated run due to the use of randomized random seeds.

Command line examples
---------------------

For small atom displacements the default tolerance is enough:

```
./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -reorder -fast
```

```
0.0506
```

but for atom displacements larger than the tolerance the assignment will fail:

```
./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -reorder -fast
```

```
Error: Assignment failed
```

Increasing the tolerance will fix the problem but the calculation will slow down:

```
./build/molalign examples/Co138_0.xyz examples/Co138_3.xyz -reorder -fast -tol 0.7
```

```
0.1973
```

Printing multiple alignments and stats can be useful to identify rotational symmetric clusters:

```
./build/molalign examples/Co138_0.xyz examples/Co138_1.xyz -reorder -fast -rec 5 -stats
```

```
 Map    Count    Steps     Angle       RMSD
-------------------------------------------
   1       10     12.1      60.5     0.0506
   2        9     12.3      66.4     0.0506
   3       15     10.1      50.2     0.0506
   4        1      9.0      41.8     0.6652
   5        1      6.0       4.1     0.6716
-------------------------------------------
Random trials = 66
Minimization steps = 595
Visited local minima > 5
```

There are three different assignments with the same RMSD indicating that the cluster
has three fold symmetry. Notice that the three optimal symmetric assignments are visited
multiple times while the suboptimal ones are visited only once.

Python examples
---------------

Minimize the RMSD between two Cobalt clusters:

```
from ase.io import read
from molalignlib import assign_atoms
mol0 = read('Co138_0.xyz')
mol1 = read('Co138_1.xyz')
# Find optimal assignment
assignment = assign_atoms(mol0, mol1, fast=True)
# Reorder mol1 with the optimal assignment
mol1 = mol1[assignment.order]
# Align mol1 to mol0 (returns RMSD)
mol1.align_to(mol0)
```

*MolAlignLib* uses *ASE* to read and write atomic coordinates.

References
----------

<a id="1">[1]</a>
J. M. Vasquez-Perez, L. A. Zarate-Hernandez, C. Z. Gomez-Castro, U. A. Nolasco-Hernandez.
A Practical Algorithm to Solve the Near-Congruence Problem for Rigid Molecules and Clusters,
Journal of Chemical Information and Modeling (2023), DOI: <https://doi.org/10.1021/acs.jcim.2c01187>
