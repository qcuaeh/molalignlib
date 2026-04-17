MolAlignLib
===========

Building
--------
To build the executables, GFortran 4.8+ (or any other Fortran 2008 compiler) is required:
```
mkdir build && cd build
cmake ..
make
```

*atormsd*
-----------
### Usage
```
atormsd file1 file2 [options]
```

### Options
`-align` Align atoms to minimize the RMSD.  
`-remap` Remap atoms to minimize the RMSD.  
`-label` Use atom labels to distinguish atom types.  
`-near` Use nearest-neighbor assignment without pruning.  
`-prune TOL` Prune assignments exceeding tolerance *TOL*.  
`-freq N` Exit if the count frequency *N* is reached.  
`-trials N` Exit if trial limit *N* is reached.  
`-records N` Record the *N* lowest RMSDs (1 by default).  
`-printmap` Print optimized atom order.  
`-aligned FILE` Write aligned coordinates to *FILE*.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect aligned coordinates.  
`-stats` Print detailed optimization stats.  
`-random` Set random seed from clock.  

*conformsd*
-----------
### Usage
```
conformsd file1 file2 [options]
```

### Options
`-align` Align atoms to minimize the RMSD.  
`-remap` Remap atoms to minimize the RMSD.  
`-bond` Set bonds from atom distances.  
`-label` Use atom labels to distinguish atom types.  
`-exhaustive` Use exhaustive search for atom remapping.  
`-stochastic` Use non-adaptive stochastic search for atom remapping.  
`-freq N` Exit if the count frequency *N* is reached.  
`-trials N` Exit if trial limit *N* is reached.  
`-records N` Record the *N* lowest RMSDs (1 by default).  
`-aligned FILE` Write aligned coordinates to *FILE*.  
`-printmap` Print optimized atom order.  
`-printtree` Print the assignment tree.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect aligned coordinates.  
`-stats` Print detailed optimization stats.  
`-random` Set random seed from clock.  
