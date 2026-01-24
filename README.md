MolAlignLib
===========

Build all programs
------------------
```
./build.sh
```

*atormsd*
-----------
### Usage
```
atormsd file1 file2 [options]
```

### Options
`-align` Align atoms to minimize the RMSD.  
`-remap` Remap atoms to minimize the RMSD.  
`-prune TOL` Prune assignments exceeding tolerance *TOL*.  
`-thres N` Exit if count threshold *N* is reached.  
`-trials N` Exit if trial limit *N* is reached.  
`-records N` Record the *N* lowest RMSDs (1 by default).  
`-atomorder` Print optimized atom order.  
`-aligned FILE` Write aligned coordinates to *FILE*.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect aligned coordinates.  
`-stats` Print detailed optimization stats.  
`-random` Set random seed from clock.  
`-test` Print output and aligned coordinates to stdout.  

*conformsd*
-----------
### Usage
```
conformsd file1 file2 [options]
```

### Options
`-align` Align atoms to minimize the RMSD.  
`-remap` Remap atoms to minimize the RMSD.  
`-bond` Set bonds from atom distances.  
`-thres N` Exit if count threshold *N* is reached.  
`-trials N` Exit if trial limit *N* is reached.  
`-records N` Record the *N* lowest RMSDs (1 by default).  
`-atomorder` Print optimized atom order.  
`-aligned FILE` Write aligned coordinates to *FILE*.  
`-tree` Print the assignment tree.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect aligned coordinates.  
`-stats` Print detailed optimization stats.  
`-random` Set random seed from clock.  
`-test` Print output and aligned coordinates to stdout.  
