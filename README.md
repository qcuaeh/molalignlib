MolAlignLib
===========

Build all programs
------------------
```
./build.sh
```

ConfRMSD program
----------------
### Usage
```
ConfRMSD file1 file2 [options]
```

### Options
`-align` Align atoms to minimize the RMSD.  
`-remap` Remap atoms to minimize the RMSD.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect the aligned coordinates.  
`-tree` Print the assignment tree.  
`-stats` Print detailed stats of the calculation.  
`-test` Use always the same seed for testing.  
`-o FILE` Write the aligned coordinates to *FILE*.  
`-n N` Find the *N* lowest RMSDs (1 by default).
