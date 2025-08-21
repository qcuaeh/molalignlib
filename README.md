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
`-rebond` Set bonds from atom distances.  
`-naive` Probe all feasible assignment combinations.  
`-heavy` Ignore hydrogen atoms.  
`-mass` Use mass weighted coordinates.  
`-mirror` Reflect aligned coordinates.  
`-mapping` Print optimal mapping.  
`-aligned FILE` Write aligned coordinates to *FILE*.  
`-stats` Print detailed stats of the calculation.  
`-tree` Print the assignment tree.  
`-test` Use always the same seed for testing.  
`-n N` Find the *N* lowest RMSDs (1 by default).
