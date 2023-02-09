#!/usr/bin/env python
# coding: utf-8

# ASE must be imported before MolAlignLib in order to patch ASE's Atoms class
# with the align_to method
from ase.io import read, write
from molalignlib import assign_atoms


# Read clusters coordinates
atoms0 = read('Co138_0.xyz', index=0)
atoms1 = read('Co138_1.xyz', index=0)


# Find the 5 best assignments between atoms0 and atoms1 using biasing and iteration
assignments = assign_atoms(atoms0, atoms1, biasing=True, iteration=True, stats=True, records=5)


# Align atoms1 to atoms0 for each assignment and write the aligned coordinates to a file
write('aligned.xyz', atoms0, comment='Reference')
for i in assignments:
    atoms2 = atoms1[i.order]
    rmsd = atoms2.align_to(atoms0)
    write('aligned.xyz', atoms2, append=True, comment='RMSD {:.4f}'.format(rmsd))

