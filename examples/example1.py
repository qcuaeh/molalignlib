#!/usr/bin/env python
# coding: utf-8

# Import ASE to read and write files
from ase.io import read, write
from molalignlib import assign_atoms


# Read clusters coordinates
atoms0 = read('Co138_0.xyz')
atoms1 = read('Co138_1.xyz')


# Find the 5 best assignments between atoms0 and atoms1
assignments = assign_atoms(atoms0, atoms1, fast=True, tol=0.35, rec=5, stats=True)


# Align atoms1 to atoms0 for each assignment and write the aligned coordinates to a file
for assignment in assignments:
    atoms2 = atoms1[assignment.order]
    rmsd = atoms2.align_to(atoms0)
    write('aligned.xyz', atoms2, append=True, comment='RMSD={:.4f}'.format(rmsd))

