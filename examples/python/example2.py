#!/usr/bin/env python
# coding: utf-8

# Import ASE to read and write files
from ase.io import read, write
from molalignlib import assign_atoms


# Read clusters coordinates
mol0 = read('../Co138_0.xyz')
mol1 = read('../Co138_1.xyz')


# Find the 5 best assignments between mol0 and mol1
assignments = assign_atoms(mol0, mol1, fast=True, tol=0.35, rec=5, stats=True)


# Align mol1 to mol0 for each assignment and write the aligned coordinates to a file
for assignment in assignments:
    auxmol = mol1[assignment.order]
    rmsd = auxmol.align_to(mol0)
    print('{:.4f}'.format(rmsd))
    write('aligned.xyz', auxmol, append=True)

