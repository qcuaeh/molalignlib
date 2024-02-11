#!/usr/bin/env python
# coding: utf-8

# Import ASE to read and write files
from ase.io import read, write
from molalignlib import assign_atoms


# Read clusters coordinates
mol0 = read('../Co138_0.xyz')
mol1 = read('../Co138_1.xyz')


# Find the optimal assignment between mol0 and mol1
assignment = assign_atoms(mol0, mol1, fast=True, tol=0.35, stats=True)


# Align mol1 to mol0 for the assignment and write the aligned coordinates to a file
mol1 = mol1[assignment.order]
rmsd = mol1.align_to(mol0)
print('{:.4f}'.format(rmsd))
write('aligned.xyz', mol1)

