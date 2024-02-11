#!/usr/bin/env python
# coding: utf-8

# Import ASE to read and write files
from ase.io import read, write
from molalignlib import assign_atoms


# Read clusters coordinates
mols = [
    read('../Co138_0.xyz'),
    read('../Co138_1.xyz'),
    read('../Co138_2.xyz'),
]


# Calculate the RMSD between all possible pairs
for mol0 in mols:
    for mol1 in mols:
        assignment = assign_atoms(mol0, mol1, fast=True, tol=0.35)
        rmsd = mol1[assignment.order].align_to(mol0)
        print('{:.4f}'.format(rmsd), end=' ')
    print()

