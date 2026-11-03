#!/usr/bin/env python3
## Built from examples.ipynb by nbtopy ##

# %%
from molalignlib import read_molecules

# %%
# Example 1 - Retrieve multiple mappings, comparing as unstructured clusters

mol0, mol1 = read_molecules('clusters.xyz', frames=(0,1))

# Ask for the top 5 lowest RMSD mappings
results = mol0.ato_rmsd(mol1, remap=True, align=True, prune=True, prune_tol=0.1, max_recs=5)

for i, result in enumerate(results, start=1):
    print(f'Mapping {i} RMSD = {result.rmsd:.4f}')

# The list is sorted best-first
best = results[0]

# %%
# Example 2 - RMSD matrix over all conformer pairs

conformers = read_molecules('conformers.sdf')   # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol0.confo_rmsd(mol1, remap=True, align=True)
        print(f'{result.rmsd:.4f}', end=2*' ')
    print()
