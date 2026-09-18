#!/usr/bin/env python3
## Built from examples.ipynb by nbtopy ##

# %%
from molalignlib import read_clusters, read_conformers

# %%
# Example 1 - Retrieve multiple mappings

mol0, mol1 = read_clusters('clusters.xyz', frames=(0,1))

# Ask for the top 5 lowest RMSD mappings
results = mol0.rmsd_to(mol1, remap=True, align=True, prune=True, prune_tol=0.1, n_records=5)

for i, result in enumerate(results, start=1):
    print(f'Mapping {i} RMSD = {result.rmsd:.4f}')

# The list is sorted best-first
best = results[0]

# %%
# Example 2a - RMSD matrix over all conformer pairs

conformers = read_conformers('conformers.sdf')   # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol1.rmsd_to(mol0, remap=True, align=True)[0]
        print(f'{result.rmsd:.4f}', end=2*' ')
    print()

# %%
# Example 2b - RMSD matrix over all conformer pairs, deriving bond connectivity
# from geometry

conformers = read_conformers('conformers.xyz')   # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol1.rmsd_to(mol0, remap=True, align=True, infer_bonds=True)[0]
        print(f'{result.rmsd:.4f}', end=2*' ')
    print()
