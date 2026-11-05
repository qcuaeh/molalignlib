#!/usr/bin/env python3
## Built from examples.ipynb by nbtopy ##

# %% [markdown]
# # MolAlignLib examples
#
# These examples use the Python bindings of MolAlignLib and two sets of data files in this folder:
#
# * `clusters/`: atom clusters of Ar₂₄Cl₂ (no bonds), one per file in `Ar24Cl2_1.xyz` … `Ar24Cl2_3.xyz`, and all together in the multi-frame file `Ar24Cl2_all.xyz`.
# * `conformers/`: conformers of C₈₄H₁₁₄Cl₂N₈O₁₇SSi₂ (with bonds), one per file in `CID168301727_1.sdf` … `CID168301727_3.sdf`, and all together in the multi-frame file `CID168301727_all.sdf`.
#
# Comparisons of two structures read them from the single-structure files with `Molecule.from_file`; the RMSD matrices read every frame of the `_all` file with `read_frames`.
#
# There are two ways to compare two structures, one per algorithm:
#
# * `mol0.ato_rmsd(mol1, ...)` compares them as **unstructured atom clusters**: bonds are ignored.
# * `mol0.confo_rmsd(mol1, ...)` compares them as **conformers**: atom assignments always respect the bond graph.
#
# In both, `mol0` is the reference: it stays fixed, and `mol1` is aligned and reordered onto it.

# %%
import numpy as np

from molalignlib import (
    Molecule, read_frames, bond_type_code, bond_type_label, BOND_TOL_DEFAULT,
)

# %% [markdown]
# ## 1. Symmetry-corrected RMSD of two conformers
#
# The most common use: read two conformers of a molecule and compare them with `confo_rmsd`. `Molecule.from_file` reads one structure from a file (the first frame, or the frame given by `frame_idx`), and names the molecule after the file. `align=True` superposes them, and `remap=True` finds the atom order that gives the lowest RMSD within the symmetry of the bond graph.

# %%
# Example 1 - RMSD between two conformers

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')
print(mol0)
print(mol1)

result = mol0.confo_rmsd(mol1, align=True, remap=True)
print(f'RMSD = {result.rmsd:.4f} Å')

# %% [markdown]
# ## 2. RMSD of two atom clusters
#
# For metal clusters, nanoparticles and other systems where bonds are absent or irrelevant, `ato_rmsd` matches atoms by element and position only.

# %%
# Example 2 - RMSD between two atom clusters

mol0 = Molecule.from_file('clusters/Ar24Cl2_1.xyz')
mol1 = Molecule.from_file('clusters/Ar24Cl2_2.xyz')

result = mol0.ato_rmsd(mol1, align=True, remap=True)
print(f'RMSD = {result.rmsd:.4f} Å')

# %% [markdown]
# ## 3. What `align` and `remap` do
#
# * Neither: the RMSD of the coordinates as they are, in the input atom order. Without `remap`, both molecules must list the same atoms in the same order (and, for conformers, the same bonds); otherwise a `ValueError` is raised.
# * `align=True`: `mol1` is rotated and translated onto `mol0`, keeping the atom order.
# * `remap=True`: the atom order is optimised (within the symmetry of the bond graph), keeping the coordinates.
# * Both: the symmetry-corrected RMSD, usually what you want.

# %%
# Example 3 - RMSD with and without alignment and remapping

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')

for align in (False, True):
    for remap in (False, True):
        try:
            result = mol0.confo_rmsd(mol1, align=align, remap=remap)
            print(f'align={align!s:5}  remap={remap!s:5}  RMSD = {result.rmsd:.4f}')
        except ValueError as error:   # without remap, the atom orders must agree
            print(f'align={align!s:5}  remap={remap!s:5}  ValueError: {error}')

# %% [markdown]
# ## 4. The best solution in detail, and the aligned structure
#
# An `RMSDResult` holds the RMSD, the atom permutation (`mapping`, 0-based: entry `j` is the atom of `mol1` placed on line `j` of `mol0`) and the 4×4 matrix that moves `mol1` onto `mol0`. `apply_to` applies both and returns a new `Molecule` whose bonds are renumbered, so it can be written with connectivity.

# %%
# Example 4 - Inspect the best solution and save the aligned structure

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')
result = mol0.confo_rmsd(mol1, align=True, remap=True)

print(result)
print('Mapping:', result.mapping)
print('Transform:')
print(np.round(result.transform, 4))

# Line j of `aligned` is atom result.mapping[j] of mol1, moved into the frame of mol0
aligned = result.apply_to(mol1)
aligned.write('CID168301727_2_aligned.sdf', comment=f'rmsd={result.rmsd:.6f}')

# The RMSD can be recomputed from the aligned coordinates
rmsd = np.sqrt(np.mean(np.sum((aligned.coords - mol0.coords)**2, axis=1)))
print(f'Recomputed RMSD = {rmsd:.6f}')

# %% [markdown]
# ## 5. RMSD matrix over all conformer pairs
#
# `read_frames` without `frames` returns every frame of a multi-frame file as a list, named `<file stem>_<frame index>` (0-based). With `frames=(i, j, ...)` it returns only those frames, as a tuple: `read_frames('conformers/CID168301727_all.sdf', frames=(0, 1))` gives the same two conformers as the `_1` and `_2` files.

# %%
# Example 5 - RMSD matrix over all conformer pairs

conformers = read_frames('conformers/CID168301727_all.sdf')   # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol0.confo_rmsd(mol1, remap=True, align=True)
        print(f'{result.rmsd:.4f}', end=2*' ')
    print()

# %% [markdown]
# ## 6. RMSD matrix with NumPy, and removing duplicate conformers
#
# With `align=True, remap=True` the RMSD is symmetric, so only half of the pairs need to be computed. A common use of the matrix is to keep one representative of each group of conformers closer than a threshold.

# %%
# Example 6 - Symmetric RMSD matrix and duplicate removal

conformers = read_frames('conformers/CID168301727_all.sdf')
n = len(conformers)

matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        rmsd = conformers[i].confo_rmsd(conformers[j], align=True, remap=True).rmsd
        matrix[i, j] = matrix[j, i] = rmsd

np.set_printoptions(precision=3, suppress=True)
print(matrix)

threshold = 0.5   # Å
unique = []
for i in range(n):
    if all(matrix[i, k] > threshold for k in unique):
        unique.append(i)
print(f'Kept {len(unique)} of {n} conformers, all more than {threshold} Å apart:',
      [conformers[k].name for k in unique])

# %% [markdown]
# ## 7. Heavy atoms only, and mass weighting
#
# `heavy_only=True` leaves hydrogens out of the RMSD (they are still paired afterwards, so the mapping covers every atom). `mass_weight=True` weights each atom by its atomic mass.

# %%
# Example 7 - Hydrogens and mass weighting

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')

for heavy_only in (False, True):
    for mass_weight in (False, True):
        result = mol0.confo_rmsd(mol1, align=True, remap=True,
                                 heavy_only=heavy_only, mass_weight=mass_weight)
        print(f'heavy_only={heavy_only!s:5}  mass_weight={mass_weight!s:5}  RMSD = {result.rmsd:.4f}')

# %% [markdown]
# ## 8. Several ranked solutions for atom clusters
#
# With `align=True, remap=True`, `max_recs=N` returns a list of up to *N* distinct solutions, best first. Pruning (`prune=True`) skips pairings between atoms with incompatible environments; it needs a tolerance `prune_tol` (Å).

# %%
# Example 8 - Retrieve multiple mappings, comparing as unstructured clusters

mol0 = Molecule.from_file('clusters/Ar24Cl2_1.xyz')
mol1 = Molecule.from_file('clusters/Ar24Cl2_2.xyz')

# Ask for the top 5 lowest RMSD mappings
results = mol0.ato_rmsd(mol1, remap=True, align=True, prune=True, prune_tol=0.1, max_recs=5)

for i, result in enumerate(results, start=1):
    print(f'Mapping {i} RMSD = {result.rmsd:.4f}')

# The list is sorted best-first
best = results[0]

# %% [markdown]
# ## 9. Molecules built in memory
#
# `Molecule.from_symbols` (or `Molecule.from_numbers`) builds a molecule from element symbols and coordinates in Å. Without `bond_data` it has no bonds, so `confo_rmsd` needs `rebond=True` to perceive them from the geometry: two atoms are bonded when closer than the sum of their covalent radii plus `bond_tol` (default `BOND_TOL_DEFAULT`).
#
# Below, a copy of ethanol with its atoms shuffled, rotated and translated is compared with the original: the RMSD is zero and `mapping` undoes the shuffle.

# %%
# Example 9 - Recover a known permutation

symbols = ['C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
coords = np.array([
    [ 1.1879, -0.3829,  0.0000],   # C (methyl)
    [ 0.0000,  0.5526,  0.0000],   # C (methylene)
    [-1.1867, -0.2472,  0.0000],   # O
    [-1.9237,  0.3850,  0.0000],   # H (hydroxyl)
    [ 2.0985,  0.2306,  0.0000],   # H (methyl)
    [ 1.1184, -1.0093,  0.8869],   # H (methyl)
    [ 1.1184, -1.0093, -0.8869],   # H (methyl)
    [ 0.0227,  1.1812,  0.8852],   # H (methylene)
    [ 0.0227,  1.1812, -0.8852],   # H (methylene)
])
ethanol = Molecule.from_symbols(symbols, coords, name='ethanol')
print(ethanol, '- has bonds:', ethanol.has_bonds)

# Atom i of the copy is atom perm[i] of ethanol
rng = np.random.default_rng(seed=1)
perm = rng.permutation(len(ethanol))
q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
rotation = q * np.sign(np.linalg.det(q))     # a proper rotation (determinant +1)
shuffled = Molecule.from_symbols([symbols[i] for i in perm],
                                 coords[perm] @ rotation.T + [1.0, -2.0, 0.5],
                                 name='ethanol_shuffled')

print('Default bond tolerance:', BOND_TOL_DEFAULT, 'Å')
result = ethanol.confo_rmsd(shuffled, align=True, remap=True, rebond=True)
print(f'RMSD = {result.rmsd:.2e}')
print('Permutation recovered:', np.array_equal(perm[result.mapping], np.arange(len(ethanol))))

# As clusters (no bonds needed) the answer is the same
result = ethanol.ato_rmsd(shuffled, align=True, remap=True)
print(f'ato_rmsd RMSD = {result.rmsd:.2e}')

# %% [markdown]
# Errors reported by the library raise `ValueError`. Without bonds every atom is a fragment of its own, which `confo_rmsd` rejects (see example 14):

# %%
# Example 9b - confo_rmsd needs bonds

try:
    ethanol.confo_rmsd(shuffled, align=True, remap=True)
except ValueError as error:
    print('ValueError:', error)

# %% [markdown]
# ## 10. Mirror images
#
# Rotations cannot superpose the two enantiomers of a chiral molecule. `mirror=True` reflects `mol1` (x → −x) before comparing, so the comparison is made with its mirror image. Here the mirror image of `mol1` is built by hand; reflecting it back with `mirror=True` gives the original RMSD, and the transform then includes the reflection (determinant −1).

# %%
# Example 10 - Compare with the mirror image

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')

mirror_image = Molecule(atom_data=mol1.atom_data, coords=mol1.coords * [-1, 1, 1],
                        bond_data=mol1.bond_data, symbols=mol1.symbols,
                        name=mol1.name + '_mirror')

plain = mol0.confo_rmsd(mol1, align=True, remap=True)
reflected = mol0.confo_rmsd(mirror_image, align=True, remap=True)
reflected_back = mol0.confo_rmsd(mirror_image, align=True, remap=True, mirror=True)

print(f'mol1                              RMSD = {plain.rmsd:.4f}')
print(f'mirror image (equal if achiral)   RMSD = {reflected.rmsd:.4f}')
print(f'mirror image with mirror=True     RMSD = {reflected_back.rmsd:.4f}')
print('det of the 3x3 block:', round(np.linalg.det(reflected_back.transform[:3, :3]), 6))

# %% [markdown]
# ## 11. Different numbers of hydrogens
#
# Molecules with different atoms are not isomers and cannot be compared. With `heavy_only=True` only the heavy atoms need to match: the smaller molecule is padded with padding atoms, so the mapping has `max(len(mol0), len(mol1))` entries, and `apply_to` puts a padding atom `X` on the line of each hydrogen without a partner.

# %%
# Example 11 - Ethoxide (ethanol without its hydroxyl hydrogen, atom 3)

ethoxide = Molecule.from_symbols(symbols[:3] + symbols[4:], np.delete(coords, 3, axis=0),
                                 name='ethoxide')

try:
    ethanol.confo_rmsd(ethoxide, align=True, remap=True, rebond=True)
except ValueError as error:
    print('ValueError:', error)

result = ethanol.confo_rmsd(ethoxide, align=True, remap=True, rebond=True, heavy_only=True)
print(f'Heavy-atom RMSD = {result.rmsd:.2e}')
print('Mapping:', result.mapping, f'(index {len(ethoxide)} is a padding atom)')

matched = result.apply_to(ethoxide)
print('Symbols of ethanol:        ', ethanol.symbols)
print('Symbols of matched ethoxide:', matched.symbols)

# %% [markdown]
# ## 12. Bond types
#
# Bond types are short labels (`'1'`, `'ar'`, `'3/2'`, ...) stored as integer codes in the third column of `bond_data`. They are compared literally, and only with `use_bond_type=True`; otherwise only connectivity matters. Bonds read from files are labelled from their bond orders (`'1'`–`'5'`, `'ar'`, `'am'`, `'dv'`).

# %%
# Example 12 - Bond type labels and codes

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')

for label in ('1', '2', 'ar', '3/2', ''):
    code = bond_type_code(label)
    print(f'{label!r:6} -> {code:5} -> {bond_type_label(code)!r}')

print('Bond types in', mol0.name, ':', sorted({bond_type_label(t) for t in mol0.bond_data[:, 2]}))
result = mol0.confo_rmsd(mol1, align=True, remap=True, use_bond_type=True)
print(f'RMSD with bond types = {result.rmsd:.4f}')

# %%
# Example 12b - Formaldehyde with explicit bonds (1-based atom indices)

single, double = bond_type_code('1'), bond_type_code('2')
hcho_symbols = ['C', 'O', 'H', 'H']
hcho_coords = np.array([[0.000, 0.000, 0.0], [0.000, 1.208, 0.0],
                        [0.935, -0.541, 0.0], [-0.935, -0.541, 0.0]])
hcho = Molecule.from_symbols(hcho_symbols, hcho_coords,
                             bond_data=[[1, 2, double], [1, 3, single], [1, 4, single]])

# Same connectivity, but the C=O bond labelled as single
relabelled = Molecule.from_symbols(hcho_symbols, hcho_coords,
                                   bond_data=[[1, 2, single], [1, 3, single], [1, 4, single]])

print('Connectivity only: RMSD =', hcho.confo_rmsd(relabelled, align=True, remap=True).rmsd)
try:
    hcho.confo_rmsd(relabelled, align=True, remap=True, use_bond_type=True)
except ValueError as error:
    print('With bond types: ValueError:', error)

# %% [markdown]
# ## 13. Atom labels
#
# Every atom carries an integer label (0 by default). With `use_atom_type=True`, atoms are only matched with atoms of the same element *and* label. Here four argon atoms sit on the corners of a 2 Å × 1 Å rectangle: the same coordinates, labelled by rows in `a` and by columns in `b`.

# %%
# Example 13 - Restrict matching with labels

rectangle = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
a = Molecule.from_symbols(['Ar'] * 4, rectangle, labels=[1, 1, 2, 2], name='rows')
b = Molecule.from_symbols(['Ar'] * 4, rectangle, labels=[1, 2, 2, 1], name='columns')
print('atom_data of a (atomic number, label):', a.atom_data.tolist())

for use_atom_type in (False, True):
    result = a.ato_rmsd(b, align=True, remap=True, use_atom_type=use_atom_type)
    print(f'use_atom_type={use_atom_type!s:5}  RMSD = {result.rmsd:.4f}  mapping = {result.mapping}')

# %% [markdown]
# ## 14. Complexes with several fragments
#
# A fragment is a connected piece of the bond graph. By default each molecule must be a single fragment; raise `max_frags` to compare complexes, salts or clusters of molecules. Here an ethanol dimer is compared with a rotated copy that lists the two molecules in the opposite order.

# %%
# Example 14 - An ethanol dimer

shift = [0.0, 0.0, 4.0]
dimer = Molecule.from_symbols(symbols * 2, np.vstack([coords, coords + shift]), name='dimer')
swapped = Molecule.from_symbols(symbols * 2, np.vstack([coords + shift, coords]) @ rotation.T,
                                name='dimer_swapped')

try:
    dimer.confo_rmsd(swapped, align=True, remap=True, rebond=True)
except ValueError as error:
    print('ValueError:', error)

result = dimer.confo_rmsd(swapped, align=True, remap=True, rebond=True, max_frags=2)
print(f'RMSD = {result.rmsd:.2e}')
print('Mapping:', result.mapping)   # the first molecule of dimer is the second of swapped

# %% [markdown]
# ## 15. Search settings and diagnostics
#
# With `align=True, remap=True` the library samples random orientations of `mol1`. The search stops when the best solution has been found `min_freq` times or after `max_trials` orientations. Results are reproducible unless `random=True` seeds the generator from the clock. For conformers, keep the default `min_freq=100`, the value validated on large benchmarks.

# %%
# Example 15 - Search settings for clusters

mol0 = Molecule.from_file('clusters/Ar24Cl2_1.xyz')
mol1 = Molecule.from_file('clusters/Ar24Cl2_2.xyz')

for settings in ({}, {'min_freq': 3}, {'max_trials': 20}, {'random': True}):
    result = mol0.ato_rmsd(mol1, align=True, remap=True, **settings)
    print(f'{str(settings):20} RMSD = {result.rmsd:.4f}')

# %% [markdown]
# `print_stats=True` prints the ranked local minima and search statistics, and `print_assignment_tree=True` (conformers only) prints the assignment tree of the HNA partition. The library prints them from compiled code, so depending on your Jupyter setup they may show up in the terminal that runs Jupyter rather than below the cell.

# %%
# Example 15b - Diagnostics

mol0 = Molecule.from_file('conformers/CID168301727_1.sdf')
mol1 = Molecule.from_file('conformers/CID168301727_2.sdf')
result = mol0.confo_rmsd(mol1, align=True, remap=True, print_stats=True, print_assignment_tree=True)
