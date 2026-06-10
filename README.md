MolAlignLib
===========

A high-performance library to compute optimal RMSD between atom clusters and symmetry-corrected RMSD between molecular conformers.

MolAlignLib uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to achieve exact, topologically valid atom assignments between conformers in milliseconds, even for highly symmetric molecules that are intractable by conventional graph-isomorphism approaches.


Table of Contents
-----------------
1. [Overview](#overview)
2. [Building from Source](#building-from-source)
3. [Command-line Programs](#command-line-programs)
   - [atormsd](#atormsd)
   - [conformsd](#conformsd)
4. [C Binding](#c-binding)
   - [atormsd_calculate](#atormsd_calculate)
   - [conformsd_calculate](#conformsd_calculate)
5. [Python API](#python-api)
   - [Installation](#installation)
   - [AtomCluster](#atomcluster)
   - [Conformer](#conformer)
   - [RMSDResult](#rmsdresult)
   - [Examples](#examples)
6. [Algorithm Notes](#algorithm-notes)


Overview
--------

MolAlignLib exposes two distinct RMSD calculation modes:

| Mode | Program / function | When to use |
|------|-------------------|-------------|
| **atoRMSD** | `atormsd` / `atormsd_calculate` | Unstructured atom clusters — no bond topology required |
| **confoRMSD** | `conformsd` / `conformsd_calculate` | Molecular conformers — bond topology guides atom matching via HNA partitioning |

Both modes support:
- **Alignment** (`-align`): optimally rotate and translate one structure onto the other.
- **Remapping** (`-remap`): find the atom permutation that minimises the RMSD.
- **Mass weighting** (`-mass`): weight each atom by its atomic mass.
- **Heavy-atom only** (`-heavy`): exclude hydrogen atoms.
- **Mirror images** (`-mirror`): reflect the second structure before comparison.


Building from Source
--------------------

### Requirements

- CMake ≥ 3.15
- GFortran ≥ 7.0 (or Intel Fortran; any Fortran 2008-compliant compiler)
- A C compiler (for the C binding)

### Steps

```bash
git clone -b devel https://github.com/qcuaeh/molalignlib.git
cd molalignlib
mkdir build && cd build
cmake ..
make
```

After a successful build the following are available inside `build/fortran`:

| Artifact | Description |
|----------|-------------|
| `atormsd` | Standalone cluster RMSD program |
| `conformsd` | Standalone conformer RMSD program |
| `libmolalignlib.a` | C binding static library |

### Python extension modules

See the [Python API](#python-api) section below.


Command-line Programs
---------------------

### atormsd

Calculate RMSD between two unstructured atom clusters. No bond information is needed; atom matching is guided by element type and spatial proximity alone.

```
atormsd file1 file2 [options]
```

#### Options

| Option | Argument | Description |
|--------|----------|-------------|
| `-align` | | Align atoms to minimise the RMSD |
| `-remap` | | Remap atoms to minimise the RMSD |
| `-label` | | Use atom labels to distinguish atom types |
| `-prune` | `TOL` | Prune assignments with pair distances exceeding *TOL* Å |
| `-freq` | `N` | Stop if the best solution is found *N* consecutive times |
| `-trials` | `N` | Stop after at most *N* optimisation trials |
| `-records` | `N` | Record the *N* lowest RMSDs found (default: 1) |
| `-printmap` | | Print the optimised atom mapping to stdout |
| `-aligned` | `FILE` | Write aligned coordinates of molecule 2 to *FILE* |
| `-heavy` | | Ignore hydrogen atoms |
| `-mass` | | Use mass-weighted coordinates |
| `-mirror` | | Reflect molecule 2 before comparison |
| `-stats` | | Print detailed optimisation statistics |
| `-random` | | Seed the random-number generator from the system clock |

#### Examples

```bash
# Plain RMSD (no alignment, no remapping)
atormsd mol1.xyz mol2.xyz

# Align and remap atoms, print the mapping
atormsd mol1.xyz mol2.xyz -align -remap -printmap

# Remap with distance pruning, keep the 3 best solutions
atormsd mol1.xyz mol2.xyz -align -remap -prune 0.5 -records 3

# Heavy-atom, mass-weighted RMSD with alignment
atormsd mol1.sdf mol2.sdf -align -remap -heavy -mass

# Write the aligned molecule 2 to a new file
atormsd mol1.xyz mol2.xyz -align -remap -aligned mol2_aligned.xyz
```


### conformsd

Calculate the symmetry-corrected RMSD between two molecular conformers. Bond topology is used to build the HNA partition, which guides atom matching and guarantees chemically valid assignments even for highly symmetric molecules. For alignment calculations, the algorithm automatically selects between stochastic orientation sampling and direct enumeration based on the structure of the assignment tree, avoiding unnecessary computation.

```
conformsd file1 file2 [options]
```

**Supported file formats:** XYZ, SDF, and Mol2.

#### Options

| Option | Argument | Description |
|--------|----------|-------------|
| `-align` | | Align atoms to minimise the RMSD |
| `-remap` | | Remap atoms to minimise the RMSD |
| `-label` | | Use atom labels to distinguish atom types |
| `-freq` | `N` | Stop if the best solution is found *N* consecutive times |
| `-trials` | `N` | Stop after at most *N* optimisation trials |
| `-records` | `N` | Record the *N* lowest RMSDs found (default: 1) |
| `-printmap` | | Print the optimised atom mapping to stdout |
| `-printtree` | | Print the internal assignment tree |
| `-aligned` | `FILE` | Write aligned coordinates of molecule 2 to *FILE* |
| `-heavy` | | Ignore hydrogen atoms |
| `-mass` | | Use mass-weighted coordinates |
| `-mirror` | | Reflect molecule 2 before comparison |
| `-bond` | | Derive bond connectivity from interatomic distances instead of the file's bond table |
| `-exhaustive` | | Force exhaustive orientation-independent search regardless of assignment tree topology |
| `-stochastic` | | Force stochastic fixed-orientation search regardless of assignment tree topology |
| `-stats` | | Print detailed optimisation statistics |
| `-random` | | Seed the random-number generator from the system clock |

`-exhaustive` and `-stochastic` are mutually exclusive. If neither is given, the strategy is chosen automatically based on the ratio of total to partial assignment combinations in the tree: stochastic fixed-orientation search is used when the ratio is high, and exhaustive orientation-independent search when it is low.

#### Examples

```bash
# Symmetry-corrected RMSD with remapping and alignment
conformsd conf1.sdf conf2.sdf -align -remap

# Derive connectivity from geometry (useful for XYZ input)
conformsd conf1.xyz conf2.xyz -align -remap -bond

# Heavy atoms only, exhaustive search
conformsd conf1.sdf conf2.sdf -align -remap -heavy -exhaustive

# Print the atom permutation that maps conf2 onto conf1
conformsd conf1.sdf conf2.sdf -align -remap -printmap

# Compute an RMSD matrix for all poses in an SDF file (shell loop)
for i in 1 2 3; do
  for j in 1 2 3; do
    conformsd pose${i}.sdf pose${j}.sdf -align -remap
  done
done
```

C Binding
---------

Compile the static library, include the appropriate header and link against `libmolalignlib`.

```c
#include "conformsd.h" /* for conformsd_calculate */
#include "atormsd.h"   /* for atormsd_calculate  */
```

### Atom data layout

Both functions receive atom information as a flat `int` array of length
`n_atoms * 2`, packed as:

```
[ elnum_0, label_0, elnum_1, label_1, ... ]
```

- `elnum` — atomic number (e.g. 6 for carbon, 8 for oxygen).
- `label` — user-defined integer label; pass `0` for unlabelled atoms.

Coordinates are passed as a flat `double` array of length `n_atoms * 3`,
in row-major (C) order: `[x_0, y_0, z_0, x_1, y_1, z_1, ...]`.

### Transform output

Both functions write a row-major 4 × 4 homogeneous transformation matrix
(16 `double` values) to the `transform` output parameter. The matrix maps
molecule-2 coordinates into the molecule-1 reference frame:

```
p_out = R * p_in + t
```

The matrix is the identity when `align_flag = false`.

### atormsd_calculate

```c
void atormsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag,
    bool stats_flag, bool random_flag,
    double prune_tol, int conv_freq, int max_trials,
    double *rmsd, int *natoms, int *atomperm,
    double *transform, int *error_code);
```

#### Parameters

| Parameter | Direction | Description |
|-----------|-----------|-------------|
| `n_atoms1` | in | Number of atoms in molecule 1 |
| `atom_data1` | in | Packed atom data for molecule 1, length `n_atoms1*2` |
| `coords1` | in | Coordinates for molecule 1, length `n_atoms1*3` |
| `n_atoms2` | in | Number of atoms in molecule 2 |
| `atom_data2` | in | Packed atom data for molecule 2, length `n_atoms2*2` |
| `coords2` | in | Coordinates for molecule 2, length `n_atoms2*3` |
| `align_flag` | in | Enable structural alignment |
| `remap_flag` | in | Enable atom remapping |
| `heavy_flag` | in | Use only heavy (non-hydrogen) atoms |
| `mass_flag` | in | Weight atoms by atomic mass |
| `mirror_flag` | in | Mirror molecule 2 before comparison |
| `label_flag` | in | Use atom labels for type matching |
| `stats_flag` | in | Print optimisation statistics to stdout |
| `random_flag` | in | Seed RNG from system clock |
| `prune_tol` | in | Pruning distance tolerance (Å); negative value disables pruning |
| `conv_freq` | in | Convergence frequency threshold |
| `max_trials` | in | Maximum number of optimisation trials |
| `rmsd` | out | Calculated RMSD (Å) |
| `natoms` | out | Number of elements written to `atomperm` |
| `atomperm` | out | Atom permutation, **0-based**; caller must allocate ≥ `natoms` elements |
| `transform` | out | 4 × 4 homogeneous transform (row-major, 16 doubles) |
| `error_code` | out | `0` = success; `1` = not isomers; `2` = atom type mismatch |

#### Minimal example

```c
#include <stdio.h>
#include <stdlib.h>
#include "atormsd.h"

int main(void)
{
    /* Water molecule 1 */
    int    ad1[] = { 8,0, 1,0, 1,0 };       /* O, H, H */
    double xy1[] = { 0.000, 0.000, 0.000,
                     0.757, 0.586, 0.000,
                    -0.757, 0.586, 0.000 };

    /* Water molecule 2 (slightly displaced) */
    int    ad2[] = { 8,0, 1,0, 1,0 };
    double xy2[] = { 0.010, 0.010, 0.000,
                     0.767, 0.596, 0.000,
                    -0.747, 0.596, 0.000 };

    double rmsd, transform[16];
    int    atomperm[3], natoms, error_code;

    atormsd_calculate(
        3, ad1, xy1,
        3, ad2, xy2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false,
        /*stats=*/false, /*random=*/false,
        /*prune_tol=*/-1.0, /*conv_freq=*/10, /*max_trials=*/10000,
        &rmsd, &natoms, atomperm, transform, &error_code);

    if (error_code != 0) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd);
    return 0;
}
```

Compile:

```bash
gcc example.c -o example -lmolalignlib -lgfortran -lm
```


### conformsd_calculate

```c
void conformsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_bonds1, const int *bond_data1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    int n_bonds2, const int *bond_data2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag, bool bond_flag,
    bool stats_flag, bool random_flag,
    int conv_freq, int max_trials,
    double *rmsd, int *natoms, int *atomperm,
    double *transform, int *error_code);
```

Bond data is passed as a flat `int` array of length `n_bonds * 3`, packed as:

```
[ atom1_0, atom2_0, type_0, atom1_1, atom2_1, type_1, ... ]
```

Atom indices are **1-based**. When `bond_flag = true` the library derives
connectivity from atomic geometry and the bond arrays may be empty
(`n_bonds = 0`, `bond_data = NULL`).

#### Parameters

| Parameter | Direction | Description |
|-----------|-----------|-------------|
| `n_atoms1` | in | Number of atoms in molecule 1 |
| `atom_data1` | in | Packed atom data for molecule 1, length `n_atoms1*2` |
| `coords1` | in | Coordinates for molecule 1, length `n_atoms1*3` |
| `n_bonds1` | in | Number of bonds in molecule 1 |
| `bond_data1` | in | Flat bond array for molecule 1, length `n_bonds1*3` |
| `n_atoms2` | in | Number of atoms in molecule 2 |
| `atom_data2` | in | Packed atom data for molecule 2, length `n_atoms2*2` |
| `coords2` | in | Coordinates for molecule 2, length `n_atoms2*3` |
| `n_bonds2` | in | Number of bonds in molecule 2 |
| `bond_data2` | in | Flat bond array for molecule 2, length `n_bonds2*3` |
| `align_flag` | in | Enable structural alignment |
| `remap_flag` | in | Enable atom remapping |
| `heavy_flag` | in | Use only heavy (non-hydrogen) atoms |
| `mass_flag` | in | Weight atoms by atomic mass |
| `mirror_flag` | in | Mirror molecule 2 before comparison |
| `label_flag` | in | Use atom labels for type matching |
| `bond_flag` | in | Derive connectivity from geometry (ignores bond arrays) |
| `stats_flag` | in | Print optimisation statistics to stdout |
| `random_flag` | in | Seed RNG from system clock |
| `conv_freq` | in | Convergence frequency threshold |
| `max_trials` | in | Maximum number of optimisation trials |
| `rmsd` | out | Calculated RMSD (Å) |
| `natoms` | out | Number of elements written to `atomperm` |
| `atomperm` | out | Atom permutation, **0-based**; caller must allocate ≥ `natoms` elements |
| `transform` | out | 4 × 4 homogeneous transform (row-major, 16 doubles) |
| `error_code` | out | `0` = success; `1` = not isomers; `2` = missing bonds; `3` = atom type mismatch |

#### Minimal example

```c
#include <stdio.h>
#include "conformsd.h"

int main(void)
{
    int    ad1[] = { 6,0, 8,0, 1,0, 1,0 };  /* C, O, H, H (formaldehyde) */
    double xy1[] = { 0.000,  0.000, 0.000,
                     0.000,  1.208, 0.000,
                     0.935, -0.541, 0.000,
                    -0.935, -0.541, 0.000 };
    int bd1[] = { 1,2,2,  1,3,1,  1,4,1 };  /* C=O, C-H, C-H (1-based) */

    int    ad2[] = { 6,0, 8,0, 1,0, 1,0 };
    double xy2[] = { 0.010,  0.010, 0.000,
                     0.010,  1.218, 0.000,
                     0.945, -0.531, 0.000,
                    -0.925, -0.531, 0.000 };
    int bd2[] = { 1,2,2,  1,3,1,  1,4,1 };

    double rmsd, transform[16];
    int    atomperm[4], natoms, error_code;

    conformsd_calculate(
        4, ad1, xy1, 3, bd1,
        4, ad2, xy2, 3, bd2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false, /*bond_flag=*/false,
        /*stats=*/false, /*random=*/false,
        /*conv_freq=*/100, /*max_trials=*/10000,
        &rmsd, &natoms, atomperm, transform, &error_code);

    if (error_code != 0) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd);
    return 0;
}
```

Python API
----------

The Python API provides a higher-level interface to MolAlignLib.

You can try the Python API without any installation on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalignlib.git/devel?urlpath=%2Fdoc%2Ftree%2Fpython%2Fexamples%2Fexamples.ipynb)

### Installation

Python ≥ 3.8, scikit-build-core, Cython, NumPy and Chemfiles are required.

```bash
# Upgrade pip (recommended)
python3 -m pip install --user --upgrade pip

# Install build tools and dependencies
python3 -m pip install --user scikit-build-core cython numpy chemfiles

# Install the package (builds the Fortran library automatically)
python3 -m pip install --user .
```

> **Note:** `pip install .` does not build the standalone executables. Use the
> plain CMake workflow above to build `conformsd` and `atormsd`.

### Quick start

```python
from molalignlib import read_conformers, read_clusters, Conformer, AtomCluster
```

### AtomCluster

Represents an unstructured set of atoms with no bond topology. Useful for comparing metal clusters, nanoparticles, or other systems where connectivity is absent or irrelevant. Wraps `atormsd_calculate`.

#### Construction

```python
# From a file (single frame)
mol = AtomCluster.from_file("cluster.xyz")

# From a file (specific frame in a multi-frame trajectory)
mol = AtomCluster.from_file("trajectory.xyz", frame_idx=3)

# From numpy arrays directly
import numpy as np
atom_data = np.array([[26, 0], [26, 0]], dtype=np.int32)  # two Fe atoms
coords    = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], dtype=np.float64)
mol = AtomCluster(atom_data=atom_data, coords=coords, name="dimer")
```

#### Reading multiple frames

```python
clusters = read_clusters("trajectory.xyz")           # all frames → list
mol0, mol1 = read_clusters("trajectory.xyz", frames=(0, 1))
```

#### Computing RMSD

```python
result = mol0.rmsd_to(mol1,
    align=True,
    remap=True,
    heavy_only=False,
    mass_weighted=False,
    mirror=False,
    use_labels=False,
    stats=False,
    random=False,
    prune_tol=-1.0,    # disable pruning (default -1.0)
    conv_freq=10,
    max_trials=10000,
)
```

#### Writing output

```python
mol.write_xyz("output.xyz", comment="my cluster")
```

### Conformer

Represents a molecule with full bond topology. The HNA partitioning is used internally to guarantee chemically valid atom assignments, handling arbitrary degrees of topological symmetry efficiently. Wraps `conformsd_calculate`.

#### Construction

```python
# From a file (single frame; bond table is parsed automatically)
conf = Conformer.from_file("molecule.sdf")

# From a file (specific frame in a multi-frame SDF)
conf = Conformer.from_file("poses.sdf", frame_idx=3)

# From numpy arrays
import numpy as np
atom_data = np.array([[6,0],[8,0],[1,0],[1,0]], dtype=np.int32)
coords    = np.zeros((4, 3), dtype=np.float64)
bond_data = np.array([[1,2,2],[1,3,1],[1,4,1]], dtype=np.int32)  # 1-based
conf = Conformer(atom_data=atom_data, coords=coords, bond_data=bond_data)
```

#### Reading multiple frames

```python
conformers = read_conformers("poses.sdf")           # all frames → list
c0, c1 = read_conformers("poses.sdf", frames=(0, 1))
```

#### Computing RMSD

```python
result = c0.rmsd_to(c1,
    align=True,
    remap=True,
    heavy_only=False,
    mass_weighted=False,
    mirror=False,
    use_labels=False,
    bond_flag=False,   # True → derive bonds from geometry
    stats=False,
    random=False,
    conv_freq=100,
    max_trials=10000,
)
print(result.rmsd)              # float, Å
print(result.atom_permutation)  # int32 array, 0-based
print(result.transform)         # 4×4 float64 array
```

### RMSDResult

Returned by every `.rmsd_to()` call.

| Attribute | Type | Description |
|-----------|------|-------------|
| `rmsd` | `float` | Root-mean-square deviation in Å |
| `atom_permutation` | `int32 ndarray (n,)` | 0-based index array mapping *other* atoms onto *self* |
| `transform` | `float64 ndarray (4, 4)` | Homogeneous rotation + translation matrix (maps *other* → *self* frame) |

#### Applying the result

```python
# Produce a new object that is aligned and reordered to match the reference
other_aligned = result.apply_to(other)
other_aligned.write_xyz("aligned.xyz")
```

### Examples

#### Example 1 — Symmetry-corrected RMSD between two conformers

```python
from molalignlib import read_conformers

c0, c1 = read_conformers("PRDCC002527_poses.sdf", frames=(0, 1))

result = c0.rmsd_to(c1, remap=True, align=True)
print(f"RMSD = {result.rmsd:.4f} Å")

c1_aligned = result.apply_to(c1)
c1_aligned.write_xyz("aligned.xyz", comment=f"RMSD={result.rmsd:.4f}")
```

#### Example 2 — RMSD matrix for a full conformer ensemble

```python
from molalignlib import read_conformers

conformers = read_conformers("PRDCC002527_poses.sdf")

for c0 in conformers:
    for c1 in conformers:
        result = c1.rmsd_to(c0, remap=True, align=True)
        print(f"{result.rmsd:.4f}", end="  ")
    print()
```

#### Example 3 — RMSD of the first frame against every other frame (clusters)

```python
from molalignlib import read_clusters

clusters = read_clusters("Co138_frames.xyz")
ref = clusters[0]

for mol in clusters[1:]:
    result = mol.rmsd_to(ref, remap=True, align=True, prune_tol=0.1)
    print(f"{result.rmsd:.4f}")
```

#### Example 4 — Align one cluster to another and save the result

```python
from molalignlib import read_clusters

mol0, mol1 = read_clusters("Co138_frames.xyz", frames=(0, 1))

result = mol0.rmsd_to(mol1, remap=True, align=True, prune_tol=0.1, stats=True)
print(f"RMSD = {result.rmsd:.4f} Å")

mol1_aligned = result.apply_to(mol1)
mol1_aligned.write_xyz("Co138_aligned.xyz", comment=f"RMSD={result.rmsd:.4f}")
```

#### Example 5 — Build AtomCluster from an ASE object

```python
import numpy as np
from molalignlib import AtomCluster

# Suppose `ase_atoms` is an ASE Atoms object
atom_data = np.column_stack([
    ase_atoms.get_atomic_numbers(),
    np.zeros(len(ase_atoms), dtype=np.int32),
]).astype(np.int32)
coords = ase_atoms.get_positions().astype(np.float64)

cluster = AtomCluster(atom_data=atom_data, coords=coords)
```


Algorithm Notes
---------------

- **atoRMSD:** uses a stochastic strategy with distance-based pruning for unstructured clusters where no bond topology is available. Algorithm described in [Vásquez-Pérez et al., *J. Chem. Inf. Model.* (2023)](https://doi.org/10.1021/acs.jcim.2c01187).
- **confoRMSD:** uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to decompose the assignment problem into independent branches, reducing the number of evaluated combinations from the product of branch possibilities to their sum. For alignment calculations, the algorithm adaptively selects between stochastic orientation sampling (efficient for highly symmetric molecules) and exhaustive enumeration (efficient for molecules with few branches), based on the ratio of total to partial combinations in the assignment tree. Benchmarks show 100% topologically correct assignments across 1.4 million molecular pairs with millisecond-scale mean execution times. Full algorithm description in [Vásquez-Pérez et al., *J. Chem. Theory Comput.* (2026)](https://doi.org/10.1021/acs.jctc.6c00545).
- **Transform output:** the 4 × 4 homogeneous transformation matrix encodes both the optimal rotation *R* and the translation *t* needed to superimpose molecule 2 on molecule 1.
- **Default maximum trials:** 10,000 random orientations for alignment. The convergence frequency threshold defaults to 10 for `atormsd` and 100 for `conformsd`; numerical experiments show that values below 100 can produce incorrect assignments for conformers.
