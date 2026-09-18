MolAlignLib
===========

A high-performance library to compute optimal RMSD between atom clusters and symmetry-corrected RMSD between molecular conformers.

MolAlignLib uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to achieve exact, topologically valid atom assignments between conformers in milliseconds, even for highly symmetric molecules that are intractable by conventional graph-isomorphism approaches.

### Try It Online

You can try the Python API without any installation on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalignlib.git/devel?urlpath=%2Fdoc%2Ftree%2Fpython%2Fexamples%2Fexamples.ipynb)


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
   - [Demo programs](#demo-programs)
5. [Python API](#python-api)
   - [Installation](#installation)
   - [Atoms](#atomcluster)
   - [Conformer](#conformer)
   - [RMSDResult](#rmsdresult)
     - [Retrieving multiple ranked solutions](#retrieving-multiple-ranked-solutions)
   - [Examples](#examples)
     - [Using Atoms/Conformer directly](#using-atomsconformer-directly)
     - [Using the read_* functions](#using-the-read_-functions)
6. [Algorithm Notes](#algorithm-notes)


Overview
--------

MolAlignLib exposes two distinct RMSD calculation modes:

| Program | Function | When to use |
|------|-------------------|-------------|
| **atormsd** | atormsd_calculate | Unstructured atom clusters (no bond topology required) |
| **conformsd** | conformsd_calculate | Molecular conformers (bond topology guides atom matching via HNA partitioning) |

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
cmake -B build
make -C build
make -C build install
```

By default, `install` places files under the system prefix (`/usr/local` on Linux/macOS), which typically requires `sudo`:

```bash
sudo make -C build install
```

### User-local install

To install without root privileges, set `CMAKE_INSTALL_PREFIX` at configure time to a directory you own, e.g. `~/.local`:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local
make -C build
make -C build install
```

This installs:

| Artifact | Destination |
|----------|-------------|
| `atormsd`, `conformsd` | `$HOME/.local/bin` |
| `libmolalign.a` | `$HOME/.local/lib` |
| `molalign.h` | `$HOME/.local/include/molalignlib` |

Make sure `$HOME/.local/bin` is on your `PATH` to run the installed executables directly:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

After a successful build (before installing), the following are also available directly inside `build/fortran`:

| Artifact | Description |
|----------|--------------|
| **atormsd** | Standalone cluster RMSD program |
| **conformsd** | Standalone conformer RMSD program |
| **libmolalign.a** | Static library containing the core algorithms and both C bindings (`atormsd_calculate`, `conformsd_calculate`) |

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
| -align | | Align atoms to minimise the RMSD |
| -remap | | Remap atoms to minimise the RMSD |
| -label | | Use atom labels to distinguish atom types |
| -prune | TOL | Prune assignments with pair distances exceeding *TOL* Å |
| -freq | N | Stop if the best solution is found *N* consecutive times |
| -trials | N | Stop after at most *N* optimisation trials |
| -records | N | Record the *N* lowest RMSDs found (default: 1) |
| -assignment | | Print the optimised atom mapping to stdout |
| -aligned | FILE | Write aligned coordinates of molecule 2 to *FILE* |
| -heavy | | Ignore hydrogen atoms |
| -mass | | Use mass-weighted coordinates |
| -mirror | | Reflect molecule 2 before comparison |
| -stats | | Print detailed optimisation statistics |
| -random | | Seed the random-number generator from the system clock |

#### Examples

```bash
# Plain RMSD (no alignment, no remapping)
atormsd mol1.xyz mol2.xyz

# Align and remap atoms, print the mapping
atormsd mol1.xyz mol2.xyz -align -remap -assignment

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
| -align | | Align atoms to minimise the RMSD |
| -remap | | Remap atoms to minimise the RMSD |
| -label | | Use atom labels to distinguish atom types |
| -freq | N | Stop if the best solution is found *N* consecutive times |
| -trials | N | Stop after at most *N* optimisation trials |
| -records | N | Record the *N* lowest RMSDs found (default: 1) |
| -assignment | | Print the optimised atom mapping to stdout |
| -assigntree | | Print the internal assignment tree |
| -aligned | FILE | Write aligned coordinates of molecule 2 to *FILE* |
| -heavy | | Ignore hydrogen atoms |
| -mass | | Use mass-weighted coordinates |
| -mirror | | Reflect molecule 2 before comparison |
| -bond | TOL | Derive bond connectivity from interatomic distances instead of the file's bond table, using detection tolerance *TOL* |
| -exhaustive | | Force exhaustive orientation-independent search regardless of assignment tree topology |
| -stochastic | | Force stochastic fixed-orientation search regardless of assignment tree topology |
| -stats | | Print detailed optimisation statistics |
| -random | | Seed the random-number generator from the system clock |

`-exhaustive` and `-stochastic` are mutually exclusive. If neither is given, the strategy is chosen automatically based on the ratio of total to partial assignment combinations in the tree: stochastic fixed-orientation search is used when the ratio is high, and exhaustive orientation-independent search when it is low.

#### Examples

```bash
# Symmetry-corrected RMSD with remapping and alignment
conformsd conf1.sdf conf2.sdf -align -remap

# Derive connectivity from geometry (useful for XYZ input)
conformsd conf1.xyz conf2.xyz -align -remap -bond 0.3

# Heavy atoms only, exhaustive search
conformsd conf1.sdf conf2.sdf -align -remap -heavy -exhaustive

# Print the atom permutation that maps conf2 onto conf1
conformsd conf1.sdf conf2.sdf -align -remap -assignment

# Compute an RMSD matrix for all poses in an SDF file (shell loop)
for i in 1 2 3; do
  for j in 1 2 3; do
    conformsd pose${i}.sdf pose${j}.sdf -align -remap
  done
done
```

C Binding
---------

Include `molalign.h` and link against `libmolalign`, which contains both `atormsd_calculate` and `conformsd_calculate`: one header, one library, no separate per-binding include or link step.

```c
#include "molalign.h" /* for atormsd_calculate and conformsd_calculate */
```

### Atom data layout

Both functions receive atom information as a flat `int` array of length
`n_atoms * 2`, packed as:

```
[ elnum_0, label_0, elnum_1, label_1, ... ]
```

- `elnum`: atomic number (e.g. 6 for carbon, 8 for oxygen).
- `label`: user-defined integer label; pass `0` for unlabelled atoms.

Coordinates are passed as a flat `double` array of length `n_atoms * 3`,
in row-major (C) order: `[x_0, y_0, z_0, x_1, y_1, z_1, ...]`.

### Transform output

Both functions write row-major 4 × 4 homogeneous transformation matrices
(16 `double` values each) to the `transform_list` output parameter, one
per returned record, flattened back-to-back. Record `i` (0-based) occupies
`transform_list[i*16 .. i*16+15]`. Each matrix maps molecule-2 coordinates
into the molecule-1 reference frame:

```
p_out = R * p_in + t
```

A matrix is the identity when `align_flag = false`.

### Multiple ranked solutions

Both functions accept an `n_records` parameter requesting up to that many
ranked candidate solutions (best RMSD first) instead of just the single
best one. `occ_records` reports how many were actually found; it can be
smaller than `n_records`, and it is always `1` unless both `align_flag`
and `remap_flag` are true. All output arrays (`rmsd_list`, `atomperm_list`,
`transform_list`) are flattened and must be pre-allocated by the caller
for `n_records` records; only the first `occ_records` entries are
meaningful. `natoms` is the same for every record and is written once.

### atormsd_calculate

```c
void atormsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag,
    bool print_stats, bool random_flag,
    double prune_tol, int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *natoms, int *atomperm_list,
    double *transform_list, int *occ_records, int *error_code);
```

#### Parameters

| Parameter | Direction | Description |
|-----------|-----------|-------------|
| **n_atoms1** | in | Number of atoms in molecule 1 |
| **atom_data1** | in | Packed atom data for molecule 1, length `n_atoms1*2` |
| **coords1** | in | Coordinates for molecule 1, length `n_atoms1*3` |
| **n_atoms2** | in | Number of atoms in molecule 2 |
| **atom_data2** | in | Packed atom data for molecule 2, length `n_atoms2*2` |
| **coords2** | in | Coordinates for molecule 2, length `n_atoms2*3` |
| **align_flag** | in | Enable structural alignment |
| **remap_flag** | in | Enable atom remapping |
| **heavy_flag** | in | Use only heavy (non-hydrogen) atoms |
| **mass_flag** | in | Weight atoms by atomic mass |
| **mirror_flag** | in | Mirror molecule 2 before comparison |
| **label_flag** | in | Use atom labels for type matching |
| **print_stats** | in | Print optimisation statistics to stdout |
| **random_flag** | in | Seed RNG from system clock |
| **prune_tol** | in | Pruning distance tolerance (Å); negative value disables pruning |
| **conv_freq** | in | Convergence frequency threshold |
| **max_trials** | in | Maximum number of optimisation trials |
| **n_records** | in | Maximum number of ranked candidate solutions to return (≥ 1). Values > 1 only take effect when `align_flag` and `remap_flag` are both true |
| **rmsd_list** | out | RMSD (Å) of each returned record, length `n_records`; caller allocates |
| **natoms** | out | Number of elements per permutation record (same for every record) |
| **atomperm_list** | out | Flattened, **0-based** atom permutations, length `n_records*natoms`. Record `i` occupies `atomperm_list[i*natoms .. i*natoms+natoms-1]`; caller allocates ≥ `n_records*natoms` elements |
| **transform_list** | out | Flattened row-major 4 × 4 homogeneous transforms, length `n_records*16`. Record `i` occupies `transform_list[i*16 .. i*16+15]`; caller allocates ≥ `n_records*16` elements |
| **occ_records** | out | Actual number of records written (≤ `n_records`) |
| **error_code** | out | `0` = success; `1` = not isomers; `2` = atom type mismatch |

#### Minimal example

```c
#include <stdio.h>
#include <stdlib.h>
#include "molalign.h"

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

    /* Request a single (best) solution */
    double rmsd_list[1], transform_list[16];
    int    atomperm_list[3], natoms, occ_records, error_code;

    atormsd_calculate(
        3, ad1, xy1,
        3, ad2, xy2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false,
        /*stats=*/false, /*random=*/false,
        /*prune_tol=*/-1.0, /*conv_freq=*/10, /*max_trials=*/10000,
        /*n_records=*/1,
        rmsd_list, &natoms, atomperm_list,
        transform_list, &occ_records, &error_code);

    if (error_code != 0) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd_list[0]);
    return 0;
}
```

To retrieve several ranked solutions, allocate the output arrays for
`n_records` records and loop over the first `occ_records` entries:

```c
#define N_RECORDS 5

double rmsd_list[N_RECORDS];
double transform_list[N_RECORDS * 16];
int    atomperm_list[N_RECORDS * 3];  /* N_RECORDS * natoms (known here to be 3) */
int    natoms, occ_records, error_code;

atormsd_calculate(
    3, ad1, xy1, 3, ad2, xy2,
    /*align=*/true, /*remap=*/true,
    /*heavy=*/false, /*mass=*/false,
    /*mirror=*/false, /*label=*/false,
    /*stats=*/false, /*random=*/false,
    /*prune_tol=*/-1.0, /*conv_freq=*/10, /*max_trials=*/10000,
    /*n_records=*/N_RECORDS,
    rmsd_list, &natoms, atomperm_list,
    transform_list, &occ_records, &error_code);

for (int i = 0; i < occ_records; i++) {
    printf("Solution %d: RMSD = %.6f Å\n", i, rmsd_list[i]);
}
```

Compile:

```bash
gcc example.c -o example -lmolalign -lgfortran -lm
```


### conformsd_calculate

```c
void conformsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_bonds1, const int *bond_data1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    int n_bonds2, const int *bond_data2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag, bool bond_flag, double bond_tol,
    bool print_stats, bool print_assigntree, bool random_flag,
    int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *natoms, int *atomperm_list,
    double *transform_list, int *occ_records, int *error_code);
```

Bond data is passed as a flat `int` array of length `n_bonds * 3`, packed as:

```
[ atom1_0, atom2_0, type_0, atom1_1, atom2_1, type_1, ... ]
```

Atom indices are **1-based**. When `bond_flag = true` the library derives
connectivity from atomic geometry and the bond arrays may be empty
(`n_bonds = 0`, `bond_data = NULL`). In that case `bond_tol` (bond detection
tolerance) has no default and must be supplied; it is ignored when
`bond_flag = false`.

#### Parameters

| Parameter | Direction | Description |
|-----------|-----------|-------------|
| **n_atoms1** | in | Number of atoms in molecule 1 |
| **atom_data1** | in | Packed atom data for molecule 1, length `n_atoms1*2` |
| **coords1** | in | Coordinates for molecule 1, length `n_atoms1*3` |
| **n_bonds1** | in | Number of bonds in molecule 1 |
| **bond_data1** | in | Flat bond array for molecule 1, length `n_bonds1*3` |
| **n_atoms2** | in | Number of atoms in molecule 2 |
| **atom_data2** | in | Packed atom data for molecule 2, length `n_atoms2*2` |
| **coords2** | in | Coordinates for molecule 2, length `n_atoms2*3` |
| **n_bonds2** | in | Number of bonds in molecule 2 |
| **bond_data2** | in | Flat bond array for molecule 2, length `n_bonds2*3` |
| **align_flag** | in | Enable structural alignment (default `false` in the Python/Cython wrapper) |
| **remap_flag** | in | Enable atom remapping (default `false` in the Python/Cython wrapper) |
| **heavy_flag** | in | Use only heavy (non-hydrogen) atoms |
| **mass_flag** | in | Weight atoms by atomic mass |
| **mirror_flag** | in | Mirror molecule 2 before comparison |
| **label_flag** | in | Use atom labels for type matching |
| **bond_flag** | in | Derive connectivity from geometry (ignores bond arrays) |
| **bond_tol** | in | Bond detection tolerance; required when `bond_flag = true`, ignored otherwise (no default) |
| **print_stats** | in | Print optimisation statistics to stdout |
| **print_assigntree** | in | Print the internal assignment tree |
| **random_flag** | in | Seed RNG from system clock |
| **conv_freq** | in | Convergence frequency threshold |
| **max_trials** | in | Maximum number of optimisation trials |
| **n_records** | in | Maximum number of ranked candidate solutions to return (≥ 1). Values > 1 only take effect when `align_flag` and `remap_flag` are both true |
| **rmsd_list** | out | RMSD (Å) of each returned record, length `n_records`; caller allocates |
| **natoms** | out | Number of elements per permutation record (same for every record) |
| **atomperm_list** | out | Flattened, **0-based** atom permutations, length `n_records*natoms`. Record `i` occupies `atomperm_list[i*natoms .. i*natoms+natoms-1]`; caller allocates ≥ `n_records*natoms` elements |
| **transform_list** | out | Flattened row-major 4 × 4 homogeneous transforms, length `n_records*16`. Record `i` occupies `transform_list[i*16 .. i*16+15]`; caller allocates ≥ `n_records*16` elements |
| **occ_records** | out | Actual number of records written (≤ `n_records`) |
| **error_code** | out | `0` = success; `1` = not isomers; `2` = atom type mismatch; `3` = missing bonds; `4` = bond connectivity mismatch (only possible when `remap_flag = false`). Numbered to match `atormsd_calculate` where applicable. |

#### Minimal example

```c
#include <stdio.h>
#include "molalign.h"

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

    /* Request a single (best) solution */
    double rmsd_list[1], transform_list[16];
    int    atomperm_list[4], natoms, occ_records, error_code;

    conformsd_calculate(
        4, ad1, xy1, 3, bd1,
        4, ad2, xy2, 3, bd2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false,
        /*bond_flag=*/false, /*bond_tol=*/0.0,
        /*stats=*/false, /*print_assigntree=*/false, /*random=*/false,
        /*conv_freq=*/100, /*max_trials=*/10000,
        /*n_records=*/1,
        rmsd_list, &natoms, atomperm_list,
        transform_list, &occ_records, &error_code);

    if (error_code != 0) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd_list[0]);
    return 0;
}
```

As with `atormsd_calculate`, pass `n_records > 1` and size the output
arrays accordingly to retrieve several ranked solutions in one call; loop
over the first `occ_records` entries.

### Demo programs

Two self-contained demo programs, [`atormsd_demo.c`](atormsd_demo.c) and
[`conformsd_demo.c`](conformsd_demo.c), show how to call `atormsd_calculate`
and `conformsd_calculate` directly from C, using flat coordinate/element
arrays read from XYZ files. Compile and run instructions are in each file's
header comment.

Python API
----------

The Python API provides a higher-level interface to MolAlignLib. (See [Try It Online](#try-it-online) above to test it on Binder without installing anything.)

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

### Atoms

Represents an unstructured set of atoms with no bond topology. Useful for comparing metal clusters, nanoparticles, or other systems where connectivity is absent or irrelevant. Wraps `atormsd_calculate`.

#### Construction

```python
# From a file (defaults to the first frame)
mol = Atoms.from_file("clusters.xyz")

# From a file (specific frame in a multi-frame trajectory)
mol = Atoms.from_file("clusters.xyz", frame_idx=2)

# From element symbols and coordinates directly
import numpy as np
symbols = ["Fe", "Fe"]
coords  = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], dtype=np.float64)
mol = Atoms.from_symbols(symbols, coords, name="dimer")

# From atomic numbers and coordinates directly
mol = Atoms.from_numbers([26, 26], coords, name="dimer")  # two Fe atoms
```

#### Reading multiple frames

```python
clusters = read_clusters("clusters.xyz")           # all frames, returns a list
mol0, mol1 = read_clusters("clusters.xyz", frames=(0, 1))
```

#### Computing RMSD

`rmsd_to()` always returns a **list** of `RMSDResult`, ranked best (lowest
RMSD) first. By default only the single best solution is computed:

```python
results = mol0.rmsd_to(mol1,
    align=True,       # default: False
    remap=True,       # default: False
    heavy_only=False,
    mass_weighted=False,
    mirror=False,
    use_labels=False,
    stats=False,
    random=False,
    prune_tol=-1.0,    # disable pruning (default -1.0)
    conv_freq=10,
    max_trials=10000,
    n_records=1,       # request up to this many ranked solutions
)
result = results[0]   # best (lowest RMSD) solution
print(result.rmsd)
```

Pass `n_records > 1` (with `align=True, remap=True`) to retrieve several
ranked candidate solutions in one call; see
[Retrieving multiple ranked solutions](#retrieving-multiple-ranked-solutions).

#### Writing output

```python
mol.write_xyz("output.xyz", comment="my cluster")
```

### Conformer

Represents a molecule with full bond topology. The HNA partitioning is used internally to guarantee chemically valid atom assignments, handling arbitrary degrees of topological symmetry efficiently. Wraps `conformsd_calculate`.

#### Construction

```python
# From a file (defaults to the first frame; bond table is parsed automatically)
conf = Conformer.from_file("conformers.sdf")

# From a file (specific frame in a multi-frame SDF)
conf = Conformer.from_file("conformers.sdf", frame_idx=2)

# From element symbols, coordinates, and bonds directly
import numpy as np
symbols   = ["C", "O", "H", "H"]
coords    = np.zeros((4, 3), dtype=np.float64)
bond_data = np.array([[1,2,2],[1,3,1],[1,4,1]], dtype=np.int32)  # 1-based
conf = Conformer.from_symbols(symbols, coords, bond_data=bond_data)

# From atomic numbers, coordinates, and bonds directly
conf = Conformer.from_numbers([6, 8, 1, 1], coords, bond_data=bond_data)
```

#### Reading multiple frames

```python
conformers = read_conformers("conformers.sdf")           # all frames, returns a list
c0, c1 = read_conformers("conformers.sdf", frames=(0, 1))
```

#### Computing RMSD

`rmsd_to()` always returns a **list** of `RMSDResult`, ranked best (lowest
RMSD) first. By default only the single best solution is computed:

```python
results = c0.rmsd_to(c1,
    align=True,        # default: False
    remap=True,         # default: False
    heavy_only=False,
    mass_weighted=False,
    mirror=False,
    use_labels=False,
    bond_flag=False,   # True means derive bonds from geometry
    bond_tol=None,     # required (float), no default, when bond_flag=True
    stats=False,
    random=False,
    conv_freq=100,
    max_trials=10000,
    n_records=1,        # request up to this many ranked solutions
)
result = results[0]             # best (lowest RMSD) solution
print(result.rmsd)              # float, Å
print(result.atom_permutation)  # int32 array, 0-based
print(result.transform)         # 4×4 float64 array
```

Pass `n_records > 1` (with `align=True, remap=True`) to retrieve several
ranked candidate solutions in one call; see
[Retrieving multiple ranked solutions](#retrieving-multiple-ranked-solutions).

### RMSDResult

`.rmsd_to()` returns a `list[RMSDResult]`, one element per requested
solution (best/lowest RMSD first). Each `RMSDResult` holds:

| Attribute | Type | Description |
|-----------|------|-------------|
| **rmsd** | `float` | Root-mean-square deviation in Å |
| **atom_permutation** | `int32 ndarray (n,)` | 0-based index array mapping *other* atoms onto *self* |
| **transform** | `float64 ndarray (4, 4)` | Homogeneous rotation + translation matrix (maps *other* to *self* frame) |

#### Applying the result

```python
# Produce a new object that is aligned and reordered to match the reference
result = mol0.rmsd_to(other, align=True, remap=True)[0]
other_aligned = result.apply_to(other)
other_aligned.write_xyz("aligned.xyz")
```

#### Retrieving multiple ranked solutions

Set `n_records` to inspect several distinct candidate mappings instead of
just the best one. This only produces more than one result when both
`align=True` and `remap=True`; the returned list is otherwise always
length 1 regardless of `n_records`, and may be shorter than `n_records`
if fewer distinct solutions were found:

```python
results = mol0.rmsd_to(mol1, align=True, remap=True, n_records=5)

for i, result in enumerate(results):
    print(f"Solution {i}: RMSD = {result.rmsd:.4f} Å")

best = results[0]
```

### Examples

The snippets below are grouped by how the structures get loaded:
build a single `Atoms`/`Conformer` object at a time with `.from_file()`
(or from in-memory data), or use `read_clusters`/`read_conformers` to pull
several frames out of a file at once. See [Atoms](#atomcluster) and
[Conformer](#conformer) above for the full constructor reference.

Examples that read a file use [`clusters.xyz`](python/examples/clusters.xyz)
(3 unstructured Co clusters) and [`conformers.sdf`](python/examples/conformers.sdf)
(3 bonded conformer poses), the same files used in
[`python/examples/examples.py`](python/examples/examples.py).

#### Using Atoms/Conformer directly

##### Single conformer pair, from a specific SDF frame

```python
from molalignlib import Conformer

c0 = Conformer.from_file("conformers.sdf", frame_idx=0)
c1 = Conformer.from_file("conformers.sdf", frame_idx=1)

result = c0.rmsd_to(c1, remap=True, align=True)[0]
print(f"RMSD = {result.rmsd:.4f} Å")
```

##### Single cluster pair, from a specific XYZ frame

```python
from molalignlib import Atoms

mol0 = Atoms.from_file("clusters.xyz", frame_idx=0)
mol1 = Atoms.from_file("clusters.xyz", frame_idx=1)

result = mol0.rmsd_to(mol1, remap=True, align=True)[0]
print(f"RMSD = {result.rmsd:.4f} Å")

mol1_aligned = result.apply_to(mol1)
mol1_aligned.write_xyz("aligned.xyz", comment=f"RMSD={result.rmsd:.4f}")
```

##### Single cluster pair, from ASE objects

The equivalent workflow for a pair of `ase.Atoms` objects, using
`Atoms` instead of `Conformer` since no bond topology is involved:

```python
import numpy as np
from molalignlib import Atoms

def cluster_from_ase(ase_atoms):
    return Atoms.from_symbols(
        ase_atoms.get_chemical_symbols(),
        ase_atoms.get_positions().astype(np.float64),
    )

# Suppose `ase_atoms0` and `ase_atoms1` are the two structures being compared
mol0 = cluster_from_ase(ase_atoms0)
mol1 = cluster_from_ase(ase_atoms1)

result = mol0.rmsd_to(mol1, remap=True, align=True)[0]
print(f"RMSD = {result.rmsd:.4f} Å")

mol1_aligned = result.apply_to(mol1)
mol1_aligned.write_xyz("aligned.xyz", comment=f"RMSD={result.rmsd:.4f}")
```

#### Using the read_* functions

##### Retrieve multiple ranked mappings for a cluster pair

```python
from molalignlib import read_clusters

mol0, mol1 = read_clusters("clusters.xyz", frames=(0, 1))

# Ask for the top 5 lowest RMSD mappings
results = mol0.rmsd_to(mol1, remap=True, align=True, prune_tol=0.1, n_records=5)

for i, result in enumerate(results, start=1):
    print(f"Mapping {i}: RMSD = {result.rmsd:.4f} Å")

best = results[0]  # the list is sorted best-first
```

##### RMSD matrix over all conformer pairs

```python
from molalignlib import read_conformers

conformers = read_conformers("conformers.sdf")  # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol1.rmsd_to(mol0, remap=True, align=True)[0]
        print(f"{result.rmsd:.4f}", end=2 * " ")
    print()
```

##### RMSD matrix, deriving bond connectivity from geometry

Useful for XYZ input, which (unlike SDF/MOL2) carries no bond table of
its own:

```python
from molalignlib import read_conformers

conformers = read_conformers("conformers.xyz")  # read all frames

for mol0 in conformers:
    for mol1 in conformers:
        result = mol1.rmsd_to(mol0, remap=True, align=True,
                               bond_flag=True, bond_tol=0.3)[0]
        print(f"{result.rmsd:.4f}", end=2 * " ")
    print()
```

For the full runnable scripts, see
[`python/examples/examples.py`](python/examples/examples.py) and the
equivalent Jupyter notebook, [`python/examples/examples.ipynb`](python/examples/examples.ipynb)
(the one launched by the [Binder link](#try-it-online) above).


Algorithm Notes
---------------

- **atoRMSD:** uses a stochastic strategy with distance-based pruning for unstructured clusters where no bond topology is available. Algorithm described in [Vásquez-Pérez et al., *J. Chem. Inf. Model.* (2023)](https://doi.org/10.1021/acs.jcim.2c01187).
- **confoRMSD:** uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to decompose the assignment problem into independent branches, reducing the number of evaluated combinations from the product of branch possibilities to their sum. For alignment calculations, the algorithm adaptively selects between stochastic orientation sampling (efficient for highly symmetric molecules) and exhaustive enumeration (efficient for molecules with few branches), based on the ratio of total to partial combinations in the assignment tree. Benchmarks show 100% topologically correct assignments across 1.4 million molecular pairs with millisecond-scale mean execution times. Full algorithm description in [Vásquez-Pérez et al., *J. Chem. Theory Comput.* (2026)](https://doi.org/10.1021/acs.jctc.6c00545).
- **Transform output:** the 4 × 4 homogeneous transformation matrix encodes both the optimal rotation *R* and the translation *t* needed to superimpose molecule 2 on molecule 1.
- **Default maximum trials:** 10,000 random orientations for alignment. The convergence frequency threshold defaults to 10 for `atormsd` and 100 for `conformsd`; numerical experiments show that values below 100 can produce incorrect assignments for conformers.
