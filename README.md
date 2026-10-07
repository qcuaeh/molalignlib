MolAlignLib
===========

A high-performance library to compute optimal RMSD between atom clusters and symmetry-corrected RMSD between molecular conformers.

MolAlignLib uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to achieve exact, topologically valid atom assignments between conformers in milliseconds, even for highly symmetric molecules that are intractable by conventional graph-isomorphism approaches.

### Try It Online

You can try the Python bindings right away on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qcuaeh/molalignlib.git/devel?urlpath=%2Fdoc%2Ftree%2Fpython%2Fexamples%2Fexamples.ipynb)


Table of Contents
-----------------
1. [Overview](#overview)
2. [Building from Source](#building-from-source)
3. [Command-line Programs](#command-line-programs)
   - [atormsd](#atormsd)
   - [conformsd](#conformsd)
4. [C Binding](#c-binding)
   - [atormsd](#atormsd-1)
   - [conformsd](#conformsd-1)
   - [Demo program](#demo-program)
5. [Python bindings](#python-api)
   - [Installation](#installation)
   - [Molecule](#molecule)
   - [atormsd_to](#atormsd_to)
   - [conformsd_to](#conformsd_to)
   - [RMSDResult](#rmsdresult)
6. [Runnable examples](#runnable-examples)
7. [Algorithm Notes](#algorithm-notes)


Overview
--------

MolAlignLib exposes two distinct RMSD calculation modes:

| Program | C function | Python method | When to use |
|------|-------------------|---------------|-------------|
| **atormsd** | atormsd_calculate | Molecule.atormsd_to | Unstructured atom clusters (no bond topology required) |
| **conformsd** | conformsd_calculate | Molecule.conformsd_to | Molecular conformers (bond topology guides atom matching via HNA partitioning) |

Both modes support:
- **Alignment** (`-align`): optimally rotate and translate one structure onto the other.
- **Remapping** (`-remap`): find the atom permutation that minimises the RMSD.
- **Mass weighting** (`-massweight`): weight each atom by its atomic mass.
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

See the [Python bindings](#python-api) section below.


Command-line Programs
---------------------

### atormsd

Calculate RMSD between two unstructured atom clusters. No bond information is needed; atom matching is guided by element type and spatial proximity alone. (With `-heavy`, bonds read from the files are only used afterwards, to pair each excluded hydrogen with the image of its bonded heavy atom.)

```
atormsd file1 file2 [options]
```

**Supported file formats:** XYZ, MOL, SDF and Mol2, selected by the lowercase file extension (`.xyz`, `.mol`, `.sdf`, `.mol2`); only the first structure of each file is read. Dummy atoms (symbol `X`, or SYBYL type `Du` in Mol2 files) are skipped, together with their bonds, and the remaining atoms are numbered consecutively. Options may appear before or after the file paths and are case-insensitive.

**Output:** one line per record with the RMSD (Å). With `-remap -mapping` it is followed by the atom permutation that maps file 2 onto file 1: a comma-separated list of 1-based indices, where entry *i* is the atom of file 2 placed on line *i* of file 1.

#### Options

| Option | Argument | Description |
|--------|----------|-------------|
| -align | | Align atoms to minimise the RMSD |
| -remap | | Remap atoms to minimise the RMSD |
| -atomtype | | Only match atoms with the same label. A label is an element symbol optionally followed by digits (e.g. `C`, `C1`, `C2`); atoms with different digit suffixes are never matched |
| -prunetol | TOL | Prune atom pairs (only with `-remap`): two atoms are never paired if their sorted distances to the atoms of some atom type differ by more than 2√3 × *TOL* Å |
| -atofreq | N | Stop once the best solution has been found *N* times (default: 10) |
| -maxtrials | N | Stop after at most *N* random orientations (default: 10000) |
| -maxrecs | N | Record the *N* lowest RMSDs found (default: 1); more than one record is only produced with `-align -remap` |
| -mapping | | Print the atom permutation after each RMSD (only with `-remap`; see *Output* above) |
| -aligned | FILE | Write molecule 2, aligned and reordered to match molecule 1, to *FILE* instead of printing the results: one structure per record, with `rmsd=` and the RMSD as its title. Only takes effect with `-align`. The format is taken from the extension (`.xyz`, `.sdf` or `.mol2`); a bare extension such as `xyz` writes to stdout |
| -heavy | | Ignore hydrogen atoms |
| -massweight | | Use mass-weighted coordinates |
| -mirror | | Reflect molecule 2 (x → −x) before comparison |
| -stats | | Print the ranked local minima found and the search statistics (only with `-align -remap`) |
| -random | | Seed the random-number generator from the system clock (otherwise results are reproducible) |

The integer arguments of `-atofreq`, `-maxtrials` and `-maxrecs` must be at least 1. Option arguments may not start with `-`.

#### Examples

```bash
# Plain RMSD (no alignment, no remapping)
atormsd mol1.xyz mol2.xyz

# Align and remap atoms, print the mapping
atormsd mol1.xyz mol2.xyz -align -remap -mapping

# Remap with distance pruning, keep the 3 best solutions
atormsd mol1.xyz mol2.xyz -align -remap -prunetol 0.5 -maxrecs 3

# Heavy-atom, mass-weighted RMSD with alignment
atormsd mol1.sdf mol2.sdf -align -remap -heavy -massweight

# Write the aligned and reordered molecule 2 to a new file (instead of printing the RMSD)
atormsd mol1.xyz mol2.xyz -align -remap -aligned mol2_aligned.xyz
```


### conformsd

Calculate the symmetry-corrected RMSD between two molecular conformers. Bond topology is used to build the HNA partition, which guides atom matching and guarantees chemically valid assignments even for highly symmetric molecules. For alignment calculations, the algorithm automatically selects between stochastic orientation sampling and direct enumeration based on the structure of the assignment tree, avoiding unnecessary computation.

```
conformsd file1 file2 [options]
```

**Supported file formats:** XYZ, MOL, SDF and Mol2, selected by the lowercase file extension (`.xyz`, `.mol`, `.sdf`, `.mol2`); only the first structure of each file is read. XYZ files have no bond table, so use `-bondtol` with them. Dummy atoms (symbol `X`, or SYBYL type `Du` in Mol2 files) are skipped, together with their bonds, and the remaining atoms are numbered consecutively. Options may appear before or after the file paths and are case-insensitive.

**Output:** the same as for [atormsd](#atormsd): one line per record with the RMSD (Å), followed with `-remap -mapping` by the comma-separated, 1-based atom permutation.

#### Options

| Option | Argument | Description |
|--------|----------|-------------|
| -align | | Align atoms to minimise the RMSD |
| -remap | | Remap atoms to minimise the RMSD |
| -atomtype | | Only match atoms with the same label. A label is an element symbol optionally followed by digits (e.g. `C`, `C1`, `C2`); atoms with different digit suffixes are never matched |
| -confofreq | N | Stop the random orientation search once the best solution has been found more than *N* times; also the strategy threshold (see below) (default: 100) |
| -maxtrials | N | Stop after at most *N* random orientations (default: 10000) |
| -maxfrags | N | Maximum number of molecular fragments allowed in each molecule (default: 1; see below) |
| -maxrecs | N | Record the *N* lowest RMSDs found (default: 1); more than one record is only produced with `-align -remap` when the stochastic search is selected (see below) |
| -mapping | | Print the atom permutation after each RMSD (only with `-remap`; see *Output* above) |
| -assigntree | | Print the assignment tree (one node per branch point, labelled element × number of equivalent atoms) and its total and partial combination counts (only with `-remap`) |
| -aligned | FILE | Write molecule 2, aligned and reordered to match molecule 1, to *FILE* instead of printing the results: one structure per record, with `rmsd=` and the RMSD as its title. Only takes effect with `-align`. The format is taken from the extension (`.xyz`, `.sdf` or `.mol2`); a bare extension such as `xyz` writes to stdout |
| -heavy | | Ignore hydrogen atoms |
| -massweight | | Use mass-weighted coordinates |
| -mirror | | Reflect molecule 2 (x → −x) before comparison |
| -bondtol | TOL | Derive bond connectivity from interatomic distances instead of the file's bond table: atoms are bonded when closer than the sum of their covalent radii plus *TOL* Å |
| -bondtype | | Use bond types from the files to guide atom matching (see *Bond types* below). No effect with `-bondtol` |
| -stats | | Print the ranked local minima found and the search statistics (only with `-align -remap`) |
| -random | | Seed the random-number generator from the system clock (otherwise results are reproducible) |

With `-align -remap`, the search strategy is chosen automatically from the assignment tree: stochastic fixed-orientation search is used when the total number of complete assignments exceeds `-confofreq` times the sum of partial combinations, and exhaustive orientation-independent search otherwise. The default of 100 reproduced the reference assignments on the full CCD and BIRD benchmarks, whereas 25 produced some incorrect ones. Either strategy can be forced at compile time with the `FORCE_STOCHASTIC` and `FORCE_EXHAUSTIVE` parameters in `fortran/parameters.f90`.

A molecular fragment is a connected component of the bond graph of the compared atoms (after `-heavy` exclusions). An atom without bonds is a fragment of its own, so with the default `-maxfrags 1` each molecule must be a single connected structure, and a molecule of two or more atoms with no bonds is rejected. Raise `-maxfrags` to compare complexes, salts or solvated systems made of several molecules. The integer arguments of `-confofreq`, `-maxtrials`, `-maxfrags` and `-maxrecs` must be at least 1. Option arguments may not start with `-`.

#### Bond types

A bond type is an opaque label: a short, case-insensitive string with no intrinsic meaning. With `-bondtype`, two bonds match only when their labels are the same string, so `3/2`, `1.5` and `6/4` are three different types. The only exception is that directed types are compared undirected (`dr`/`dl` as `dv`, `up`/`dn` as `1`). A label is one of:

- a digit `1`–`9`;
- a letter followed by a letter or a digit (e.g. `ar`, `h3`);
- a digit, one of the separators `/ : - . ,`, and a digit (e.g. `3/2`, `3:2`).

The readers give bonds conventional labels, so that bond types read from different formats can be compared: Mol2 bond types are kept as they are (`1`, `2`, `3`, `ar`, `am`, or any other valid label as an extension), and MOL/SDF bond numbers 1–10 other than 8 become `1`, `2`, `3`, `ar`, `sd`, `sa`, `da`, `co` and `hb`. Mol2 bonds of type `un` and MOL/SDF bonds of type 8 (any) are of undefined type, which has no label. Entries that are not bonds (Mol2 `du` and `nc`, MOL/SDF type 0) are skipped. The full list of conventional labels is in `fortran/parameters.f90`.

Bonds of undefined type (see above, and all bonds derived with `-bondtol`) cannot be compared, so with `-bondtype` a bond of undefined type between compared atoms is an error. The same molecule written with Kekulé bonds in one file and aromatic bonds in the other does not match either.

#### Examples

```bash
# Symmetry-corrected RMSD with remapping and alignment
conformsd conf1.sdf conf2.sdf -align -remap

# Derive connectivity from geometry (useful for XYZ input)
conformsd conf1.xyz conf2.xyz -align -remap -bondtol 0.3

# Also distinguish bond types (labels are compared literally, see Bond types)
conformsd conf1.sdf conf2.mol2 -align -remap -bondtype

# Heavy atoms only
conformsd conf1.sdf conf2.sdf -align -remap -heavy

# Allow up to two fragments per molecule (e.g. a host-guest complex)
conformsd complex1.sdf complex2.sdf -align -remap -maxfrags 2

# Print the atom permutation that maps conf2 onto conf1
conformsd conf1.sdf conf2.sdf -align -remap -mapping

# Compute an RMSD matrix for a set of pose files (shell loop; each
# program call reads only the first structure of each file)
for i in 1 2 3; do
  for j in 1 2 3; do
    conformsd pose${i}.sdf pose${j}.sdf -align -remap
  done
done
```

C Binding
---------

Include `molalign.h` and link against `libmolalign`, which contains both `atormsd_calculate` and `conformsd_calculate`: one header and one library, no separate per-binding include or link step. `molalign.h` also defines the `MOLALIGN_SUCCESS` and `MOLALIGN_ERROR_*` constants returned through `error_code`, and documents the bond type codes.

```c
#include "molalign.h"
```

The header is installed in the `molalignlib` subdirectory of the include directory, which is not on the compiler's default search path, so compile with `-I<prefix>/include/molalignlib`. Since `libmolalign.a` is a static Fortran library, list it before the Fortran runtime and math libraries when linking: `-lmolalign -lgfortran -lm`.

### Atom data layout

Both functions receive atom information as a flat `int` array of length
`n_atoms * 2`, packed as:

```
[ elnum_0, label_0, elnum_1, label_1, ... ]
```

- `elnum`: atomic number, 1 to 104 (e.g. 6 for carbon, 8 for oxygen); 104
  is the Lennard-Jones pseudo-element "LJ". There is no dummy element, so
  leave dummy atoms out of the arrays; other values give
  `MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER`.
- `label`: user-defined integer label; pass `0` for unlabelled atoms. With
  `useatomtype_flag = true`, only atoms with the same label are matched.

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

With `mirror_flag = true`, *R* includes the reflection of molecule 2
(x → −x) as well as the rotation, so its determinant is −1. When
`align_flag = false`, *R* is the identity (or just the reflection when
`mirror_flag = true`) and *t* is zero.

### Multiple ranked solutions

Both functions accept a `max_records` parameter requesting up to that many
ranked candidate solutions (best RMSD first) instead of just the single
best one. `n_records` reports how many were actually found; it can be
smaller than `max_records`, and it is always `1` unless both `align_flag`
and `remap_flag` are true. With `conformsd_calculate` it is also `1` when the
exhaustive search is selected (see [conformsd_calculate](#conformsd-1)). All output arrays (`rmsd_list`, `mapping_list`,
`transform_list`) are flattened and must be pre-allocated by the caller
for `max_records` records; only the first `n_records` entries are
meaningful.

### Atom permutations and padding

Each permutation has `n_padding = max(n_atoms1, n_atoms2)` entries, so
`mapping_list` needs `max_records*n_padding` elements, and record `i`
occupies `mapping_list[i*n_padding .. i*n_padding+n_padding-1]`. Entry `j`
(0-based) is the atom of molecule 2 placed on line `j` of molecule 1. The
sizes can only differ when `heavy_flag = true`, since only the compared
atoms must match: the smaller molecule is then padded with padding atoms
appended after its real atoms, so values `>= n_atoms2` denote padding atoms
of molecule 2 and entries `j >= n_atoms1` hold the extra atoms of molecule
2. Atoms that are not part of the RMSD (hydrogens with `heavy_flag = true`)
are paired afterwards: with an atom of the same element, following their bonded heavy
atom in `conformsd_calculate` and by distance otherwise, and then with
whatever atoms are left.

### atormsd

```c
void atormsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool massweight_flag,
    bool mirror_flag, bool useatomtype_flag,
    bool printstats_flag, bool random_flag,
    bool pruning_flag, double prune_tol, int ato_freq, int max_trials,
    int max_records,
    double *rmsd_list, int *mapping_list,
    double *transform_list, int *n_records, int *error_code);
```

With alignment and remapping, random orientations of molecule 2 are each
refined by alternating atom assignment (within atom types) and optimal
superposition, and the distinct local minima are ranked.

Pruning speeds up the search by never pairing two atoms whose environments
are incompatible: their sorted distances to the atoms of some atom type
differ by more than 2√3 × `prune_tol`. When `pruning_flag = true`,
`prune_tol` (Å) has no default and must be supplied; it is ignored when
`pruning_flag = false`.

#### Parameters

| Parameter | Direction | Description |
|-----------|-----------|-------------|
| **n_atoms1** | in | Number of atoms in molecule 1 |
| **atom_data1** | in | Packed atom data for molecule 1, length `n_atoms1*2` |
| **coords1** | in | Coordinates for molecule 1, length `n_atoms1*3` |
| **n_atoms2** | in | Number of atoms in molecule 2 |
| **atom_data2** | in | Packed atom data for molecule 2, length `n_atoms2*2` |
| **coords2** | in | Coordinates for molecule 2, length `n_atoms2*3` |
| **align_flag** | in | Optimally rotate and translate molecule 2 onto molecule 1 |
| **remap_flag** | in | Find the atom permutation that minimises the RMSD (otherwise the input order is kept) |
| **heavy_flag** | in | Exclude hydrogen atoms from the RMSD |
| **massweight_flag** | in | Weight atoms by atomic mass |
| **mirror_flag** | in | Reflect molecule 2 (x → −x) before comparison |
| **useatomtype_flag** | in | Only match atoms with the same label |
| **printstats_flag** | in | Print optimisation statistics to stdout |
| **random_flag** | in | Seed the random-number generator from the system clock (otherwise results are reproducible) |
| **pruning_flag** | in | Enable pruning of atom pairs (see above) |
| **prune_tol** | in | Pruning tolerance (Å); required when `pruning_flag = true`, ignored otherwise (no default) |
| **ato_freq** | in | Stop once the best solution has been found this many times (≥ 1) |
| **max_trials** | in | Maximum number of random orientations (≥ 1) |
| **max_records** | in | Maximum number of ranked candidate solutions to return (≥ 1). Values > 1 only take effect when `align_flag` and `remap_flag` are both true |
| **rmsd_list** | out | RMSD (Å) of each returned record, length `max_records`; caller allocates |
| **mapping_list** | out | Flattened, **0-based** atom permutations, length `max_records*n_padding` with `n_padding = max(n_atoms1, n_atoms2)` (see [Atom permutations and padding](#atom-permutations-and-padding)); caller allocates ≥ `max_records*n_padding` elements |
| **transform_list** | out | Flattened row-major 4 × 4 homogeneous transforms, length `max_records*16`. Record `i` occupies `transform_list[i*16 .. i*16+15]`; caller allocates ≥ `max_records*16` elements |
| **n_records** | out | Actual number of records written (≤ `max_records`) |
| **error_code** | out | A constant from `enum molalign_error_code` in `molalign.h`: `MOLALIGN_SUCCESS`; `MOLALIGN_ERROR_INVALID_BOUND` (`ato_freq`, `max_trials` or `max_records` less than 1; checked before any output is written); `MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER` (an `elnum` outside 1–104); `MOLALIGN_ERROR_NOT_ISOMERS` (molecules are not isomers); `MOLALIGN_ERROR_ATOM_TYPE_MISMATCH` (atom types do not match; only when `remap_flag = false`); `MOLALIGN_ERROR_ASSIGNMENT_FAILED` (assignment failed, e.g. `prune_tol` too tight; only when `remap_flag = true`) |

#### Minimal example

```c
#include <stdio.h>
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
    int    mapping_list[3], n_records, error_code;   /* n_padding = 3 */

    atormsd_calculate(
        3, ad1, xy1,
        3, ad2, xy2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false,
        /*stats=*/false, /*random=*/false,
        /*pruning_flag=*/false, /*prune_tol=*/0.0, /*ato_freq=*/10, /*max_trials=*/10000,
        /*max_records=*/1,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd_list[0]);
    return 0;
}
```

To retrieve several ranked solutions, allocate the output arrays for
`max_records` records and loop over the first `n_records` entries:

```c
#define N_RECORDS 5

double rmsd_list[N_RECORDS];
double transform_list[N_RECORDS * 16];
int    mapping_list[N_RECORDS * 3];  /* N_RECORDS * n_padding */
int    n_records, error_code;

atormsd_calculate(
    3, ad1, xy1, 3, ad2, xy2,
    /*align=*/true, /*remap=*/true,
    /*heavy=*/false, /*mass=*/false,
    /*mirror=*/false, /*label=*/false,
    /*stats=*/false, /*random=*/false,
    /*pruning_flag=*/false, /*prune_tol=*/0.0, /*ato_freq=*/10, /*max_trials=*/10000,
    /*max_records=*/N_RECORDS,
    rmsd_list, mapping_list,
    transform_list, &n_records, &error_code);

for (int i = 0; i < n_records; i++) {
    printf("Solution %d: RMSD = %.6f Å\n", i, rmsd_list[i]);
}
```

Compile:

```bash
gcc example.c -o example -I/usr/local/include/molalignlib -L/usr/local/lib \
    -lmolalign -lgfortran -lm
```

Replace `/usr/local` with the install prefix (e.g. `$HOME/.local` for a
[user-local install](#user-local-install)), or with the build tree when
using the library before installing it.


### conformsd

```c
void conformsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_bonds1, const int *bond_data1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    int n_bonds2, const int *bond_data2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool massweight_flag,
    bool mirror_flag, bool useatomtype_flag, bool bonding_flag, double bond_tol,
    bool usebondtype_flag,
    bool printstats_flag, bool printassigntree_flag, bool random_flag,
    int confo_freq, int max_trials, int max_fragments,
    int max_records,
    double *rmsd_list, int *mapping_list,
    double *transform_list, int *n_records, int *error_code);
```

Atom assignments always respect the bond topology: they are searched over
the assignment tree of the Hierarchical Neighborhood of Atoms (HNA)
partition, so symmetry-equivalent atoms are permuted without ever breaking
a bond. With alignment and remapping, the search strategy is chosen from
the assignment tree: random orientation sampling when the total number of
complete assignments exceeds `confo_freq` times the sum of partial
combinations, and exhaustive enumeration otherwise (see
[Algorithm Notes](#algorithm-notes)). The exhaustive search returns a single
record, whatever `max_records` is.

Bond data is passed as a flat `int` array of length `n_bonds * 3`, packed as:

```
[ atom1_0, atom2_0, type_0, atom1_1, atom2_1, type_1, ... ]
```

Atom indices are **1-based**. When `bonding_flag = true` the library derives
connectivity from atomic geometry and the bond arrays may be empty
(`n_bonds = 0`, `bond_data = NULL`): two atoms are bonded when closer than
the sum of their covalent radii plus `bond_tol` (Å). `bond_tol` then has no
default and must be supplied; it is ignored when `bonding_flag = false`.

Each molecule may have at most `max_fragments` molecular fragments: connected
components of the bond graph of the compared atoms (after `heavy_flag`
exclusions). An atom without bonds is a fragment of its own, so with
`max_fragments = 1` a molecule of two or more atoms without bonds gives
`MOLALIGN_ERROR_TOO_MANY_FRAGMENTS`.

Every entry of the bond array is a bond. The bond `type` is the integer
code of a bond type label (see [Bond types](#bond-types)). Types are only
used when `usebondtype_flag = true`, and then the labels are compared
literally, never interpreted, so the same kind of bond must have the same
label in both molecules. Without it the type is ignored and any value is
accepted. Use the conventional labels that the file
readers use:

| Label | Meaning | Code |
|-------|---------|------|
| `1`, `2`, `3` | single, double, triple | 1, 2, 3 |
| `ar` | aromatic | 37 |
| `am` | amide | 32 |
| `3/2` | fractional order 3/2 | 978 |
| `3:2` | three-center two-electron bond | 1078 |

In general, a digit `d` has code `d`; a letter and a letter-or-digit `ab`
have code `10 + 36*a + b` (`a` = 0–25 for a–z, `b` = 0–35 for 0–9, a–z);
and a digit, separator and digit `nsd` have code `946 + 100*s + 10*n + d`
(`s` = 0–4 for `/ : - . ,`), so valid codes are 1–1445. Code 1446 is not
the code of any label: it marks a bond of undefined type, such as a bond
perceived from geometry or read as Mol2 `un` or MOL/SDF type 8. With
`usebondtype_flag = true`, a bond of undefined type between compared atoms
gives `MOLALIGN_ERROR_UNDEFINED_BOND_TYPE`, because it cannot be compared,
and any other invalid code gives `MOLALIGN_ERROR_INVALID_BOND_TYPE`. Note that the same molecule encoded with different
Kekulé or aromatic bond types is reported as
`MOLALIGN_ERROR_NOT_CONFORMERS` (or `MOLALIGN_ERROR_BOND_MISMATCH` when
`remap_flag = false`).

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
| **align_flag** | in | Optimally rotate and translate molecule 2 onto molecule 1 |
| **remap_flag** | in | Find the atom permutation that minimises the RMSD (otherwise the input order is kept) |
| **heavy_flag** | in | Exclude hydrogen atoms from the RMSD |
| **massweight_flag** | in | Weight atoms by atomic mass |
| **mirror_flag** | in | Reflect molecule 2 (x → −x) before comparison |
| **useatomtype_flag** | in | Only match atoms with the same label |
| **bonding_flag** | in | Derive connectivity from geometry (ignores bond arrays) |
| **bond_tol** | in | Bond detection tolerance (Å), see above; required when `bonding_flag = true`, ignored otherwise (no default) |
| **usebondtype_flag** | in | Use bond types to guide atom matching (labels compared literally, never interpreted; see above). Ignored when `bonding_flag = true` |
| **printstats_flag** | in | Print optimisation statistics to stdout |
| **printassigntree_flag** | in | Print the assignment tree and its combination counts to stdout (only when `remap_flag = true`) |
| **random_flag** | in | Seed the random-number generator from the system clock (otherwise results are reproducible) |
| **confo_freq** | in | Stop the random orientation search once the best solution has been found more than this many times; also the strategy threshold (see above). 100 is the validated value (≥ 1) |
| **max_trials** | in | Maximum number of random orientations (≥ 1) |
| **max_fragments** | in | Maximum number of molecular fragments allowed in each molecule (≥ 1; see above). 1 is the usual choice |
| **max_records** | in | Maximum number of ranked candidate solutions to return (≥ 1). Values > 1 only take effect when `align_flag` and `remap_flag` are both true and the stochastic search is selected |
| **rmsd_list** | out | RMSD (Å) of each returned record, length `max_records`; caller allocates |
| **mapping_list** | out | Flattened, **0-based** atom permutations, length `max_records*n_padding` with `n_padding = max(n_atoms1, n_atoms2)` (see [Atom permutations and padding](#atom-permutations-and-padding)); caller allocates ≥ `max_records*n_padding` elements |
| **transform_list** | out | Flattened row-major 4 × 4 homogeneous transforms, length `max_records*16`. Record `i` occupies `transform_list[i*16 .. i*16+15]`; caller allocates ≥ `max_records*16` elements |
| **n_records** | out | Actual number of records written (≤ `max_records`) |
| **error_code** | out | A constant from `enum molalign_error_code` in `molalign.h`: `MOLALIGN_SUCCESS`; `MOLALIGN_ERROR_INVALID_BOUND` (`confo_freq`, `max_trials`, `max_fragments` or `max_records` less than 1; checked before any output is written); `MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER` (an `elnum` outside 1–104); `MOLALIGN_ERROR_NOT_ISOMERS` (molecules are not isomers); `MOLALIGN_ERROR_INVALID_BOND_TYPE` (a bond type code that is not the code of any label) and `MOLALIGN_ERROR_UNDEFINED_BOND_TYPE` (a bond of undefined type), both only when `usebondtype_flag = true` and only for bonds between compared atoms; `MOLALIGN_ERROR_TOO_MANY_FRAGMENTS` (a molecule has more than `max_fragments` fragments); `MOLALIGN_ERROR_ATOM_TYPE_MISMATCH` (atom types do not match; only when `remap_flag = false`); `MOLALIGN_ERROR_BOND_MISMATCH` (bond connectivity does not match; only when `remap_flag = false`); `MOLALIGN_ERROR_NOT_CONFORMERS` (same composition but different connectivity; only when `remap_flag = true`) |

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
    int    mapping_list[4], n_records, error_code;   /* n_padding = 4 */

    conformsd_calculate(
        4, ad1, xy1, 3, bd1,
        4, ad2, xy2, 3, bd2,
        /*align=*/true, /*remap=*/true,
        /*heavy=*/false, /*mass=*/false,
        /*mirror=*/false, /*label=*/false,
        /*bonding_flag=*/false, /*bond_tol=*/0.0,
        /*usebondtype_flag=*/false,
        /*stats=*/false, /*printassigntree_flag=*/false, /*random=*/false,
        /*confo_freq=*/100, /*max_trials=*/10000, /*max_fragments=*/1,
        /*max_records=*/1,
        rmsd_list, mapping_list,
        transform_list, &n_records, &error_code);

    if (error_code != MOLALIGN_SUCCESS) { fprintf(stderr, "Error %d\n", error_code); return 1; }
    printf("RMSD = %.6f Å\n", rmsd_list[0]);
    return 0;
}
```

As with `atormsd_calculate`, pass `max_records > 1` and size the output
arrays accordingly to retrieve several ranked solutions in one call; loop
over the first `n_records` entries. Compile it the same way.

### Demo program

The demo program `bindings/atormsd_demo.c` is a C version of the
[atormsd](#atormsd) program built on `atormsd_calculate`. It takes the same
options, with the same defaults, and prints the same output (it also skips
dummy atoms), but only reads and writes XYZ files; it also accepts `-help`. It shows how to pack the atom
data and coordinates read from the files, size the output buffers for
several records and padded permutations, and apply the returned transforms
(`-aligned`). Compile and run instructions are in its header comment.

Python bindings
----------

The Python bindings provide a higher-level interface to MolAlignLib. (See [Try It Online](#try-it-online) above to test it on Binder without installing anything.)

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

### Molecule

The Python wrapper has a single structure class, `Molecule`: a set of atoms with
3-D coordinates and an optional bond table. The same object can be compared
in two ways, one method per algorithm of the library:

| Method | Wraps | Bonds |
|--------|-------|-------|
| [`atormsd_to`](#atormsd_to) | `atormsd_calculate` | Ignored; atoms are matched within atom types only. Useful for metal clusters, nanoparticles, or other systems where connectivity is absent or irrelevant |
| [`conformsd_to`](#conformsd_to) | `conformsd_calculate` | Required (from the bond table, or inferred with `bond_tol`); by default each molecule must be a single connected fragment (see `max_fragments`). HNA partitioning guarantees chemically valid assignments, handling arbitrary degrees of topological symmetry efficiently |

```python
from molalignlib import Molecule, read_molecules
```

#### Construction

```python
# From a file (defaults to the first frame; bonds are read when the format has
# them, and dummy atoms "X" are skipped)
mol = Molecule.from_file("conformers.sdf")

# From a file (specific frame in a multi-frame file)
mol = Molecule.from_file("conformers.sdf", frame_idx=2)

# From element symbols and coordinates directly (no bonds)
import numpy as np
symbols = ["Fe", "Fe"]
coords  = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], dtype=np.float64)
mol = Molecule.from_symbols(symbols, coords, name="dimer")

# From element symbols, coordinates, and bonds: [atom1, atom2, type], with
# 1-based atom indices and bond type codes (see Bond types below)
symbols   = ["C", "O", "H", "H"]
coords    = np.zeros((4, 3), dtype=np.float64)
bond_data = np.array([[1,2,2],[1,3,1],[1,4,1]], dtype=np.int32)  # C=O, C-H, C-H
mol = Molecule.from_symbols(symbols, coords, bond_data=bond_data)

# From atomic numbers, coordinates, and (optionally) bonds
mol = Molecule.from_numbers([6, 8, 1, 1], coords, bond_data=bond_data)
```

`from_symbols` and `from_numbers` accept `labels=` (one integer per atom,
default 0) for use with `use_atom_type=True`, and `name=` (default
`"molecule"`). Molecules read from files are unlabelled and named after the
file stem.

#### Bond types

Bond types are opaque labels (`"1"`, `"ar"`, `"3/2"`, ...; see
[Bond types](#bond-types) in the conformsd section), stored in the third
column of `bond_data` as integer codes. `bond_type_code` and
`bond_type_label` convert between the two:

```python
from molalignlib import bond_type_code, bond_type_label

bond_type_code("ar")    # 37
bond_type_label(978)    # '3/2'
```

Bonds read from files get the conventional labels of their chemfiles bond
orders (`"1"` to `"5"`, `"am"`, `"ar"`; bonds without an order are of
undefined type), so molecules read from different formats can be compared
with `use_bond_type=True`. Molecules read from files also record the file format
as their `bond_source`; it is informational only.

#### Reading multiple frames

```python
molecules = read_molecules("conformers.sdf")            # all frames, returns a list
mol0, mol1 = read_molecules("conformers.sdf", frames=(0, 1))
```

#### Properties

| Property | Description |
|----------|-------------|
| **name** | Name (file stem plus frame index for molecules read with `read_molecules`) |
| **n_atoms** | Number of atoms (also `len(mol)`) |
| **n_bonds**, **has_bonds** | Number of bonds, and whether there are any |
| **symbols** | Element symbols |
| **atom_data** | `int32 (n_atoms, 2)` array of `[atomic_number, label]` (atomic numbers 1–104) |
| **coords** | `float64 (n_atoms, 3)` coordinates in Å |
| **bond_data** | `int32 (n_bonds, 3)` array of `[atom1, atom2, type]`, 1-based, with bond type codes (empty when there are no bonds) |
| **bond_source** | Where the bonds came from (file format, or `None`); informational |

#### Writing output

`write()` uses chemfiles, so it supports the same range of formats as
reading (XYZ, PDB, SDF, MOL2, ...); the output format is inferred from
the file extension. Bond connectivity is written too when the molecule has
bonds and the target format supports it (e.g. SDF/MOL2); bond *types* are
not currently preserved, only which atoms are bonded. Writing to a format
with no bond table (e.g. XYZ) simply omits connectivity:

```python
mol.write("output.sdf", comment="my molecule")
mol.write("output.xyz")  # coordinates and symbols only
```

### atormsd_to

Compares the two molecules as unstructured atom clusters. Bonds are
ignored, so it works the same for molecules with or without a bond table.

```python
mol0, mol1 = read_molecules("clusters.xyz", frames=(0, 1))

results = mol0.atormsd_to(mol1, align=True, remap=True)
result = results[0]   # best (lowest RMSD) solution
print(result.rmsd)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| **other** | `Molecule` | *required* | The structure to compare against `self` |
| **align** | `bool` | `False` | Optimally rotate and translate `other` onto `self` |
| **remap** | `bool` | `False` | Find the atom permutation that minimises the RMSD |
| **heavy_only** | `bool` | `False` | Exclude hydrogen atoms from the calculation |
| **mass_weight** | `bool` | `False` | Weight each atom by its atomic mass |
| **mirror** | `bool` | `False` | Reflect `other` before comparison |
| **use_atom_type** | `bool` | `False` | Only match atoms with the same label (see `labels` in the constructors) |
| **print_stats** | `bool` | `False` | Print detailed optimisation statistics |
| **random** | `bool` | `False` | Seed the random-number generator from the system clock (otherwise results are reproducible) |
| **prune_tol** | `float` | `None` | Pruning tolerance (Å): two atoms are never paired if their sorted distances to the atoms of some atom type differ by more than 2√3 × `prune_tol`. `None` disables pruning |
| **ato_freq** | `int` | `10` | Stop searching once the best solution has been found this many times (≥ 1) |
| **max_trials** | `int` | `10000` | Stop after at most this many random orientations (≥ 1) |
| **max_records** | `int` | `1` | Request up to this many ranked solutions (≥ 1; see [Retrieving multiple ranked solutions](#retrieving-multiple-ranked-solutions)) |

Values below 1 for `ato_freq`, `max_trials` or `max_records` raise `ValueError`.

### conformsd_to

Compares the two molecules as conformers. Atom assignments always respect
the bond topology, so both molecules must have the same bond graph. Bonds
come from each molecule's bond table, or are inferred from geometry when
`bond_tol` is given. Each molecule may have at most `max_fragments` molecular
fragments (connected components of the bond graph, after `heavy_only`
exclusions), 1 by default. An atom without bonds is a fragment of its own,
so a molecule without bonds (e.g. read from XYZ) needs `bond_tol`;
otherwise it exceeds the default `max_fragments` and a `ValueError` is raised.

```python
c0, c1 = read_molecules("conformers.sdf", frames=(0, 1))

results = c0.conformsd_to(c1, align=True, remap=True)
result = results[0]             # best (lowest RMSD) solution
print(result.rmsd)              # float, Å
print(result.mapping)           # int32 array, 0-based
print(result.transform)         # 4×4 float64 array

# XYZ input has no bond table: infer connectivity from geometry
x0, x1 = read_molecules("conformers.xyz", frames=(0, 1))
result = x0.conformsd_to(x1, align=True, remap=True, bond_tol=0.3)[0]
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| **other** | `Molecule` | *required* | The structure to compare against `self` |
| **align** | `bool` | `False` | Optimally rotate and translate `other` onto `self` |
| **remap** | `bool` | `False` | Find the atom permutation that minimises the RMSD |
| **heavy_only** | `bool` | `False` | Exclude hydrogen atoms from the calculation |
| **mass_weight** | `bool` | `False` | Weight each atom by its atomic mass |
| **mirror** | `bool` | `False` | Reflect `other` before comparison |
| **use_atom_type** | `bool` | `False` | Only match atoms with the same label (see `labels` in the constructors) |
| **bond_tol** | `float` | `None` | Bond-detection tolerance (Å): infer bond connectivity from geometry instead of using each molecule's bond table. `None` uses the bond tables |
| **use_bond_type** | `bool` | `False` | Use bond types to guide atom matching. Labels are compared literally (see [Bond types](#bond-types-1)); bonds of undefined type between compared atoms raise `ValueError`. No effect when `bond_tol` is given |
| **print_stats** | `bool` | `False` | Print detailed optimisation statistics |
| **print_assignment_tree** | `bool` | `False` | Print the assignment tree and its combination counts |
| **random** | `bool` | `False` | Seed the random-number generator from the system clock (otherwise results are reproducible) |
| **confo_freq** | `int` | `100` | Stop the random orientation search once the best solution has been found more than this many times; also the threshold that selects it over exhaustive enumeration (see [Algorithm Notes](#algorithm-notes)). 100 is the validated value (≥ 1) |
| **max_trials** | `int` | `10000` | Stop after at most this many random orientations (≥ 1) |
| **max_fragments** | `int` | `1` | Maximum number of molecular fragments allowed in each molecule (≥ 1); raise it to compare complexes or other multi-molecule systems |
| **max_records** | `int` | `1` | Request up to this many ranked solutions (≥ 1; see [Retrieving multiple ranked solutions](#retrieving-multiple-ranked-solutions)) |

Values below 1 for `confo_freq`, `max_trials`, `max_fragments` or `max_records`
raise `ValueError`.

### RMSDResult

`atormsd_to()` and `conformsd_to()` both return a `list[RMSDResult]`, one
element per solution found (best/lowest RMSD first). Each `RMSDResult`
holds:

| Attribute | Type | Description |
|-----------|------|-------------|
| **rmsd** | `float` | Root-mean-square deviation in Å |
| **mapping** | `int32 ndarray (n_padding,)` | 0-based index array: entry `j` is the atom of *other* placed on line `j` of *self*. `n_padding = max(len(self), len(other))`; the lengths can only differ with `heavy_only=True`, and then values `>= len(other)` denote padding atoms of *other* and entries `j >= len(self)` hold its extra atoms |
| **transform** | `float64 ndarray (4, 4)` | Homogeneous rotation + translation matrix (maps *other* to *self* frame); includes the reflection when `mirror=True` |

#### Applying the result

Use the `apply_to` method to produce a new `Molecule` that is aligned and
reordered to match the reference. Its bonds, if any, are renumbered to the
new atom order, so the result can be written with connectivity. If the
reference has more atoms (only possible with `heavy_only=True`), the missing
lines are filled with padding atoms (`"X"`), which can be written but not
compared:

```python
result = mol0.conformsd_to(mol1, align=True, remap=True)[0]
mol1_aligned = result.apply_to(mol1)
mol1_aligned.write("aligned.sdf")
```

#### Retrieving multiple ranked solutions

Set `max_records` to inspect several distinct candidate mappings instead of
just the best one. This only produces more than one result when both
`align=True` and `remap=True` (and, for `conformsd_to`, when the stochastic
search is selected); the returned list is otherwise always length 1
regardless of `max_records`, and may be shorter than `max_records` if fewer
distinct solutions were found:

```python
results = mol0.atormsd_to(mol1, align=True, remap=True, max_records=5)

for i, result in enumerate(results):
    print(f"Solution {i}: RMSD = {result.rmsd:.4f} Å")

best = results[0]
```


Runnable examples
-----------------

For fully runnable examples see
[`python/examples/examples.py`](python/examples/examples.py)
and the equivalent Jupyter notebook
[`python/examples/examples.ipynb`](python/examples/examples.ipynb)
(the one launched by the [Binder link](#try-it-online) above).


Algorithm Notes
---------------

- **atoRMSD:** uses a stochastic strategy with distance-based pruning for unstructured clusters where no bond topology is available. Algorithm described in [Vásquez-Pérez et al., *J. Chem. Inf. Model.* (2023)](https://doi.org/10.1021/acs.jcim.2c01187).
- **confoRMSD:** uses the Hierarchical Neighborhood of Atoms (HNA) partitioning to decompose the assignment problem into independent branches, reducing the number of evaluated combinations from the product of branch possibilities to their sum. For alignment calculations, the algorithm adaptively selects between stochastic orientation sampling (efficient for highly symmetric molecules) and exhaustive enumeration (efficient for molecules with few branches), based on the ratio of total to partial combinations in the assignment tree. Benchmarks show 100% topologically correct assignments across 1.4 million molecular pairs with millisecond-scale mean execution times. Full algorithm description in [Vásquez-Pérez et al., *J. Chem. Theory Comput.* (2026)](https://doi.org/10.1021/acs.jctc.6c00545).
- **Transform output:** the 4 × 4 homogeneous transformation matrix encodes both the optimal rotation *R* and the translation *t* needed to superimpose molecule 2 on molecule 1. When mirroring, *R* also includes the reflection, so the matrix maps the original (unmirrored) coordinates of molecule 2.
- **Default maximum trials:** 10,000 random orientations for alignment. The convergence frequency threshold defaults to 10 for `atormsd` and 100 for `conformsd`; numerical experiments show that values below 100 can produce incorrect assignments for conformers.
