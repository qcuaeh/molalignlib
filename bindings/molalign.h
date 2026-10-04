#ifndef MOLALIGN_C_BINDING_H
#define MOLALIGN_C_BINDING_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Error codes returned by atormsd_calculate and conformsd_calculate, and by
 * the Fortran isormsd_calculate, through their error_code argument.
 *
 * Hand-maintained mirror of the Fortran error_codes module
 * (error_codes.f90): keep the values of both in sync.
 */
enum molalign_error_code {
    /* No error. */
    MOLALIGN_SUCCESS                        = 0,
    /* Different atom counts or compositions. */
    MOLALIGN_ERROR_NOT_ISOMERS              = 1,
    /* Atom types differ in input order. */
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH       = 2,
    /* More molecular fragments than max_fragments. */
    MOLALIGN_ERROR_TOO_MANY_FRAGMENTS       = 3,
    /* Bonds differ in input order. */
    MOLALIGN_ERROR_BOND_MISMATCH            = 4,
    /* Same composition but non-isomorphic bond graphs. */
    MOLALIGN_ERROR_NOT_CONFORMERS           = 5,
    /* No valid assignment under the pruning constraints. */
    MOLALIGN_ERROR_ASSIGNMENT_FAILED = 6,
    /* Atomic number outside the element tables. */
    MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER    = 7,
    /* Count parameter less than 1. */
    MOLALIGN_ERROR_INVALID_BOUND        = 8
};

/**
 * Calculate the RMSD between two atom clusters (no bond topology), optionally
 * returning several ranked candidate solutions instead of just the best one.
 * With alignment and remapping, random orientations are each refined by
 * alternating atom assignment (within atom types) and superposition, and the
 * distinct local minima are ranked.
 *
 * Molecule data is passed as flat arrays; the library performs no file I/O.
 *
 * @param n_atoms1      Number of atoms in molecule 1
 * @param atom_data1    Packed atom data, length n_atoms1*2:
 *                       [elnum0, label0, elnum1, label1, ...]
 *                       elnum is the atomic number, 0 to 104 (0 = dummy
 *                       atom "X", 104 = Lennard-Jones "LJ"); dummy atoms
 *                       are always excluded from the comparison.
 *                       label = 0 means unlabelled; labels only matter
 *                       with useatomtype_flag=true.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 *
 * @param align_flag    Optimally rotate and translate molecule 2 onto molecule 1
 * @param remap_flag    Find the atom permutation that minimises the RMSD
 *                      (otherwise the input order is kept)
 * @param heavy_flag    Exclude hydrogen atoms from the RMSD
 * @param massweight_flag     Weight atoms by their atomic masses
 * @param mirror_flag   Mirror second molecule (reflect it on the yz plane,
 *                      x -> -x) before comparing; transform_list includes
 *                      the reflection
 * @param useatomtype_flag Only match atoms with the same label
 * @param printstats_flag   Print optimisation statistics to stdout
 * @param random_flag   Seed the random number generator from the clock
 *                      (otherwise results are reproducible)
 * @param pruning_flag Enable distance-based pruning of atom pairs
 * @param prune_tol     Pruning tolerance (Angstrom): two atoms are never
 *                      paired if their sorted distances to the atoms of some
 *                      atom type differ by more than 2*sqrt(3)*prune_tol.
 *                      Required (no default) when pruning_flag=true;
 *                      unused otherwise.
 * @param ato_freq     Stop the random search once the best solution has
 *                      been found this many times (>=1)
 * @param max_trials    Maximum number of random orientations (>=1)
 *
 * @param max_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see n_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length max_records.
 *                            Only the first *n_records entries are valid.
 * @param mapping_list [out] Flattened, 0-based atom permutations, length
 *                            max_records*n_padding, n_padding = max(n_atoms1, n_atoms2).
 *                            Record i (0-based) occupies
 *                            mapping_list[i*n_padding .. i*n_padding + n_padding - 1].
 *                            Same padding convention as conformsd_calculate:
 *                            values >= n_atoms2 are dummy atoms of cluster 2.
 *                            Sizes can only differ when heavy_flag=true.
 *                            Caller allocates >= max_records*n_padding.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length max_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            the input molecule-2 coords to the molecule-1
 *                            frame: p_out = R*p_in + t. With mirror_flag=true,
 *                            R includes the reflection (det R = -1). With
 *                            align_flag=false, R is the identity, or just the
 *                            reflection when mirror_flag=true, and t = 0.
 *                            Caller allocates >= max_records*16.
 * @param n_records   [out] Actual number of records written (<= max_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (above): MOLALIGN_SUCCESS,
 *                            MOLALIGN_ERROR_INVALID_BOUND (ato_freq,
 *                            max_trials or max_records less than 1),
 *                            MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER,
 *                            MOLALIGN_ERROR_NOT_ISOMERS,
 *                            MOLALIGN_ERROR_ATOM_TYPE_MISMATCH (only when
 *                            remap_flag=false), or
 *                            MOLALIGN_ERROR_ASSIGNMENT_FAILED (only when
 *                            remap_flag=true).
 */
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

/**
 * Calculate the symmetry-corrected RMSD between two molecular conformers,
 * optionally returning several ranked candidate solutions instead of just
 * the best one. Atom assignments always respect the bond topology: they are
 * searched over the assignment tree of the Hierarchical Neighborhood of Atoms
 * (HNA) partition (J. Chem. Theory Comput., doi:10.1021/acs.jctc.6c00545).
 * With alignment, the strategy is chosen from the shape of that tree:
 * stochastic orientation sampling, or exhaustive enumeration when the number
 * of complete assignments is small.
 *
 * Molecule data is passed as flat arrays; the library performs no file I/O.
 *
 * @param n_atoms1      Number of atoms in molecule 1
 * @param atom_data1    Packed atom data, length n_atoms1*2:
 *                       [elnum0, label0, elnum1, label1, ...]
 *                       elnum is the atomic number, 0 to 104 (0 = dummy
 *                       atom "X", 104 = Lennard-Jones "LJ"); dummy atoms
 *                       are always excluded from the comparison.
 *                       label = 0 means unlabelled; labels only matter
 *                       with useatomtype_flag=true.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_bonds1      Number of bonds in molecule 1 (may be 0 when bonding_flag=true)
 * @param bond_data1    Flat bond array [a1,a2,type,...], 1-based, length n_bonds1*3.
 *                      type is any integer bond-type code; it is only used
 *                      when usebondtype_flag=true, and then only compared for
 *                      equality, never interpreted.
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 * @param n_bonds2      Number of bonds in molecule 2 (may be 0 when bonding_flag=true)
 * @param bond_data2    Flat bond array [a1,a2,type,...], 1-based, length n_bonds2*3
 *
 * @param align_flag    Optimally rotate and translate molecule 2 onto molecule 1
 * @param remap_flag    Find the atom permutation that minimises the RMSD
 *                      (otherwise the input order is kept)
 * @param heavy_flag    Exclude hydrogen atoms from the RMSD
 * @param massweight_flag     Weight atoms by their atomic masses
 * @param mirror_flag   Mirror second molecule (reflect it on the yz plane,
 *                      x -> -x) before comparing; transform_list includes
 *                      the reflection
 * @param useatomtype_flag Only match atoms with the same label
 * @param bonding_flag  Derive connectivity from geometry, not bond table
 * @param bond_tol      Bond detection tolerance (Angstrom): atoms are bonded
 *                      when closer than the sum of their covalent radii plus
 *                      bond_tol. Required (no default) when bonding_flag=true;
 *                      unused otherwise.
 * @param usebondtype_flag Use bond types to guide atom matching: bonds of
 *                      different type are distinguished in the HNA
 *                      partition. Types are compared, not interpreted, so
 *                      both molecules must use the same bond-type
 *                      convention (e.g. read from the same file format with
 *                      the same parser). Differing Kekule/aromatic encodings
 *                      of the same molecule yield
 *                      MOLALIGN_ERROR_NOT_CONFORMERS, or
 *                      MOLALIGN_ERROR_BOND_MISMATCH when remap_flag=false.
 *                      Ignored when bonding_flag=true.
 * @param printstats_flag   Print optimisation statistics to stdout
 * @param printassigntree_flag Print the assignment tree and its combination
 *                      counts to stdout
 * @param random_flag   Seed the random number generator from the clock
 *                      (otherwise results are reproducible)
 * @param confo_freq     Stop the random orientation search once the best
 *                      solution has been found more than this many times;
 *                      also the threshold on the ratio of total to partial
 *                      assignment combinations above which that search is
 *                      used instead of exhaustive enumeration. 100 is the
 *                      value validated on the CCD and BIRD benchmarks.
 *                      Must be >=1.
 * @param max_trials    Maximum number of random orientations (>=1)
 * @param max_fragments     Maximum number of molecular fragments (connected
 *                      components of the bond graph of the compared atoms,
 *                      i.e. after heavy_flag exclusions) allowed in each
 *                      molecule (>=1; 1 is the usual choice). An atom
 *                      without bonds is a fragment of its own. Exceeding it
 *                      yields MOLALIGN_ERROR_TOO_MANY_FRAGMENTS.
 *
 * @param max_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see n_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length max_records.
 *                            Only the first *n_records entries are valid.
 * @param mapping_list [out] Flattened, 0-based atom permutations, length
 *                            max_records*n_padding, where
 *                            n_padding = max(n_atoms1, n_atoms2). Record i (0-based)
 *                            occupies mapping_list[i*n_padding .. i*n_padding + n_padding - 1]
 *                            and is a permutation of 0..n_padding-1: entry j is
 *                            the atom of molecule 2 placed on line j of
 *                            molecule 1. The smaller molecule is padded with
 *                            dummy atoms appended after its real atoms, so
 *                            values >= n_atoms2 denote dummy atoms of
 *                            molecule 2, and entries j >= n_atoms1 hold the
 *                            extra atoms of molecule 2. The sizes can only
 *                            differ when heavy_flag=true; otherwise
 *                            n_padding == n_atoms1.
 *                            Caller allocates >= max_records*n_padding.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length max_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            the input molecule-2 coords to the molecule-1
 *                            frame: p_out = R*p_in + t. With mirror_flag=true,
 *                            R includes the reflection (det R = -1). With
 *                            align_flag=false, R is the identity, or just the
 *                            reflection when mirror_flag=true, and t = 0.
 *                            Caller allocates >= max_records*16.
 * @param n_records   [out] Actual number of records written (<= max_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (above): MOLALIGN_SUCCESS,
 *                            MOLALIGN_ERROR_INVALID_BOUND (confo_freq,
 *                            max_trials, max_fragments or max_records less
 *                            than 1),
 *                            MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER,
 *                            MOLALIGN_ERROR_NOT_ISOMERS,
 *                            MOLALIGN_ERROR_TOO_MANY_FRAGMENTS,
 *                            MOLALIGN_ERROR_ATOM_TYPE_MISMATCH or
 *                            MOLALIGN_ERROR_BOND_MISMATCH (both only when
 *                            remap_flag=false), or
 *                            MOLALIGN_ERROR_NOT_CONFORMERS (only when
 *                            remap_flag=true).
 */
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

#ifdef __cplusplus
}
#endif

#endif /* MOLALIGN_C_BINDING_H */
