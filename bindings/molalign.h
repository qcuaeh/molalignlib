#ifndef MOLALIGN_C_BINDING_H
#define MOLALIGN_C_BINDING_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 *                       with atomlabel_flag=true.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 *
 * @param align_flag    Optimally rotate and translate molecule 2 onto molecule 1
 * @param remap_flag    Find the atom permutation that minimises the RMSD
 *                      (otherwise the input order is kept)
 * @param heavy_flag    Exclude hydrogen atoms from the RMSD
 * @param mass_flag     Weight atoms by their atomic masses
 * @param mirror_flag   Mirror second molecule (reflect it on the yz plane,
 *                      x -> -x) before comparing; transform_list includes
 *                      the reflection
 * @param atomlabel_flag Only match atoms with the same label
 * @param print_stats   Print optimisation statistics to stdout
 * @param random_flag   Seed the random number generator from the clock
 *                      (otherwise results are reproducible)
 * @param prunetol_flag Enable distance-based pruning of atom pairs
 * @param prune_tol     Pruning tolerance (Angstrom): two atoms are never
 *                      paired if their sorted distances to the atoms of some
 *                      atom type differ by more than 2*sqrt(3)*prune_tol.
 *                      Required (no default) when prunetol_flag=true;
 *                      unused otherwise.
 * @param conv_freq     Stop the random search once the best solution has
 *                      been found this many times
 * @param max_trials    Maximum number of random orientations
 *
 * @param n_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see occ_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length n_records.
 *                            Only the first *occ_records entries are valid.
 * @param mapping_list [out] Flattened, 0-based atom permutations, length
 *                            n_records*n_padding, n_padding = max(n_atoms1, n_atoms2).
 *                            Record i (0-based) occupies
 *                            mapping_list[i*n_padding .. i*n_padding + n_padding - 1].
 *                            Same padding convention as conformsd_calculate:
 *                            values >= n_atoms2 are dummy atoms of cluster 2.
 *                            Sizes can only differ when heavy_flag=true.
 *                            Caller allocates >= n_records*n_padding.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length n_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            the input molecule-2 coords to the molecule-1
 *                            frame: p_out = R*p_in + t. With mirror_flag=true,
 *                            R includes the reflection (det R = -1). With
 *                            align_flag=false, R is the identity, or just the
 *                            reflection when mirror_flag=true, and t = 0.
 *                            Caller allocates >= n_records*16.
 * @param occ_records   [out] Actual number of records written (<= n_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (error_codes.h): MOLALIGN_SUCCESS,
 *                            MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER,
 *                            MOLALIGN_ERROR_NOT_ISOMERS,
 *                            MOLALIGN_ERROR_ATOM_TYPE_MISMATCH (only when
 *                            remap_flag=false), or
 *                            MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED (only when
 *                            remap_flag=true).
 */
void atormsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool atomlabel_flag,
    bool print_stats, bool random_flag,
    bool prunetol_flag, double prune_tol, int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *mapping_list,
    double *transform_list, int *occ_records, int *error_code);

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
 *                       with atomlabel_flag=true.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_bonds1      Number of bonds in molecule 1 (may be 0 when bondtol_flag=true)
 * @param bond_data1    Flat bond array [a1,a2,type,...], 1-based, length n_bonds1*3.
 *                      type is any integer bond-type code; it is only used
 *                      when bondtype_flag=true, and then only compared for
 *                      equality, never interpreted.
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 * @param n_bonds2      Number of bonds in molecule 2 (may be 0 when bondtol_flag=true)
 * @param bond_data2    Flat bond array [a1,a2,type,...], 1-based, length n_bonds2*3
 *
 * @param align_flag    Optimally rotate and translate molecule 2 onto molecule 1
 * @param remap_flag    Find the atom permutation that minimises the RMSD
 *                      (otherwise the input order is kept)
 * @param heavy_flag    Exclude hydrogen atoms from the RMSD
 * @param mass_flag     Weight atoms by their atomic masses
 * @param mirror_flag   Mirror second molecule (reflect it on the yz plane,
 *                      x -> -x) before comparing; transform_list includes
 *                      the reflection
 * @param atomlabel_flag Only match atoms with the same label
 * @param bondtol_flag  Derive connectivity from geometry, not bond table
 * @param bond_tol      Bond detection tolerance (Angstrom): atoms are bonded
 *                      when closer than the sum of their covalent radii plus
 *                      bond_tol. Required (no default) when bondtol_flag=true;
 *                      unused otherwise.
 * @param bondtype_flag Use bond types to guide atom matching: bonds of
 *                      different type are distinguished in the HNA
 *                      partition. Types are compared, not interpreted, so
 *                      both molecules must use the same bond-type
 *                      convention (e.g. read from the same file format with
 *                      the same parser). Differing Kekule/aromatic encodings
 *                      of the same molecule yield
 *                      MOLALIGN_ERROR_NOT_CONFORMERS, or
 *                      MOLALIGN_ERROR_BOND_MISMATCH when remap_flag=false.
 *                      Ignored when bondtol_flag=true.
 * @param print_stats   Print optimisation statistics to stdout
 * @param print_assigntree Print the assignment tree and its combination
 *                      counts to stdout
 * @param random_flag   Seed the random number generator from the clock
 *                      (otherwise results are reproducible)
 * @param conv_freq     Stop the random orientation search once the best
 *                      solution has been found more than this many times;
 *                      also the threshold on the ratio of total to partial
 *                      assignment combinations above which that search is
 *                      used instead of exhaustive enumeration. 100 is the
 *                      value validated on the CCD and BIRD benchmarks.
 * @param max_trials    Maximum number of random orientations
 *
 * @param n_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see occ_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length n_records.
 *                            Only the first *occ_records entries are valid.
 * @param mapping_list [out] Flattened, 0-based atom permutations, length
 *                            n_records*n_padding, where
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
 *                            Caller allocates >= n_records*n_padding.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length n_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            the input molecule-2 coords to the molecule-1
 *                            frame: p_out = R*p_in + t. With mirror_flag=true,
 *                            R includes the reflection (det R = -1). With
 *                            align_flag=false, R is the identity, or just the
 *                            reflection when mirror_flag=true, and t = 0.
 *                            Caller allocates >= n_records*16.
 * @param occ_records   [out] Actual number of records written (<= n_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (error_codes.h): MOLALIGN_SUCCESS,
 *                            MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER,
 *                            MOLALIGN_ERROR_NOT_ISOMERS,
 *                            MOLALIGN_ERROR_MISSING_BONDS,
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
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool atomlabel_flag, bool bondtol_flag, double bond_tol,
    bool bondtype_flag,
    bool print_stats, bool print_assigntree, bool random_flag,
    int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *mapping_list,
    double *transform_list, int *occ_records, int *error_code);

#ifdef __cplusplus
}
#endif

#endif /* MOLALIGN_C_BINDING_H */