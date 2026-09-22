#ifndef MOLALIGN_C_BINDING_H
#define MOLALIGN_C_BINDING_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Calculate RMSD between two atom clusters, optionally returning several
 * ranked candidate solutions instead of just the best one.
 *
 * Molecule data is passed as pre-read flat arrays; no file I/O is performed
 * inside the Fortran library.
 *
 * @param n_atoms1      Number of atoms in molecule 1
 * @param atom_data1    Packed atom data, length n_atoms1*2:
 *                       [elnum0, label0, elnum1, label1, ...]
 *                       label = 0 means unlabelled.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 *
 * @param align_flag    Enable structural alignment
 * @param remap_flag    Enable atom remapping
 * @param heavy_flag    Use only heavy atoms
 * @param mass_flag     Use atomic masses as weights
 * @param mirror_flag   Mirror second molecule
 * @param label_flag    Use atom type labels for matching
 * @param print_stats    Print optimisation statistics
 * @param random_flag   Use random algorithm
 * @param prune_flag    Enable pruning
 * @param prune_tol     Pruning tolerance. Required (no default) when
 *                      prune_flag=true; unused otherwise.
 * @param conv_freq     Convergence frequency
 * @param max_trials    Maximum number of trials
 *
 * @param n_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see occ_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length n_records.
 *                            Only the first *occ_records entries are valid.
 * @param atomperm_list [out] Flattened, 0-based atom permutations, length
 *                            n_records*n_atoms1 (each permutation has one
 *                            entry per atom of molecule 1). Record i (0-based)
 *                            occupies atomperm_list[i*n_atoms1 .. i*n_atoms1 + n_atoms1 - 1].
 *                            Caller allocates >= n_records*n_atoms1.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length n_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            molecule-2 coords to molecule-1 frame:
 *                            p_out = R*p_in + t. Identity when align_flag=false.
 *                            Caller allocates >= n_records*16.
 * @param occ_records   [out] Actual number of records written (<= n_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (error_codes.h): MOLALIGN_SUCCESS,
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
    bool mirror_flag, bool label_flag,
    bool print_stats, bool random_flag,
    bool prune_flag, double prune_tol, int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *atomperm_list,
    double *transform_list, int *occ_records, int *error_code);

/**
 * Calculate RMSD between two molecular conformers, optionally returning
 * several ranked candidate solutions instead of just the best one.
 *
 * Molecule data is passed as pre-read flat arrays; no file I/O is performed
 * inside the Fortran library.
 *
 * @param n_atoms1      Number of atoms in molecule 1
 * @param atom_data1    Packed atom data, length n_atoms1*2:
 *                       [elnum0, label0, elnum1, label1, ...]
 *                       label = 0 means unlabelled.
 * @param coords1       Coordinates, row-major [x0,y0,z0,...], length n_atoms1*3
 * @param n_bonds1      Number of bonds in molecule 1 (may be 0 when bond_flag=true)
 * @param bond_data1    Flat bond array [a1,a2,type,...], 1-based, length n_bonds1*3
 * @param n_atoms2      Number of atoms in molecule 2
 * @param atom_data2    Packed atom data, length n_atoms2*2 (same layout)
 * @param coords2       Coordinates, row-major [x0,y0,z0,...], length n_atoms2*3
 * @param n_bonds2      Number of bonds in molecule 2 (may be 0 when bond_flag=true)
 * @param bond_data2    Flat bond array [a1,a2,type,...], 1-based, length n_bonds2*3
 *
 * @param align_flag    Enable structural alignment
 * @param remap_flag    Enable atom remapping
 * @param heavy_flag    Use only heavy atoms
 * @param mass_flag     Use atomic masses as weights
 * @param mirror_flag   Mirror second molecule
 * @param label_flag    Use atom type labels for matching
 * @param bond_flag     Derive connectivity from geometry, not bond table
 * @param bond_tol      Bond detection tolerance for geometry-based connectivity.
 *                      Required (no default) when bond_flag=true; unused otherwise.
 * @param print_stats    Print optimisation statistics
 * @param print_assigntree Print the atom-assignment search tree
 * @param random_flag   Use random algorithm
 * @param conv_freq   Convergence frequency
 * @param max_trials    Maximum number of trials
 *
 * @param n_records     Maximum number of ranked candidate solutions to
 *                       return (>=1). Multiple records are only ever
 *                       produced when both align_flag and remap_flag are
 *                       true; otherwise exactly one record is written
 *                       regardless of this value (see occ_records).
 *
 * @param rmsd_list     [out] RMSD of each returned record, length n_records.
 *                            Only the first *occ_records entries are valid.
 * @param atomperm_list [out] Flattened, 0-based atom permutations, length
 *                            n_records*n_atoms1 (each permutation has one
 *                            entry per atom of molecule 1). Record i (0-based)
 *                            occupies atomperm_list[i*n_atoms1 .. i*n_atoms1 + n_atoms1 - 1].
 *                            Caller allocates >= n_records*n_atoms1.
 * @param transform_list [out] Flattened row-major 4x4 homogeneous transforms,
 *                            length n_records*16. Record i occupies
 *                            transform_list[i*16 .. i*16 + 15]. Maps
 *                            molecule-2 coords to molecule-1 frame:
 *                            p_out = R*p_in + t. Identity when align_flag=false.
 *                            Caller allocates >= n_records*16.
 * @param occ_records   [out] Actual number of records written (<= n_records)
 * @param error_code    [out] A value from enum molalign_error_code
 *                            (error_codes.h): MOLALIGN_SUCCESS,
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
    bool mirror_flag, bool label_flag, bool bond_flag, double bond_tol,
    bool print_stats, bool print_assigntree, bool random_flag,
    int conv_freq, int max_trials,
    int n_records,
    double *rmsd_list, int *atomperm_list,
    double *transform_list, int *occ_records, int *error_code);

#ifdef __cplusplus
}
#endif

#endif /* MOLALIGN_C_BINDING_H */