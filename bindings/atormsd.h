#ifndef ATORMSD_C_BINDING_H
#define ATORMSD_C_BINDING_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Calculate RMSD between two atom clusters.
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
 * @param stats_flag    Print optimisation statistics
 * @param random_flag   Use random algorithm
 * @param prune_tol     Pruning tolerance (negative = disable pruning)
 * @param ato_thres     Atom threshold
 * @param max_trials    Maximum number of trials
 *
 * @param rmsd          [out] Calculated RMSD value
 * @param natoms        [out] Number of atoms in the permutation array
 * @param atomperm      [out] Atom permutation, 0-based, caller allocates >= natoms
 * @param transform     [out] Row-major 4x4 homogeneous transform (16 doubles).
 *                            Maps molecule-2 coords to molecule-1 frame:
 *                            p_out = R*p_in + t. Identity when align_flag=false.
 * @param error_code    [out] 0=success, 1=not isomers, 2=atom types mismatch
 */
void atormsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag,
    bool stats_flag, bool random_flag,
    double prune_tol, int ato_thres, int max_trials,
    double *rmsd, int *natoms, int *atomperm,
    double *transform, int *error_code);

#ifdef __cplusplus
}
#endif

#endif /* ATORMSD_C_BINDING_H */