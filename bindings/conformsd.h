#ifndef CONFORMSD_C_BINDING_H
#define CONFORMSD_C_BINDING_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Calculate RMSD between two molecular conformers.
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
 * @param print_stats    Print optimisation statistics
 * @param print_assigntree Print the atom-assignment search tree
 * @param random_flag   Use random algorithm
 * @param conv_freq   Convergence frequency
 * @param max_trials    Maximum number of trials
 *
 * @param rmsd          [out] Calculated RMSD value
 * @param natoms        [out] Number of atoms in the permutation array
 * @param atomperm      [out] Atom permutation, 0-based, caller allocates >= natoms
 * @param transform     [out] Row-major 4x4 homogeneous transform (16 doubles).
 *                            Maps molecule-2 coords to molecule-1 frame:
 *                            p_out = R*p_in + t. Identity when align_flag=false.
 * @param error_code    [out] 0=success, 1=not isomers, 2=missing bonds,
 *                            3=atom types mismatch
 */
void conformsd_calculate(
    int n_atoms1, const int *atom_data1, const double *coords1,
    int n_bonds1, const int *bond_data1,
    int n_atoms2, const int *atom_data2, const double *coords2,
    int n_bonds2, const int *bond_data2,
    bool align_flag, bool remap_flag, bool heavy_flag, bool mass_flag,
    bool mirror_flag, bool label_flag, bool bond_flag,
    bool print_stats, bool print_assigntree, bool random_flag,
    int conv_freq, int max_trials,
    double *rmsd, int *natoms, int *atomperm,
    double *transform, int *error_code);

#ifdef __cplusplus
}
#endif

#endif /* CONFORMSD_C_BINDING_H */