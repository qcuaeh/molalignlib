#ifndef ATORMSD_C_BINDING_H
#define ATORMSD_C_BINDING_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Calculate RMSD between two atom clusters
 * 
 * @param file1 Path to first molecule file
 * @param file2 Path to second molecule file
 * @param align Enable alignment (0=false, 1=true)
 * @param remap Enable atom remapping (0=false, 1=true)
 * @param heavy Use only heavy atoms (0=false, 1=true)
 * @param mass Use atomic masses as weights (0=false, 1=true)
 * @param mirror Mirror second molecule (0=false, 1=true)
 * @param label Use atom labels (0=false, 1=true)
 * @param mapping Print atom mapping (0=false, 1=true)
 * @param stats Print statistics (0=false, 1=true)
 * @param random Use random algorithm (0=false, 1=true)
 * @param prune_tol Pruning tolerance (use negative value for -near option)
 * @param ato_thres Atom threshold
 * @param max_trials Maximum number of trials
 * @param num_records Number of records to keep
 * @param rmsd Output: calculated RMSD value
 * @param natoms Output: number of atoms in permutation
 * @param atomperm Output: atom permutation array (0-based indexing, size must be >= natoms)
 * @param error_code Output: error code (0=success, 1=not isomers, 2=atom types mismatch)
 */
void atormsd_calculate(const char *file1, const char *file2,
                       int align, int remap, int heavy, int mass,
                       int mirror, int label, int mapping,
                       int stats, int random,
                       double prune_tol, int ato_thres, 
                       int max_trials, int num_records,
                       double *rmsd, int *natoms, int *atomperm, int *error_code);

#ifdef __cplusplus
}
#endif

#endif /* ATORMSD_C_BINDING_H */
