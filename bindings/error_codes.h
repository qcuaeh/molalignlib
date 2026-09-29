#ifndef ERROR_CODES_H
#define ERROR_CODES_H

/**
 * Error codes returned by atormsd_calculate and conformsd_calculate (see
 * molalign.h) through their error_code argument.
 *
 * Hand-maintained mirror of the Fortran error_codes module
 * (error_codes.f90): keep the values of both in sync.
 *
 * Each value has one meaning; the notes say which functions can return it.
 */

#ifdef __cplusplus
extern "C" {
#endif

enum molalign_error_code {
    /* No error. Returned by: both. */
    MOLALIGN_SUCCESS                  = 0,
    /* Different atom counts or compositions. Returned by: both. */
    MOLALIGN_ERROR_NOT_ISOMERS        = 1,
    /* Atom types differ in input order. Returned by: both
     * (remap_flag=false only). */
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH = 2,
    /* One or both molecules have no bonds. Returned by: conformsd. */
    MOLALIGN_ERROR_MISSING_BONDS      = 3,
    /* Bonds differ in input order. Returned by: conformsd
     * (remap_flag=false only). */
    MOLALIGN_ERROR_BOND_MISMATCH      = 4,
    /* Same composition but non-isomorphic bond graphs. Returned by:
     * conformsd (remap_flag=true only). */
    MOLALIGN_ERROR_NOT_CONFORMERS     = 5,
    /* No valid assignment under the pruning constraints (pruning tolerance
     * might be too tight). Returned by: atormsd (remap_flag=true only). */
    MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED  = 6,
    /* An atomic number is outside the element tables (0 to 104, where 0
     * is the dummy atom "X"). Returned by: atormsd, conformsd. */
    MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER     = 7
};

#ifdef __cplusplus
}
#endif

#endif /* ERROR_CODES_H */
