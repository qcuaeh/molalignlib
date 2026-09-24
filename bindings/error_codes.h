#ifndef ERROR_CODES_H
#define ERROR_CODES_H

/**
 * Error codes returned by atormsd_calculate, conformsd_calculate and
 * isormsd_calculate (see molalign.h) via their `error_code` output
 * parameter.
 *
 * This header is the C-facing mirror of the Fortran `error_codes` module
 * (error_codes.f90); the two declare the same values by hand and must be
 * kept numerically in sync - there is no generation step tying them
 * together. If you change a value here, change it there too, and vice
 * versa.
 *
 * Every value has exactly one meaning, whichever function returns it. Not
 * every function can return every code; the "Returned by" notes say which
 * can.
 */

#ifdef __cplusplus
extern "C" {
#endif

enum molalign_error_code {
    /* No error. Returned by: all. */
    MOLALIGN_SUCCESS                  = 0,
    /* Different atom counts or compositions. Returned by: all. */
    MOLALIGN_ERROR_NOT_ISOMERS        = 1,
    /* Atom types differ in input order. Returned by: all
     * (remap_flag=false only). */
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH = 2,
    /* One or both molecules have no bonds. Returned by: conformsd, isormsd. */
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
