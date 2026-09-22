# cython: language_level=3
"""
Cython wrapper around the atormsd_calculate / conformsd_calculate C functions.

Both wrappers share nearly identical boilerplate (input validation,
output-buffer allocation, result packing); that logic lives once, in the
private helpers below. The two public entry points, atormsd_calculate(...)
and conformsd_calculate(...), match the underlying C functions declared in
molalign.h; the extern declarations below are aliased to
c_atormsd_calculate / c_conformsd_calculate internally to avoid colliding
with the Python-level wrapper names.
"""

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

np.import_array()

cdef extern from "error_codes.h":
    # Error codes returned via each function's error_code output argument.
    # Declared once in error_codes.h (itself kept in sync by hand with
    # the Fortran error_codes module, error_codes.f90) and reused here
    # instead of hardcoding the numbers a third time.
    enum: MOLALIGN_SUCCESS
    enum: MOLALIGN_ERROR_NOT_ISOMERS
    enum: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
    enum: MOLALIGN_ERROR_MISSING_BONDS
    enum: MOLALIGN_ERROR_BOND_MISMATCH
    enum: MOLALIGN_ERROR_NOT_CONFORMERS
    enum: MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED

cdef extern from "molalign.h":
    void c_atormsd_calculate "atormsd_calculate"(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint mass_flag,
        bint mirror_flag, bint label_flag,
        bint print_stats, bint random_flag,
        bint prune_flag, double prune_tol, int conv_freq, int max_trials,
        int n_records,
        double *rmsd_list, int *atomperm_list,
        double *transform_list, int *occ_records, int *error_code)

    void c_conformsd_calculate "conformsd_calculate"(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_bonds1, const int *bond_data1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        int n_bonds2, const int *bond_data2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint mass_flag,
        bint mirror_flag, bint label_flag, bint bond_flag, double bond_tol,
        bint print_stats, bint print_assigntree, bint random_flag,
        int conv_freq, int max_trials,
        int n_records,
        double *rmsd_list, int *atomperm_list,
        double *transform_list, int *occ_records, int *error_code)


_ATORMSD_ERROR_MESSAGES = {
    MOLALIGN_ERROR_NOT_ISOMERS: "Clusters are not isomers.",
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH: "Atom types mismatch between the two clusters.",
    MOLALIGN_ERROR_PRUNED_ASSIGNMENT_FAILED: "Assignment failed (pruning tolerance might be too tight).",
}

_CONFORMSD_ERROR_MESSAGES = {
    MOLALIGN_ERROR_NOT_ISOMERS: "Molecules are not isomers (different atom counts or compositions).",
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH: "Atom type mismatch between the two conformers.",
    MOLALIGN_ERROR_MISSING_BONDS: "Missing bond data for one or both conformers.",
    MOLALIGN_ERROR_BOND_MISMATCH: "Bond connectivity mismatch between the two conformers (only raised "
       "when remap_flag=False).",
    MOLALIGN_ERROR_NOT_CONFORMERS: "Molecules are not conformers (same composition but different "
       "bond graphs; only raised when remap_flag=True).",
}

# Sentinel used by rmsd.py when bond_flag=True (geometry-derived connectivity).
_EMPTY_BONDS = np.empty((0, 3), dtype=np.int32)


# --------------------------------------------------------------------------- #
# Shared helpers                                                              #
# --------------------------------------------------------------------------- #

cdef _validate_atoms(
    np.ndarray[np.int32_t,   ndim=2] atom_data1,
    np.ndarray[np.float64_t, ndim=2] coords1,
    np.ndarray[np.int32_t,   ndim=2] atom_data2,
    np.ndarray[np.float64_t, ndim=2] coords2,
    int n_records,
):
    """Validate atom/coord shapes and n_records; raise ValueError on failure."""
    if atom_data1.ndim != 2 or atom_data1.shape[1] != 2:
        raise ValueError("atom_data1 must have shape (n, 2)")
    if atom_data2.ndim != 2 or atom_data2.shape[1] != 2:
        raise ValueError("atom_data2 must have shape (n, 2)")
    if coords1.ndim != 2 or coords1.shape[1] != 3:
        raise ValueError("coords1 must have shape (n, 3)")
    if coords2.ndim != 2 or coords2.shape[1] != 3:
        raise ValueError("coords2 must have shape (n, 3)")
    if n_records < 1:
        raise ValueError("n_records must be >= 1")


cdef _validate_bonds(np.ndarray[np.int32_t, ndim=2] bd1,
                      np.ndarray[np.int32_t, ndim=2] bd2):
    """Validate bond-array shapes; raise ValueError on failure."""
    if bd1.ndim != 2 or bd1.shape[1] != 3:
        raise ValueError("bond_data1 must have shape (b, 3)")
    if bd2.ndim != 2 or bd2.shape[1] != 3:
        raise ValueError("bond_data2 must have shape (b, 3)")


cdef int _alloc_buffers(int n_records, int n_atoms,
                         double **rmsd_list, int **atomperm_list,
                         double **transform_list) except -1:
    """
    Malloc the three output buffers. Returns 0 on success; raises
    MemoryError (after freeing any partially-allocated buffers) on failure.
    """
    rmsd_list[0]      = <double *>malloc(n_records * sizeof(double))
    atomperm_list[0]  = <int *>malloc(n_records * n_atoms * sizeof(int))
    transform_list[0] = <double *>malloc(n_records * 16 * sizeof(double))

    if rmsd_list[0] == NULL or atomperm_list[0] == NULL or transform_list[0] == NULL:
        if rmsd_list[0] != NULL: free(rmsd_list[0])
        if atomperm_list[0] != NULL: free(atomperm_list[0])
        if transform_list[0] != NULL: free(transform_list[0])
        raise MemoryError("Could not allocate output buffers.")
    return 0


cdef _pack_outputs(double *rmsd_list, int *atomperm_list, double *transform_list,
                    int n_atoms, int occ_records):
    """Copy the raw C output buffers into numpy arrays."""
    rmsd = np.array([rmsd_list[i] for i in range(occ_records)], dtype=np.float64)

    perm = np.empty((occ_records, n_atoms), dtype=np.int32)
    for r in range(occ_records):
        for a in range(n_atoms):
            perm[r, a] = atomperm_list[r * n_atoms + a]

    transform = np.empty((occ_records, 4, 4), dtype=np.float64)
    for r in range(occ_records):
        for a in range(16):
            transform[r, a // 4, a % 4] = transform_list[r * 16 + a]

    return rmsd, perm, transform


# --------------------------------------------------------------------------- #
# atormsd                                                                     #
# --------------------------------------------------------------------------- #

def atormsd_calculate(
    np.ndarray[np.int32_t,   ndim=2] atom_data1,
    np.ndarray[np.float64_t, ndim=2] coords1,
    np.ndarray[np.int32_t,   ndim=2] atom_data2,
    np.ndarray[np.float64_t, ndim=2] coords2,
    bint   align_flag,
    bint   remap_flag,
    bint   heavy_flag,
    bint   mass_flag,
    bint   mirror_flag,
    bint   label_flag,
    bint   print_stats,
    bint   random_flag,
    bint   prune_flag,
    prune_tol,
    int    conv_freq,
    int    max_trials,
    int    n_records = 1,
):
    """
    Thin wrapper around the C ``atormsd_calculate`` function.

    Parameters
    ----------
    atom_data1, atom_data2 : int32 array, shape (n, 2)
        ``[[atomic_number, label], ...]`` for each molecule.
    coords1, coords2 : float64 array, shape (n, 3)
        Cartesian coordinates in Angstrom.
    n_records : int, default 1
        Maximum number of ranked candidate solutions to return. Records
        beyond the first are only ever produced when both ``align_flag``
        and ``remap_flag`` are true; otherwise exactly one solution is
        returned regardless of this value.
    prune_tol : float or None
        Pruning tolerance. Required (no default) when ``prune_flag`` is
        true; unused otherwise.
    (remaining keyword arguments map 1-to-1 onto the C flags)

    Returns
    -------
    rmsd : float64 ndarray, shape (occ_records,)
        RMSD of each returned candidate solution, best first.
    atom_permutation : int32 ndarray, shape (occ_records, n_atoms1) — 0-based
    transform : float64 ndarray, shape (occ_records, 4, 4)

    ``occ_records`` (<= n_records) is however many distinct solutions the
    library actually found; it may be smaller than requested.
    """
    _validate_atoms(atom_data1, coords1, atom_data2, coords2, n_records)

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # prune_tol has no default: it is required (and only used) when
    # prune_flag=True, since that is what enables pruning inside the
    # library.
    cdef double c_prune_tol = 0.0
    if prune_flag:
        if prune_tol is None:
            raise ValueError("prune_tol is required when prune_flag=True")
        c_prune_tol = <double>prune_tol

    cdef int    occ_records = 0
    cdef int    err         = 0

    cdef double *rmsd_list
    cdef int    *atomperm_list
    cdef double *transform_list
    _alloc_buffers(n_records, n1, &rmsd_list, &atomperm_list, &transform_list)

    try:
        c_atormsd_calculate(
            n1,
            <const int    *>atom_data1.data,
            <const double *>coords1.data,
            n2,
            <const int    *>atom_data2.data,
            <const double *>coords2.data,
            align_flag, remap_flag, heavy_flag, mass_flag,
            mirror_flag, label_flag,
            print_stats, random_flag,
            prune_flag, c_prune_tol, conv_freq, max_trials,
            n_records,
            rmsd_list, atomperm_list,
            transform_list, &occ_records, &err,
        )

        if err != 0:
            raise ValueError(
                _ATORMSD_ERROR_MESSAGES.get(err, "atormsd_calculate error code {}.".format(err))
            )

        rmsd, perm, transform = _pack_outputs(
            rmsd_list, atomperm_list, transform_list, n1, occ_records
        )
    finally:
        free(rmsd_list)
        free(atomperm_list)
        free(transform_list)

    return rmsd, perm, transform


# --------------------------------------------------------------------------- #
# conformsd                                                                   #
# --------------------------------------------------------------------------- #

def conformsd_calculate(
    np.ndarray[np.int32_t,   ndim=2] atom_data1,
    np.ndarray[np.float64_t, ndim=2] coords1,
    np.ndarray[np.int32_t,   ndim=2] atom_data2,
    np.ndarray[np.float64_t, ndim=2] coords2,
    bond_data1                        = None,
    bond_data2                        = None,
    bint align_flag    = False,
    bint remap_flag    = False,
    bint heavy_flag    = False,
    bint mass_flag     = False,
    bint mirror_flag   = False,
    bint label_flag    = False,
    bint bond_flag     = False,
    bond_tol                          = None,
    bint print_stats    = False,
    bint print_assigntree = False,
    bint random_flag   = False,
    int  conv_freq      = 100,
    int  max_trials    = 10000,
    int  n_records     = 1,
):
    """
    Thin wrapper around the C ``conformsd_calculate`` function.

    Parameters
    ----------
    atom_data1, atom_data2 : int32 array, shape (n, 2)
        ``[[atomic_number, label], ...]`` for each molecule.
    coords1, coords2 : float64 array, shape (n, 3)
        Cartesian coordinates in Angstrom.
    bond_data1, bond_data2 : int32 array, shape (b, 3), or None
        ``[[a1, a2, bond_type], ...]`` (1-based atom indices).
        Pass ``None`` (or omit) when ``bond_flag=True``; in that case the
        library derives connectivity from geometry.
    bond_tol : float, required when bond_flag=True
        Bond detection tolerance for geometry-based connectivity. Has no
        default and is ignored when bond_flag=False.
    align_flag, remap_flag : bool
        Both default to ``False``: by default no structural alignment or
        atom remapping is performed.
    print_assigntree : bool
        Print the atom-assignment search tree during optimisation.
    n_records : int, default 1
        Maximum number of ranked candidate solutions to return. Records
        beyond the first are only ever produced when both ``align_flag``
        and ``remap_flag`` are true; otherwise exactly one solution is
        returned regardless of this value.
    (remaining keyword arguments map 1-to-1 onto the C flags)

    Returns
    -------
    rmsd : float64 ndarray, shape (occ_records,)
        RMSD of each returned candidate solution, best first.
    atom_permutation : int32 ndarray, shape (occ_records, n_atoms1) — 0-based
    transform : float64 ndarray, shape (occ_records, 4, 4)

    ``occ_records`` (<= n_records) is however many distinct solutions the
    library actually found; it may be smaller than requested.
    """
    _validate_atoms(atom_data1, coords1, atom_data2, coords2, n_records)

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # Resolve bond arrays (None → empty placeholder).
    cdef np.ndarray[np.int32_t, ndim=2] bd1, bd2
    bd1 = np.ascontiguousarray(bond_data1 if bond_data1 is not None else _EMPTY_BONDS,
                                dtype=np.int32)
    bd2 = np.ascontiguousarray(bond_data2 if bond_data2 is not None else _EMPTY_BONDS,
                                dtype=np.int32)
    _validate_bonds(bd1, bd2)

    cdef int nb1 = bd1.shape[0]
    cdef int nb2 = bd2.shape[0]

    # bond_tol has no default: it is required (and only used) when
    # bond_flag=True, since that is what triggers geometry-based bond
    # perception inside the library.
    cdef double c_bond_tol = 0.0
    if bond_flag:
        if bond_tol is None:
            raise ValueError("bond_tol is required when bond_flag=True")
        c_bond_tol = <double>bond_tol

    cdef int    occ_records = 0
    cdef int    err         = 0

    cdef double *rmsd_list
    cdef int    *atomperm_list
    cdef double *transform_list
    _alloc_buffers(n_records, n1, &rmsd_list, &atomperm_list, &transform_list)

    try:
        c_conformsd_calculate(
            n1,
            <const int    *>atom_data1.data,
            <const double *>coords1.data,
            nb1,
            <const int    *>bd1.data,
            n2,
            <const int    *>atom_data2.data,
            <const double *>coords2.data,
            nb2,
            <const int    *>bd2.data,
            align_flag, remap_flag, heavy_flag, mass_flag,
            mirror_flag, label_flag, bond_flag, c_bond_tol,
            print_stats, print_assigntree, random_flag,
            conv_freq, max_trials,
            n_records,
            rmsd_list, atomperm_list,
            transform_list, &occ_records, &err,
        )

        if err != 0:
            raise ValueError(
                _CONFORMSD_ERROR_MESSAGES.get(err, "conformsd_calculate error code {}.".format(err))
            )

        rmsd, perm, transform = _pack_outputs(
            rmsd_list, atomperm_list, transform_list, n1, occ_records
        )
    finally:
        free(rmsd_list)
        free(atomperm_list)
        free(transform_list)

    return rmsd, perm, transform
