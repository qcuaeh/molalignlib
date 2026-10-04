# cython: language_level=3
"""
Cython wrappers of the C functions atormsd_calculate and conformsd_calculate
(molalign.h). The C functions are imported as c_atormsd_calculate and
c_conformsd_calculate so that the Python wrappers can keep their names;
validation, buffer allocation, error handling and result packing are shared
helpers.
"""

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

np.import_array()

cdef extern from "molalign.h":
    # Error codes of the error_code argument, taken from molalign.h so
    # their values are not repeated here
    enum: MOLALIGN_SUCCESS
    enum: MOLALIGN_ERROR_NOT_ISOMERS
    enum: MOLALIGN_ERROR_ATOM_TYPE_MISMATCH
    enum: MOLALIGN_ERROR_TOO_MANY_FRAGMENTS
    enum: MOLALIGN_ERROR_BOND_MISMATCH
    enum: MOLALIGN_ERROR_NOT_CONFORMERS
    enum: MOLALIGN_ERROR_ASSIGNMENT_FAILED
    enum: MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER
    enum: MOLALIGN_ERROR_INVALID_BOUND

    void c_atormsd_calculate "atormsd_calculate"(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint massweight_flag,
        bint mirror_flag, bint useatomtype_flag,
        bint printstats_flag, bint random_flag,
        bint pruning_flag, double prune_tol, int ato_freq, int max_trials,
        int max_records,
        double *rmsd_list, int *mapping_list,
        double *transform_list, int *n_records, int *error_code)

    void c_conformsd_calculate "conformsd_calculate"(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_bonds1, const int *bond_data1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        int n_bonds2, const int *bond_data2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint massweight_flag,
        bint mirror_flag, bint useatomtype_flag, bint bonding_flag, double bond_tol,
        bint usebondtype_flag,
        bint printstats_flag, bint printassigntree_flag, bint random_flag,
        int confo_freq, int max_trials, int max_fragments,
        int max_records,
        double *rmsd_list, int *mapping_list,
        double *transform_list, int *n_records, int *error_code)


# Messages for every error code of both C functions. Each function can only
# return a subset of them (see molalign.h); codes that a function never
# returns are simply never looked up for it.
_ERROR_MESSAGES = {
    MOLALIGN_ERROR_INVALID_ATOMIC_NUMBER: "Atomic number out of range (valid: 0 to 104, where 0 is the dummy atom 'X').",
    MOLALIGN_ERROR_NOT_ISOMERS: "Molecules are not isomers (different atom counts or compositions).",
    MOLALIGN_ERROR_ATOM_TYPE_MISMATCH: "Atom type mismatch between the two molecules.",
    MOLALIGN_ERROR_TOO_MANY_FRAGMENTS: "One or both molecules have more molecular fragments "
       "(connected components of the bond graph) than max_fragments.",
    MOLALIGN_ERROR_BOND_MISMATCH: "Bond connectivity mismatch between the two conformers, or bond "
       "type mismatch when usebondtype_flag=True (only raised when remap_flag=False).",
    MOLALIGN_ERROR_NOT_CONFORMERS: "Molecules are not conformers (same composition but different "
       "bond graphs, or different bond types when usebondtype_flag=True; only raised when "
       "remap_flag=True).",
    MOLALIGN_ERROR_ASSIGNMENT_FAILED: "Assignment failed (pruning tolerance might be too tight).",
    MOLALIGN_ERROR_INVALID_BOUND: "A bound parameter is out of range.",
}

# Bond array passed when no bonds are given (bond_data=None)
_EMPTY_BONDS = np.empty((0, 3), dtype=np.int32)


# --------------------------------------------------------------------------- #
# Shared helpers                                                              #
# --------------------------------------------------------------------------- #

cdef _validate_atoms(
    np.ndarray[np.int32_t,   ndim=2] atom_data1,
    np.ndarray[np.float64_t, ndim=2] coords1,
    np.ndarray[np.int32_t,   ndim=2] atom_data2,
    np.ndarray[np.float64_t, ndim=2] coords2,
):
    """Validate atom/coord shapes; raise ValueError on failure."""
    if atom_data1.ndim != 2 or atom_data1.shape[1] != 2:
        raise ValueError("atom_data1 must have shape (n, 2)")
    if atom_data2.ndim != 2 or atom_data2.shape[1] != 2:
        raise ValueError("atom_data2 must have shape (n, 2)")
    if coords1.ndim != 2 or coords1.shape[1] != 3:
        raise ValueError("coords1 must have shape (n, 3)")
    if coords2.ndim != 2 or coords2.shape[1] != 3:
        raise ValueError("coords2 must have shape (n, 3)")


def _validate_counts(**counts):
    """
    Validate count parameters (*_freq, max_*):
    each must be >= 1. Raise ValueError naming the first offending one.
    The C functions check them as well (MOLALIGN_ERROR_INVALID_BOUND),
    but checking here gives a precise message before any buffer is
    allocated.
    """
    for name, value in counts.items():
        if value < 1:
            raise ValueError("{} must be >= 1 (got {})".format(name, value))


cdef _validate_bonds(np.ndarray[np.int32_t, ndim=2] bd1,
                      np.ndarray[np.int32_t, ndim=2] bd2):
    """Validate bond-array shapes; raise ValueError on failure."""
    if bd1.ndim != 2 or bd1.shape[1] != 3:
        raise ValueError("bond_data1 must have shape (b, 3)")
    if bd2.ndim != 2 or bd2.shape[1] != 3:
        raise ValueError("bond_data2 must have shape (b, 3)")


cdef _raise_on_error(int err, str func_name):
    """Raise ValueError for a nonzero error code from the C library."""
    if err != MOLALIGN_SUCCESS:
        raise ValueError(
            _ERROR_MESSAGES.get(err, "{} error code {}.".format(func_name, err))
        )


cdef int _alloc_buffers(int max_records, int n_atoms,
                         double **rmsd_list, int **mapping_list,
                         double **transform_list) except -1:
    """
    Allocate the three output buffers for max_records records of n_atoms
    mapping entries. Returns 0; on failure frees any allocated buffer and
    raises MemoryError.
    """
    rmsd_list[0]      = <double *>malloc(max_records * sizeof(double))
    mapping_list[0]  = <int *>malloc(max_records * n_atoms * sizeof(int))
    transform_list[0] = <double *>malloc(max_records * 16 * sizeof(double))

    if rmsd_list[0] == NULL or mapping_list[0] == NULL or transform_list[0] == NULL:
        if rmsd_list[0] != NULL: free(rmsd_list[0])
        if mapping_list[0] != NULL: free(mapping_list[0])
        if transform_list[0] != NULL: free(transform_list[0])
        raise MemoryError("Could not allocate output buffers.")
    return 0


cdef _pack_outputs(double *rmsd_list, int *mapping_list, double *transform_list,
                    int n_atoms, int n_records):
    """Copy the first n_records records of the C buffers into numpy arrays."""
    rmsd = np.array([rmsd_list[i] for i in range(n_records)], dtype=np.float64)

    mapping = np.empty((n_records, n_atoms), dtype=np.int32)
    for r in range(n_records):
        for a in range(n_atoms):
            mapping[r, a] = mapping_list[r * n_atoms + a]

    transform = np.empty((n_records, 4, 4), dtype=np.float64)
    for r in range(n_records):
        for a in range(16):
            transform[r, a // 4, a % 4] = transform_list[r * 16 + a]

    return rmsd, mapping, transform


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
    bint   massweight_flag,
    bint   mirror_flag,
    bint   useatomtype_flag,
    bint   printstats_flag,
    bint   random_flag,
    bint   pruning_flag,
    prune_tol,
    int    ato_freq,
    int    max_trials,
    int    max_records = 1,
):
    """
    Thin wrapper around the C ``atormsd_calculate`` function (see molalign.h
    for the meaning of the flags).

    Parameters
    ----------
    atom_data1, atom_data2 : int32 array, shape (n, 2)
        ``[[atomic_number, label], ...]`` for each molecule.
    coords1, coords2 : float64 array, shape (n, 3)
        Cartesian coordinates in Angstrom.
    ato_freq, max_trials : int
        Convergence frequency and maximum number of random orientations of
        the search; both must be >= 1.
    max_records : int, default 1
        Maximum number of ranked candidate solutions to return (>= 1).
        Records beyond the first are only ever produced when both
        ``align_flag`` and ``remap_flag`` are true; otherwise exactly one
        solution is returned regardless of this value.
    prune_tol : float or None
        Pruning tolerance (Angstrom). Required (no default) when
        ``pruning_flag`` is true; unused otherwise.
    (the remaining arguments map one-to-one onto the C arguments)

    Returns
    -------
    rmsd : float64 ndarray, shape (n_records,)
        RMSD of each returned candidate solution, best first.
    mapping : int32 ndarray, shape (n_records, n_padding) — 0-based
        ``n_padding = max(n_atoms1, n_atoms2)``; same padding convention as
        ``conformsd_calculate``. Sizes can only differ when ``heavy_flag=True``.
    transform : float64 ndarray, shape (n_records, 4, 4)
        Homogeneous matrices mapping the input coordinates of molecule 2 to
        the molecule-1 frame. With ``mirror_flag=True`` the 3x3 block
        includes the reflection (determinant -1). With ``align_flag=False``
        it is the identity, or just the reflection when mirroring.

    ``n_records`` (<= max_records) is however many distinct solutions the
    library actually found; it may be smaller than requested.
    """
    _validate_atoms(atom_data1, coords1, atom_data2, coords2)
    _validate_counts(ato_freq=ato_freq, max_trials=max_trials, max_records=max_records)

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]
    cdef int n_padding = max(n1, n2)

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # prune_tol is required, and only used, when pruning_flag=True
    cdef double c_prune_tol = 0.0
    if pruning_flag:
        if prune_tol is None:
            raise ValueError("prune_tol is required when pruning_flag=True")
        c_prune_tol = <double>prune_tol

    cdef int    n_records = 0
    cdef int    err         = 0

    cdef double *rmsd_list
    cdef int    *mapping_list
    cdef double *transform_list
    _alloc_buffers(max_records, n_padding, &rmsd_list, &mapping_list, &transform_list)

    try:
        c_atormsd_calculate(
            n1,
            <const int    *>atom_data1.data,
            <const double *>coords1.data,
            n2,
            <const int    *>atom_data2.data,
            <const double *>coords2.data,
            align_flag, remap_flag, heavy_flag, massweight_flag,
            mirror_flag, useatomtype_flag,
            printstats_flag, random_flag,
            pruning_flag, c_prune_tol, ato_freq, max_trials,
            max_records,
            rmsd_list, mapping_list,
            transform_list, &n_records, &err,
        )

        _raise_on_error(err, "atormsd_calculate")

        rmsd, mapping, transform = _pack_outputs(
            rmsd_list, mapping_list, transform_list, n_padding, n_records
        )
    finally:
        free(rmsd_list)
        free(mapping_list)
        free(transform_list)

    return rmsd, mapping, transform


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
    bint massweight_flag     = False,
    bint mirror_flag   = False,
    bint useatomtype_flag = False,
    bint bonding_flag  = False,
    bond_tol                          = None,
    bint usebondtype_flag = False,
    bint printstats_flag    = False,
    bint printassigntree_flag = False,
    bint random_flag   = False,
    int  confo_freq      = 100,
    int  max_trials    = 10000,
    int  max_fragments     = 1,
    int  max_records     = 1,
):
    """
    Thin wrapper around the C ``conformsd_calculate`` function (see
    molalign.h for the meaning of the flags).

    Parameters
    ----------
    atom_data1, atom_data2 : int32 array, shape (n, 2)
        ``[[atomic_number, label], ...]`` for each molecule.
    coords1, coords2 : float64 array, shape (n, 3)
        Cartesian coordinates in Angstrom.
    bond_data1, bond_data2 : int32 array, shape (b, 3), or None
        ``[[a1, a2, bond_type], ...]`` (1-based atom indices).
        Pass ``None`` (or omit) when ``bonding_flag=True``; in that case the
        library derives connectivity from geometry.
    bond_tol : float, required when bonding_flag=True
        Bond detection tolerance (Angstrom) for geometry-based connectivity.
        Has no default and is ignored when bonding_flag=False.
    usebondtype_flag : bool, default False
        Use the bond types (third column of ``bond_data``) to guide atom
        matching. Types are compared for equality, never interpreted, so
        both molecules must use the same bond-type convention (read from
        the same file format with the same parser). Ignored when
        ``bonding_flag=True``.
    align_flag, remap_flag : bool, default False
        Optimally superpose molecule 2, and search the atom permutation
        that minimises the RMSD.
    printassigntree_flag : bool
        Print the assignment tree and its combination counts to stdout.
    confo_freq : int, default 100
        Stop the random orientation search once the best solution has been
        found more than this many times; also the threshold on the ratio of
        total to partial assignment combinations above which that search is
        used instead of exhaustive enumeration. Must be >= 1.
    max_trials : int, default 10000
        Maximum number of random orientations (>= 1).
    max_fragments : int, default 1
        Maximum number of molecular fragments (connected components of the
        bond graph of the compared atoms, i.e. after ``heavy_flag``
        exclusions) allowed in each molecule (>= 1). An atom without bonds
        is a fragment of its own, so a molecule without bonds has as many
        fragments as atoms. Exceeding it raises ValueError.
    max_records : int, default 1
        Maximum number of ranked candidate solutions to return (>= 1).
        Records beyond the first are only ever produced when both
        ``align_flag`` and ``remap_flag`` are true; otherwise exactly one
        solution is returned regardless of this value.
    (the remaining keyword arguments map one-to-one onto the C arguments)

    Returns
    -------
    rmsd : float64 ndarray, shape (n_records,)
        RMSD of each returned candidate solution, best first.
    mapping : int32 ndarray, shape (n_records, n_padding) — 0-based
        ``n_padding = max(n_atoms1, n_atoms2)``. Each row is a permutation of
        ``0..n_padding-1``; entry ``j`` is the atom of molecule 2 placed on line
        ``j`` of molecule 1. The smaller molecule is padded with dummy atoms
        appended after its real atoms: values ``>= n_atoms2`` denote dummy
        atoms of molecule 2, and entries ``j >= n_atoms1`` hold the extra
        atoms of molecule 2. Sizes can only differ when ``heavy_flag=True``.
    transform : float64 ndarray, shape (n_records, 4, 4)
        Homogeneous matrices mapping the input coordinates of molecule 2 to
        the molecule-1 frame. With ``mirror_flag=True`` the 3x3 block
        includes the reflection (determinant -1). With ``align_flag=False``
        it is the identity, or just the reflection when mirroring.

    ``n_records`` (<= max_records) is however many distinct solutions the
    library actually found; it may be smaller than requested.
    """
    _validate_atoms(atom_data1, coords1, atom_data2, coords2)
    _validate_counts(confo_freq=confo_freq, max_trials=max_trials,
                     max_fragments=max_fragments, max_records=max_records)

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]
    cdef int n_padding = max(n1, n2)

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # None means no bonds
    cdef np.ndarray[np.int32_t, ndim=2] bd1, bd2
    bd1 = np.ascontiguousarray(bond_data1 if bond_data1 is not None else _EMPTY_BONDS,
                                dtype=np.int32)
    bd2 = np.ascontiguousarray(bond_data2 if bond_data2 is not None else _EMPTY_BONDS,
                                dtype=np.int32)
    _validate_bonds(bd1, bd2)

    cdef int nb1 = bd1.shape[0]
    cdef int nb2 = bd2.shape[0]

    # bond_tol is required, and only used, when bonding_flag=True
    cdef double c_bond_tol = 0.0
    if bonding_flag:
        if bond_tol is None:
            raise ValueError("bond_tol is required when bonding_flag=True")
        c_bond_tol = <double>bond_tol

    cdef int    n_records = 0
    cdef int    err         = 0

    cdef double *rmsd_list
    cdef int    *mapping_list
    cdef double *transform_list
    _alloc_buffers(max_records, n_padding, &rmsd_list, &mapping_list, &transform_list)

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
            align_flag, remap_flag, heavy_flag, massweight_flag,
            mirror_flag, useatomtype_flag, bonding_flag, c_bond_tol,
            usebondtype_flag,
            printstats_flag, printassigntree_flag, random_flag,
            confo_freq, max_trials, max_fragments,
            max_records,
            rmsd_list, mapping_list,
            transform_list, &n_records, &err,
        )

        _raise_on_error(err, "conformsd_calculate")

        rmsd, mapping, transform = _pack_outputs(
            rmsd_list, mapping_list, transform_list, n_padding, n_records
        )
    finally:
        free(rmsd_list)
        free(mapping_list)
        free(transform_list)

    return rmsd, mapping, transform
