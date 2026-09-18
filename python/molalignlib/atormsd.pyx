# cython: language_level=3
"""Cython wrapper around the atormsd_calculate C function."""

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

np.import_array()

cdef extern from "atormsd.h":
    void atormsd_calculate(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint mass_flag,
        bint mirror_flag, bint label_flag,
        bint print_stats, bint random_flag,
        double prune_tol, int conv_freq, int max_trials,
        int n_records,
        double *rmsd_list, int *natoms, int *atomperm_list,
        double *transform_list, int *occ_records, int *error_code)


_ERROR_MESSAGES = {
    1: "Clusters are not isomers.",
    2: "Atom types mismatch between the two clusters.",
    3: "Assignment failed (pruning tolerance is too tight).",
}


def calculate(
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
    double prune_tol,
    int    conv_freq,
    int    max_trials,
    int    n_records = 1,
):
    """
    Thin wrapper around ``atormsd_calculate``.

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
    (remaining keyword arguments map 1-to-1 onto the C flags)

    Returns
    -------
    rmsd : float64 ndarray, shape (occ_records,)
        RMSD of each returned candidate solution, best first.
    atom_permutation : int32 ndarray, shape (occ_records, n_atoms) — 0-based
    transform : float64 ndarray, shape (occ_records, 4, 4)

    ``occ_records`` (<= n_records) is however many distinct solutions the
    library actually found; it may be smaller than requested.
    """
    # ------------------------------------------------------------------ #
    # Validate & ensure C-contiguous layout                               #
    # ------------------------------------------------------------------ #
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

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]
    cdef int max_atoms = n1 if n1 > n2 else n2

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # ------------------------------------------------------------------ #
    # Output buffers - sized for up to n_records candidate solutions       #
    # ------------------------------------------------------------------ #
    cdef int    natoms      = 0
    cdef int    occ_records = 0
    cdef int    err         = 0

    cdef double *rmsd_list      = <double *>malloc(n_records * sizeof(double))
    cdef int    *atomperm_list  = <int *>malloc(n_records * max_atoms * sizeof(int))
    cdef double *transform_list = <double *>malloc(n_records * 16 * sizeof(double))

    if rmsd_list == NULL or atomperm_list == NULL or transform_list == NULL:
        if rmsd_list != NULL: free(rmsd_list)
        if atomperm_list != NULL: free(atomperm_list)
        if transform_list != NULL: free(transform_list)
        raise MemoryError("Could not allocate output buffers.")

    # ------------------------------------------------------------------ #
    # Call the Fortran/C library                                           #
    # ------------------------------------------------------------------ #
    try:
        atormsd_calculate(
            n1,
            <const int    *>atom_data1.data,
            <const double *>coords1.data,
            n2,
            <const int    *>atom_data2.data,
            <const double *>coords2.data,
            align_flag, remap_flag, heavy_flag, mass_flag,
            mirror_flag, label_flag,
            print_stats, random_flag,
            prune_tol, conv_freq, max_trials,
            n_records,
            rmsd_list, &natoms, atomperm_list,
            transform_list, &occ_records, &err,
        )

        if err != 0:
            raise ValueError(
                _ERROR_MESSAGES.get(err, "atormsd_calculate error code {}.".format(err))
            )

        rmsd = np.array([rmsd_list[i] for i in range(occ_records)], dtype=np.float64)

        perm = np.empty((occ_records, natoms), dtype=np.int32)
        for r in range(occ_records):
            for a in range(natoms):
                perm[r, a] = atomperm_list[r * natoms + a]

        transform = np.empty((occ_records, 4, 4), dtype=np.float64)
        for r in range(occ_records):
            for a in range(16):
                transform[r, a // 4, a % 4] = transform_list[r * 16 + a]
    finally:
        free(rmsd_list)
        free(atomperm_list)
        free(transform_list)

    return rmsd, perm, transform
