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
        bint stats_flag, bint random_flag,
        double prune_tol, int conv_freq, int max_trials,
        double *rmsd, int *natoms, int *atomperm,
        double *transform, int *error_code)


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
    bint   stats_flag,
    bint   random_flag,
    double prune_tol,
    int    conv_freq,
    int    max_trials,
):
    """
    Thin wrapper around ``atormsd_calculate``.

    Parameters
    ----------
    atom_data1, atom_data2 : int32 array, shape (n, 2)
        ``[[atomic_number, label], ...]`` for each molecule.
    coords1, coords2 : float64 array, shape (n, 3)
        Cartesian coordinates in Angstrom.
    (remaining keyword arguments map 1-to-1 onto the C flags)

    Returns
    -------
    rmsd : float
    atom_permutation : int32 ndarray, shape (n_atoms,)  — 0-based
    transform : float64 ndarray, shape (4, 4)
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

    cdef int n1 = atom_data1.shape[0]
    cdef int n2 = atom_data2.shape[0]

    atom_data1 = np.ascontiguousarray(atom_data1)
    atom_data2 = np.ascontiguousarray(atom_data2)
    coords1    = np.ascontiguousarray(coords1)
    coords2    = np.ascontiguousarray(coords2)

    # ------------------------------------------------------------------ #
    # Output buffers                                                       #
    # ------------------------------------------------------------------ #
    cdef double rmsd_val = 0.0
    cdef int    natoms   = 0
    cdef int    err      = 0

    cdef int    *atomperm = <int *>malloc(max(n1, n2) * sizeof(int))
    cdef double  tf[16]

    if atomperm == NULL:
        raise MemoryError("Could not allocate atom permutation buffer.")

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
            stats_flag, random_flag,
            prune_tol, conv_freq, max_trials,
            &rmsd_val, &natoms, atomperm,
            tf, &err,
        )

        if err != 0:
            raise ValueError(
                _ERROR_MESSAGES.get(err, "atormsd_calculate error code {}.".format(err))
            )

        perm = np.array([atomperm[i] for i in range(natoms)], dtype=np.int32)
    finally:
        free(atomperm)

    transform = np.array(tf, dtype=np.float64).reshape(4, 4)
    return rmsd_val, perm, transform
