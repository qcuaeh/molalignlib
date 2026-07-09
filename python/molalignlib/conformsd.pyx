# cython: language_level=3
"""Cython wrapper around the conformsd_calculate C function."""

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

np.import_array()

cdef extern from "conformsd.h":
    void conformsd_calculate(
        int n_atoms1, const int *atom_data1, const double *coords1,
        int n_bonds1, const int *bond_data1,
        int n_atoms2, const int *atom_data2, const double *coords2,
        int n_bonds2, const int *bond_data2,
        bint align_flag, bint remap_flag, bint heavy_flag, bint mass_flag,
        bint mirror_flag, bint label_flag, bint bond_flag,
        bint print_stats, bint print_assigntree, bint random_flag,
        int conv_freq, int max_trials,
        double *rmsd, int *natoms, int *atomperm,
        double *transform, int *error_code)


_ERROR_MESSAGES = {
    1: "Molecules are not isomers (different atom counts or compositions).",
    2: "Missing bond data for one or both conformers.",
    3: "Atom types mismatch between the two conformers.",
}

# Sentinel used by rmsd.py when bond_flag=True (geometry-derived connectivity).
_EMPTY_BONDS = np.empty((0, 3), dtype=np.int32)


def calculate(
    np.ndarray[np.int32_t,   ndim=2] atom_data1,
    np.ndarray[np.float64_t, ndim=2] coords1,
    np.ndarray[np.int32_t,   ndim=2] atom_data2,
    np.ndarray[np.float64_t, ndim=2] coords2,
    bond_data1                        = None,
    bond_data2                        = None,
    bint align_flag    = True,
    bint remap_flag    = True,
    bint heavy_flag    = False,
    bint mass_flag     = False,
    bint mirror_flag   = False,
    bint label_flag    = False,
    bint bond_flag     = False,
    bint print_stats    = False,
    bint print_assigntree = False,
    bint random_flag   = False,
    int  conv_freq      = 100,
    int  max_trials    = 10000,
):
    """
    Thin wrapper around ``conformsd_calculate``.

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
    print_assigntree : bool
        Print the atom-assignment search tree during optimisation.
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

    # Resolve bond arrays (None → empty placeholder).
    cdef np.ndarray[np.int32_t, ndim=2] bd1, bd2
    bd1 = np.ascontiguousarray(bond_data1 if bond_data1 is not None else _EMPTY_BONDS,
                                dtype=np.int32)
    bd2 = np.ascontiguousarray(bond_data2 if bond_data2 is not None else _EMPTY_BONDS,
                                dtype=np.int32)

    if bd1.ndim != 2 or bd1.shape[1] != 3:
        raise ValueError("bond_data1 must have shape (b, 3)")
    if bd2.ndim != 2 or bd2.shape[1] != 3:
        raise ValueError("bond_data2 must have shape (b, 3)")

    cdef int nb1 = bd1.shape[0]
    cdef int nb2 = bd2.shape[0]

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
        conformsd_calculate(
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
            mirror_flag, label_flag, bond_flag,
            print_stats, print_assigntree, random_flag,
            conv_freq, max_trials,
            &rmsd_val, &natoms, atomperm,
            tf, &err,
        )

        if err != 0:
            raise ValueError(
                _ERROR_MESSAGES.get(err, "conformsd_calculate error code {}.".format(err))
            )

        perm = np.array([atomperm[i] for i in range(natoms)], dtype=np.int32)
    finally:
        free(atomperm)

    transform = np.array(tf, dtype=np.float64).reshape(4, 4)
    return rmsd_val, perm, transform