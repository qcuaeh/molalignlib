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
        bint mirror_flag, bint label_flag, bint bond_flag, double bond_tol,
        bint print_stats, bint print_assigntree, bint random_flag,
        int conv_freq, int max_trials,
        int n_records,
        double *rmsd_list, int *natoms, int *atomperm_list,
        double *transform_list, int *occ_records, int *error_code)


_ERROR_MESSAGES = {
    1: "Molecules are not isomers (different atom counts or compositions).",
    2: "Atom type mismatch between the two conformers.",
    3: "Missing bond data for one or both conformers.",
    4: "Bond connectivity mismatch between the two conformers (only raised "
       "when remap_flag=False).",
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

    # bond_tol has no default: it is required (and only used) when
    # bond_flag=True, since that is what triggers geometry-based bond
    # perception inside the library.
    cdef double c_bond_tol = 0.0
    if bond_flag:
        if bond_tol is None:
            raise ValueError("bond_tol is required when bond_flag=True")
        c_bond_tol = <double>bond_tol

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
            mirror_flag, label_flag, bond_flag, c_bond_tol,
            print_stats, print_assigntree, random_flag,
            conv_freq, max_trials,
            n_records,
            rmsd_list, &natoms, atomperm_list,
            transform_list, &occ_records, &err,
        )

        if err != 0:
            raise ValueError(
                _ERROR_MESSAGES.get(err, "conformsd_calculate error code {}.".format(err))
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
