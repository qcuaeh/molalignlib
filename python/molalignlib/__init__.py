"""
High-level Python interface of MolAlignLib, built on the molalign Cython
extension.

Supported file formats (via chemfiles):
  - XYZ        (.xyz)
  - SDF / MOL  (.sdf, .mol)
  - MOL2       (.mol2)
  - and many others

Classes
-------
Molecule
    A set of atoms with 3-D coordinates and an optional bond table. The same
    object can be compared in two ways, one method per algorithm of the
    library:

    atormsd_to(other, ...)
        Unstructured comparison (atormsd_calculate). Bonds are ignored and
        atoms are matched only within atom types.

    conformsd_to(other, ...)
        Conformer comparison (conformsd_calculate). Atom assignments always
        respect the bond topology (HNA partitioning); both molecules must
        have the same bond graph.

RMSDResult
    Return value of both methods. Each method returns a list of
    RMSDResult objects, best (lowest RMSD) first. By default only the single
    best solution is computed (a list of length 1); pass a larger max_records
    to get several ranked candidate solutions at once. The returned list may
    be shorter than max_records if the library did not find that many distinct
    solutions.
"""

from pathlib import Path
import numpy as np
import chemfiles

from . import molalign as _molalign


__all__ = [
    "ATOMIC_SYMBOLS",
    "DUMMY_ATOMIC_NUMBER",
    "symbol_to_atomic_number",
    "atomic_number_to_symbol",
    "RMSDResult",
    "Molecule",
    "read_molecules",
]


# ---------------------------------------------------------------------------
# Element symbol <-> atomic number table
#
# Mirrors atomic_symbols(0:n_elems) in chemdata.f90, so index i is element
# number i of the Fortran core. Entry 0 ("X") is the dummy atom and entry
# 104 ("LJ") the Lennard-Jones placeholder; every real element's number
# equals its atomic number.
# ---------------------------------------------------------------------------

ATOMIC_SYMBOLS = (
    "X",
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "LJ",
)

# Element number of the dummy atom ("X"). Used to pad the smaller molecule;
# dummy atoms are always excluded from the comparison.
DUMMY_ATOMIC_NUMBER = 0

# Case-insensitive symbol -> element number lookup (0 = dummy atom).
_SYMBOL_TO_NUMBER = {sym.upper(): i for i, sym in enumerate(ATOMIC_SYMBOLS)}

# Highest element number that is a real element (Lr). Above it, chemfiles'
# numbering (104 = Rf, ...) and the library's (104 = LJ) disagree.
_LAST_REAL_ELEMENT = _SYMBOL_TO_NUMBER["LR"]


def symbol_to_atomic_number(symbol):
    """Look up the atomic number for an element symbol (case-insensitive)."""
    key = symbol.strip().upper()
    try:
        return _SYMBOL_TO_NUMBER[key]
    except KeyError:
        raise ValueError("Unrecognized element symbol: {!r}".format(symbol))


def atomic_number_to_symbol(number):
    """Element symbol of an element number, or the number as a string if unknown."""
    number = int(number)
    if 0 <= number < len(ATOMIC_SYMBOLS):
        return ATOMIC_SYMBOLS[number]
    return str(number)


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

class RMSDResult(object):
    """Return value of Molecule.atormsd_to / conformsd_to."""

    def __init__(self, rmsd, mapping, transform):
        self.rmsd = rmsd
        """Root-mean-square deviation (Angstrom)."""
        self.mapping = mapping
        """0-based index array mapping other's atoms onto self's atoms.

        Entry j is the atom of `other` placed on line j of `self`. Its length
        is max(len(self), len(other)); the two differ only with
        heavy_only=True. Values >= len(other) denote dummy atoms ("X") that
        pad `other`; entries j >= len(self) hold the extra atoms of `other`.
        """
        self.transform = transform
        """4x4 homogeneous matrix mapping other's input coordinates to self's frame.

        Rotation plus translation; with mirror=True the 3x3 block also
        includes the reflection (determinant -1), so apply_to returns the
        mirrored, aligned structure.
        """

    def apply_to(self, molecule):
        """
        Apply the stored transform and atom permutation to a Molecule,
        returning a new Molecule aligned and reordered to match the reference.

        The result has len(mapping) atoms, so line j of it corresponds
        to line j of the reference. If the reference has more atoms (possible
        only with heavy_only=True), the missing lines are filled with dummy
        atoms (element "X", atomic number 0) whose coordinates are
        placeholders. If it has fewer, the extra atoms come last. Bonds are
        renumbered to the new atom order.
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("Expected Molecule, got {}".format(type(molecule).__name__))

        R = self.transform[:3, :3]
        t = self.transform[:3, 3]
        idx = np.asarray(self.mapping)

        # When the permutation is longer than the molecule, pad it with dummy
        # atoms ("X") appended after the real atoms, as the Fortran core does,
        # so that every index is valid.
        coords = molecule.coords
        atom_data = molecule.atom_data
        symbols = list(molecule.symbols)
        n_dummy = len(idx) - len(coords)
        if n_dummy > 0:
            coords = np.vstack([coords, np.zeros((n_dummy, 3), dtype=coords.dtype)])
            dummy_row = np.array([[DUMMY_ATOMIC_NUMBER, 0]], dtype=atom_data.dtype)
            atom_data = np.vstack([atom_data, np.repeat(dummy_row, n_dummy, axis=0)])
            symbols += [ATOMIC_SYMBOLS[DUMMY_ATOMIC_NUMBER]] * n_dummy

        new_coords = coords[idx] @ R.T + t
        new_atom_data = atom_data[idx]
        new_symbols = [symbols[i] for i in idx]

        # Renumber the (1-based) bond atoms to the new atom order
        new_bond_data = molecule.bond_data.copy()
        if new_bond_data.shape[0] > 0:
            inv_idx = np.argsort(idx)
            new_bond_data[:, :2] = inv_idx[new_bond_data[:, :2] - 1] + 1

        return Molecule(
            atom_data=new_atom_data,
            coords=new_coords,
            bond_data=new_bond_data,
            name=molecule.name,
            symbols=new_symbols,
            bond_source=molecule.bond_source,
        )

    def __repr__(self):
        return "RMSDResult(rmsd={:.6f})".format(self.rmsd)


# ---------------------------------------------------------------------------
# Chemfiles helpers
# ---------------------------------------------------------------------------

def _file_format(path):
    """
    Format tag of a file, taken from its extension (the same way chemfiles
    infers the format). Molecules read from files with the same tag were
    parsed by the same chemfiles reader, so their bond types follow the same
    convention.
    """
    return Path(path).suffix.lower().lstrip(".")


def _extract_frame_data(frame):
    """
    Extract atom_data, coords, bond_data and symbols from a chemfiles Frame.
    Atoms are unlabelled; bonds without a known order get type 1.
    """
    n_atoms = len(frame.atoms)
    atom_data = np.zeros((n_atoms, 2), dtype=np.int32)
    symbols = []

    for i, atom in enumerate(frame.atoms):
        # chemfiles reports 0 (or None) for types it doesn't know, such as
        # "X" and "LJ", and uses its own numbering beyond Lr. In those cases
        # fall back to the library's table, which raises ValueError for
        # symbols it doesn't know either.
        elnum = atom.atomic_number
        if not elnum or elnum > _LAST_REAL_ELEMENT:
            elnum = symbol_to_atomic_number(atom.type)
        atom_data[i, 0] = elnum
        atom_data[i, 1] = 0
        symbols.append(atom.type)

    coords = np.array(frame.positions, dtype=np.float64)

    bonds = frame.topology.bonds
    try:
        orders = frame.topology.bond_orders
    except AttributeError:
        orders = []

    n_bonds = len(bonds)
    bond_data = np.empty((n_bonds, 3), dtype=np.int32)
    for i in range(n_bonds):
        bond_data[i, 0] = bonds[i][0] + 1
        bond_data[i, 1] = bonds[i][1] + 1
        if orders and i < len(orders):
            bond_data[i, 2] = int(orders[i])
        else:
            bond_data[i, 2] = 1

    return atom_data, coords, bond_data, symbols


def _build_frame(coords, symbols, comment=None, bond_data=None):
    """
    Build a chemfiles.Frame from coordinates, symbols, and (optionally)
    1-based bond_data (as stored on Molecule), for writing out through
    chemfiles.Trajectory.
    """
    frame = chemfiles.Frame()
    for sym, (x, y, z) in zip(symbols, coords):
        frame.add_atom(chemfiles.Atom(sym), [x, y, z])

    if bond_data is not None:
        for row in bond_data:
            # bond_data atom indices are 1-based; chemfiles is 0-based.
            frame.add_bond(int(row[0]) - 1, int(row[1]) - 1)

    if comment:
        try:
            frame.comment = comment
        except AttributeError:
            # Some chemfiles versions/formats don't expose a settable comment.
            pass

    return frame


def _write_frame(path, frame):
    """
    Write a single-frame chemfiles.Frame to path. The output format is
    inferred from the file extension (.xyz, .pdb, .sdf, .mol2, ...), the
    same way chemfiles infers the format when reading.
    """
    with chemfiles.Trajectory(str(path), "w") as trajectory:
        trajectory.write(frame)


def _molecule_from_frame(frame, name, bond_source):
    atom_data, coords, bond_data, symbols = _extract_frame_data(frame)
    return Molecule(
        atom_data=atom_data, coords=coords, bond_data=bond_data,
        name=name, symbols=symbols, bond_source=bond_source,
    )


# ---------------------------------------------------------------------------
# Public frame reader
# ---------------------------------------------------------------------------

def read_molecules(path, frames=None):
    """
    Read the frames of a file as Molecule objects: a tuple with the frames
    of the given indices, or a list of all frames if frames is None. Bonds
    are read when the format provides them, and bond_source is set to the
    file format.
    """
    stem = Path(path).stem
    source = _file_format(path)
    molecules = []

    with chemfiles.Trajectory(str(path)) as trajectory:
        if frames is not None:
            for idx in frames:
                frame = trajectory.read_step(idx)
                name = "{}_{}".format(stem, idx)
                molecules.append(_molecule_from_frame(frame, name, source))
            return tuple(molecules)

        for idx, frame in enumerate(trajectory):
            name = "{}_{}".format(stem, idx)
            molecules.append(_molecule_from_frame(frame, name, source))
        return molecules


# ---------------------------------------------------------------------------
# Molecule
# ---------------------------------------------------------------------------

class Molecule(object):
    """
    A set of atoms with 3-D coordinates and an optional bond table.

    Use atormsd_to to compare it as an unstructured cluster (bonds ignored)
    and conformsd_to to compare it as a conformer (same bond graph
    required).
    """

    def __init__(self, atom_data=None, coords=None, bond_data=None, name="molecule",
                 symbols=None, labels=None, bond_source=None):
        """
        Either ``atom_data`` (an (n, 2) array of element numbers and labels)
        or ``symbols`` (element symbols, e.g. ["C", "H", "H", "H"]) must be
        provided; the other is derived. If both are given, ``atom_data`` is
        used for the calculations and ``symbols`` only for output.

        bond_data : array of shape (b, 3), optional
            ``[[a1, a2, bond_type], ...]`` with 1-based atom indices.
            Defaults to no bonds. Required by conformsd_to unless
            connectivity is inferred with ``bond_tol``.

        labels : sequence of int, optional
            Per-atom labels, which restrict matching to atoms with the same
            label when ``use_atom_type=True``. Only used when building
            ``atom_data`` from ``symbols``; defaults to all zero
            (unlabelled).

        bond_source : str, optional
            Tag naming the convention of the bond types in ``bond_data``
            (the file format for molecules read from files). Bond types are
            only compared when both molecules share the same tag; see
            ``use_bond_type``. Defaults to None, which is the tag of all
            molecules built directly from arrays.
        """
        self._coords = np.asarray(coords, dtype=np.float64)

        if atom_data is None:
            if symbols is None:
                raise ValueError("Either atom_data or symbols must be provided")
            elnums = [symbol_to_atomic_number(s) for s in symbols]
            if labels is None:
                labels = [0] * len(elnums)
            elif len(labels) != len(elnums):
                raise ValueError("labels must have the same length as symbols")
            self._atom_data = np.array(list(zip(elnums, labels)), dtype=np.int32).reshape(-1, 2)
            self._symbols = list(symbols)
        else:
            self._atom_data = np.asarray(atom_data, dtype=np.int32)
            if symbols is None:
                self._symbols = [atomic_number_to_symbol(a[0]) for a in self._atom_data]
            else:
                self._symbols = list(symbols)

        if self._coords.shape[0] != self._atom_data.shape[0]:
            raise ValueError("coords and atom_data/symbols must have the same length")

        if bond_data is not None:
            self._bond_data = np.asarray(bond_data, dtype=np.int32).reshape(-1, 3)
        else:
            self._bond_data = np.empty((0, 3), dtype=np.int32)

        self._name = name
        self._bond_source = bond_source

    # -- constructors -------------------------------------------------------

    @classmethod
    def from_file(cls, path, frame_idx=0):
        """Read a single frame from a file as a Molecule."""
        with chemfiles.Trajectory(str(path)) as trajectory:
            frame = trajectory.read_step(frame_idx)
            return _molecule_from_frame(frame, Path(path).stem, _file_format(path))

    @classmethod
    def from_symbols(cls, symbols, coords, bond_data=None, name="molecule", labels=None,
                     bond_source=None):
        """Build a Molecule directly from element symbols, coordinates, and bonds."""
        return cls(coords=coords, bond_data=bond_data, name=name, symbols=symbols,
                   labels=labels, bond_source=bond_source)

    @classmethod
    def from_numbers(cls, atomic_numbers, coords, bond_data=None, name="molecule", labels=None,
                     bond_source=None):
        """Build a Molecule directly from atomic numbers, coordinates, and bonds."""
        atomic_numbers = list(atomic_numbers)
        if labels is None:
            labels = [0] * len(atomic_numbers)
        elif len(labels) != len(atomic_numbers):
            raise ValueError("labels must have the same length as atomic_numbers")
        atom_data = np.array(list(zip(atomic_numbers, labels)), dtype=np.int32).reshape(-1, 2)
        return cls(atom_data=atom_data, coords=coords, bond_data=bond_data, name=name,
                   bond_source=bond_source)

    # -- properties ---------------------------------------------------------

    @property
    def name(self):
        return self._name

    @property
    def n_atoms(self):
        return self._atom_data.shape[0]

    @property
    def n_bonds(self):
        return self._bond_data.shape[0]

    @property
    def has_bonds(self):
        return self.n_bonds > 0

    @property
    def atom_data(self):
        return self._atom_data

    @property
    def coords(self):
        return self._coords

    @property
    def bond_data(self):
        return self._bond_data

    @property
    def symbols(self):
        return self._symbols

    @property
    def bond_source(self):
        """Tag of the bond-type convention (file format, or None)."""
        return self._bond_source

    def __len__(self):
        return self.n_atoms

    # -- output -------------------------------------------------------------

    def write(self, path, comment=None):
        """
        Write this molecule, with its bond connectivity, using chemfiles.
        The format is inferred from the file extension (e.g. .xyz, .sdf,
        .mol2, .pdb), as when reading. Bond types are not written.
        """
        frame = _build_frame(
            self._coords, self._symbols, comment=comment,
            bond_data=self._bond_data if self.has_bonds else None,
        )
        _write_frame(path, frame)

    # -- comparison helpers -------------------------------------------------

    def _check_other(self, other):
        if not isinstance(other, Molecule):
            raise TypeError("Expected Molecule, got {}".format(type(other).__name__))

    def _check_bond_types(self, other, use_bond_type, bonding_flag):
        """Bond types are only comparable within the same bond_source."""
        if use_bond_type and not bonding_flag and self._bond_source != other._bond_source:
            raise ValueError(
                "use_bond_type=True requires both molecules to have the same bond "
                "source (got {!r} and {!r})".format(self._bond_source, other._bond_source)
            )

    # -- atormsd ------------------------------------------------------------

    def atormsd_to(
        self,
        other,
        align=False,
        remap=False,
        heavy_only=False,
        mass_weight=False,
        mirror=False,
        use_atom_type=False,
        print_stats=False,
        random=False,
        prune_tol=None,
        ato_freq=10,
        max_trials=10000,
        max_records=1,
    ):
        """
        Compare as unstructured atom clusters (atormsd_calculate); bonds are
        ignored. Returns up to max_records RMSDResult, best (lowest RMSD) first.

        Records beyond the first are only ever produced when both align=True
        and remap=True; otherwise the returned list always has length 1
        regardless of max_records. The list may be shorter than max_records
        if the library did not find that many distinct solutions.

        With heavy_only=True the two molecules may differ in their number of
        hydrogens. The smaller one is then padded with dummy atoms, and each
        mapping has max(len(self), len(other)) entries (see
        RMSDResult.mapping). Hydrogens are not part of the RMSD; they are
        paired afterwards by distance.

        prune_tol (Å) enables pruning: two atoms are never paired if their
        sorted distances to the atoms of some atom type differ by more than
        2*sqrt(3)*prune_tol. Defaults to None, which disables pruning.

        The search over random orientations stops once the best solution
        has been found ato_freq times, or after max_trials orientations
        (ato_freq, max_trials and max_records must be >= 1).
        random=True seeds it from the clock (otherwise results are
        reproducible), and print_stats=True prints its statistics.
        """

        self._check_other(other)
        pruning_flag = prune_tol is not None

        rmsd_vals, maps, tfs = _molalign.atormsd_calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            massweight_flag=mass_weight,
            mirror_flag=mirror,
            useatomtype_flag=use_atom_type,
            printstats_flag=print_stats,
            random_flag=random,
            prune_tol=prune_tol,
            pruning_flag=pruning_flag,
            ato_freq=ato_freq,
            max_trials=max_trials,
            max_records=max_records,
        )
        return [
            RMSDResult(rmsd=rmsd_vals[i], mapping=maps[i], transform=tfs[i])
            for i in range(len(rmsd_vals))
        ]

    # -- conformsd ----------------------------------------------------------

    def conformsd_to(
        self,
        other,
        align=False,
        remap=False,
        heavy_only=False,
        mass_weight=False,
        mirror=False,
        use_atom_type=False,
        bond_tol=None,
        use_bond_type=False,
        print_stats=False,
        print_assignment_tree=False,
        random=False,
        confo_freq=100,
        max_trials=10000,
        max_fragments=1,
        max_records=1,
    ):
        """
        Compare as conformers (conformsd_calculate): atom assignments always
        respect the bond topology, so both molecules must have the same bond
        graph. Returns up to max_records RMSDResult, best (lowest RMSD) first.

        Records beyond the first are only ever produced when both align=True
        and remap=True; otherwise the returned list always has length 1
        regardless of max_records. The list may be shorter than max_records
        if the library did not find that many distinct solutions.

        With heavy_only=True the two molecules may differ in their number of
        hydrogens. The smaller one is then padded with dummy atoms, and each
        mapping has max(len(self), len(other)) entries (see
        RMSDResult.mapping). Hydrogens are not part of the RMSD;
        they are paired afterwards, following their heavy neighbour where
        bonds are known and by distance otherwise.

        bond_tol (Å) enables bond detection: connectivity is inferred from
        geometry with this tolerance instead of using each molecule's bond
        table. Defaults to None, which uses the bond tables.

        max_fragments (>= 1, default 1) is the maximum number of molecular
        fragments (connected components of the bond graph of the compared
        atoms, i.e. after heavy_only exclusions) allowed in each molecule;
        a ValueError is raised if either molecule has more. An atom without
        bonds is a fragment of its own, so with the default a molecule of
        two or more atoms must have bonds connecting all of them.

        use_bond_type=True also uses the bond types to guide atom
        matching. Types are compared, never interpreted, so both molecules
        must have the same bond_source (e.g. both read from files of the
        same format); otherwise a ValueError is raised. It has no effect
        when bond_tol is given, since inferred bonds are untyped.

        With align=True, the search strategy is chosen from the assignment
        tree: random orientations when the total number of assignments
        exceeds confo_freq times the sum of partial combinations, exhaustive
        enumeration otherwise. The random search stops once the best
        solution has been found more than confo_freq times, or after
        max_trials orientations. confo_freq, max_trials and max_records must
        be >= 1. The default confo_freq=100 is the value
        validated on the CCD and BIRD benchmarks. random=True seeds the
        search from the clock (otherwise results are reproducible),
        print_stats=True prints its statistics and print_assignment_tree=True prints the
        assignment tree and its combination counts.
        """

        self._check_other(other)
        bonding_flag = bond_tol is not None
        self._check_bond_types(other, use_bond_type, bonding_flag)
        bond_data1 = None if bonding_flag else self._bond_data
        bond_data2 = None if bonding_flag else other._bond_data

        rmsd_vals, maps, tfs = _molalign.conformsd_calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            bond_data1=bond_data1,
            bond_data2=bond_data2,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            massweight_flag=mass_weight,
            mirror_flag=mirror,
            useatomtype_flag=use_atom_type,
            usebondtype_flag=use_bond_type,
            bond_tol=bond_tol,
            bonding_flag=bonding_flag, 
            printstats_flag=print_stats,
            printassigntree_flag=print_assignment_tree,
            random_flag=random,
            confo_freq=confo_freq,
            max_trials=max_trials,
            max_fragments=max_fragments,
            max_records=max_records,
        )
        return [
            RMSDResult(rmsd=rmsd_vals[i], mapping=maps[i], transform=tfs[i])
            for i in range(len(rmsd_vals))
        ]

    def __repr__(self):
        return "Molecule(name={!r}, n_atoms={}, n_bonds={})".format(
            self._name, self.n_atoms, self.n_bonds
        )
