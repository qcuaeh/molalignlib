"""
High-level Python wrappers around the molalign Cython extension.

Supported file formats (via chemfiles):
  - XYZ        (.xyz)
  - SDF / MOL  (.sdf, .mol)
  - MOL2       (.mol2)
  - and many others

Classes
-------
Atoms
    Represents an unstructured set of atoms (no bond topology).
    Uses atormsd_calculate under the hood.

Conformer
    Represents a molecule with bond topology (a conformer).
    Uses conformsd_calculate under the hood.

Both expose a .rmsd_to(other, n_records=1, **kwargs) method that returns a
list of RMSDResult objects, best (lowest RMSD) first. By default only the
single best solution is computed (a list of length 1); pass a larger
n_records to get several ranked candidate solutions at once. The returned
list may be shorter than n_records if the library did not find that many
distinct solutions.
"""

from pathlib import Path
import numpy as np
import chemfiles

from . import molalign as _molalign


# ---------------------------------------------------------------------------
# Element symbol <-> atomic number table
#
# Mirrors atomic_symbols(:) in chemdata.f90 (index i -> atomic number i+1),
# so this stays in lock-step with the Fortran core's own numbering. Entries
# 104/105 ("X", "LJ") are the library's dummy / Lennard-Jones placeholder
# "elements".
# ---------------------------------------------------------------------------

ATOMIC_SYMBOLS = (
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
    "Md", "No", "Lr", "X", "LJ",
)

# Case-insensitive symbol -> atomic number lookup (1-based).
_SYMBOL_TO_NUMBER = {sym.upper(): i + 1 for i, sym in enumerate(ATOMIC_SYMBOLS)}


def symbol_to_atomic_number(symbol):
    """Look up the atomic number for an element symbol (case-insensitive)."""
    key = symbol.strip().upper()
    try:
        return _SYMBOL_TO_NUMBER[key]
    except KeyError:
        raise ValueError("Unrecognized element symbol: {!r}".format(symbol))


def atomic_number_to_symbol(number):
    """Look up the element symbol for an atomic number, if known."""
    number = int(number)
    if 1 <= number <= len(ATOMIC_SYMBOLS):
        return ATOMIC_SYMBOLS[number - 1]
    return str(number)


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

class RMSDResult(object):
    """Return value of .rmsd_to()."""

    def __init__(self, rmsd, atom_permutation, transform):
        self.rmsd = rmsd
        """Root-mean-square deviation (Angstrom)."""
        self.atom_permutation = atom_permutation
        """0-based index array mapping other's atoms onto self's atoms."""
        self.transform = transform
        """4x4 homogeneous rotation+translation matrix (maps other -> self frame)."""

    def apply_to(self, cluster_or_conformer):
        """
        Apply the stored transform and atom permutation to a cluster or conformer,
        returning a new object aligned and reordered to match the reference.
        """
        R = self.transform[:3, :3]
        t = self.transform[:3, 3]
        idx = self.atom_permutation

        new_coords = cluster_or_conformer.coords[idx] @ R.T + t
        new_atomdata = cluster_or_conformer.atom_data[idx]
        new_symbols = [cluster_or_conformer.symbols[i] for i in idx]

        if isinstance(cluster_or_conformer, Conformer):
            inv_idx = np.argsort(idx)
            new_bonddata = cluster_or_conformer.bond_data.copy()
            
            for i in range(new_bonddata.shape[0]):
                new_bonddata[i, 0] = inv_idx[new_bonddata[i, 0] - 1] + 1
                new_bonddata[i, 1] = inv_idx[new_bonddata[i, 1] - 1] + 1
                
            c = Conformer(
                atom_data=new_atomdata,
                coords=new_coords,
                bond_data=new_bonddata,
                name=cluster_or_conformer.name,
                symbols=new_symbols
            )
        else:
            c = Atoms(
                atom_data=new_atomdata,
                coords=new_coords,
                name=cluster_or_conformer.name,
                symbols=new_symbols
            )
        return c

    def __repr__(self):
        return "RMSDResult(rmsd={:.6f})".format(self.rmsd)


# ---------------------------------------------------------------------------
# Chemfiles helpers
# ---------------------------------------------------------------------------

def _extract_frame_data(frame):
    """
    Extract atom_data, coords, bond_data, and symbols from a chemfiles Frame.
    """
    n_atoms = len(frame.atoms)
    atom_data = np.zeros((n_atoms, 2), dtype=np.int32)
    symbols = []
    
    for i, atom in enumerate(frame.atoms):
        atom_data[i, 0] = atom.atomic_number
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
    1-based bond_data (as stored on Conformer), for writing out through
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


# ---------------------------------------------------------------------------
# Public frame readers
# ---------------------------------------------------------------------------

def read_clusters(path, frames=None):
    """
    Read frames from a molecular file and return Atoms objects.
    If 'frames' is provided (e.g., a tuple of indices), it returns a tuple 
    of only those specific frames. Otherwise, it returns a list of all frames.
    Bond information is ignored.
    """
    stem = Path(path).stem
    clusters = []
    
    with chemfiles.Trajectory(str(path)) as trajectory:
        if frames is not None:
            for idx in frames:
                frame = trajectory.read_step(idx)
                atom_data, coords, _, symbols = _extract_frame_data(frame)
                name = "{}_{}".format(stem, idx)
                clusters.append(Atoms(
                    atom_data=atom_data, coords=coords, name=name, symbols=symbols
                ))
            return tuple(clusters)
        else:
            for idx, frame in enumerate(trajectory):
                atom_data, coords, _, symbols = _extract_frame_data(frame)
                name = "{}_{}".format(stem, idx)
                clusters.append(Atoms(
                    atom_data=atom_data, coords=coords, name=name, symbols=symbols
                ))
            return clusters

def read_conformers(path, frames=None):
    """
    Read frames from a molecular file and return Conformer objects.
    If 'frames' is provided (e.g., a tuple of indices), it returns a tuple 
    of only those specific frames. Otherwise, it returns a list of all frames.
    Bond information is parsed and used.
    """
    stem = Path(path).stem
    conformers = []
    
    with chemfiles.Trajectory(str(path)) as trajectory:
        if frames is not None:
            for idx in frames:
                frame = trajectory.read_step(idx)
                atom_data, coords, bond_data, symbols = _extract_frame_data(frame)
                name = "{}_{}".format(stem, idx)
                conformers.append(Conformer(
                    atom_data=atom_data, coords=coords, bond_data=bond_data, 
                    name=name, symbols=symbols
                ))
            return tuple(conformers)
        else:
            for idx, frame in enumerate(trajectory):
                atom_data, coords, bond_data, symbols = _extract_frame_data(frame)
                name = "{}_{}".format(stem, idx)
                conformers.append(Conformer(
                    atom_data=atom_data, coords=coords, bond_data=bond_data, 
                    name=name, symbols=symbols
                ))
            return conformers


# ---------------------------------------------------------------------------
# Atoms
# ---------------------------------------------------------------------------

class Atoms(object):
    """
    An unstructured cluster of atoms with 3-D coordinates but no bond topology.
    """

    def __init__(self, atom_data=None, coords=None, name="cluster", symbols=None, labels=None):
        """
        Either ``atom_data`` (atomic numbers, as before) or ``symbols``
        (element symbols, e.g. ["C", "H", "H", "H"]) must be provided; the
        other is derived automatically. Providing both is fine too -- in
        that case ``atom_data`` wins and ``symbols`` is only used for
        display/labelling.

        labels : sequence of int, optional
            Per-atom labels (see molalign.h); only used when building
            ``atom_data`` from ``symbols``. Defaults to all-zero
            (unlabelled). Ignored if ``atom_data`` is given directly.
        """
        self._coords = np.asarray(coords, dtype=np.float64)
        self._name = name

        if atom_data is None:
            if symbols is None:
                raise ValueError("Either atom_data or symbols must be provided")
            elnums = [symbol_to_atomic_number(s) for s in symbols]
            if labels is None:
                labels = [0] * len(elnums)
            elif len(labels) != len(elnums):
                raise ValueError("labels must have the same length as symbols")
            self._atom_data = np.array(list(zip(elnums, labels)), dtype=np.int32)
            self._symbols = list(symbols)
        else:
            self._atom_data = np.asarray(atom_data, dtype=np.int32)
            if symbols is None:
                self._symbols = [atomic_number_to_symbol(a[0]) for a in self._atom_data]
            else:
                self._symbols = list(symbols)

        if self._coords.shape[0] != self._atom_data.shape[0]:
            raise ValueError("coords and atom_data/symbols must have the same length")

    @classmethod
    def from_file(cls, path, frame_idx=0):
        """Read a single frame from a file as an Atoms."""
        with chemfiles.Trajectory(str(path)) as trajectory:
            frame = trajectory.read_step(frame_idx)
            atom_data, coords, _, symbols = _extract_frame_data(frame)
            name = Path(path).stem
            return cls(atom_data=atom_data, coords=coords, name=name, symbols=symbols)

    @classmethod
    def from_symbols(cls, symbols, coords, name="cluster", labels=None):
        """Build an Atoms directly from element symbols and coordinates."""
        return cls(coords=coords, name=name, symbols=symbols, labels=labels)

    @classmethod
    def from_numbers(cls, atomic_numbers, coords, name="cluster", labels=None):
        """Build an Atoms directly from atomic numbers and coordinates."""
        atomic_numbers = list(atomic_numbers)
        if labels is None:
            labels = [0] * len(atomic_numbers)
        elif len(labels) != len(atomic_numbers):
            raise ValueError("labels must have the same length as atomic_numbers")
        atom_data = np.array(list(zip(atomic_numbers, labels)), dtype=np.int32)
        return cls(atom_data=atom_data, coords=coords, name=name)

    @property
    def name(self):
        return self._name

    @property
    def n_atoms(self):
        return self._atom_data.shape[0]

    @property
    def atom_data(self):
        return self._atom_data

    @property
    def coords(self):
        return self._coords

    @property
    def symbols(self):
        return self._symbols

    def write(self, path, comment=None):
        """
        Write this cluster using chemfiles. The output format is inferred
        from the file extension (e.g. .xyz, .pdb, .sdf, .mol2, ...), the
        same way read_clusters()/from_file() infer the input format.
        """
        frame = _build_frame(self._coords, self._symbols, comment=comment)
        _write_frame(path, frame)

    def rmsd_to(
        self,
        other,
        align=False,
        remap=False,
        heavy_only=False,
        mass_weighted=False,
        mirror=False,
        use_labels=False,
        stats=False,
        random=False,
        prune=False,
        prune_tol=None,
        conv_freq=10,
        max_trials=10000,
        n_records=1,
    ):
        """
        Compute up to n_records ranked candidate solutions.

        Returns a list of RMSDResult, best (lowest RMSD) first. Records
        beyond the first are only ever produced when both align=True and
        remap=True; otherwise the returned list always has length 1
        regardless of n_records. The list may be shorter than n_records
        if the library did not find that many distinct solutions.
        """
        if not isinstance(other, Atoms):
            raise TypeError("Expected Atoms, got {}".format(type(other).__name__))

        # prune_tol has no default: it is required whenever prune=True.
        if prune and prune_tol is None:
            raise ValueError("prune_tol is required when prune=True")

        rmsd_vals, perms, tfs = _molalign.atormsd_calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            mass_flag=mass_weighted,
            mirror_flag=mirror,
            label_flag=use_labels,
            print_stats=stats,
            random_flag=random,
            prune_flag=prune,
            prune_tol=prune_tol,
            conv_freq=conv_freq,
            max_trials=max_trials,
            n_records=n_records,
        )
        return [
            RMSDResult(rmsd=rmsd_vals[i], atom_permutation=perms[i], transform=tfs[i])
            for i in range(len(rmsd_vals))
        ]

    def __repr__(self):
        return "Atoms(name={!r}, n_atoms={})".format(self._name, self.n_atoms)


# ---------------------------------------------------------------------------
# Conformer
# ---------------------------------------------------------------------------

class Conformer(object):
    """
    A molecule with full bond topology and 3-D coordinates (a conformer).
    """

    def __init__(self, atom_data=None, coords=None, bond_data=None, name="conformer",
                 symbols=None, labels=None):
        """
        Either ``atom_data`` (atomic numbers, as before) or ``symbols``
        (element symbols, e.g. ["C", "H", "H", "H"]) must be provided; the
        other is derived automatically. Providing both is fine too -- in
        that case ``atom_data`` wins and ``symbols`` is only used for
        display/labelling.

        labels : sequence of int, optional
            Per-atom labels (see molalign.h); only used when building
            ``atom_data`` from ``symbols``. Defaults to all-zero
            (unlabelled). Ignored if ``atom_data`` is given directly.
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
            self._atom_data = np.array(list(zip(elnums, labels)), dtype=np.int32)
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
            self._bond_data = np.asarray(bond_data, dtype=np.int32)
        else:
            self._bond_data = np.empty((0, 3), dtype=np.int32)

        self._name = name

    @classmethod
    def from_file(cls, path, frame_idx=0):
        """Read a single frame from a file as a Conformer."""
        with chemfiles.Trajectory(str(path)) as trajectory:
            frame = trajectory.read_step(frame_idx)
            atom_data, coords, bond_data, symbols = _extract_frame_data(frame)
            name = Path(path).stem
            return cls(atom_data=atom_data, coords=coords, bond_data=bond_data, name=name, symbols=symbols)

    @classmethod
    def from_symbols(cls, symbols, coords, bond_data=None, name="conformer", labels=None):
        """Build a Conformer directly from element symbols, coordinates, and bonds."""
        return cls(coords=coords, bond_data=bond_data, name=name, symbols=symbols, labels=labels)

    @classmethod
    def from_numbers(cls, atomic_numbers, coords, bond_data=None, name="conformer", labels=None):
        """Build a Conformer directly from atomic numbers, coordinates, and bonds."""
        atomic_numbers = list(atomic_numbers)
        if labels is None:
            labels = [0] * len(atomic_numbers)
        elif len(labels) != len(atomic_numbers):
            raise ValueError("labels must have the same length as atomic_numbers")
        atom_data = np.array(list(zip(atomic_numbers, labels)), dtype=np.int32)
        return cls(atom_data=atom_data, coords=coords, bond_data=bond_data, name=name)

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

    def write(self, path, comment=None):
        """
        Write this conformer, including bond connectivity, using chemfiles.
        The output format is inferred from the file extension (e.g. .sdf,
        .mol2, .pdb, .xyz, ...), the same way read_conformers()/from_file()
        infer the input format. Bond orders are not currently preserved on
        write, only connectivity.
        """
        frame = _build_frame(
            self._coords, self._symbols, comment=comment, bond_data=self._bond_data
        )
        _write_frame(path, frame)

    def rmsd_to(
        self,
        other,
        align=False,
        remap=False,
        heavy_only=False,
        mass_weighted=False,
        mirror=False,
        use_labels=False,
        infer_bonds=False,
        bond_tol=0.3,
        stats=False,
        random=False,
        conv_freq=100,
        max_trials=10000,
        n_records=1,
    ):
        """
        Compute up to n_records ranked candidate solutions.

        Returns a list of RMSDResult, best (lowest RMSD) first. Records
        beyond the first are only ever produced when both align=True and
        remap=True; otherwise the returned list always has length 1
        regardless of n_records. The list may be shorter than n_records
        if the library did not find that many distinct solutions.

        bond_tol (default 0.3 Å) is the bond-detection tolerance used when
        infer_bonds=True; it is ignored otherwise.
        """
        if not isinstance(other, Conformer):
            raise TypeError("Expected Conformer, got {}".format(type(other).__name__))

        rmsd_vals, perms, tfs = _molalign.conformsd_calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            bond_data1=self._bond_data if not infer_bonds else None,
            bond_data2=other._bond_data if not infer_bonds else None,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            mass_flag=mass_weighted,
            mirror_flag=mirror,
            label_flag=use_labels,
            bond_flag=infer_bonds,
            bond_tol=bond_tol,
            print_stats=stats,
            random_flag=random,
            conv_freq=conv_freq,
            max_trials=max_trials,
            n_records=n_records,
        )
        return [
            RMSDResult(rmsd=rmsd_vals[i], atom_permutation=perms[i], transform=tfs[i])
            for i in range(len(rmsd_vals))
        ]

    def __repr__(self):
        return "Conformer(name={!r}, n_atoms={}, n_bonds={})".format(
            self._name, self.n_atoms, self.n_bonds
        )