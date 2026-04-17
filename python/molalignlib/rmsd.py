"""
High-level Python wrappers around the atormsd / conformsd Cython extensions.

Supported file formats (via chemfiles):
  - XYZ        (.xyz)
  - SDF / MOL  (.sdf, .mol)
  - MOL2       (.mol2)
  - and many others

Classes
-------
AtomCluster
    Represents an unstructured set of atoms (no bond topology).
    Uses atormsd_calculate under the hood.

Conformer
    Represents a molecule with bond topology (a conformer).
    Uses conformsd_calculate under the hood.

Both expose a .rmsd_to(other, **kwargs) method that returns an RMSDResult object.
"""

from pathlib import Path
import numpy as np
import chemfiles

from . import atormsd
from . import conformsd


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
            c = AtomCluster(
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


# ---------------------------------------------------------------------------
# Public frame readers
# ---------------------------------------------------------------------------

def read_clusters(path, frames=None):
    """
    Read frames from a molecular file and return AtomCluster objects.
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
                clusters.append(AtomCluster(
                    atom_data=atom_data, coords=coords, name=name, symbols=symbols
                ))
            return tuple(clusters)
        else:
            for idx, frame in enumerate(trajectory):
                atom_data, coords, _, symbols = _extract_frame_data(frame)
                name = "{}_{}".format(stem, idx)
                clusters.append(AtomCluster(
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
# AtomCluster
# ---------------------------------------------------------------------------

class AtomCluster(object):
    """
    An unstructured cluster of atoms with 3-D coordinates but no bond topology.
    """

    def __init__(self, atom_data, coords, name="cluster", symbols=None):
        self._atom_data = np.asarray(atom_data, dtype=np.int32)
        self._coords = np.asarray(coords, dtype=np.float64)
        self._name = name
        
        if symbols is None:
            self._symbols = [str(a[0]) for a in self._atom_data]
        else:
            self._symbols = symbols

    @classmethod
    def from_file(cls, path, frame_idx=0):
        """Read a single frame from a file as an AtomCluster."""
        with chemfiles.Trajectory(str(path)) as trajectory:
            frame = trajectory.read_step(frame_idx)
            atom_data, coords, _, symbols = _extract_frame_data(frame)
            name = Path(path).stem
            return cls(atom_data=atom_data, coords=coords, name=name, symbols=symbols)

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

    def write_xyz(self, path, comment=""):
        """Write the cluster to an XYZ file using chemfiles."""
        frame = chemfiles.Frame()
        for sym, (x, y, z) in zip(self._symbols, self._coords):
            frame.add_atom(chemfiles.Atom(sym), [x, y, z])
        
        with chemfiles.Trajectory(str(path), "w") as trajectory:
            trajectory.write(frame)

    def rmsd_to(
        self,
        other,
        align=True,
        remap=True,
        heavy_only=False,
        mass_weighted=False,
        mirror=False,
        use_labels=False,
        stats=False,
        random=False,
        prune_tol=-1.0,
        conv_freq=10,
        max_trials=10000,
    ):
        if not isinstance(other, AtomCluster):
            raise TypeError("Expected AtomCluster, got {}".format(type(other).__name__))

        rmsd_val, perm, tf = atormsd.calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            mass_flag=mass_weighted,
            mirror_flag=mirror,
            label_flag=use_labels,
            stats_flag=stats,
            random_flag=random,
            prune_tol=prune_tol,
            conv_freq=conv_freq,
            max_trials=max_trials,
        )
        return RMSDResult(rmsd=rmsd_val, atom_permutation=perm, transform=tf)

    def __repr__(self):
        return "AtomCluster(name={!r}, n_atoms={})".format(self._name, self.n_atoms)


# ---------------------------------------------------------------------------
# Conformer
# ---------------------------------------------------------------------------

class Conformer(object):
    """
    A molecule with full bond topology and 3-D coordinates (a conformer).
    """

    def __init__(self, atom_data, coords, bond_data=None, name="conformer", symbols=None):
        self._atom_data = np.asarray(atom_data, dtype=np.int32)
        self._coords = np.asarray(coords, dtype=np.float64)
        
        if bond_data is not None:
            self._bond_data = np.asarray(bond_data, dtype=np.int32)
        else:
            self._bond_data = np.empty((0, 3), dtype=np.int32)
            
        self._name = name
        
        if symbols is None:
            self._symbols = [str(a[0]) for a in self._atom_data]
        else:
            self._symbols = symbols

    @classmethod
    def from_file(cls, path, frame_idx=0):
        """Read a single frame from a file as a Conformer."""
        with chemfiles.Trajectory(str(path)) as trajectory:
            frame = trajectory.read_step(frame_idx)
            atom_data, coords, bond_data, symbols = _extract_frame_data(frame)
            name = Path(path).stem
            return cls(atom_data=atom_data, coords=coords, bond_data=bond_data, name=name, symbols=symbols)

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

    def rmsd_to(
        self,
        other,
        align=True,
        remap=True,
        heavy_only=False,
        mass_weighted=False,
        mirror=False,
        use_labels=False,
        bond_flag=False,
        stats=False,
        random=False,
        conv_freq=100,
        max_trials=10000,
    ):
        if not isinstance(other, Conformer):
            raise TypeError("Expected Conformer, got {}".format(type(other).__name__))

        rmsd_val, perm, tf = conformsd.calculate(
            self._atom_data, self._coords,
            other._atom_data, other._coords,
            bond_data1=self._bond_data if not bond_flag else None,
            bond_data2=other._bond_data if not bond_flag else None,
            align_flag=align,
            remap_flag=remap,
            heavy_flag=heavy_only,
            mass_flag=mass_weighted,
            mirror_flag=mirror,
            label_flag=use_labels,
            bond_flag=bond_flag,
            stats_flag=stats,
            random_flag=random,
            conv_freq=conv_freq,
            max_trials=max_trials,
        )
        return RMSDResult(rmsd=rmsd_val, atom_permutation=perm, transform=tf)

    def __repr__(self):
        return "Conformer(name={!r}, n_atoms={}, n_bonds={})".format(
            self._name, self.n_atoms, self.n_bonds
        )
