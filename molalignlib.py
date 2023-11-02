# MolAlignLib
# Copyright (C) 2022 José M. Vásquez

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import molalignlibext
from ase import Atoms

#print(molalignlibext.library.align_atoms.__doc__)
#print(molalignlibext.library.assign_atoms.__doc__)

class Assignment:
    def __init__(self, count, order):
        self.count = count
        self.order = order

def assign_atoms(
    mol0,
    mol1,
    fast = False,
    test = False,
    stats = False,
    rec = None,
    tol = 0.35,
    count = 10,
    trials = None,
    # dummy arguments
    file = None,
    reorder = None,
    mass = None,
    mirror = None,
    out = None,
):
    if not isinstance(mol0, Atoms):
        raise TypeError('An Atoms object was expected')
    if not isinstance(mol1, Atoms):
        raise TypeError('An Atoms object was expected as argument')
    if fast is False:
        molalignlibext.flags.bias_flag = False
        molalignlibext.flags.iter_flag = False
    elif fast is True:
        molalignlibext.flags.bias_flag = True
        molalignlibext.flags.iter_flag = True
    else:
        raise TypeError('"fast" must be boolean')
    if type(test) is bool:
        molalignlibext.flags.test_flag = test
    else:
        raise TypeError('"test" must be boolean')
    if type(stats) is bool:
        molalignlibext.flags.stats_flag = stats
    else:
        raise TypeError('"stats" must be boolean')
    if rec is None:
        maxrec = 1
        multirec = False
    elif type(rec) is int:
        maxrec = rec
        multirec = True
    else:
        raise TypeError('"rec" must be an integer')
    if type(tol) is float:
        molalignlibext.biasing.bias_tol = tol
    else:
        raise TypeError('"tol" must be a float number')
    if type(count) is int:
        molalignlibext.bounds.maxcount = count
    else:
        raise TypeError('"count" must be an integer')
    if trials is None:
        molalignlibext.flags.trial_flag = False
    elif type(trials) is int:
        molalignlibext.flags.trial_flag = True
        molalignlibext.bounds.maxtrials = trials
    else:
        raise TypeError('"trials" must be an integer')
    natom0 = len(mol0)
    natom1 = len(mol1)
    molalignlibext.bounds.maxrec = maxrec
    molalignlibext.biasing.bias_scale = 1000.
    znums0 = mol0.get_atomic_numbers()
    types0 = np.ones(natom0, dtype=np.int32)
    coords0 = mol0.positions.T # Convert to column-major order
    znums1 = mol1.get_atomic_numbers()
    types1 = np.ones(natom1, dtype=np.int32)
    coords1 = mol1.positions.T # Convert to column-major order
#    adjmat0 = np.zeros((natom0, natom0), dtype=np.int8, order='F')
#    adjmat1 = np.zeros((natom1, natom1), dtype=np.int8, order='F')
    weights0 = mol0.get_weights()
    weights1 = mol1.get_weights()
#    permlist = np.empty((maxrec, natom0), dtype=np.int32).T
    permlist = np.empty((natom0, maxrec), dtype=np.int32, order='F')
    countlist = np.empty(maxrec, dtype=np.int32)
    nrec, error = \
        molalignlibext.library.assign_atoms(
            natom0,
            znums0,
            types0,
            weights0,
            coords0,
            natom1,
            znums1,
            types1,
            weights1,
            coords1,
            permlist,
            countlist,
        )
    if error:
        raise RuntimeError('Assignment failed')
    if multirec:
        return [Assignment(countlist[i], permlist[:, i] - 1) for i in range(nrec)]
    else:
        return Assignment(countlist[0], permlist[:, 0] - 1)

def align_to(self, other):
    if not isinstance(other, Atoms):
        raise TypeError('An Atoms object was expected')
    natom0 = len(other)
    natom1 = len(self)
    znums0 = other.get_atomic_numbers()
    types0 = np.ones(natom0, dtype=np.int32)
    coords0 = other.positions.T # Convert to column-major order
    znums1 = self.get_atomic_numbers()
    types1 = np.ones(natom1, dtype=np.int32)
    coords1 = self.positions.T # Convert to column-major order
    weights0 = other.get_weights()
    weights1 = self.get_weights()
    molalignlibext.bounds.natom0 = natom0
    molalignlibext.bounds.natom1 = natom1
    travec, rotmat, error = \
        molalignlibext.library.align_atoms(
            natom0,
            znums0,
            types0,
            weights0,
            coords0,
            natom1,
            znums1,
            types1,
            weights1,
            coords1,
        )
    if error:
        raise RuntimeError('Alignment failed')
    aligned_positions = np.dot(self.positions, rotmat.T) + travec
    rmsd = (np.sum(weights0*np.sum((other.positions - aligned_positions)**2, axis=1))/np.sum(weights0))**0.5
    self.positions = aligned_positions
    return rmsd

def get_ones(self):
    return np.ones(len(self), dtype=np.float64)

def mirror_x(self):
    self.positions[:, 0] = -self.positions[:, 0]

# Monkey patch ASE's Atoms class
Atoms.align_to = align_to
Atoms.get_weights = get_ones
Atoms.mirror_x = mirror_x
