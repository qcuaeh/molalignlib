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
    atoms0,
    atoms1,
    biasing = None,
    testing = None,
    bonding = None,
    iteration = None,
    massweighted = None,
    stats = None,
    records = None,
    count = None,
    trials = None,
    tolerance = None,
    scale = None,
):
    if not isinstance(atoms0, Atoms):
        raise TypeError('An Atoms object was expected')
    if not isinstance(atoms1, Atoms):
        raise TypeError('An Atoms object was expected as argument')
    if biasing is None:
        biasing = False
    elif type(biasing) is not bool:
        raise TypeError('"biasing" must be a boolean')
    if testing is None:
        testing = False
    elif type(testing) is not bool:
        raise TypeError('"testing" must be a boolean')
    if bonding is None:
        bonding = False
    elif type(bonding) is not bool:
        raise TypeError('"bonding" must be a boolean')
    if iteration is None:
        iteration = False
    elif type(iteration) is not bool:
        raise TypeError('"iteration" must be a boolean')
    if massweighted is None:
        massweighted = False
    elif type(massweighted) is not bool:
        raise TypeError('"massweighted" must be a boolean')
    if stats is None:
        stats = False
    elif type(stats) is not bool:
        raise TypeError('"stats" must be a boolean')
    if records is None:
        records = 1
    elif type(records) is not int:
        raise TypeError('"records" must be an integer')
    if count is None:
        count = 10
    elif type(count) is not int:
        raise TypeError('"count" must be an integer')
    if tolerance is None:
        tolerance = 0.35
    elif type(tolerance) is not float:
        raise TypeError('"tolerance" must be a float number')
    if scale is None:
        scale = 1000.
    elif type(scale) is not float:
        raise TypeError('"scale" must be a float number')
    if trials is None:
        molalignlibext.flags.trial_flag = False
    else:
        if isinstance(trials, int):
            molalignlibext.flags.trial_flag = True
            molalignlibext.bounds.maxtrials = trials
        else:
            raise TypeError('"trials" must be an integer')
    natom0 = len(atoms0)
    natom1 = len(atoms1)
    molalignlibext.flags.bias_flag = biasing
    molalignlibext.flags.test_flag = testing
    molalignlibext.flags.bond_flag = bonding
    molalignlibext.flags.iter_flag = iteration
    molalignlibext.flags.stats_flag = stats
    molalignlibext.bounds.maxrec = records
    molalignlibext.bounds.maxcount = count
    molalignlibext.biasing.bias_tol = tolerance
    molalignlibext.biasing.bias_scale = scale
    znums0 = atoms0.get_atomic_numbers()
    types0 = np.ones(natom0, dtype=np.int32)
    coords0 = atoms0.positions.T # Convert to column-major order
    znums1 = atoms1.get_atomic_numbers()
    types1 = np.ones(natom1, dtype=np.int32)
    coords1 = atoms1.positions.T # Convert to column-major order
    adjmat0 = np.zeros((natom0, natom0), dtype=np.int8, order='F')
    adjmat1 = np.zeros((natom1, natom1), dtype=np.int8, order='F')
    if massweighted:
        weights0 = atoms0.get_masses()
        weights1 = atoms1.get_masses()
    else:
        weights0 = np.ones(natom0, dtype=np.float64)
        weights1 = np.ones(natom1, dtype=np.float64)
#    permlist = np.empty((records, natom0), dtype=np.int32).T
    permlist = np.empty((natom0, records), dtype=np.int32, order='F')
    countlist = np.empty(records, dtype=np.int32)
    nrec, error = \
        molalignlibext.library.assign_atoms(
            natom0,
            znums0,
            types0,
            weights0,
            coords0,
            adjmat0,
            natom1,
            znums1,
            types1,
            weights1,
            coords1,
            adjmat1,
            permlist,
            countlist,
        )
    if error:
        raise RuntimeError('Assignment failed')
    return [Assignment(countlist[i], permlist[:, i] - 1) for i in range(nrec)]

def align_to(self, other, massweighted=None):
    if not isinstance(other, Atoms):
        raise TypeError('An Atoms object was expected')
    if massweighted is None:
        massweighted = False
    elif type(massweighted) is not bool:
        raise TypeError('"massweighted" must be a boolean')
    natom0 = len(other)
    natom1 = len(self)
    znums0 = other.get_atomic_numbers()
    types0 = np.ones(natom0, dtype=np.int32)
    coords0 = other.positions.T # Convert to column-major order
    znums1 = self.get_atomic_numbers()
    types1 = np.ones(natom1, dtype=np.int32)
    coords1 = self.positions.T # Convert to column-major order
    if massweighted:
        weights0 = other.get_masses()
        weights1 = self.get_masses()
    else:
        weights0 = np.ones(natom0, dtype=np.float64)
        weights1 = np.ones(natom1, dtype=np.float64)
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

# Monkey patch Atoms class
Atoms.align_to = align_to
