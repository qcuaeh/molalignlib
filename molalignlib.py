# MolAlign
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
from ase import Atoms
from molalignlibext import molalignlib

class Alignable(Atoms):
    def __init__(
        self,
        atoms,
        mass_weighted = False,
    ):
        if not isinstance(atoms, Atoms):
            raise TypeError('An Atoms object was expected')
        self.__dict__.update(atoms.__dict__)
        self.mass_weighted = mass_weighted
    def alignto(self, atoms):
        if not isinstance(atoms, Atoms):
            raise TypeError('An Atoms object was expected')
        znums0 = atoms.get_atomic_numbers()
        types0 = np.ones(len(atoms), dtype=np.int32)
        coords0 = atoms.positions.T # Convert to column-major order
        znums1 = self.get_atomic_numbers()
        types1 = np.ones(len(self), dtype=np.int32)
        coords1 = self.positions.T # Convert to column-major order
        if self.mass_weighted:
            weights0 = atoms.get_masses()
            weights1 = self.get_masses()
        else:
            weights0 = np.ones(len(atoms), dtype=np.float64)
            weights1 = np.ones(len(self), dtype=np.float64)
        travec, rotmat, error = \
            molalignlib.align_atoms(
                znums0,
                types0,
                coords0,
                weights0,
                znums1,
                types1,
                coords1,
                weights1,
            )
        if error:
            raise RuntimeError('Alignment failed')
        self.positions = np.dot(self.positions, rotmat.T) + travec
    def disto(self, atoms):
        if self.mass_weighted:
            weights = self.get_masses()/self.get_masses().sum()
        else:
            weights = 1./len(self)
        return ((weights*((self.positions - atoms.positions)**2).sum(axis=1)).sum())**0.5

class Assignment:
    def __init__(
        self,
        atoms0,
        atoms1,
        biasing = False,
        bias_tol = 0.2,
        bias_scale = 1.e3,
        iteration = False,
        testing = False,
        records = 1,
        count = 10,
        trials = None,
        mass_weighted = False,
    ):
        if not isinstance(atoms0, Atoms):
            raise TypeError('An Atoms object was expected')
        if not isinstance(atoms1, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if type(biasing) is not bool:
            raise TypeError('"biasing" must be a boolean')
        if type(iteration) is not bool:
            raise TypeError('"iteration" must be a boolean')
        if type(testing) is not bool:
            raise TypeError('"testing" must be a boolean')
        if type(records) is not int:
            raise TypeError('"records" must be an integer')
        if type(count) is not int:
            raise TypeError('"count" must be an integer')
        if type(bias_tol) is not float:
            raise TypeError('"bias_tol" must be a real')
        if type(bias_scale) is not float:
            raise TypeError('"bias_scale" must be a real')
        if trials is None:
            molalignlib.set_free_flag(True)
        elif type(trials) is int:
            molalignlib.set_free_flag(False)
            molalignlib.set_max_trials(trials)
        else:
            raise TypeError('"trials" must be an integer')
        molalignlib.set_bias_flag(biasing)
        molalignlib.set_bias_scale(bias_scale)
        molalignlib.set_bias_tol(bias_tol)
        molalignlib.set_conv_flag(iteration)
        molalignlib.set_test_flag(testing)
        molalignlib.set_max_count(count)
        znums0 = atoms0.get_atomic_numbers()
        types0 = np.ones(len(atoms0), dtype=int)
        coords0 = atoms0.positions.T # Convert to column-major order
        znums1 = atoms1.get_atomic_numbers()
        types1 = np.ones(len(atoms1), dtype=int)
        coords1 = atoms1.positions.T # Convert to column-major order
        if mass_weighted:
            weights0 = atoms0.get_masses()
            weights1 = atoms1.get_masses()
        else:
            weights0 = np.ones(len(atoms0), dtype=np.float64)
            weights1 = np.ones(len(atoms1), dtype=np.float64)
        nmap, mapind, mapcount, mapdist, error = \
            molalignlib.assign_atoms(
                znums0,
                types0,
                coords0,
                weights0,
                znums1,
                types1,
                coords1,
                weights1,
                records,
            )
        if error:
            raise RuntimeError('Assignment failed')
        self.dist = list(mapdist[:nmap])
        self.count = list(mapcount[:nmap])
        self.mapind = [mapind[:, i] - 1 for i in range(nmap)]
