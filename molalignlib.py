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
from ase import Atoms
from molalignlibext import settings, library

#print(library.align_atoms.__doc__)
#print(library.assign_atoms.__doc__)

class Alignment:
    def __init__(self, atoms0, atoms1, mw=False):
        if not isinstance(atoms0, Atoms):
            raise TypeError('An Atoms object was expected')
        if not isinstance(atoms1, Atoms):
            raise TypeError('An Atoms object was expected')
        znums0 = atoms0.get_atomic_numbers()
        types0 = np.ones(len(atoms0), dtype=np.int32)
        coords0 = atoms0.positions.T # Convert to column-major order
        znums1 = atoms1.get_atomic_numbers()
        types1 = np.ones(len(atoms1), dtype=np.int32)
        coords1 = atoms1.positions.T # Convert to column-major order
        if mw:
            weights0 = atoms0.get_masses()
            weights1 = atoms1.get_masses()
        else:
            weights0 = np.ones(len(atoms0), dtype=np.float64)
            weights1 = np.ones(len(atoms1), dtype=np.float64)
        travec, rotmat, dist2, error = \
            library.align_atoms(
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
        self.travec = travec
        self.rotmat = rotmat.T
        self.rmsd = dist2**0.5
    def align(self, atoms):
        return Atoms(
            numbers = atoms.numbers,
            positions = np.dot(atoms.positions, self.rotmat) + self.travec,
        )

class Assignment:
    def __init__(
        self,
        atoms0,
        atoms1,
        biasing = False,
        iteration = False,
        reproducible = False,
        stats = False,
        rec = 1,
        count = 10,
        trials = None,
        tol = 0.2,
        mw = False,
    ):
        if not isinstance(atoms0, Atoms):
            raise TypeError('An Atoms object was expected')
        if not isinstance(atoms1, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if not isinstance(biasing, bool):
            raise TypeError('"biasing" must be a boolean')
        if not isinstance(iteration, bool):
            raise TypeError('"iteration" must be a boolean')
        if not isinstance(stats, bool):
            raise TypeError('"stats" must be a boolean')
        if not isinstance(rec, int):
            raise TypeError('"rec" must be an integer')
        if not isinstance(count, int):
            raise TypeError('"count" must be an integer')
        if not isinstance(tol, float):
            raise TypeError('"tol" must be a real')
        if mw:
            weights0 = atoms0.get_masses()
            weights1 = atoms1.get_masses()
        else:
            weights0 = np.ones(len(atoms0), dtype=np.float64)
            weights1 = np.ones(len(atoms1), dtype=np.float64)
        if trials is None:
            settings.trial_flag = False
        else:
            if isinstance(trials, int):
                settings.trial_flag = True
                settings.maxtrials = trials
            else:
                raise TypeError('"trials" must be an integer')
        settings.bias_flag = biasing
        settings.iter_flag = iteration
        settings.repro_flag = reproducible
        settings.stats_flag = stats
        settings.maxcount = count
        settings.biastol = tol
        znums0 = atoms0.get_atomic_numbers()
        types0 = np.ones(len(atoms0), dtype=np.int32)
        coords0 = atoms0.positions.T # Convert to column-major order
        znums1 = atoms1.get_atomic_numbers()
        types1 = np.ones(len(atoms1), dtype=np.int32)
        coords1 = atoms1.positions.T # Convert to column-major order
        nmap, maplist, countlist, dist2list, error = \
            library.assign_atoms(
                znums0,
                types0,
                coords0,
                weights0,
                znums1,
                types1,
                coords1,
                weights1,
                rec,
            )
        if error:
            raise RuntimeError('Assignment failed')
        maplist = maplist - 1
        self.maplist = [maplist[:, i] for i in range(nmap)]
        self.rmsdlist = [dist2list[i]**0.5 for i in range(nmap)]
        self.countlist = [countlist[i] for i in range(nmap)]
    def __iter__(self):
        return iter(self.maplist)
