#!/usr/bin/env python
# coding: utf-8

import molalign
import numpy as np
from ase import Atoms

#print(molalign.library.remap.__doc__)
#print(molalign.library.align.__doc__)

class Aligner(Atoms):
    def __init__(self, atoms, records=1, scale=1000., count=None, trial=None, bias=None, \
                 iteration=False, weighted=False, testing=False):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            print('An Atoms object was expected as argument')
            raise SystemExit
        if type(records) is int:
            self.records = records
        else:
            print('"records" must be an integer')
            raise SystemExit
        if type(scale) is float:
            molalign.options.lenscale = scale
        else:
            print('"scale" must be a real')
            raise SystemExit
        if count is None:
            molalign.options.converged = False
        else:
            if type(count) is int:
                molalign.options.converged = True
                molalign.options.mincount = count
            else:
                print('"count" must be an integer')
                raise SystemExit
        if trial is None:
            molalign.options.bounded = False
        else:
            if type(trial) is int:
                molalign.options.bounded = True
                molalign.options.maxtrial = trial
            else:
                print('"trial" must be an integer')
                raise SystemExit
        if type(iteration) is bool:
            molalign.options.iteration = iteration
        else:
            print('"iteration" must be a boolean')
            raise SystemExit
        if type(testing) is bool:
            molalign.options.testing = testing
        else:
            print('"testing" must be a boolean')
            raise SystemExit
        if bias is None:
            molalign.options.biased = False
        else:
            if type(bias) is float:
                molalign.options.biased = True
                molalign.options.tolerance = bias
            else:
                print('"bias" must be a real')
                raise SystemExit
        if type(weighted) is bool:
            if weighted:
                self.weights = atoms.get_masses()/sum(atoms.get_masses())
            else:
                self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        else:
            print('"weighted" must be a boolean')
            raise SystemExit
    def remapping(self, other):
        if not isinstance(other, Atoms):
            print('An Atoms object was expected as argument')
            raise SystemExit
        if len(other) != len(self):
            print('Argument does no have the right length')
            raise SystemExit
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()
        n, maplist, mapcount = molalign.library.remap(znums0, znums1, types0, types1, \
            coords0, coords1, self.weights, self.records)
        return [i - 1 for i in maplist.transpose()[:n]], mapcount[:n]
    def aligned(self, other, mapping):
        if not isinstance(other, Atoms):
            print('An Atoms object was expected as first argument')
            raise SystemExit
        if len(other) != len(self):
            print('First argument does no have the right length')
            raise SystemExit
        if type(mapping) is not np.ndarray or mapping.dtype is not np.dtype('int32'):
            print('An integer numpy array was expected as second argument')
            raise SystemExit
        if len(mapping) != len(self):
            print('Second argument does no have the right length')
            raise SystemExit
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()[mapping]
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()[:, mapping]
        travec, rotmat = molalign.library.align(znums0, znums1, types0, types1, coords0, \
             coords1, self.weights)
        coords1 = np.matmul(rotmat, coords1).transpose() + travec
        return Atoms(numbers=znums1, positions=coords1)

