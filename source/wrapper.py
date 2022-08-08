#!/usr/bin/env python
# coding: utf-8

import molalign
import numpy as np
from ase import Atoms

#print(molalign.library.remap.__doc__)
#print(molalign.library.align.__doc__)

class Aligner(Atoms):
    def __init__(self, atoms, rec=1, count=None, max=None, scale=1000., bias=None, \
                 iteration=False, mass_weighted=False, testing=False):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            raise TypeError('An Atoms object was expected as argument')
        if type(rec) is int:
            self.records = rec
        else:
            raise TypeError('"rec" must be an integer')
        if type(scale) is float:
            molalign.options.lenscale = scale
        else:
            raise TypeError('"scale" must be a real')
        if count is None:
            molalign.options.converge = False
        else:
            if type(count) is int:
                molalign.options.converge = True
                molalign.options.mincount = count
            else:
                raise TypeError('"count" must be an integer')
        if max is None:
            molalign.options.complete = False
        else:
            if type(max) is int:
                molalign.options.complete = True
                molalign.options.maxtrial = max
            else:
                raise TypeError('"max" must be an integer')
        if not molalign.options.converge and not molalign.options.complete:
            raise ValueError('either "count" or "max" must be defined')
        if type(iteration) is bool:
            molalign.options.iteration = iteration
        else:
            raise TypeError('"iteration" must be a boolean')
        if type(testing) is bool:
            molalign.options.testing = testing
        else:
            raise TypeError('"testing" must be a boolean')
        if bias is None:
            molalign.options.biased = False
        else:
            if type(bias) is float:
                molalign.options.biased = True
                molalign.options.tolerance = bias
            else:
                raise TypeError('"bias" must be a real')
        if type(mass_weighted) is bool:
            if mass_weighted:
                self.weights = atoms.get_masses()/sum(atoms.get_masses())
            else:
                self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        else:
            raise TypeError('"mass_weighted" must be a boolean')
    def remapping(self, other):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if len(other) != len(self):
            raise ValueError('Argument does no have the right length')
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()
        n, maplist, mapcount, mindist = molalign.library.remap(znums0, znums1, types0, \
            types1, coords0, coords1, self.weights, self.records)
        return [i - 1 for i in maplist.transpose()[:n]], mapcount[:n], mindist[:n]
    def aligned(self, other, mapping):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as first argument')
        if len(other) != len(self):
            raise ValueError('First argument does no have the right length')
        if type(mapping) is not np.ndarray or mapping.dtype is not np.dtype('int32'):
            raise TypeError('An integer numpy array was expected as second argument')
        if len(mapping) != len(self):
            raise ValueError('Second argument does no have the right length')
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

