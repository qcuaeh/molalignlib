#!/usr/bin/env python
# coding: utf-8

import molalign
import numpy as np
from ase import Atoms

#print(molalign.library.remap.__doc__)
#print(molalign.library.align.__doc__)

class Aligner(Atoms):
    def __init__(self, atoms, records=1, count=None, trials=None, biasing=False, \
                 bias_scale=1000., bias_tol=None, weighted=False, testing=False):
        if not isinstance(atoms, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if type(records) is not int:
            raise TypeError('"records" must be integer')
        if type(bias_scale) is not float:
            raise TypeError('"bias_scale" must be real')
        if type(weighted) is not bool:
            raise TypeError('"weighted" must be boolean')
        if type(testing) is not bool:
            raise TypeError('"testing" must be boolean')
        if type(biasing) is not bool:
            raise TypeError('"biasing" must be boolean')
        if count is not None and type(count) is not int:
            raise TypeError('"count" must be integer')
        if trials is not None and type(trials) is not int:
            raise TypeError('"trials" must be integer')
        if bias_tol is not None and type(bias_tol) is not float:
            raise TypeError('"bias_tol" must be real')
        self.__dict__.update(atoms.__dict__)
        self.records = records
        molalign.options.test_flag = testing
        molalign.options.bias_flag = biasing
        molalign.options.bias_scale = bias_scale
        if count is None and trials is None:
            raise ValueError('either "count" or "trials" must be set')
        if count is None:
            molalign.options.converge = False
        else:
            molalign.options.converge = True
            molalign.options.maxcount = count
        if trials is None:
            molalign.options.halt_flag = False
        else:
            molalign.options.halt_flag = True
            molalign.options.maxtrials = trials
        if biasing:
            if bias_tol is None:
                raise ValueError('"biasing" is True but "bias_tol" is not set')
            else:
                molalign.options.bias_tol = bias_tol
        if weighted:
            self.weights = atoms.get_masses()/sum(atoms.get_masses())
        else:
            self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
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

