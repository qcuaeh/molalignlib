#!/usr/bin/env python
# coding: utf-8

import numpy as np
from ase import Atoms
from molalignlib import library, options

#print(library.remap.__doc__)
#print(library.align.__doc__)

class Alignment(Atoms):
    def __init__(self, atoms, biased=False, weighted=False, testing=False, records=1,
                 count=10, trials=None, bias_scale=1000., bias_tol=0.2):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            raise TypeError('An Atoms object was expected as argument')
        if type(biased) is not bool:
            raise TypeError('"biased" must be boolean')
        if type(weighted) is not bool:
            raise TypeError('"weighted" must be boolean')
        if type(testing) is not bool:
            raise TypeError('"testing" must be boolean')
        if type(records) is not int:
            raise TypeError('"records" must be integer')
        if type(count) is not int:
            raise TypeError('"count" must be integer')
        if type(bias_tol) is not float:
            raise TypeError('"bias_tol" must be real')
        if type(bias_scale) is not float:
            raise TypeError('"bias_scale" must be real')
        if trials is None:
            options.halt_flag = False
        elif type(trials) is int:
            options.halt_flag = True
            options.maxtrials = trials
        else:
            raise TypeError('"trials" must be integer')
        if weighted:
            self.weights = atoms.get_masses()/sum(atoms.get_masses())
        else:
            self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        self.records = records
        options.test_flag = testing
        options.bias_flag = biased
        options.bias_scale = bias_scale
        options.bias_tol = bias_tol
        options.maxcount = count
    def sort(self, other):
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
        n, maplist, mapcount, mindist = library.remap(znums0, znums1, types0, \
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
        travec, rotmat = library.align(znums0, znums1, types0, types1, coords0, \
             coords1, self.weights)
        coords1 = np.matmul(rotmat, coords1).transpose() + travec
        return Atoms(numbers=znums1, positions=coords1)

