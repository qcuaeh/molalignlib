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
            print('"atoms" must be an Atoms object')
            raise SystemExit
        if isinstance(records, int):
            self.records = records
        else:
            print('"records" must be an integer')
            raise SystemExit
        if isinstance(scale, float):
            molalign.options.lenscale = scale
        else:
            print('"scale" must be a real')
            raise SystemExit
        if count is None:
            molalign.options.converged = False
        else:
            if isinstance(count, int):
                molalign.options.converged = True
                molalign.options.mincount = count
            else:
                print('"count" must be an integer')
                raise SystemExit
        if trial is None:
            molalign.options.bounded = False
        else:
            if isinstance(trial, int):
                molalign.options.bounded = True
                molalign.options.maxtrial = trial
            else:
                print('"trial" must be an integer')
                raise SystemExit
        if isinstance(iteration, bool):
            molalign.options.iteration = iteration
        else:
            print('"iteration" must be a boolean')
            raise SystemExit
        if isinstance(testing, bool):
            molalign.options.testing = testing
        else:
            print('"testing" must be a boolean')
            raise SystemExit
        if bias is None:
            molalign.options.biased = False
        else:
            if isinstance(bias, float):
                molalign.options.biased = True
                molalign.options.tolerance = bias
            else:
                print('"bias" must be a real')
                raise SystemExit
        if isinstance(weighted, bool):
            if weighted:
                self.weights = atoms.get_masses()/sum(atoms.get_masses())
            else:
                self.weights = np.ones(len(atoms), dtype=float)/len(atoms)
        else:
            print('"weighted" must be a boolean')
            raise SystemExit
    def remap(self, other):
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()
        n, maplist, mapcount = molalign.library.remap(znums0, znums1, types0, types1, \
            coords0, coords1, self.weights, self.records)
        return [i - 1 for i in maplist.transpose()[:n]], mapcount[:n]
    def align(self, other, mapping):
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()[mapping]
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()[:, mapping]
        travec, rotmat = molalign.library.align(znums0, znums1, types0, types1, coords0, \
             coords1, self.weights)
        return Atoms(numbers=znums1, positions=np.matmul(rotmat, coords1).transpose()+travec)

