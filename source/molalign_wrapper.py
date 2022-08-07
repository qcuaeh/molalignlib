#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
from ase import Atoms
import molalign

#print(molalign.library.remap.__doc__)
#print(molalign.library.align.__doc__)

def remap(atoms0, atoms1, weights, iteration=False, test=False, bias=None, count=None, trial=None, scale=1000., records=1):
    if isinstance(iteration, bool):
        molalign.options.iteration = iteration
    else:
        print('"iter" must be boolean')
        sys.exit()
    if isinstance(test, bool):
        molalign.options.testing = test
    else:
        print('"test" must be boolean')
        sys.exit()
    if isinstance(scale, float):
        molalign.options.lenscale = scale
    else:
        print('"scale" must be a real number')
        sys.exit()
    if bias is None:
        molalign.options.biased = False
    elif isinstance(bias, float):
        molalign.options.biased = True
        molalign.options.tolerance = bias
    else:
        print('"bias" must be a real number')
        sys.exit()
    if trial is None:
        molalign.options.bounded = False
    elif isinstance(trial, int):
        molalign.options.bounded = True
        molalign.options.maxtrial = trial
    else:
        print('"trial" must be an integer')
        sys.exit()
    if count is None:
        molalign.options.converged = False
    elif isinstance(count, int):
        molalign.options.converged = True
        molalign.options.mincount = count
    else:
        print('"count" must be an integer')
        sys.exit()
    znums0 = atoms0.get_atomic_numbers()
    znums1 = atoms1.get_atomic_numbers()
    types0 = np.ones(len(atoms0), dtype=int)
    types1 = np.ones(len(atoms1), dtype=int)
    coords0 = atoms0.get_positions().transpose()
    coords1 = atoms1.get_positions().transpose()
    n, maplist, mapcount = molalign.library.remap(znums0, znums1, types0, types1, coords0, coords1, weights/sum(weights), records)
    return [i - 1 for i in maplist.transpose()[:n]], mapcount[:n]

def align(atoms0, atoms1, weights, mapping):
    znums0 = atoms0.get_atomic_numbers()
    znums1 = atoms1.get_atomic_numbers()[mapping]
    types0 = np.ones(len(atoms0), dtype=int)
    types1 = np.ones(len(atoms1), dtype=int)
    coords0 = atoms0.get_positions().transpose()
    coords1 = atoms1.get_positions().transpose()[:, mapping]
    travec, rotmat = molalign.library.align(znums0, znums1, types0, types1, coords0, coords1, weights/sum(weights))
    return Atoms(numbers=znums1, positions=np.matmul(rotmat, coords1).transpose()+travec)

