import numpy as np
from ase import io
from ralign import options, remap, align

def optstr(string):
    return string + ' '*(options.optlen - len(string))

atoms0 = io.read('r005/100cobalt_j5.xyz', index=0)
atoms1 = io.read('r005/100cobalt_j5.xyz', index=1)

znums0 = atoms0.get_atomic_numbers()
znums1 = atoms1.get_atomic_numbers()
types0 = np.ones(len(atoms0), dtype=int)
types1 = np.ones(len(atoms1), dtype=int)
weights0 = np.ones(len(atoms0), dtype=float)
weights1 = np.ones(len(atoms1), dtype=float)
coords0 = np.transpose(atoms0.get_positions())
coords1 = np.transpose(atoms1.get_positions())

# Normalize weights
weights0 = weights0/sum(weights0)
weights1 = weights1/sum(weights1)

options.biased = True
options.testing = True
options.counting = True
options.iterative = True
options.maxcount = 10
options.lenscale = 1000.0
options.tolerance = 0.17

remap(znums0, znums1, types0, types1, weights0, weights1, coords0, coords1, 10)
