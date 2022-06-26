import numpy as np
from ase import io
from ralign import options, superpose

def optstr(string):
    return string + ' '*(options.optlen - len(string))

atoms0 = io.read('r005/100cobalt_j5.xyz', index=0)
atoms1 = io.read('r005/100cobalt_j5.xyz', index=1)

znums0 = atoms0.get_atomic_numbers()
znums1 = atoms1.get_atomic_numbers()
types0 = np.zeros(len(atoms0), dtype=int)
types1 = np.zeros(len(atoms1), dtype=int)
coords0 = np.transpose(atoms0.get_positions())
coords1 = np.transpose(atoms1.get_positions())

options.remap = True
options.biased = True
options.testing = True
options.matching = True
options.iterative = True
options.maxmatch = 10
options.biasscale = 1000.0
options.tolerance = 0.17
options.weighter = optstr('none')

superpose(znums0, znums1, types0, types1, coords0, coords1, 10)
