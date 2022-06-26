import numpy as np
from ase import io
from ralign import options, routines

atoms0 = io.read('r005/100cobalt_j5.xyz', index=0)
atoms1 = io.read('r005/100cobalt_j5.xyz', index=1)

def lblstr(string):
    return string + (options.lbllen - len(string))*' '

def optstr(string):
    return string + (options.optlen - len(string))*' '

labels0 = np.empty((len(atoms0), options.lbllen), dtype='c')
for i, symbol in enumerate(atoms0.get_chemical_symbols()):
    labels0[i] = lblstr(symbol)

labels1 = np.empty((len(atoms1), options.lbllen), dtype='c')
for i, symbol in enumerate(atoms1.get_chemical_symbols()):
    labels1[i] = lblstr(symbol)

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

routines.superpose(labels0, labels1, coords0, coords1, 10)

