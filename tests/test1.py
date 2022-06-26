import numpy as np
from ase import io
from ralign import options, superpose

atoms0 = io.read('r005/100cobalt_j5.xyz', index=0)
atoms1 = io.read('r005/100cobalt_j5.xyz', index=1)

labels0 = np.empty((len(atoms0), 32), dtype='c')
for i, symbol in enumerate(atoms0.get_chemical_symbols()):
    labels0[i] = symbol + (32 - len(symbol))*' '

labels1 = np.empty((len(atoms1), 32), dtype='c')
for i, symbol in enumerate(atoms1.get_chemical_symbols()):
    labels1[i] = symbol + (32 - len(symbol))*' '

coords0 = np.transpose(atoms0.get_positions())
coords1 = np.transpose(atoms1.get_positions())

maxrecord = 10
options.remap = True
options.biased = True
options.testing = True
options.matching = True
options.iterative = True
options.maxmatch = 10
options.scale = 1000.0
options.tolerance = 0.17
options.weighter = 'none                            '

superpose(maxrecord, labels0, labels1, coords0, coords1)

