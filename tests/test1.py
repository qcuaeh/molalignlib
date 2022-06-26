import numpy as np
from ase import io
from ralign import options, superpose

atoms0 = io.read('r005/100cobalt_j5.xyz', index=0)
atoms1 = io.read('r005/100cobalt_j5.xyz', index=1)

natom = len(atoms0)

labels0 = np.empty((natom, 32), dtype='c')
for i, symbol in enumerate(atoms0.get_chemical_symbols()):
    labels0[i] = symbol + (32 - len(symbol))*' '

labels1 = np.empty((natom, 32), dtype='c')
for i, symbol in enumerate(atoms1.get_chemical_symbols()):
    labels1[i] = symbol + (32 - len(symbol))*' '

coords0 = np.transpose(atoms0.get_positions())
coords1 = np.transpose(atoms1.get_positions())
center0 = np.empty(3, dtype=float)
center1 = np.empty(3, dtype=float)
nrecord = np.empty((), dtype=int)
atomaplist = np.empty((natom, 10), dtype=int)
rotmatlist = np.empty((3, 3, 10), dtype=float)

maxrecord = 10
options.remap = False
options.biased = True
options.testing = True
options.matching = True
options.iterative = True
options.maxmatch = 10
options.tolerance = 0.017
options.weighter = 'none                            '

superpose(maxrecord, labels0, labels1, coords0, coords1)

