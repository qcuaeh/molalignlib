#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sys
import numpy as np
from ase import io

sys.path.insert(0, '../lib')
from ralign import options, remap, align

print(remap.__doc__)
print(align.__doc__)


# In[4]:


def optstr(string):
    return string + ' '*(options.optlen - len(string))

atoms0 = io.read('r005/Co100.xyz', index=0)
atoms1 = io.read('r005/Co100.xyz', index=1)

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
options.converged = True
options.iterative = True
options.testing = True
options.mincount = 10
options.lenscale = 1000.0
options.tolerance = 0.17


# In[7]:


n, maplist, mapcount = remap(znums0, znums1, types0, types1, weights0, weights1, coords0, coords1, 10)


# In[55]:


for i, m in enumerate([i - 1 for i in maplist.transpose()[:n]]):
    travec, rotmat = align(znums0, znums1[m], types0, types1[m], weights0, weights1[m], coords0, coords1[:, m])
    atoms1.set_positions(np.matmul(rotmat, coords1).transpose() + travec)
    io.write('aligned_{}.xyz'.format(i + 1), atoms0,)
    io.write('aligned_{}.xyz'.format(i + 1), atoms1, append = True)


# In[ ]:




