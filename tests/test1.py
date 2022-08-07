#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import numpy as np
from ase import io
from molalign_wrapper import remap, align


# In[2]:


atoms0 = io.read('r005/Co100.xyz', index=0)
atoms1 = io.read('r005/Co100.xyz', index=1)
weights = np.ones(len(atoms0), dtype=float) # Unweighted
#weights = atoms0.get_masses() # Mass weighted
maplist, mapcount = remap(atoms0, atoms1, weights, count=10, bias=0.17, iteration=True, test=True, records=10)


# In[3]:


for i, mapping in enumerate(maplist, start=1):
    atoms1_aligned = align(atoms0, atoms1, weights, mapping)
    io.write('aligned_{}.xyz'.format(i), atoms0,)
    io.write('aligned_{}.xyz'.format(i), atoms1_aligned, append = True)

