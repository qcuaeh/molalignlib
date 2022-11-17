#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import modules
from ase import io
from molalignlib import Alignment


# In[ ]:


# Read clusters coordinates
atoms0 = io.read('tests/r005/Co100.xyz', index=0)
atoms1 = io.read('tests/r005/Co100.xyz', index=1)


# In[ ]:


# Create an alignment object with atoms0 as reference
alignment0 = Alignment(atoms0, biasing=True, iteration=True)


# In[ ]:


# Sort atoms1 to minimize the RMSD respect to atoms0
maplist, mapcount, mindist = alignment0.sorted(atoms1)


# In[ ]:


# Align atoms1 to atoms0 for each calculated mapping and write coordinates to file
for i, mapping in enumerate(maplist, start=1):
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), alignment0.aligned(atoms1[mapping]), append=True)

