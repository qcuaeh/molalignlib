#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import modules
from ase import io
from molalignlib import Align, Assign


# In[ ]:


# Read clusters coordinates
atoms0 = io.read('Co138_0.xyz', index=0)
atoms1 = io.read('Co138_1.xyz', index=0)


# In[ ]:


# Create an "assign" object with atoms0 as reference,
# biasing and iteration enabled and up to 5 recorded mappings
assign0 = Assign(atoms0, biasing=True, iteration=True, records=5)


# In[ ]:


# Assign atoms in atoms1 and atoms0 to minimize the RMSD
assignments = assign0(atoms1)


# In[ ]:


# Create an "align" object with atoms0 as reference
align0 = Align(atoms0)


# In[ ]:


# Align atoms in atoms1 and atoms0 for each calculated mapping and
# write coordinates to a file
for i, a in enumerate(assignments, start=1):
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), align0(atoms1[a.map]).atoms, append=True)

