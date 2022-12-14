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
assign_to_atoms0 = Assign(atoms0, biasing=True, iteration=True, records=5)


# In[ ]:


# Assign atoms in atoms1 and atoms0 to minimize the RMSD
assignments = assign_to_atoms0(atoms1)


# In[ ]:


# Create an "align" object with atoms0 as reference
align_to_atoms0 = Align(atoms0)


# In[ ]:


# Align atoms1 to atoms0 for each assignment and write aligned coordinates to file
for i, (order, count, rmsd) in enumerate(assignments, start=1):
    print(i, count, rmsd)
    atoms1_aligned = align_to_atoms0(atoms1[order])
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), atoms1_aligned, append=True)

