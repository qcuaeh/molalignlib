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


# Align atoms1 to atoms0 for each assignment,
# print count, RMSD from assignment, RMSD from alignment
# and write aligned coordinates to file
for i, assignment in enumerate(assignments, start=1):
    aligned1 = align0(atoms1[assignment.map])
    print(i, assignment.count, assignment.rmsd, aligned1.rmsd)
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), aligned1.atoms, append=True)

