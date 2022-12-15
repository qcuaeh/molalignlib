#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import modules
from ase import io
from molalignlib import Alignable, Assignment


# In[ ]:


# Read clusters coordinates
atoms0 = io.read('examples/Co138_0.xyz', index=0)
atoms1 = io.read('examples/Co138_1.xyz', index=0)


# In[ ]:


# Create an assignment object between atoms0 and atoms1,
# enable biasing and iteration and record up to 5 assignments
assignments = Assignment(atoms0, atoms1, biasing=True, iteration=True, records=5)


# In[ ]:


# Align atoms1 to atoms0 for each assignment and
# write the aligned coordinates to a file
for i, order in enumerate(assignments.mapind, start=1):
    alignable1 = Alignable(atoms1[order])
    alignable1.alignto(atoms0)
    print(i, alignable1.disto(atoms0))
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), alignable1, append=True)

