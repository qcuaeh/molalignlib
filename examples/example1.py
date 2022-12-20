#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import modules
from ase.io import read, write
from molalignlib import Assignment, Alignment


# In[ ]:


# Read clusters coordinates
atoms0 = read('Co138_0.xyz', index=0)
atoms1 = read('Co138_1.xyz', index=0)


# In[ ]:


# Create an assignment object between atoms0 and atoms1,
# enable biasing and iteration and record up to 5 assignments
assignments = Assignment(atoms0, atoms1, biasing=True, iteration=True, stats=True, rec=5)


# In[ ]:


# Align atoms1 to atoms0 for each assignment and
# write the aligned coordinates to a file
for i, mapping in enumerate(assignments, start=1):
    alignment = Alignment(atoms0, atoms1[mapping])
    write('aligned_{}.xyz'.format(i), atoms0)
    write('aligned_{}.xyz'.format(i), alignment.align(atoms1[mapping]), append=True)

