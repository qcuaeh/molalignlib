#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import modules
import sys
sys.path.append('/home/jmvasquez/repos/molalign')
from ase import io
from molalignlib import Assignment, Alignment


# In[ ]:


# Read clusters coordinates
atoms0 = io.read('Co138_0.xyz', index=0)
atoms1 = io.read('Co138_1.xyz', index=0)


# In[ ]:


# Create an assignment object between atoms0 and atoms1,
# enable biasing and iteration and record up to 5 assignments
assignments = Assignment(atoms0, atoms1, biasing=True, iteration=True, records=5)


# In[ ]:


# Align atoms1 to atoms0 for each assignment and
# write the aligned coordinates to a file
for i, mapping in enumerate(assignments, start=1):
    alignment = Alignment(atoms0, atoms1[mapping])
    print(i, alignment.rmsd)
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), alignment.align(atoms1[mapping]), append=True)

