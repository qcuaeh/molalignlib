#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from ase import io
from molalign_wrapper import Aligner


# In[ ]:


atoms0 = io.read('r005/Co100.xyz', index=0)
atoms1 = io.read('r005/Co100.xyz', index=1)
aligner = Aligner(atoms0, records=10, count=10, bias=0.17, iteration=True, weighted=True, testing=True)


# In[ ]:


maplist, mapcount = aligner.remap(atoms1)


# In[ ]:


for i, mapping in enumerate(maplist, start=1):
    atoms1_aligned = aligner.align(atoms1, mapping)
    io.write('aligned_{}.xyz'.format(i), atoms0,)
    io.write('aligned_{}.xyz'.format(i), atoms1_aligned, append = True)

