#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from ase import io
from molalign import Alignment


# In[ ]:


atoms0 = io.read('r005/Co100.xyz', index=0)
atoms1 = io.read('r005/Co100.xyz', index=1)
alignment0 = Alignment(atoms0, biased=True)


# In[ ]:


maplist, mapcount, mindist = alignment0.sort(atoms1)


# In[ ]:


for i, mapping in enumerate(maplist, start=1):
    io.write('aligned_{}.xyz'.format(i), atoms0)
    io.write('aligned_{}.xyz'.format(i), alignment0.aligned(atoms1, mapping), append=True)

