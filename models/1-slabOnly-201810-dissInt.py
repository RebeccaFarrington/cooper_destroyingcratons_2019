
# coding: utf-8

# **Integrates the dissipation field**
# 
# Using z=0.4 to distinguish between the slab and craton regions

# In[1]:

import underworld as uw
import math
from underworld import function as fn
import glucifer
import os
import numpy as np
import h5py
from shutil import copyfile


# In[2]:

inputPath = os.path.join(
    "/short/m18/rjf565/CratonStability-201711/cratonSlab-dryOlivene-3D-201710-proc128-res128-moreMem-noCraton/")


# In[3]:

#length in [km]
dim = 3

boxLength = 5.0  # x 1e6 m = 5,000 km
boxHeight = 0.66 # x 1e6 m =   660 km
boxWidth  = 3.0  # x 1e6 m = 3,000 km
minCoord = [-boxLength/2., -boxHeight, -boxWidth/2.]  # x centre at zero
maxCoord = [ boxLength/2., 0.,          boxWidth/2.]  # y = 0 at surface
# assume boxHeight = 1
aspectRatioX = boxLength/ 1.
aspectRatioZ = boxWidth / 1.

unitRes     = 128 
resolution  = [int(unitRes*aspectRatioX),unitRes,int(unitRes*aspectRatioZ)]


# In[4]:

mesh = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
                                 elementRes  = resolution, 
                                 minCoord    = minCoord, 
                                 maxCoord    = maxCoord,
                               ) 


# In[5]:

dissipationField          = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 ) 


# In[6]:

dissipationCratonField = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
dissipationSlabField   = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 


# In[7]:

dissipationCratonField.data[:] = 0.
dissipationSlabField.data[:]   = 0.


# In[8]:

dissipationField.load(inputPath+'dissipationField.00001.h5')


# In[9]:

dissipationIntegral = mesh.integrate(dissipationField)


# In[10]:

diss_max = dissipationField.data.max()


# In[11]:

diss_min = dissipationField.data.min()


# In[12]:

# dissipation RMS
diss_rms = dissipationIntegral[0] / ( boxLength * boxHeight * boxWidth)


# In[55]:

print 'Dissipation, Total = {0:.3e}; min = {1:.3e}; max = {2:.3e}; v_rms = {3:.3e}'.format(dissipationIntegral[0], diss_min, diss_max, diss_rms)


# In[ ]:

# define region with slab, z > 0.4
# define the volume
fn_conditional_slab_vol   = fn.branching.conditional(  ( (fn.coord()[2] >= 0.4, 1.),
                                                         (                True, 0.) ) )
fn_conditional_craton_vol = fn.branching.conditional(  ( (fn.coord()[2] <  0.4, 1.),
                                                (                     True, 0.) ) )


# In[49]:

slab_vol   = mesh.integrate(fn_conditional_slab_vol)[0]


# In[50]:

craton_vol = mesh.integrate(fn_conditional_craton_vol)[0]


# In[51]:

print slab_vol, craton_vol, slab_vol+craton_vol


# In[61]:

print 'Volume, slab region = {0:.3e}; craton region = {1:.3e}; total = {2:.3e}'.format(slab_vol, craton_vol, slab_vol+craton_vol)


# In[47]:

# integrate to get volume
fn_conditional_slab_diss   = fn.branching.conditional(  ( (fn.coord()[2] >= 0.4, dissipationField),
                                                (                          True, 0.) ) )
fn_conditional_craton_diss = fn.branching.conditional(  ( (fn.coord()[2] <  0.4, dissipationField),
                                                          (                True, 0.) ) )


# In[52]:

slab_diss   = mesh.integrate(fn_conditional_slab_diss)[0]


# In[53]:

craton_diss = mesh.integrate(fn_conditional_craton_diss)[0]


# In[60]:

print 'Dissipation, slab region = {0:.3e}; craton region = {1:.3e}; total = {2:.3e}'.format(slab_diss, craton_diss, slab_diss+craton_diss)


# In[72]:

# figsize=(1050.,350.)
# figDiss = glucifer.Figure(figsize=figsize, axis=True, )
# figDiss.append( glucifer.objects.Surface(mesh, dissipationField, onMesh=True) )
# figDiss.append( glucifer.objects.IsoSurface(mesh, dissipationField, isovalue=1e8, colourBar=False))
# # figDiss.append( glucifer.objects.Contours(mesh, dissipationField, interval=5e9, 
# #                                            colours='Black', colourBar=False) )
# # figDiss.append( glucifer.objects.Mesh(mesh))
# figDiss.show()


# In[73]:

#figDiss.open_viewer()


# **Quantifying dissipation**
# 
# Craton and slab regions defined by z = 0.4 plane
# 
# | Model | Description             | Total | Min   | Max   | Craton | Slab |
# | :---: |:-----------             | -----:| -----:| -----:| -----: | ----:| 
# | 1     | Step, slab only         |       |       |       |        |      |
# | 2     | Step, craton only       | | | | | | 
# | 3     | Step, craton & slab     | | | | | | 
# | 4     | Straight, craton only   | | | | | |  
# | 5     | Straight, craton & slab | | | | | |  
# | 6     | Diagonal, craton only   | | | | | |  
# | 7     | Diagonal, craton & slab | | | | | |  |
# 

# In[ ]:




# In[ ]:



