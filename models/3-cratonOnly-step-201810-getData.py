
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

#import matplotlib.pyplot as pyplot
#get_ipython().magic(u'matplotlib inline')


# In[2]:


inputPath = os.path.join(
			"/short/m18/rjf565/CratonStability-201810/3-cratonOnly-step-201810/"
			)


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


# define field names of interest
# define lines of interest
# 
# for each name
#     for each line
#         evaluate the field along the line
#         save swarm data to file

# In[5]:


# create fields
dissipationField      = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
pressureField         = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 ) 
viscosityField        = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
vorticityXField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
vorticityYField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
vorticityZField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
temperatureField      = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
velocityField         = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=3 ) 
strainRateField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=6 )
stressField           = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=6 )
velocityGradientField = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=9 ) 


# In[6]:


# create field dictionary
fieldDict = {
            'dissipationField'     :dissipationField, 
            'pressureField'        :pressureField,
            'strainRateField'      :strainRateField,
            'stressField'          :stressField,
            'temperatureField'     :temperatureField,
            'velocityField'        :velocityField,
            'velocityGradientField':velocityGradientField,
            'viscosityField'       :viscosityField,
            'vorticityXField'      :vorticityXField,
            'vorticityYField'      :vorticityYField,
            'vorticityZField'      :vorticityZField,
            }


# In[7]:


# load field variable data
for field in sorted(fieldDict):
    print field, fieldDict[field]
    path = inputPath+str(field)+'.00001.h5'
    fieldDict[field].load(path)


# In[8]:

print '*** Get transect data! ***'


fn_stressField = 2. * viscosityField * strainRateField
stressField.data[:] = fn_stressField.evaluate(mesh)


# In[9]:


# create tracer lines
tracerSwarm01 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm02 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm03 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm04 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm05 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm06 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm07 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm08 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm09 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm10 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm11 = uw.swarm.Swarm( mesh=mesh )
tracerSwarm12 = uw.swarm.Swarm( mesh=mesh )


# In[10]:


# create swarm dictionary
swarmDict = {
           'y130+x660':tracerSwarm01, 
            'y130+x825':tracerSwarm02, 
            'y130+x990':tracerSwarm03, 
            'y200+x660':tracerSwarm04, 
            'y200+x825':tracerSwarm05, 
            'y200+x990':tracerSwarm06, 
            'y130-z165':tracerSwarm07, 
            'y130+z000':tracerSwarm08, 
            'y130+z165':tracerSwarm09, 
            'y200-z165':tracerSwarm10, 
            'y200+z000':tracerSwarm11, 
            'y200+z165':tracerSwarm12, 
            }


# In[11]:


tracerRes = int(unitRes*boxLength*10+1)
# y = -130km, -0.197
# x = -660km, 1.0
particleCoord = np.zeros((tracerRes,3))
particleCoord[:,0] = -1.0
particleCoord[:,1] = -0.197
particleCoord[:,2] = np.linspace(minCoord[2], maxCoord[2], tracerRes)
temp = tracerSwarm01.add_particles_with_coordinates(particleCoord)

# x = -825km, 1.25
particleCoord[:,0] = -1.25
temp = tracerSwarm02.add_particles_with_coordinates(particleCoord)

# x = -990km, 1.5
particleCoord[:,0] = -1.50
temp = tracerSwarm03.add_particles_with_coordinates(particleCoord)

# y = -200km, -0.303
# x = -660km, 1.0
particleCoord[:,0] = -1.0
particleCoord[:,1] = -0.303
temp = tracerSwarm04.add_particles_with_coordinates(particleCoord)

# x = -825km, 1.0
particleCoord[:,0] = -1.25
temp = tracerSwarm05.add_particles_with_coordinates(particleCoord)

# x = -990km, 1.0
particleCoord[:,0] = -1.5
temp = tracerSwarm06.add_particles_with_coordinates(particleCoord)

# y = -130km, -0.197
# z = -165km, -0.25
particleCoord[:,0] = np.linspace(minCoord[0], maxCoord[0], tracerRes)
particleCoord[:,1] = -0.197
particleCoord[:,2] = -0.25+0.125
temp = tracerSwarm07.add_particles_with_coordinates(particleCoord)

# z = -km, 0.0
particleCoord[:,2] = 0.0+0.125
temp = tracerSwarm08.add_particles_with_coordinates(particleCoord)

# z = 165km, 0.25
particleCoord[:,2] = 0.25+0.125
temp = tracerSwarm09.add_particles_with_coordinates(particleCoord)

# y = -200km, -0.197
# z = -165km, -0.25
particleCoord[:,1] = -0.303
particleCoord[:,2] = -0.25+0.125
temp = tracerSwarm10.add_particles_with_coordinates(particleCoord)

# z = -km, 0.0
particleCoord[:,2] = 0.0+0.125
temp = tracerSwarm11.add_particles_with_coordinates(particleCoord)

# z = 165km, 0.25
particleCoord[:,2] = 0.25+0.125
temp = tracerSwarm12.add_particles_with_coordinates(particleCoord)


# In[12]:


# # have a look
# figsize=(1050.,350.)
# figtemp = glucifer.Figure(figsize=figsize, axis=True)
# figtemp.append( glucifer.objects.Surface(mesh, temperatureField) )
# figtemp.append( glucifer.objects.IsoSurface(mesh, viscosityField, isovalue=1e4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm01, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm02, temperatureField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm03, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm04, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm05, temperatureField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm06, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm07, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm08, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm09, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm10, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm11, viscosityField, fn_size=4, colourBar=False))
# figtemp.append( glucifer.objects.Points(tracerSwarm12, viscosityField, fn_size=4, colourBar=False))
# figtemp.show()


# In[13]:


# figtemp.open_viewer()


# In[14]:


for swarm in sorted(swarmDict):
    swarmObject = swarmDict[swarm]
    swarmObject.myVariables = {}
    swarmObject.save( inputPath+'/'+swarm+'.h5')
    print swarm, swarmObject
    for field in sorted(fieldDict):
        fieldObject = fieldDict[field]
        var = swarmObject.add_variable('double',fieldObject.nodeDofCount)
        swarmObject.myVariables[field] = var
        swarmObject.myVariables[field].data[:] = fieldObject.evaluate(swarmObject)
        swarmObject.myVariables[field].save( inputPath+'/'+swarm+'_'+field+'.h5')

print '*** get dissipation integration ***'

dissipationCratonField = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 
dissipationSlabField   = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 ) 


# In[7]:

dissipationCratonField.data[:] = 0.
dissipationSlabField.data[:]   = 0.

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
fn_conditional_slab_vol   = fn.branching.conditional(  ( (fn.coord()[2] >= 0.5, 1.),
                                                         (                True, 0.) ) )
fn_conditional_craton_vol = fn.branching.conditional(  ( (fn.coord()[2] <  0.5, 1.),
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
fn_conditional_slab_diss   = fn.branching.conditional(  ( (fn.coord()[2] >= 0.5, dissipationField),
                                                (                          True, 0.) ) )
fn_conditional_craton_diss = fn.branching.conditional(  ( (fn.coord()[2] <  0.5, dissipationField),
                                                          (                True, 0.) ) )


# In[52]:

slab_diss   = mesh.integrate(fn_conditional_slab_diss)[0]


# In[53]:

craton_diss = mesh.integrate(fn_conditional_craton_diss)[0]


# In[60]:

print 'Dissipation, slab region = {0:.3e}; craton region = {1:.3e}; total = {2:.3e}'.format(slab_diss, craton_diss, slab_diss+craton_diss)




print '*** Finished ***'


