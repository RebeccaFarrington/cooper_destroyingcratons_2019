
# coding: utf-8

# Craton stability - 201710
# ======
# 
# This model solves the craton destruction cartoon plate/craton/subduction system geometry sketched out by Katie in January.

# **We need to consider carefully the viscosity/strength profile within the craton**
# 
# I have used a half-space cooling temperature profile, I don't think this is as linear within the craton as a continental geotherm would be. It does produce a strong craton. However, as we are testing the strength of this craton in the face of a fast moving slab/paddle, I'd like to confirm it gives the desired strength/viscosity profile. See 'Half-space cooling model' notebook for details on the temperature profile.

# **Need to fix**
# 1. We should model the whole mantle
# 2. What resultion do we want - highest possible.
# 3. need to include transisition regions in the material layout
#     b. as per KT's model
# 4. slab temperature profile is not correct

# In[1]:

import numpy as np
import underworld as uw
import math
from underworld import function as fn
#import glucifer
import os


# In[2]:

#try:
#    workdir
#except NameError:
#    workdir = os.path.abspath(".")

#outputPath = os.path.join(workdir,"cratonSlab-201710-output/Test01")
outputPath = './'

if uw.rank()==0:
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
uw.barrier()


# Set simulation parameters.

# In[3]:

# physical parameters
g        = 9.8      # [m/(s.s)],   gravity 
alpha    = 3*1e-5   # [K^-1],      thermal expansivity coefficient
kappa    = 1e-6     # [m.m/s],     thermal diffusivity
rho0     = 3300.    # [kg/m^3],    reference density
Temp_Min = 300.0    # [K],         surface temperature, 26.85C, 0C = 273.15K
Temp_Max = 1573.0   # [K],         mantle temperature
R        = 8.3145   # [J/(K.mol)], gas constant

deltaTemp = Temp_Max - Temp_Min


# In[4]:

ageD3 =  780.*1e6*(60.*60.*24*365) # in seconds, 100 Myr, thermal age of cratonic lithosphere, 
ageD2 =  280.*1e6*(60.*60.*24*365)
ageD1 =   30.*1e6*(60.*60.*24*365)
ageSL =   20.*1e6*(60.*60.*24*365) 
ageOL =   20.*1e6*(60.*60.*24*365)

transLength = 200.*1e3           # transition length, 200 km

#length in [km]
dim = 3
# minCoord = [-2.5,0.0,-1.5]  # min coordinates in (x,y)
# maxCoord = [ 2.5,1.0, 1.5]  # max coordinates in (x,y)

boxLength = 5.0  # x 1e3 m = 5,000 km
boxHeight = 0.66 # x 1e3 m =   660 km
boxWidth  = 3.0  # x 1e3 m = 3,000 km
minCoord = [-boxLength/2., -boxHeight, -boxWidth/2.]  # x centre at zero
maxCoord = [ boxLength/2., 0.,          boxWidth/2.]           # y = 0 at surface
# assume boxHeight = 1
aspectRatioX = boxLength/ 1.
aspectRatioZ = boxWidth / 1.

unitRes     = 128 
resolution  = [int(unitRes*aspectRatioX),unitRes,int(unitRes*aspectRatioZ)]


# In[5]:

# scale vectors for 1 cm/yr
# v = (kappa / depth) * v' [m/s]
# [cm/yr] = (1e2cm/m).(60*60*24*365s/yr)[m/s]
velScalingFactor = (kappa/(boxHeight*1e6))*(1e2*(60*60*24*365))
#print velScalingFactor


# Create mesh and finite element variables
# ------

# In[6]:
print '*** Making mesh ***'
mesh = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
                                 elementRes  = resolution, 
                                 minCoord    = minCoord, 
                                 maxCoord    = maxCoord,
                               ) 

velocityField       = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=dim )
pressureField       = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 )
temperatureField    = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )
temperatureDotField = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )
viscosityField      = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )
deformationField    = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )


# In[7]:

#print minCoord, maxCoord


# In[8]:
print '*** Initialise Fields ***'
# initialise fields
temperatureField.data[:]    = 0.
temperatureDotField.data[:] = 0.
velocityField.data[:]       = [0.,0.,0.]
pressureField.data[:]       = 0.
viscosityField.data[:]      = 0.
deformationField.data[:]    = 0.


# In[9]:

# Analysis fields
velocityGradientField = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=dim*dim )
strainRateField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=dim*dim-dim )
vorticityXField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
vorticityYField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
vorticityZField       = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
stressField           = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=dim*dim-dim  )
dissipationField      = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )

velocityGradientField.data[:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
strainRateField.data[:]       = [0.,0.,0.,0.,0.,0.]
vorticityXField.data[:]       = [0.]
vorticityYField.data[:]       = [0.]
vorticityZField.data[:]       = [0.]
stressField.data[:]           = [0.,0.,0.,0.,0.,0.]
dissipationField.data[:]      = [0.]


# In[10]:

# mesh.reset()


# In[11]:

# with mesh.deform_mesh():
#     mesh.data[:,0] = mesh.data[:,0] * np.exp(mesh.data[:,0]*mesh.data[:,0]) / np.exp(maxCoord[0]*maxCoord[0])
#     mesh.data[:,1] = mesh.data[:,1] * np.exp(mesh.data[:,1]*mesh.data[:,1]*2) / np.exp(boxHeight*boxHeight*2)


# In[12]:

#figsize=(1050.,350.)


# In[13]:

#figMesh = glucifer.Figure( figsize=figsize )
#figMesh.append( glucifer.objects.Mesh(mesh)) 
#figMesh.show()
#figMesh.open_viewer()


# **Initialise & apply boundary conditions on the temperature Field**

# In[14]:
print '*** initial temperature ***'
# box margin shape
craton_d3_Shape = np.array([ (2.6,-1.6), (2.6,-0.1), (-1.5,-0.1), (-1.5, -1.6) ])
craton_d2_Shape = np.array([ (2.6,-0.1), (2.6, 0.2), (-1.5, 0.2), (-1.5, -0.1) ])
craton_d1_Shape = np.array([ (2.6, 0.2), (2.6, 0.5), (-1.5, 0.5), (-1.5,  0.2) ])

slab_xMin       = -1.5
slab_xMax       = -1.3
slab_xMid       = slab_xMin+0.5*(slab_xMax-slab_xMin)
slab_yMin       = -0.495
slab_yMax       = 0.0
slab_zMin       = 0.5
slab_zMax       = 1.6

craton_d3 = fn.shape.Polygon( craton_d3_Shape )
craton_d2 = fn.shape.Polygon( craton_d2_Shape )
craton_d1 = fn.shape.Polygon( craton_d1_Shape )


# In[15]:

print 'element width, unit res 128 = ', 1./64*1e3, '[km]'
print 0.02*1e3
print 2./unitRes


# In[16]:

#dTrans = 1./unitRes+0.01
dTrans = 0.02 # 40 km transition length
# create transition regions between steps in z direction
craton_trans_d1dOL_Shape = np.array([(2.6, 0.5-dTrans),(2.6, 0.5+dTrans),(-1.5, 0.5+dTrans),(-1.5, 0.5-dTrans)])
craton_trans_d1d2_Shape  = np.array([(2.6, 0.2-dTrans),(2.6, 0.2+dTrans),(-1.5, 0.2+dTrans),(-1.5, 0.2-dTrans)])
craton_trans_d2d3_Shape  = np.array([(2.6,-0.1-dTrans),(2.6,-0.1+dTrans),(-1.5,-0.1+dTrans),(-1.5,-0.1-dTrans)])

craton_trans_d1dOL = fn.shape.Polygon( craton_trans_d1dOL_Shape )
craton_trans_d1d2  = fn.shape.Polygon( craton_trans_d1d2_Shape ) 
craton_trans_d2d3  = fn.shape.Polygon( craton_trans_d2d3_Shape ) 

# create transition regions in X direction
craton_d1dOL_Shape = np.array([ (-1.5-dTrans, 0.2), (-1.5-dTrans, 0.5), (-1.5+dTrans, 0.5), (-1.5+dTrans,  0.2) ])
craton_d2dOL_Shape = np.array([ (-1.5-dTrans,-0.1), (-1.5-dTrans, 0.2), (-1.5+dTrans, 0.2), (-1.5+dTrans, -0.1) ])
craton_d3dOL_Shape = np.array([ (-1.5-dTrans,-1.6), (-1.5-dTrans,-0.1), (-1.5+dTrans,-0.1), (-1.5+dTrans, -1.6) ])

craton_trans_d1dOL  = fn.shape.Polygon( craton_d1dOL_Shape ) 
craton_trans_d2dOL  = fn.shape.Polygon( craton_d2dOL_Shape ) 
craton_trans_d3dOL  = fn.shape.Polygon( craton_d3dOL_Shape ) 


# In[17]:

# calculate half-space cooling model
for index, coord in enumerate(mesh.data):
    depth  = -1.0*coord[1]*1e6  # in [m]
    xCoord = coord[0]*1e6
    slabXcoord = coord[1]*np.tan(-np.pi/6)
    # if within slab shape
    if ( (coord[0]>slab_xMin+slabXcoord) & (coord[0]<slab_xMax+slabXcoord) 
          & (coord[1]>slab_yMin) & (coord[2] > slab_zMin)):
        # 0. define midplane of slab
        midPlane = slab_xMid+coord[1]*np.tan(-np.pi/6)
        # 1. distant from mid plane
        distance = math.fabs(coord[0]-midPlane)
        distance_depth = distance*1e6  # distance to midplane in [m]
        # 2. temperature at distace from mid plane
        tempSlab = deltaTemp*math.erf(distance_depth/(2.*math.sqrt(kappa*ageSL)))/deltaTemp
        # take the greater temperature, slab or surface plate
        tempOL   = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*ageOL)))/deltaTemp# if within current surface plate
        temp = min(tempSlab, tempOL) 

    # if within craton
    elif craton_d3.evaluate((coord[0],coord[2])):
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*ageD3)))/deltaTemp
    elif craton_d2.evaluate((coord[0],coord[2])):
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*ageD2)))/deltaTemp
    elif craton_d1.evaluate((coord[0],coord[2])):
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*ageD1)))/deltaTemp
    else:
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*ageOL)))/deltaTemp

    # if within craton perimeter - Z direction
    if craton_trans_d1dOL.evaluate((coord[0],coord[2])):
        ZtransStart =  0.5-dTrans
        phi  = (coord[2]-ZtransStart)/(2.*dTrans)
        age = (1.-phi) * ageD1 + phi * ageOL
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp
    elif craton_trans_d1d2.evaluate((coord[0],coord[2])):
        ZtransStart =  0.2-dTrans
        phi  = (coord[2]-ZtransStart)/(2.*dTrans)
        age = (1.-phi) * ageD2 + phi * ageD1
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp
    elif craton_trans_d2d3.evaluate((coord[0],coord[2])):
        ZtransStart = -0.1-dTrans
        phi  = (coord[2]-ZtransStart)/(2.*dTrans)
        age = (1.-phi) * ageD3 + phi * ageD2
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp        

    # if within craton perimeter - X direction
    if craton_trans_d1dOL.evaluate((coord[0],coord[2])):
        XtransStart = -1.5-dTrans
        phi  = (coord[0]-XtransStart)/(2.*dTrans)
        age = (1.-phi) * ageOL + phi * ageD1
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp
    elif craton_trans_d2dOL.evaluate((coord[0],coord[2])):
        XtransStart = -1.5-dTrans
        phi  = (coord[0]-XtransStart)/(2.*dTrans)
        age = (1.-phi) * ageOL + phi * ageD2
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp
    elif craton_trans_d3dOL.evaluate((coord[0],coord[2])):
        XtransStart = -1.5-dTrans
        phi  = (coord[0]-XtransStart)/(2.*dTrans)
        age = (1.-phi) * ageOL + phi * ageD3
        temp = deltaTemp*math.erf(depth/(2.*math.sqrt(kappa*age)))/deltaTemp
        
    temperatureField.data[index] = temp
    


# **Plot initial temperature**

# In[18]:

#figtemp = glucifer.Figure(figsize=figsize, axis=True, )
#figtemp.append( glucifer.objects.Surface(mesh, temperatureField, onMesh=True) )
#figtemp.append( glucifer.objects.Contours(mesh, temperatureField, interval=(200/deltaTemp), 
#                                           limits=(0.,(1400-300)/deltaTemp), 
#                                           colours='Black', colourBar=False))
#figtemp.append( glucifer.objects.IsoSurface(mesh, temperatureField, isovalue=0.9, colourBar=False))
#figtemp.append( glucifer.objects.Mesh(mesh))
#figtemp.show()
#figtemp.open_viewer()


# In[19]:

# for index in mesh.specialSets["MinJ_VertexSet"]:
#     temperatureField.data[index] = 1.
# for index in mesh.specialSets["MaxJ_VertexSet"]:
#     temperatureField.data[index] = 0.


# In[20]:
print '*** set boundary conditions***'
# set boundary conditions of velocity and temperature field
iWalls  = mesh.specialSets["MinI_VertexSet"] + mesh.specialSets["MaxI_VertexSet"]
jWalls  = mesh.specialSets["MinJ_VertexSet"] + mesh.specialSets["MaxJ_VertexSet"]
kWalls  = mesh.specialSets["MinK_VertexSet"] + mesh.specialSets["MaxK_VertexSet"]
topWall = mesh.specialSets["MaxJ_VertexSet"]

freeslipBC = uw.conditions.DirichletCondition( variable      = velocityField, 
                                               indexSetsPerDof = ( iWalls, jWalls, kWalls) )
tempBC     = uw.conditions.DirichletCondition( variable      = temperatureField, 
                                               indexSetsPerDof = ( jWalls, ) )


# In[ ]:




# **Define the rheology**

# In[21]:
print '*** define rheology ***'
# strain rate invariant
strainRateFn              = fn.tensor.symmetric( velocityField.fn_gradient )
strainRate_2ndInvariantFn = fn.tensor.second_invariant(strainRateFn)
strainRate_2ndInvariantFn_scaled = strainRate_2ndInvariantFn * (kappa / ((boxHeight*1e6)**2))


# In[22]:

# hydrostatice pressure 
yCoord = fn.input()[1]*1e6
z_hat  = -1.0*(yCoord)
P_stat = rho0 * g * z_hat


# In[23]:

# adiabatic gradiet to temperature solution, 0.5 K/km
Tm = (z_hat/1e3)*0.5+(temperatureField*deltaTemp)+Temp_Min


# In[24]:

# limiters
eta_max = 1e24 # Pa.s, maximum viscosity
eta_min = 1e19 # Pa.s, minimum viscosity


# In[25]:

# Diffusion creep (dry olivine)
E_diff = 300.*1e3     # J/mol, activation energy
V_diff = 4.5/(1e2**3) # m^3/mol, activation volume
A_diff = 4.0e-10      # Pa^-1.s^-1, prefactor
n_diff = 1.


# In[26]:

diffusionCreep  = 0.5*fn.math.pow(A_diff,(-1./n_diff))
diffusionCreep *= fn.math.exp((E_diff+P_stat*V_diff)/(n_diff*R*Tm))
diffusionCreep *= fn.math.pow((strainRate_2ndInvariantFn_scaled+1.0e-24),(1-n_diff)/n_diff)


# In[27]:

# Dislocation creep (dry olivine)
E_disl = 540.*1e3     # J/mol, activation energy
V_disl = 10./(1e2**3) # m^3/mol, activation volume
A_disl = 1.0e-15      # Pa^-n.s^-1, prefactor
n_disl = 3.5


# In[28]:

dislocationCreep  = 0.5*fn.math.pow(A_disl,(-1./n_disl))
dislocationCreep *= fn.math.exp((E_disl+P_stat*V_disl)/(n_disl*R*Tm))
dislocationCreep *= fn.math.pow((strainRate_2ndInvariantFn_scaled+1.0e-24),(1-n_disl)/n_disl)


# In[29]:

creep_harmonic = fn.math.pow((1./diffusionCreep + 1./dislocationCreep)/2.,-1)


# In[30]:

fn_viscosity = fn.misc.max(fn.misc.min(creep_harmonic, eta_max), eta_min )/eta_min


# In[31]:

#viscosity = glucifer.objects.Surface(mesh, fn_viscosity*eta_min, logScale=True, valueRange=[1e18, 1e24], )
#viscosity.colourBar["tickvalues"] = [1e18, 1e19, 1e20, 1e21, 1e22, 1e23, 1e24]

#figEta = glucifer.Figure(figsize=figsize, title='Viscosity field, Temperature contours')
#figEta.append( viscosity)
#figEta.append( glucifer.objects.Contours(mesh, temperatureField, interval=(100/deltaTemp), 
#                                         limits=(0.,(1400-300)/deltaTemp),                                          
#                                         colours='Black', colourBar=False))

#figEta.show()
#figEta.open_viewer()


# **Define stress, density, gravity and buoyancy force**

# In[32]:
print '*** set analysis functions***'
# velocity gradients function
velocityGradientFn = velocityField.fn_gradient

# symmetric and antisymmetric components of velocity gradient
strainRateFn = fn.tensor.symmetric(velocityField.fn_gradient)

# w = \nambla v, vorticity = curl of velocity
vorticityXFn  = velocityGradientFn[7]-velocityGradientFn[5]
vorticityYFn  = velocityGradientFn[2]-velocityGradientFn[6]
vorticityZFn  = velocityGradientFn[3]-velocityGradientFn[1]

#stressFn = fn_viscosity * strainRateFn
stressFn = fn_viscosity * strainRateFn
dissipationFn = (  strainRateFn[0]*stressFn[0] + 2.*strainRateFn[3]*stressFn[3] + 2.*strainRateFn[4]*stressFn[4]
                                               +    strainRateFn[1]*stressFn[1] + 2.*strainRateFn[5]*stressFn[5] 
                                                                                +    strainRateFn[2]*stressFn[2]  
                )


# In[33]:

eta0 = eta_min
Ra   = (alpha*rho0*g*deltaTemp*((boxHeight*1e6)**3)) / (eta0*kappa)


# In[34]:

densityFn = Ra*temperatureField
gravity = ( 0.0, 1.0, 0.0 )
buoyancyFn = gravity*densityFn


# **Setup a Stokes system**

# In[35]:
print '*** set up systems - stokes ***'
stokes = uw.systems.Stokes( velocityField = velocityField, 
                            pressureField = pressureField,
                               conditions = [freeslipBC,],
                             fn_viscosity = fn_viscosity, 
                             fn_bodyforce = buoyancyFn,
                          )


# In[36]:
print '*** set up systems - solver ***'
solver = uw.systems.Solver(stokes)


# In[37]:

#if uw.nProcs()==0:
#    solver.set_inner_method("lu")


# **Checkpointing**

# In[38]:s
print '*** checkpointing 00000 ***'
#dissipationFn.evaluate(mesh)


# In[39]

# evaluate functions
viscosityField.data[:]  = fn_viscosity.evaluate( mesh )

#velocityGradientField.data[:] = velocityGradientFn.evaluate(mesh)
#strainRateField.data[:]  = strainRateFn.evaluate(mesh)
#vorticityXField.data[:]  = vorticityXFn.evaluate(mesh)
#vorticityYField.data[:]  = vorticityYFn.evaluate(mesh)
#vorticityZField.data[:]  = vorticityZFn.evaluate(mesh)
#stressField.data[:]      = strainRateFn.evaluate(mesh)
#dissipationField.data[:] = dissipationFn.evaluate(mesh)


# In[40]:

# save variables
meshHnd         = mesh.save(outputPath+'mesh.h5')
#velocityHnd     = velocityField.save(   outputPath+   'velocityField.00000.h5',meshHnd)
#pressureHnd     = pressureField.save(   outputPath+   'pressureField.00000.h5',meshHnd)
temperatureHnd  = temperatureField.save(outputPath+'temperatureField.00000.h5',meshHnd)
viscosityHnd    = viscosityField.save(  outputPath+  'viscosityField.00000.h5',meshHnd)

#velocityGradientHnd = velocityGradientField.save(outputPath+'velocityGradientField.00000.h5',meshHnd)
#strainRateHnd       = strainRateField.save(      outputPath+      'strainRateField.00000.h5',meshHnd)
#vorticityXHnd       = vorticityXField.save(      outputPath+      'vorticityXField.00000.h5',meshHnd)
#vorticityYHnd       = vorticityYField.save(      outputPath+      'vorticityYField.00000.h5',meshHnd)
#vorticityZHnd       = vorticityZField.save(      outputPath+      'vorticityZField.00000.h5',meshHnd)
#stressHnd           = stressField.save(          outputPath+          'stressField.00000.h5',meshHnd)
#dissipationHnd      = dissipationField.save(     outputPath+     'dissipationField.00000.h5',meshHnd)


# In[ ]:

# # and the xdmf files
#velocityField.xdmf(   outputPath+'velocityField.00000.xdmf',   velocityHnd,   "velocity",   meshHnd,"mesh",modeltime=0.)
#pressureField.xdmf(   outputPath+'pressureField.00000.xdmf',   pressureHnd,   "pressure",   meshHnd,"mesh",modeltime=0.)
temperatureField.xdmf(outputPath+'temperatureField.00000.xdmf',temperatureHnd,"temperature",meshHnd,"mesh",modeltime=0.)
viscosityField.xdmf(  outputPath+'viscosityField.00000.xdmf',  viscosityHnd,  "viscosity",  meshHnd,"mesh",modeltime=0.)

#velocityGradientField.xdmf(outputPath+'velocityGradientField.00000.xdmf',velocityGradientHnd,"velocityGradient",meshHnd,"mesh",modeltime=0.)
#strainRateField.xdmf( outputPath+ 'strainRateField.00000.xdmf', strainRateHnd,"strainRate", meshHnd,"mesh",modeltime=0.)
#vorticityXField.xdmf( outputPath+ 'vorticityXField.00000.xdmf', vorticityXHnd,"vorticityX", meshHnd,"mesh",modeltime=0.)
#vorticityYField.xdmf( outputPath+ 'vorticityYField.00000.xdmf', vorticityYHnd,"vorticityY", meshHnd,"mesh",modeltime=0.)
#vorticityZField.xdmf( outputPath+ 'vorticityZField.00000.xdmf', vorticityZHnd,"vorticityZ", meshHnd,"mesh",modeltime=0.)
#stressField.xdmf(     outputPath+     'stressField.00000.xdmf',     stressHnd,"stress",     meshHnd,"mesh",modeltime=0.)
#dissipationField.xdmf(outputPath+'dissipationField.00000.xdmf',dissipationHnd,"dissipation",meshHnd,"mesh",modeltime=0.)


# **Set up and solve the Stokes system**

# In[ ]:
print '*** SOLVE ***'
solver.solve(nonLinearIterate=True)


# In[ ]:

#figEta.append( glucifer.objects.VectorArrows(mesh, velocityField*1e-4 ) )
#figEta.show()


# In[ ]:
print '*** deformation mechanism calculation***'
minViscosityIndex     = 0
maxViscosityIndex     = 1
diffusionCreepIndex   = 2
dislocationCreepIndex = 3

conditions = [ 
               (fn_viscosity*eta_min <= eta_min,          minViscosityIndex),
               (fn_viscosity*eta_min >= eta_max,          maxViscosityIndex),
               (diffusionCreep       <  dislocationCreep, diffusionCreepIndex), 
               (dislocationCreep     <= diffusionCreep,   dislocationCreepIndex),
              ]


# In[ ]:

# plot deformation mechanism using the swarm 
# deformationSwarmVar.data[:]  = uw.function.branching.conditional( conditions ).evaluate(swarm)
# plot deformation mechanism using the mesh
deformationField.data[:]  = uw.function.branching.conditional( conditions ).evaluate(mesh)


# In[ ]:

#figMech = glucifer.Figure(figsize=figsize, title='Deformation Mechanism')
#figMech.append( glucifer.objects.Contours(mesh, temperatureField, interval=(100/deltaTemp), 
#                                          limits=(0.,(1400-300)/deltaTemp), 
#                                          colours='Black',colourBar=False ))

# plot deformation mechanism using the mesh
#figMech.append( glucifer.objects.Surface(mesh, deformationField, colours = 'Red Blue Grey', discrete=True))


# In[ ]:

#figMech.show()
#figMech.open_viewer()


# In[ ]:
print '*** calculate analysis functions***'
# evaluate functions
viscosityField.data[:]  = fn_viscosity.evaluate( mesh )

velocityGradientField.data[:] = velocityGradientFn.evaluate(mesh)
strainRateField.data[:]  = strainRateFn.evaluate(mesh)
vorticityXField.data[:]  = vorticityXFn.evaluate(mesh)
vorticityYField.data[:]  = vorticityYFn.evaluate(mesh)
vorticityZField.data[:]  = vorticityZFn.evaluate(mesh)
stressField.data[:]      = strainRateFn.evaluate(mesh)
dissipationField.data[:] = dissipationFn.evaluate(mesh)


# In[ ]:
print '*** Save fields***'
# save variables
velocityHnd     = velocityField.save(   outputPath+   'velocityField.00001.h5',meshHnd)
pressureHnd     = pressureField.save(   outputPath+   'pressureField.00001.h5',meshHnd)
temperatureHnd  = temperatureField.save(outputPath+'temperatureField.00001.h5',meshHnd)
viscosityHnd    = viscosityField.save(  outputPath+  'viscosityField.00001.h5',meshHnd)

velocityGradientHnd = velocityGradientField.save(outputPath+'velocityGradientField.00001.h5',meshHnd)
strainRateHnd       = strainRateField.save(      outputPath+      'strainRateField.00001.h5',meshHnd)
vorticityXHnd       = vorticityXField.save(      outputPath+      'vorticityXField.00001.h5',meshHnd)
vorticityYHnd       = vorticityYField.save(      outputPath+      'vorticityYField.00001.h5',meshHnd)
vorticityZHnd       = vorticityZField.save(      outputPath+      'vorticityZField.00001.h5',meshHnd)
stressHnd           = stressField.save(          outputPath+          'stressField.00001.h5',meshHnd)
dissipationHnd      = dissipationField.save(     outputPath+     'dissipationField.00001.h5',meshHnd)


# In[ ]:

# and the xdmf files
velocityField.xdmf(   outputPath+   'velocityField.00001.xdmf',   velocityHnd,   "velocity",meshHnd,"mesh",modeltime=1.)
pressureField.xdmf(   outputPath+   'pressureField.00001.xdmf',   pressureHnd,   "pressure",meshHnd,"mesh",modeltime=1.)
temperatureField.xdmf(outputPath+'temperatureField.00001.xdmf',temperatureHnd,"temperature",meshHnd,"mesh",modeltime=1.)
viscosityField.xdmf(  outputPath+  'viscosityField.00001.xdmf',  viscosityHnd,  "viscosity",meshHnd,"mesh",modeltime=1.)

velocityGradientField.xdmf(outputPath+'velocityGradientField.00001.xdmf',velocityGradientHnd,"velocityGradient",meshHnd,"mesh",modeltime=1.)
strainRateField.xdmf( outputPath+ 'strainRateField.00001.xdmf', strainRateHnd,"strainRate", meshHnd,"mesh",modeltime=1.)
vorticityXField.xdmf( outputPath+ 'vorticityXField.00001.xdmf', vorticityXHnd,"vorticityX", meshHnd,"mesh",modeltime=1.)
vorticityYField.xdmf( outputPath+ 'vorticityYField.00001.xdmf', vorticityYHnd,"vorticityY", meshHnd,"mesh",modeltime=1.)
vorticityZField.xdmf( outputPath+ 'vorticityZField.00001.xdmf', vorticityZHnd,"vorticityZ", meshHnd,"mesh",modeltime=1.)
stressField.xdmf(     outputPath+     'stressField.00001.xdmf',     stressHnd,"stress",     meshHnd,"mesh",modeltime=1.)
dissipationField.xdmf(outputPath+'dissipationField.00001.xdmf',dissipationHnd,"dissipation",meshHnd,"mesh",modeltime=1.)


# In[ ]:

print '*** Finished ***'

