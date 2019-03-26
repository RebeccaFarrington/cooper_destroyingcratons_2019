# I need to put pvpython in my $PYTHONPATH, but what path I do not know...

from paraview.simple import *
import os

inputPath = os.path.join("/Users/rfarrington/OneDrive - The University of Melbourne/Research/08_CratonStability/")
outputPathName = 'paraview-output-20181106'

# make figures
model = "3 CratonSlab-step"
path = inputPath+model+'/'
outputPath  = str(inputPath)+str(outputPathName)+'/'
print outputPath

# 316.888 v' = 1 cm/yr. Use this for colour bar dimless scale
velocity_range = [-317., 317.]
viewSize = [4000, 4000]


temperatureField00001xdmf = XDMFReader(FileNames=[str(path)+'/temperatureField.00000.xdmf'])
temperatureField00001xdmf.GridStatus = ['FEM_Mesh_mesh']

# create renderview to modify
renderView1 = GetActiveViewOrCreate('RenderView')
temperatureField00001xdmfDisplay = Show(temperatureField00001xdmf, renderView1)
temperatureField00001xdmfDisplay.Representation = 'Outline'

# create a contour
contour1 = Contour(Input=temperatureField00001xdmf)
contour1.Isosurfaces = [0.864]

# dispolay contour
contour1Display = Show(contour1, renderView1)
contour1Display.DiffuseColor = [0.9, 0.9, 0.9]
contour1Display.Opacity = 0.5

# slice 1
slice1 = Slice(Input=temperatureField00001xdmf)
slice1.SliceType = 'Plane'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [2.5, -0.33, 0.0]
slice1.SliceType.Normal = [1.0, 0.0, 0.0]

slice1Display = Show(slice1, renderView1)
slice1Display.Representation = 'Surface'

# slice 2
slice2 = Slice(Input=temperatureField00001xdmf)
slice2.SliceType = 'Plane'

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.0, -0.33, -1.49]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

slice2Display = Show(slice2, renderView1)
slice2Display.Representation = 'Surface'
ColorBy(slice2Display, ('POINTS', 'temperature'))

# set presepctive camera (bird's eye)
camera=GetActiveCamera()
camera.SetFocalPoint(0,0,0)
camera.SetPosition(0,0,-10)
camera.SetViewUp(0,1,0)

camera.SetFocalPoint(0,-0.66,-0.5)
camera.SetPosition(-6,2,3)

renderView1.ViewSize = viewSize
#renderView1.AxesGrid.Visibility = 1

# save to model specific directory
#SaveScreenshot(outputPath+'perspective_tempIC.png')


