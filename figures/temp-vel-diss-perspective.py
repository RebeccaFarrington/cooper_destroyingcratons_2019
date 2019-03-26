# I need to put pvpython in my $PYTHONPATH, but what path I do not know...

from paraview.simple import *
import os

inputPath = os.path.join("/Users/rfarrington/OneDrive - The University of Melbourne/Research/08_CratonStability/")
outputPathName = 'paraview-output-20181026'
# 316.888 v' = 1 cm/yr. Use this for colour bar dimless scale
velocity_range = [-317., 317.]
viewSize = [4000, 4000]

# for directories in dictionary
slabOnly        = []
cratonOnlyStep  = []
cratonSlabStep  = []
cratonOnlyDiag  = []
cratonSlabDiag  = []
cratonOnlyStra  = []
cratonSlabStra  = []
cratonOnlyStra2 = []
cratonSlabDtra2 = []

# make figures
modelRunDict = {
#                    "1 SlabOnly"             : slabOnly,
#                    "2 CratonOnly-step"      : cratonOnlyStep,
#                    "3 CratonSlab-step"      : cratonSlabStep,
                    "4 CratonOnly-diagonal"  : cratonOnlyDiag,
                    "5 CratonSlab-diagonal"  : cratonSlabDiag,
#                    "6 CratonOnly-straight"  : cratonOnlyStra,
#                    "7 CratonSlab-straight"  : cratonSlabStra,
#                    "8 CratonOnly-straight2" : cratonOnlyStra2,
#                    "9 CratonSlab-straight2" : cratonSlabDtra2,
                }

def perspective_temp_vel_diss(path, outputPath, outputPath2, model):
    # load temperature file
    temperatureField00001xdmf = XDMFReader(FileNames=[str(path)+'/temperatureField.00001.xdmf'])
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

    # load dissipation file
    dissipationField00001xdmf = XDMFReader(FileNames=[str(path)+'/dissipationField.00001.xdmf'])
    dissipationField00001xdmf.PointArrayStatus = ['dissipation']
    dissipationField00001xdmf.GridStatus = ['FEM_Mesh_mesh']

    renderView1 = GetActiveViewOrCreate('RenderView')
    dissipationField00001xdmfDisplay = Show(dissipationField00001xdmf, renderView1)
    dissipationField00001xdmfDisplay.Representation = 'Outline'

    # create a contour
    contour2 = Contour(Input=dissipationField00001xdmf)
    contour2.Isosurfaces = [2.0e8]

    # dispolay contour
    contour2Display = Show(contour2, renderView1)
    contour2Display.DiffuseColor = [0.6667, 0.0, 0.0]
    contour2Display.Opacity = 0.5

    # creat velocity glyph
    velocityField00001xdmf = XDMFReader(FileNames=[str(path)+'/velocityField.00001.xdmf'])
    velocityField00001xdmf.GridStatus = ['FEM_Mesh_mesh']
    glyph1 = Glyph(Input=velocityField00001xdmf, GlyphType='Arrow')
    glyph1.ScaleMode = 'vector'
    if model[2:12] == 'CratonOnly' :
        print 'CratonOnly vector factor x4'
        glyph1.ScaleFactor = 0.00025*4.
    else:
        glyph1.ScaleFactor = 0.00025
#        glyph1.ScaleFactor = 0.0000025
    glyph1.Vectors = ['POINTS', 'velocity']
    glyph1.MaximumNumberOfSamplePoints = 2000
#    glyph1.MaximumNumberOfSamplePoints = 10000

    glyph1Display = Show(glyph1, renderView1)
    ColorBy(glyph1Display, ('POINTS', 'GlyphVector', 'Y'))
    glyph1Display.SetScalarBarVisibility(renderView1, False)

    glyphVectorLUT = GetColorTransferFunction('GlyphVector')#'GlyphVector', 'velocity', "RTData"
    glyphVectorLUT.ApplyPreset('Rainbow Desaturated', True)
    glyphVectorLUT.RescaleTransferFunction(velocity_range[0], velocity_range[1])
    glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')
    glyphVectorPWF.RescaleTransferFunction(velocity_range[0], velocity_range[1])

    # set presepctive camera (bird's eye)
    camera=GetActiveCamera()
    camera.SetFocalPoint(0,0,0)
    camera.SetPosition(0,0,-10)
    camera.SetViewUp(0,1,0)

    camera.SetFocalPoint(0,-0.66,-0.5)
    camera.SetPosition(-6,2,3)

    renderView1.ViewSize = viewSize

    # save to model specific directory
    SaveScreenshot(outputPath+'perspective_temp_vel_diss.png')
    # save to collective viz directory
    SaveScreenshot(outputPath2+'/'+str(model)+'_perspective_temp_vel_diss.png')

    # set zoom camera (devil's eye)
    camera=GetActiveCamera()
    camera.SetFocalPoint(0,0,0)
    camera.SetPosition(0,0,-10)
    camera.SetViewUp(0,1,0)

    camera.Azimuth(180)
    camera.Elevation(-90)
    camera.Dolly(1.5)

    # save to model specific directory
    SaveScreenshot(outputPath+'perspective_temp_vel_diss_zoom.png')
    # save to collective viz directory
    SaveScreenshot(outputPath2+'/'+str(model)+'_perspective_temp_vel_diss_zoom.png')

    Delete(temperatureField00001xdmf)
    Delete(temperatureField00001xdmfDisplay)
    Delete(contour1)
    Delete(contour1Display)
    Delete(dissipationField00001xdmf)
#    Delete(dissipationField00001xdmfDisplay)
#    Delete(threshold1Display)
#    Delete(dissipationLUT)
    Delete(contour2)
    Delete(contour2Display)
    Delete(velocityField00001xdmf)
    Delete(glyph1)
    Delete(glyph1Display)
    Delete(glyphVectorLUT)
    Delete(glyphVectorPWF)
	
for model in sorted(modelRunDict):
    print model
    path = inputPath+model+'/'
    outputPath  = str(path)+str(outputPathName)+'/'
    outputPath2 = str(inputPath)+'/'+str(outputPathName)+'/'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    if not os.path.exists(outputPath2):
        os.makedirs(outputPath2)
    perspective_temp_vel_diss(path, outputPath, outputPath2, model)


