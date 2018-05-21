#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

cur_path = '/Users/rahul/Desktop/Project/LSMLIB/lsmpqs-1.0/examples/finney_pack/212_3/imbibe/'
mask_name = cur_path + 'mask.raw'
# create a new 'Image Reader'
maskraw = ImageReader(FilePrefix=mask_name)

# Properties modified on maskraw
maskraw.DataScalarType = 'float'
maskraw.DataByteOrder = 'LittleEndian'
maskraw.DataExtent = [0, 217, 0, 217, 0, 217]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [920, 470]

# show data in view
maskrawDisplay = Show(maskraw, renderView1)
# trace defaults for the display properties.
maskrawDisplay.Representation = 'Outline'
maskrawDisplay.ColorArrayName = ['POINTS', '']
maskrawDisplay.ScalarOpacityUnitDistance = 1.7320508075688776
maskrawDisplay.Slice = 108

renderView1.ResetCamera()

# current camera placement for renderView1
#renderView1.CameraPosition = [508.42138099408646, -713.8426071332491, 557.8301576975391]
#renderView1.CameraFocalPoint = [152.5, 152.5, 152.5]
#renderView1.CameraViewUp = [-0.2582205645805154, 0.32033960932990924, 0.91143001636005]
#renderView1.CameraParallelScale = 264.1377481542538

# Camera placement for Finney packs
renderView1.CameraPosition = [592.6368477026884, -359.2562773141341, 380.5867493329233]
renderView1.CameraFocalPoint = [108.49999999999972, 108.50000000000007, 108.49999999999991]
renderView1.CameraViewUp = [-0.31681453488445915, 0.21010508804263228, 0.9249239982098755]
renderView1.CameraParallelScale = 187.92751262122317

# create a new 'Contour'
contour_mask = Contour(Input=maskraw)
contour_mask.ContourBy = ['POINTS', 'ImageFile']
contour_mask.Isosurfaces = [-0.23000073432922363]
contour_mask.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour_mask.Isosurfaces = [0.0]

# show data in view
contour_maskDisplay = Show(contour_mask, renderView1)
# trace defaults for the display properties.
contour_maskDisplay.ColorArrayName = [None, '']

# hide data in view
Hide(maskraw, renderView1)

# Properties modified on contour1Display
contour_maskDisplay.Opacity = 0.2
# flag for resetting data to current data
flag = 0
for i in range(1,46):
    filename = cur_path + 'data_step_' + str(i) + '_nw.raw'
    # create a new 'Image Reader'
    cur_image = ImageReader(FilePrefix=filename)

    # Properties modified on data_step_28_nwraw
    cur_image.DataScalarType = 'float'
    cur_image.DataByteOrder = 'LittleEndian'
    cur_image.DataExtent = [0, 217, 0, 217, 0, 217]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [920, 470]

    # show data in view
    cur_imageDisplay = Show(cur_image, renderView1)
    # trace defaults for the display properties.
    cur_imageDisplay.Representation = 'Outline'
    cur_imageDisplay.ColorArrayName = ['POINTS', '']
    cur_imageDisplay.ScalarOpacityUnitDistance = 1.732050807568878
    cur_imageDisplay.Slice = 152

    # reset view to fit data
    
    #if flag == 0:
     #   renderView1.ResetCamera()
      #  flag = 1

    #### saving camera placements for all active views

    # find source
    #imageReader1 = FindSource('ImageReader1')

    # create a new 'Contour'
    contour1 = Contour(Input=cur_image)
    contour1.ContourBy = ['POINTS', 'ImageFile']
    contour1.Isosurfaces = [0.1971762478351593]
    contour1.PointMergeMethod = 'Uniform Binning'

    # Properties modified on contour1
    contour1.Isosurfaces = [0.0]

    # get active view
    #renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [920, 470]

    # show data in view
    contour1Display = Show(contour1, renderView1)
    # trace defaults for the display properties.
    contour1Display.ColorArrayName = [None, '']

    # Properties modified on renderView1
    renderView1.Background = [1.0, 1.0, 1.0]

    # hide data in view
    Hide(cur_image, renderView1)

    # change solid color
    contour1Display.DiffuseColor = [0.6666666666666666, 0.0, 0.0]

    save_name = cur_path + 'Screenshots/step_' + str(i) + '.png'
    # save screenshot
    SaveScreenshot(save_name, magnification=1, quality=100, view=renderView1)

    # set active source
    #SetActiveSource(imageReader1)

    # hide data in view
    Hide(contour1, renderView1)

    # show data in view
    #imageReader1Display = Show(imageReader1, renderView1)
    # trace defaults for the display properties.
    #imageReader1Display.Representation = 'Outline'
    #imageReader1Display.ColorArrayName = ['POINTS', '']
    #imageReader1Display.ScalarOpacityUnitDistance = 1.732050807568878
    #imageReader1Display.Slice = 152

    # destroy contour1
    Delete(contour1)
    del contour1

    # destroy imageReader1
    Delete(cur_image)
    del cur_image


Delete(contour_mask)
del contour_mask

# destroy imageReader1
Delete(maskraw)
del maskraw

