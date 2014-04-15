try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
try:
  from paraview import data_exploration as wx
except:
  wx = None

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  # MANUAL: Change entries here to change the dataset
  #         to export and its destination. Available datasets:
  #           X_Y_NLAYER-primal
  #           X_Y_NLAYER-dual
  #           X_Y_Z_1LAYER-primal
  #           X_Y_Z_1LAYER-dual
  #           X_Y_Z_NLAYER-primal
  #           X_Y_Z_NLAYER-dual
  #           LON_LAT_1LAYER-primal
  #           LON_LAT_1LAYER-dual
  #           LON_LAT_NLAYER-primal
  #           LON_LAT_NLAYER-dual
  dataset = 'LON_LAT_NLAYER-primal'
  image_freq = 100
  image_pattern = 'lonlatn_%t.png'
  grid_freq = 50
  grid_pattern = 'lonlatn_%t.pvtu'
  grid_field = 'temperature'

  # MANUAL: Web outputs
  resolution = (1000, 500)
  nslices = 30
  colorByArray = {
    "temperature": {
      "lut": GetLookupTableForArray( "temperature", 1,
                RGBPoints=[0, 0.0, 0.0, 1.0, 28.658824920654297, 1.0, 0.0, 0.0],
                VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ColorSpace='HSV',
                ScalarRangeInitialized=1.0 ),
      "type": 'CELL_DATA' },
    "pressure": {
      "lut": GetLookupTableForArray( "pressure", 1,
                RGBPoints=[0, 0.0, 0.0, 1.0, 54817300, 1.0, 0.0, 0.0],
                VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ColorSpace='HSV',
                ScalarRangeInitialized=1.0 ),
      "type": 'CELL_DATA' },
    "density": {
      "lut": GetLookupTableForArray( "density", 1,
                NanColor = [0.498039, 0.498039, 0.498039],
                RGBPoints = [1006.2609252929688, 0.0, 0.0, 1.0, 1027.6470947265625, 1.0, 0.0, 0.0], ColorSpace='HSV',
                ScalarRangeInitialized=1.0 ),
      "type": 'CELL_DATA' },
    "salinity": {
      "lut": GetLookupTableForArray( "salinity", 1,
                RGBPoints=[0, 0.0, 0.0, 1.0, 39, 1.0, 0.0, 0.0],
                VectorMode='Magnitude', NanColor=[0.498039, 0.498039, 0.498039], ColorSpace='HSV',
                ScalarRangeInitialized=1.0 ),
      "type": 'CELL_DATA' }
  }

  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      RenderView1 = coprocessor.CreateView( CreateRenderView, image_pattern, image_freq, 1, 1, 1547, 844 )
      RenderView1.CacheKey = 0.0
      RenderView1.StereoType = 0
      RenderView1.UseLight = 1
      RenderView1.StereoRender = 0
      RenderView1.CameraPosition = [15000.0, 99592.9248046875, 647834.944767225]
      RenderView1.StereoCapableWindow = 0
      RenderView1.CameraClippingRange = [441356.5953195527, 1361552.4689387335]
      RenderView1.LightSwitch = 0
      RenderView1.ViewTime = 0.0
      RenderView1.Background = [0, 0, 0]
      RenderView1.CameraFocalPoint = [15000.0, 99592.9248046875, -200000.0]
      RenderView1.CameraParallelScale = 219435.83080920158
      RenderView1.CenterOfRotation = [15000.0, 99592.9248046875, -200000.0]

      RenderView1.ViewSize = resolution
      RenderView1.OrientationAxesVisibility = 0
      RenderView1.CenterAxesVisibility = 0

      # MANUAL: More datasets may be output with one run by:
      #   - making a new producer with a different 'dataset' value;
      #   - creating slice outputs involves making a new filename generator for
      #     it (to avoid overwriting other images) and creating a new writer
      #     using it and the other processor;
      #   - if a grid writer is wanted, don't forget to call SetActiveSource
      #     with the other dataset before creating it (sliceExplorer takes the
      #     dataset in, so the active source is ignored there).

      # Take an image of the dataset.
      mpas_ocean_00 = coprocessor.CreateProducer( datadescription, dataset )

      if wx is not None:
        fileGenerator = wx.FileNameGenerator('slices', '{time}_{sliceColor}_{slicePosition}.png')
        def mkproxy(filegen, view, data):
          def sliceFactory():
            return wx.SliceExplorer(filegen, view, data, colorByArray, nslices, [0, 0, 1], [0, 1, 0], [0.2, 0.99], 2)
          return sliceFactory
        # Writes slice files for web viewing.
        sliceExplorer = coprocessor.CreateWriter(mkproxy(fileGenerator, RenderView1, mpas_ocean_00), '', grid_freq)

      # Writes .pvtu files for viewing in ParaView.
      SetActiveSource(mpas_ocean_00)
      ParallelUnstructuredGridWriter1 = coprocessor.CreateWriter( XMLPUnstructuredGridWriter, grid_pattern, grid_freq )

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {dataset: [image_freq, grid_freq]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
