
try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing


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
  datasetname = "X_Y_NLAYER-dual"
  writefrequency = 1
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      junk_vtu = coprocessor.CreateProducer( datadescription, datasetname )

      # when I know more about the fields I can add in the PassArrays filter
      #PassArrays1 = PassArrays( guiName="PassArrays1", CellDataArrays=['salinity', 'temperature'] )

      outputname = datasetname+"_%t.pvtu"
      ParallelUnstructuredGridWriter1 = coprocessor.CreateWriter( XMLPUnstructuredGridWriter, outputname, writefrequency )

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {datasetname: [writefrequency]}
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
    print 'PYTHONSCRIPT::CATALYST time step and time are ', datadescription.GetTimeStep(), datadescription.GetTime()

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

    print 'PYTHONSCRIPT::CATALYST running script -- time step and time are ', datadescription.GetTimeStep(), datadescription.GetTime()

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
