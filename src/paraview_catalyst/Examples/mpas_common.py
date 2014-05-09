try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing

DEFAULT_VIEW_PROPERTIES = {
    'CacheKey': 0.0,
    'StereoType': 0,
    'UseLight': 1,
    'StereoRender': 0,
    'CameraPosition': [15000.0, 99592.9248046875, 647834.944767225],
    'StereoCapableWindow': 0,
    'CameraClippingRange': [441356.5953195527, 1361552.4689387335],
    'LightSwitch': 0,
    'ViewTime': 0.0,
    'Background': [0, 0, 0],
    'CameraFocalPoint': [15000.0, 99592.9248046875, -200000.0],
    'CameraParallelScale': 219435.83080920158,
    'CenterOfRotation': [15000.0, 99592.9248046875, -200000.0],
    'ViewSize': [500, 500],
    'OrientationAxesVisibility': 0,
    'CenterAxesVisibility': 0
}

def MPASCreateCoProcessor(datasets, options={}):
    freqs = {}
    for dataset in datasets:
        grid = dataset['grid']

        # Prepare the frequency map.
        if grid not in freqs:
            freqs[grid] = []

        # Any frequencies must be known here.
        for writer in dataset['writers']:
            freqs[grid].append(writer['frequency'])

    def mpas_create_pipeline(coprocessor, datadescription):
        class MPASPipeline(object):
            for dataset in datasets:
                grid = dataset['grid']

                view_props = DEFAULT_VIEW_PROPERTIES.copy()
                if 'view_properties' in options:
                    view_props.update(options['view_properties'])

                view = CreateRenderView()
                for (k, v) in view_props.items():
                    setattr(view, k, v)

                filters = {}
                producer = coprocessor.CreateProducer(datadescription, grid)
                SetActiveSource(producer)
                filters['simulation'] = producer

                fields = dataset.get('fields', [])
                if fields:
                    filters['simulation_orig'] = producer
                    filters['simulation'] = PassArrays(CellDataArrays=fields)

                if 'filters' in dataset:
                    for filter_desc in dataset['filters']:
                        args = filter_desc.get('args', [])
                        kwargs = filter_desc.get('kwargs', {})
                        filt = filter_desc['function'](*args, **kwargs)
                        if 'source' in filter_desc:
                            source = filter_desc['source']
                            if source in filters:
                                SetActiveSource(filters[source])
                            else:
                                raise RuntimeError('Unknown source for filter: %s' % source)
                        if 'name' in filter_desc:
                            filters[filter_desc['name']] = filt
                        if 'show' in filter_desc:
                            if filter_desc['show']:
                                Show(filt)
                            else:
                                Hide(filt)
                        if 'properties' in filter_desc:
                            for (k, v) in filter_desc['properties'].items():
                                setattr(filt, k, v)

                for writer in dataset['writers']:
                    if 'source' in writer:
                        source = writer['source']
                        if source in filters:
                            SetActiveSource(filters[source])
                        else:
                            raise RuntimeError('Unknown source for web view: %s' % source)
                    pattern = writer.get('pattern', '')
                    writer_obj = coprocessor.CreateWriter(writer['function'], pattern, writer['frequency'])
                    if 'properties' in writer:
                        for (k, v) in writer['properties'].items():
                            setattr(writer_obj, k, v)

        return MPASPipeline()

    class MPASCoProcessor(coprocessing.CoProcessor):
        def CreatePipeline(self, datadescription):
            self.Pipeline = mpas_create_pipeline(self, datadescription)

    coprocessor = MPASCoProcessor()
    coprocessor.SetUpdateFrequencies(freqs)

    # Enable Live-Visualizaton with ParaView
    coprocessor.EnableLiveVisualization(False)

    return coprocessor

# ---------------------- Data Selection method ----------------------

def MPASRequestDataDescription(coprocessor, datadescription):
    "Callback to populate the request for current timestep"
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

def MPASDoCoProcessing(coprocessor, datadescription):
    "Callback to do co-processing for current timestep"
    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
