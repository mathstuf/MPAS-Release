try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
from paraview import data_exploration as wx

import math

DEFAULT_VIEW_OPTIONS = {
    'fit_to_screen': 1,
    'magnification': 1.0,
    'width': 960,
    'height': 540,
}
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
    'ViewSize': [1000, 500],
    'OrientationAxesVisibility': 0,
    'CenterAxesVisibility': 0
}

def MPASCreateCoProcessor(datasets, options={}):
    freqs = {}
    for (name, dataset) in datasets.items():
        grid = dataset['grid']

        # Prepare the frequency map.
        if grid not in freqs:
            freqs[grid] = []

        # Any frequencies must be known here.
        if 'image_frequency' in dataset:
            freqs[grid].append(dataset['image_frequency'])
        if 'grid_frequency' in dataset:
            freqs[grid].append(dataset['grid_frequency'])
        if 'web_view_frequency' in dataset:
            freqs[grid].append(dataset['web_view_frequency'])

    def mpas_create_pipeline(coprocessor, datadescription):
        class MPASPipeline(object):
            for (name, dataset) in datasets.items():
                grid = dataset['grid']

                image_pattern = dataset.get('image_pattern', '%s_%%t.png' % name)
                # Use pi if no images are wanted since it will never have a
                # zero modulo with an integer. Using zero causes
                # ZeroDivisionError exceptions in coprocessing.
                image_frequency = dataset.get('image_frequency', math.pi)

                view_options = DEFAULT_VIEW_OPTIONS.copy()
                if 'view' in options:
                    view_options.update(options['view'])

                view_props = DEFAULT_VIEW_PROPERTIES.copy()
                if 'view_properties' in options:
                    view_props.update(options['view_properties'])

                view = coprocessor.CreateView(CreateRenderView,
                        image_pattern,
                        image_frequency,
                        view_options['fit_to_screen'],
                        view_options['magnification'],
                        view_options['width'],
                        view_options['height'])
                for (k, v) in view_props.items():
                    setattr(view, k, v)

                producer = coprocessor.CreateProducer(datadescription, grid)
                SetActiveSource(producer)

                fields = dataset.get('fields', [])
                if fields:
                    PassArrays(CellDataArrays=fields)

                if 'filters' in dataset:
                    for filter_desc in dataset['filters']:
                        args = filter_desc.get('args', [])
                        kwargs = filter_desc.get('kwargs', {})
                        filt = filter_desc['function'](*args, **kwargs)
                        if 'properties' in filter_desc:
                            for (k, v) in filter_desc['properties'].items():
                                setattr(filt, k, v)

                # Support for in-situ image exporting
                if 'web_view_slice_frequency' in dataset:
                    def slice_wrapper(view, producer, dataset):
                        def create_slice_explorer():
                            file_generator = wx.FileNameGenerator(dataset['web_view_slice_dir'], dataset['web_view_slice_pattern'])
                            return wx.SliceExplorer(file_generator, view, producer,
                                    dataset['web_view_slice_colors'],
                                    dataset.get('web_view_slice_slices', 10),
                                    dataset.get('web_view_slice_normal', [0, 0, 1]),
                                    dataset.get('web_view_slice_viewup', [0, 1, 0]),
                                    dataset.get('web_view_slice_bound_range', [0, 1]),
                                    dataset.get('web_view_slice_scale_ratio', 2))
                        return create_slice_explorer
                    web_view_exporter = coprocessor.CreateWriter(slice_wrapper(view, producer, dataset),
                            '', dataset['web_view_slice_frequency'])

                # Support for in-situ image exporting
                # XXX: Disabled for now. Will this work?
                if False and 'web_view_rotate_frequency' in dataset:
                    def rotate_wrapper(view, producer, dataset):
                        def create_rotate_explorer():
                            file_generator = wx.FileNameGenerator(dataset['web_view_rotate_dir'], dataset['web_view_rotate_pattern'])
                            return wx.ThreeSixtyImageExporter(file_generator, view,
                                    dataset.get('web_view_rotate_focal_point', [0, 0, 0]),
                                    dataset.get('web_view_rotate_distance', 100),
                                    dataset.get('web_view_rotate_axis', [0, 0, 1]),
                                    dataset.get('web_view_rotate_step', [10, 15]))
                        return create_rotate_explorer
                    web_view_exporter = coprocessor.CreateWriter(rotate_wrapper(view, producer, dataset),
                            '', dataset['web_view_rotate_frequency'])

                if 'grid_frequency' in dataset:
                    grid_class = dataset.get('grid_class', XMLPUnstructuredGridWriter)
                    grid_pattern = dataset.get('grid_pattern', '%s_%%t.pvtu' % name)
                    grid_frequency = dataset['grid_frequency']

                    coprocessor.CreateWriter(grid_class, grid_pattern, grid_frequency)

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
