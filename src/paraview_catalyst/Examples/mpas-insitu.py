from paraview.simple import *
from mpas_common import *
from mpas_exploration import *

datasets = []

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_1LAYER-primal',
    'exporter': 'grid',
    'fields': [], # empty list == all fields
    'frequency': 5,
    'configuration': {
        'pattern': 'xyz1layer_%t.pvtu'
    }
})

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_1LAYER-primal',
    'exporter': 'colorby3d',
    'fields': ('temperature', 'salinity'),
    'frequency': 5,
    'configuration': {
        'view_properties': {
            'ViewSize': (500, 500)
        },
        'contour_arrays': {
            'temperature': {
                'range': [-1.6428141593933105, 28.691740036010742],
                'nlines': 30,
                'nsurfaces': 10,
                'colorBy': ('POINT_DATA', 'temperature'),
                'colors': 'red_to_blue',
                'vector_mode': 'Component',

                # Information for the isolines explorer.
                'isoLinesArray': 'salinity'
            },
            'salinity': {
                'range': [33.391498565673828, 36.110965728759766],
                'nlines': 30,
                'nsurfaces': 10,
                'colorBy': ('POINT_DATA', 'salinity'),
                'colors': 'rainbow',
                'vector_mode': 'Component',

                # Information for the isolines explorer.
                'isoLinesArray': 'temperature'
            }
        },
        # The directory to use for images from the explorer.
        'output': 'path/to/output/into',
        # The layers to use (not applicable to contour3d).
        'layers': range(15),

        'title': 'colorby3d',
        'description': 'color by 3d simulation',

        # Options for the rotation of the camera to dump images; all are
        # optional.
        'rotate_options': {
            'distance': 25000000,
            'focal_point': (0, 0, 0),
            'axis': (0, 0, 1),
            'step': (15, 15)
        },

        'earth_core': '/.../earth-core-smooth-no-data.vtk'
    }
})

coprocessor = MPASCreateCoProcessor(datasets)

# Add paths to any other Catalyst scripts here. They will be run with the
# current simulation.
scripts = []
modules = []
for script in scripts:
    modules.append(importlib.import_module(script))

# These functions is required and is called from Catalyst without arguments.
# Instead, pass the datasets we want to export to MPASCreateCoProcessor.
def RequestDataDescription(datadescription):
    MPASRequestDataDescription(coprocessor, datadescription)
    for module in modules:
        module.RequestDataDescription(datadescription)

def DoCoProcessing(datadescription):
    MPASDoCoProcessing(coprocessor, datadescription)
    for module in modules:
        module.DoCoProcessing(datadescription)
