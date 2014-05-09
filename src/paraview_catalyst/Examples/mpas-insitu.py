from paraview.simple import *
from mpas_common import *
from mpas_exploration import *

datasets = []

earth_core_mesh = mpas_earth_core('/.../earth-core-smooth-no-data.vtk')

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_1LAYER-primal',
    'fields': ('salinity', 'temperature')

    'image_pattern': 'xyz1_%t.jpg',
    'grid_pattern': 'my_grid_%t.pvtu',
    'write_grid': True,
    'grid_frequency': 5,

    'view_size': (500, 500),

    # Available: countour3d, colorBy3d, isolines3d
    'explorer': 'isolines',
    'explorer_frequency': 5,
    'explorer_contours': {
        'temperature': {
            'range': [-1.6428141593933105, 28.691740036010742],
            'nlines': 30,
            'nsurfaces': 10,
            'colorBy': ('POINT_DATA', 'temperature'),
            'colors': 'red_to_blue',

            # Information for the isolines explorer.
            'isoLinesArray': 'salinity'
        },
        'salinity': {
            'range': [33.391498565673828, 36.110965728759766],
            'nlines': 30,
            'nsurfaces': 10,
            'colorBy': ('POINT_DATA', 'salinity'),
            'colors': 'rainbow',

            # Information for the isolines explorer.
            'isoLinesArray': 'temperature'
        }
    }

    # Explorer-specific options.
    'explorer_options': {
        # The directory to use for images from the explorer.
        'output': 'path/to/output/into',
        # The layers to use (not applicable to contour3d).
        'layers': range(15),

        'title': 'temp/salinity',
        'description': 'my simulation',

        # Options for the rotation of the camera to dump images; all are
        # optional.
        'rotate_options': {
            'distance': 25000000,
            'focal_point': (0, 0, 0),
            'axis': (0, 0, 1),
            'step': (15, 15)
        }
    },

    'earth_core': earth_core_mesh
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
