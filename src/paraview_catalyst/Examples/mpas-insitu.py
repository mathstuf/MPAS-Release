from paraview.simple import *
from mpas_common import *
from mpas_exploration import *

datasets = []

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_NLAYER-primal',
    'exporter': 'grid',
    'fields': (), # empty list == all fields
    'frequency': 5,
    'configuration': {
        'pattern': 'xyznlayer_%t.pvtu'
    }
})

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_NLAYER-primal',
    'exporter': 'colorby3d',
    'fields': ('temperature', 'salinity'),
    'frequency': 5,
    'configuration': {
        'contour_arrays': {
            'temperature': {
                'range': (-1.6428141593933105, 28.691740036010742),
                'colorBy': ('POINT_DATA', 'temperature'),
                'colors': 'red_to_blue',
                'vector_mode': 'Component'
            },
            'salinity': {
                'range': (33.391498565673828, 36.110965728759766),
                'colorBy': ('POINT_DATA', 'salinity'),
                'colors': 'rainbow',
                'vector_mode': 'Component'
            }
        },
        # The directory to use for images from the explorer.
        'output': 'colorby3d',
        # The layers to use.
        'layers': range(15),

        'title': 'colorby3d',
        'description': 'color by 3d simulation',

        'earth_core': '/.../earth-core-smooth-no-data.vtk'
    }
})

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_NLAYER-primal',
    'exporter': 'isolines3d',
    'fields': ('temperature', 'salinity'),
    'frequency': 5,
    'configuration': {
        'contour_arrays': {
            'temperature': {
                'range': (-1.6428141593933105, 28.691740036010742),
                'nlines': 30,
                'colorBy': ('POINT_DATA', 'temperature'),
                'colors': 'red_to_blue',
                'vector_mode': 'Component',

                # Information for the isolines explorer.
                'isoLinesArray': 'salinity'
            },
            'salinity': {
                'range': (33.391498565673828, 36.110965728759766),
                'nlines': 30,
                'colorBy': ('POINT_DATA', 'salinity'),
                'colors': 'rainbow',
                'vector_mode': 'Component',

                # Information for the isolines explorer.
                'isoLinesArray': 'temperature'
            }
        },
        # The directory to use for images from the explorer.
        'output': 'isolines3d',
        # The layers to use.
        'layers': range(15),

        'title': 'isolines3d',
        'description': 'isolines 3d simulation',

        'earth_core': '/.../earth-core-smooth-no-data.vtk'
    }
})

mpas_add_pipeline(datasets, {
    'grid': 'X_Y_Z_NLAYER-primal',
    'exporter': 'contour3d',
    'fields': ('temperature', 'salinity'),
    'frequency': 5,
    'configuration': {
        'contour_arrays': {
            'temperature': {
                'range': (-1.6428141593933105, 28.691740036010742),
                'colorBy': ('POINT_DATA', 'temperature'),
                'colors': 'red_to_blue',
                'vector_mode': 'Magnitude',

                # Information for the contour explorer.
                'nsurfaces': 10,
                'isoLinesArray': 'salinity'
            },
            'salinity': {
                'range': (33.391498565673828, 36.110965728759766),
                'colorBy': ('POINT_DATA', 'salinity'),
                'colors': 'rainbow',
                'vector_mode': 'Magnitude',

                # Information for the contour explorer.
                'nsurfaces': 10,
                'isoLinesArray': 'temperature'
            }
        },
        # The directory to use for images from the explorer.
        'output': 'contour3d',

        'title': 'contour3d',
        'description': 'contour 3d simulation',

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

def Finalize():
    MPASFinalize(coprocessor)
    for module in modules:
        module.Finalize(datadescription)
