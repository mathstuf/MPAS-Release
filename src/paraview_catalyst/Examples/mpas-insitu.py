from paraview.simple import *
from mpas_common import *

################################################################################
# Simulation settings
################################################################################
# Available: countour3d, colorBy3d, isolines3d
simulation = 'countour3d'

################################################################################
# View settings
################################################################################
# The view size
view_width = 600
view_height = 600

################################################################################
# Grid output settings
################################################################################
# The grid to write.
grid = 'X_Y_Z_1LAYER-primal'
# The path to write the MPAS output to
grid_pattern = '%s/%s_%t.pvtu' % (grid, grid)
# How often the grid should be written
grid_frequency = 5

################################################################################
# Spherical output settings
################################################################################
# The location to write in-situ spherical outputs to
output_directory = 'output/%s' % simulation
# How many layers to use
layer_count = 15
# How far from the center the camera should be
distance = 25000000
# The steps for theta and phi the camera should use to take snapshots
angle_steps = (15, 30)
# How often to write the sphere output
sphere_frequency = 5
# The number of contour surfaces to use
number_of_contour_surface = 10

################################################################################
# Base model settings
################################################################################
# The path to the earth core input file.
earth_core = '../earth-core-smooth-no-data.vtk'
# The VTK class to use for reading the earth core file.
earth_core_reader = LegacyVTKReader

################################################################################
################################################################################
# Editing below here should be unnecessary.
################################################################################
################################################################################

################################################################################
# IsoLines3d
################################################################################

class IsoLines3dWriter(object):
    def __init__(self):
        self.fng = wx.FileNameGenerator(output_directory, '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg')

        self.layers = tuple(range(0, layer_count, 1))
        self.function_pattern = '%s_%d'

        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.scalar = Calculator(ResultArrayName='contour',
                            Function=self.function_pattern % ('temperature', self.layers[0]))
        self.iso_lines = Contour(PointMergeMethod='Uniform Binning',
                           #ContourBy=('POINTS', 'salinity'),
                           ComputeScalars=0)
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface With Edges'
        self.thresh_rep.ColorArrayName = ('CELL_DATA', 'temperature')
        self.iso_lines_rep = Show(self.iso_lines)
        self.iso_lines_rep.ColorArrayName = ('CELL_DATA', 'temperature')
        self.iso_lines_rep.LineWidth = 3.0
        self.explorer = rotate_writer({
            'distance': distance,
            'step': angle_steps
        }, self.fng)()

    def add_attribute(self, name, value):
        setattr(self, name, value)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(colorBy=field)
            self.fng.update_label_arguments(colorBy="Color")

            self.thresh_rep.LookupTable = opts['lut']
            self.thresh_rep.ColorArrayName = ('CELL_DATA', field)

            linefield = opts['isoLinesArray']
            lineopts = self.luts[linefield]
            self.iso_lines.Isosurfaces = lineopts['isoLines']
            self.iso_lines_rep.LookupTable = lineopts['lut']
            self.iso_lines_rep.ColorArrayName = lineopts['colorBy']

            for layer in self.layers:
                self.thresh_rep.LookupTable.VectorComponent = layer

                self.fng.update_active_arguments(layer=layer)
                self.scalar.Function = self.function_pattern % (linefield, layer)

                self.explorer.UpdatePipeline(time)

################################################################################
# Contour3d
################################################################################

class Contour3dWriter(object):
    def __init__(self):
        self.fng = wx.FileNameGenerator(output_directory, '{time}/{surfaceContour}/{contourIdx}/{theta}_{phi}.jpg')

        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.surfcont = Contour(PointMergeMethod='Uniform Binning',
                           #ContourBy=('POINTS', 'temperature'),
                           ComputeNormals=0,
                           Isosurfaces=[0])
        self.linecont = Contour(PointMergeMethod='Uniform Binning',
                           #ContourBy=('POINTS', 'salinity'),
                           ComputeNormals=0,
                           Isosurfaces=[0])
        self.surfcont_rep = Show(self.surfcont)
        self.linecont_rep = Show(self.linecont)
        self.explorer = rotate_writer({
            'distance': distance,
            'step': angle_steps
        }, self.fng)()

    def add_attribute(self, name, value):
        setattr(self, name, value)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(surfaceContour=field)
            self.fng.update_label_arguments(surfaceContour="Surface contour")

            self.surfcont.ContourBy = field
            self.surfcont_rep.LookupTable = opts['lut']
            self.surfcont_rep.ColorArrayName = opts['colorBy']

            linefield = opts['isoLinesArray']
            lineopts = self.luts[linefield]
            self.linecont.ContourBy = linefield
            self.linecont.Isosurfaces = lineopts['isoLines']
            self.linecont_rep.LookupTable = lineopts['lut']
            self.linecont_rep.ColorArrayName = lineopts['colorBy']

            for idx, value in itertools.izip(itertools.count(0), opts['isoSurfaces']):
                self.surfcont.Isosurfaces = [value]

                self.fng.update_active_arguments(contourIdx=idx)

                self.explorer.UpdatePipeline(time)

################################################################################
# ColorBy3d
################################################################################

class ColorBy3dWriter(object):
    def __init__(self):
        self.fng = wx.FileNameGenerator(output_directory, '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg')

        self.layers = tuple(range(0, layer_count, 1))

        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface With Edges'
        self.thresh_rep.ColorArrayName = ('CELL_DATA', 'temperature')
        self.explorer = rotate_writer({
            'distance': distance,
            'step': angle_steps
        }, self.fng)()

    def add_attribute(self, name, value):
        setattr(self, name, value)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(colorBy=field)
            self.fng.update_label_arguments(colorBy="Color")

            self.thresh_rep.LookupTable = opts['lut']
            self.thresh_rep.ColorArrayName = ('CELL_DATA', field)

            for layer in self.layers:
                self.thresh_rep.LookupTable.VectorComponent = layer

                self.fng.update_active_arguments(layer=layer)

                self.explorer.UpdatePipeline(time)

def buildIsoValues(rangeValues, nbContour):
    inc = float(rangeValues[1] - rangeValues[0]) / float(nbContour)
    values = []
    for i in range(nbContour+1):
        values.append(rangeValues[0] + (i * inc))
    return values

vector_modes = {
    'colorby3d': 'Component',
    'countour3d': 'Magnitude'
    'isolines3d': 'Component'
}

def buildLookupTables(luts, nbSurfaces):
    for key, lut in luts.items():
        dataRange = lut['range']
        if key == 'temperature':
            lut['lut'] = GetLookupTableForArray(
                    key, 1,
                    RGBPoints=[dataRange[0], 0.231373, 0.298039,
                               0.752941, (dataRange[0]+dataRange[1])/2, 0.865003,
                               0.865003, 0.865003, dataRange[1],
                               0.705882, 0.0156863, 0.14902],
                    VectorMode=vector_modes[simulation],
                    NanColor=[0.0, 0.0, 0.0],
                    ColorSpace='Diverging',
                    ScalarRangeInitialized=1.0,
                    LockScalarRange=1)
        else:
            lut['lut'] = GetLookupTableForArray(
                    key, 1,
                    RGBPoints=[dataRange[0], 0.0, 0.0,
                               1.0, dataRange[1], 1.0,
                               0.0, 0.0, 0.0],
                    VectorMode=vector_modes[simulation],
                    NanColor=[0.0, 0.0, 0.0],
                    ColorSpace='HSV',
                    ScalarRangeInitialized=1.0,
                    LockScalarRange=1)

        lut['isoLines'] = buildIsoValues(dataRange, lut['nbLines'])
        lut['isoSurfaces'] = buildIsoValues(dataRange, nbSurfaces)

luts = {
    "temperature" : {
        "range": [-1.6428141593933105, 28.691740036010742],
        "nbLines": 30,
        "isoLinesArray": "salinity",
        "colorBy": ('POINT_DATA', 'temperature')
    },
    "salinity": {
        "range": [33.391498565673828, 36.110965728759766],
        "nbLines": 10,
        "isoLinesArray": "temperature",
        "colorBy": ('POINT_DATA', 'salinity')
    }
}
buildLookupTables(luts, number_of_contour_surface)

writers = {
    'colorby3d': ColorBy3dWriter,
    'countour3d': Contour3dWriter
    'isolines3d': IsoLines3dWriter
}

datasets = {
    simulation: {
        'grid': grid,
        'fields': ['salinity', 'temperature'],

        'view': {
            'width': view_width,
            'height': view_height
        },
        'view_properties': {
            'Background': (1, 1, 1)
        },

        'filters': [
            {
                'name': 'earth',
                'function': earth_core_reader,
                'show': True,
                'kwargs': {
                    'FileNames': [earth_core]
                }
            }
        ],

        'writers': [
            {
                'source': 'simulation',
                'function': XMLUnstructuredGridWriter,
                'pattern': grid_pattern,
                'frequency': grid_frequency
            },
            {
                'source': 'simulation',
                'function': writers[simulation],
                'frequency': sphere_frequency,
                'properties': {
                    'luts': luts
                }
            }
        ]
    }
}

coprocessor = MPASCreateCoProcessor(datasets)

# To use other scripts, you may import them here:
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
