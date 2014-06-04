from paraview.simple import *
from mpas_common import *
import itertools
try:
    from paraview import data_exploration as wx
except:
    wx = None

################################################################################
# Explorer pipeline elements
################################################################################

def slice_writer(options, fng):
    def create_slice_explorer():
        data = GetActiveSource()
        view = GetActiveView()
        return wx.SliceExplorer(fng, view, data,
                options['colors'],
                options.get('slices', 10),
                options.get('normal', (0, 0, 1)),
                options.get('viewup', (0, 1, 0)),
                options.get('bound_range', (0, 1)),
                options.get('scale_ratio', 2))
    return create_slice_explorer

DEFAULT_ROTATE_OPTIONS = {
    'distance': 25000000,
    'step': (15, 15)
}

def rotate_writer(options, fng):
    opts = DEFAULT_ROTATE_OPTIONS.copy()
    opts.update(options)

    def create_rotate_explorer():
        view = GetActiveView()
        return wx.ThreeSixtyImageStackExporter(fng, view,
                opts.get('focal_point', (0, 0, 0)),
                opts.get('distance', 100),
                opts.get('axis', (0, 0, 1)),
                opts.get('step', (10, 15)))
    return create_rotate_explorer

def contour_explorer(writer, options, fng):
    def create_contour_explorer():
        data = GetActiveSource()
        view = GetActiveView()
        explorer = wx.ContourExplorer(fng, data, options['contours'], options['range'], options['step'])
        proxy = explorer.getContour()
        rep = Show(proxy)

        rep.LookupTable = options['lut']
        rep.ColorArrayName = options['contours']

        internal = writer(options['inner'], fng)

        class ContourProxy(object):
            def __init__(self, explorer, loopee):
                self.explorer = explorer
                self.loopee = loopee

            def add_attribute(self, name, value):
                setattr(self, name, value)

            def UpdatePipeline(self, time=0):
                for steps in self.explorer:
                    self.loopee.UpdatePipeline(time)
                self.explorer.reset()

        return ContourProxy(explorer, internal)
    return create_contour_explorer

################################################################################
# Base explorer class
################################################################################

class MPASExplorer(object):
    def __init__(self, defaults, options):
        self.options = defaults.copy()
        self.options.update(options)

        self.sim = GetActiveSource()

        self.analysis = None
        if 'title' in self.options and 'description' in self.options:
            self.analysis = wx.AnalysisManager(
                    self.options['output'],
                    self.options['title'],
                    self.options['description'])

    def register_analysis(self, key, pattern, data_type):
        if self.analysis is None:
            self.fng = wx.FileNameGenerator(self.options['output'], pattern)
        else:
            self.analysis.register_analysis(key,
                    self.options['title'],
                    self.options['description'],
                    pattern,
                    data_type)
            self.fng = self.analysis.get_file_name_generator(key)

    def set_analysis(self, explorer):
        if self.analysis is None:
            return
        self.analysis.begin()
        explorer.set_analysis(self.analysis)

    def add_attribute(self, name, value):
        setattr(self, name, value)

    def update_layers(self, calc, field, lut, layers, explorer, time):
        if layers is not None and self.sim.CellData[field].GetNumberOfComponents() > 1:
            lut.VectorMode = 'Component'
            for layer in layers:
                lut.VectorComponent = layer
                self.fng.update_active_arguments(layer=layer)
                calc.Function = '%s_%d' % (field, layer)

                explorer.UpdatePipeline(time)
        else:
                lut.VectorMode = 'Magnitude'
                calc.Function = field
                explorer.UpdatePipeline(time)

    def Finalize(self):
        self.analysis.end()

################################################################################
# IsoLines3d
################################################################################

DEFAULT_ISOLINES_OPTIONS = {
    'output': 'output/isolines',
    'rotate_options': {},
    'layers': (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35)
}

class IsoLines3dWriter(MPASExplorer):
    def __init__(self, options):
        super(IsoLines3dWriter, self).__init__(DEFAULT_ISOLINES_OPTIONS, options)

        self.register_analysis('isolines', '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        self.layers = self.options['layers']
        self.function_pattern = '%s_%d'

        self.field = self.options['field']

        self.calc = Calculator(AttributeMode='Cell Data')
        self.thresh = Threshold(Scalars=('CELLS', self.field), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.iso_lines = Contour(PointMergeMethod='Uniform Binning',
                           ComputeScalars=0)
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface'
        self.iso_lines_rep = Show(self.iso_lines)
        self.iso_lines_rep.LineWidth = 3.0
        self.explorer = rotate_writer(self.options['rotate_options'], self.fng)()

        self.set_analysis(self.explorer)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(colorBy=field)
            self.fng.update_label_arguments(colorBy="Color")

            self.thresh_rep.LookupTable = opts['lut']
            self.thresh_rep.ColorArrayName = ('CELL_DATA', field)

            linefield = opts['isoLinesArray']
            lineopts = self.luts[linefield]
            self.iso_lines.Isosurfaces = lineopts['isoSurfaces']
            self.iso_lines_rep.LookupTable = lineopts['lut']
            self.iso_lines_rep.ColorArrayName = ('POINT_DATA', linefield)

            self.update_layers(self.calc, field, self.thresh_rep.LookupTable, self.layers, self.explorer, time)

################################################################################
# Contour3d
################################################################################

DEFAULT_CONTOUR3D_OPTIONS = {
    'output': 'output/contour3d',
    'rotate_options': {}
}

class Contour3dWriter(MPASExplorer):
    def __init__(self, options):
        super(Contour3dWriter, self).__init__(DEFAULT_CONTOUR3D_OPTIONS, options)

        self.register_analysis('contour3d', '{time}/{surfaceContour}/{contourIdx}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        self.field = self.options['field']

        self.calc = Calculator(AttributeMode='Cell Data', Function='')
        self.thresh = Threshold(Scalars=('CELLS', self.field), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.surfcont = Contour(PointMergeMethod='Uniform Binning',
                           ComputeNormals=0,
                           Isosurfaces=(0))
        self.linecont = Contour(PointMergeMethod='Uniform Binning',
                           ComputeNormals=0,
                           Isosurfaces=(0))
        self.surfcont_rep = Show(self.surfcont)
        self.linecont_rep = Show(self.linecont)
        self.explorer = rotate_writer(self.options['rotate_options'], self.fng)()

        self.set_analysis(self.explorer)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(surfaceContour=field)
            self.fng.update_label_arguments(surfaceContour="Surface contour")

            self.surfcont.ContourBy = field
            self.surfcont_rep.LookupTable = opts['lut']
            self.surfcont_rep.ColorArrayName = ('POINT_DATA', field)

            linefield = opts['isoLinesArray']
            lineopts = self.luts[linefield]
            self.linecont.ContourBy = linefield
            self.linecont.Isosurfaces = lineopts['isoSurfaces']
            self.linecont_rep.LookupTable = lineopts['lut']
            self.linecont_rep.ColorArrayName = ('POINT_DATA', linefield)

            for idx, value in itertools.izip(itertools.count(0), opts['isoSurfaces']):
                self.surfcont.Isosurfaces = [value]

                self.fng.update_active_arguments(contourIdx=idx)

                self.update_layers(self.calc, field, self.surfcont_rep.LookupTable, None, self.explorer, time)

################################################################################
# ColorBy3d
################################################################################

DEFAULT_COLORBY3D_OPTIONS = {
    'output': 'output/colorby3d',
    'rotate_options': {},
    'layers': (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35)
}

class ColorBy3dWriter(MPASExplorer):
    def __init__(self, options):
        super(ColorBy3dWriter, self).__init__(DEFAULT_COLORBY3D_OPTIONS, options)

        self.register_analysis('colorby3d', '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        self.layers = self.options['layers']
        self.field = self.options['field']

        self.calc = Calculator(AttributeMode='Cell Data')
        self.thresh = Threshold(Scalars=('CELLS', self.field), ThresholdRange=(-1000.0, 1000.0))
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface'
        self.explorer = rotate_writer(self.options['rotate_options'], self.fng)()

        self.set_analysis(self.explorer)

    def UpdatePipeline(self, time):
        for field, opts in self.luts.items():
            self.fng.update_active_arguments(colorBy=field)
            self.fng.update_label_arguments(colorBy="Color")

            self.thresh_rep.LookupTable = opts['lut']
            self.thresh_rep.ColorArrayName = ('CELL_DATA', field)

            self.update_layers(self.calc, field, self.thresh_rep.LookupTable, self.layers, self.explorer, time)

def buildIsoValues(rangeValues, nbContour):
    inc = float(rangeValues[1] - rangeValues[0]) / float(nbContour)
    values = []
    for i in range(nbContour+1):
        values.append(rangeValues[0] + (i * inc))
    return values

def buildLookupTables(luts):
    for key, lut in luts.items():
        dataRange = lut['range']
        colors = lut.get('colors', 'rainbow')
        if colors == 'red_to_blue':
            lut['lut'] = GetLookupTableForArray(
                    key, 1,
                    RGBPoints=[dataRange[0], 0.231373, 0.298039,
                               0.752941, (dataRange[0]+dataRange[1])/2, 0.865003,
                               0.865003, 0.865003, dataRange[1],
                               0.705882, 0.0156863, 0.14902],
                    NanColor=[0.0, 0.0, 0.0],
                    ColorSpace='Diverging',
                    ScalarRangeInitialized=1.0,
                    LockScalarRange=1)
        elif colors == 'rainbow':
            lut['lut'] = GetLookupTableForArray(
                    key, 1,
                    RGBPoints=[dataRange[0], 0.0, 0.0,
                               1.0, dataRange[1], 1.0,
                               0.0, 0.0, 0.0],
                    NanColor=[0.0, 0.0, 0.0],
                    ColorSpace='HSV',
                    ScalarRangeInitialized=1.0,
                    LockScalarRange=1)

        if 'nlines' in lut:
            lut['isoSurfaces'] = buildIsoValues(dataRange, lut['nlines'])
        if 'nsurfaces' in lut:
            lut['isoSurfaces'] = buildIsoValues(dataRange, lut['nsurfaces'])

    return luts


def writers(frequency, name, opts):
    clsmap = {
        'contour3d': Contour3dWriter,
        'colorby3d': ColorBy3dWriter,
        'isolines3d': IsoLines3dWriter,
        'grid': XMLUnstructuredGridWriter
    }

    def pipeline_element():
        cls = clsmap[name]
        if isinstance(cls, type(MPASExplorer)) and issubclass(cls, MPASExplorer):
            return cls(opts)
        return cls()

    writer = {
        'source': 'simulation',
        'function': pipeline_element,
        'frequency': frequency,
        'properties': {}
    }

    if 'pattern' in opts:
        writer['pattern'] = opts['pattern']

    if 'contour_arrays' in opts:
        writer['properties']['luts'] = buildLookupTables(opts['contour_arrays'])

    return writer

def mpas_add_pipeline(datasets, desc, **kwargs):
    pipe = kwargs.copy()

    if 'filters' not in pipe:
        pipe['filters'] = []
    if 'writers' not in pipe:
        pipe['writers'] = []

    # Grid setup
    pipe['grid'] = desc['grid']
    if 'fields' in desc:
        pipe['fields'] = desc['fields']

    opts = desc['configuration']

    if 'view_properties' in opts:
        pipe['view_properties'] = opts['view_properties']

    if 'earth_core' in opts:
        pipe['filters'].append({
            'name': 'earth',
            'function': LegacyVTKReader,
            'show': True,
            'kwargs': {
                'FileNames': [opts['earth_core']]
            }
        })

    pipe['writers'].append(writers(desc['frequency'], desc['exporter'], desc['configuration']))

    datasets.append(pipe)
