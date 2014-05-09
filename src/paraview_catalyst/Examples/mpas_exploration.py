from paraview.simple import *
from mpas_common import *
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
                options.get('normal', [0, 0, 1]),
                options.get('viewup', [0, 1, 0]),
                options.get('bound_range', [0, 1]),
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
                opts.get('focal_point', [0, 0, 0]),
                opts.get('distance', 100),
                opts.get('axis', [0, 0, 1]),
                opts.get('step', [10, 15]))
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
        explorer.set_analysis(self.analysis)

    def add_attribute(self, name, value):
        setattr(self, name, value)

################################################################################
# IsoLines3d
################################################################################

DEFAULT_ISOLINES_OPTIONS = {
    'output': 'output/isolines',
    'rotate_options': {},
    'layers': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35]
}

class IsoLines3dWriter(MPASExplorer):
    def __init__(self, options):
        super(IsoLines3dWriter, self).__init__(DEFAULT_ISOLINES_OPTIONS, options)

        self.register_analysis('isolines', '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        self.layers = self.options['layers']
        self.function_pattern = '%s_%d'

        # TODO: Un-hardcode temperature here.
        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.scalar = Calculator(ResultArrayName='contour',
                            Function=self.function_pattern % ('temperature', self.layers[0]))
        self.iso_lines = Contour(PointMergeMethod='Uniform Binning',
                           ComputeScalars=0)
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface With Edges'
        self.thresh_rep.ColorArrayName = ('CELL_DATA', 'temperature')
        self.iso_lines_rep = Show(self.iso_lines)
        self.iso_lines_rep.ColorArrayName = ('CELL_DATA', 'temperature')
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

DEFAULT_CONTOUR3D_OPTIONS = {
    'output': 'output/contour3d',
    'rotate_options': {}
}

class Contour3dWriter(MPASExplorer):
    def __init__(self, options):
        super(Contour3dWriter, self).__init__(DEFAULT_CONTOUR3D_OPTIONS, options)

        self.register_analysis('contour3d', '{time}/{surfaceContour}/{contourIdx}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        # TODO: Un-hardcode temperature here.
        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.dataconv = CellDatatoPointData()
        self.surfcont = Contour(PointMergeMethod='Uniform Binning',
                           ComputeNormals=0,
                           Isosurfaces=[0])
        self.linecont = Contour(PointMergeMethod='Uniform Binning',
                           ComputeNormals=0,
                           Isosurfaces=[0])
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

DEFAULT_COLORBY3D_OPTIONS = {
    'output': 'output/colorby3d',
    'rotate_options': {},
    'layers': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35]
}

class ColorBy3dWriter(MPASExplorer):
    def __init__(self, options):
        super(ColorBy3dWriter, self).__init__(DEFAULT_COLORBY3D_OPTIONS, options)

        self.register_analysis('colorby3d', '{time}/{colorBy}/{layer}/{theta}_{phi}.jpg',
                wx.ThreeSixtyImageStackExporter.get_data_type())

        self.layers = self.options['layers']

        self.thresh = Threshold(Scalars=('CELLS', 'temperature'), ThresholdRange=(-1000.0, 1000.0))
        self.thresh_rep = Show(self.thresh)
        self.thresh_rep.EdgeColor = (0.0, 0.0, 0.0)
        self.thresh_rep.Representation = 'Surface With Edges'
        self.thresh_rep.ColorArrayName = ('CELL_DATA', 'temperature')
        self.explorer = rotate_writer(self.options['rotate_options'], self.fng)()

        self.set_analysis(self.explorer)

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
                    VectorMode=lut['vector_mode'],
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
                    VectorMode=lut['vector_mode'],
                    NanColor=[0.0, 0.0, 0.0],
                    ColorSpace='HSV',
                    ScalarRangeInitialized=1.0,
                    LockScalarRange=1)

        if 'nlines' in lut:
            lut['isoLines'] = buildIsoValues(dataRange, lut['nlines'])
        if 'nsurfaces' in lut:
            lut['isoSurfaces'] = buildIsoValues(dataRange, lut['nsurfaces'])

    return luts


def mpas_earth_core(path, reader=LegacyVTKReader):
    return {
        'name': 'earth',
        'function': reader,
        'show': True,
        'kwargs': {
            'FileNames': [path]
        }
    }

def writers(name, opts):
    clsmap = {
        'contour3d': Contour3dWriter,
        'colorby3d': ColorBy3dWriter,
        'isolines3d': IsoLines3dWriter
    }

    def pipeline_element():
        cls = clsmap[name]
        if isinstance(cls, type(MPASExplorer)) and issubclass(cls, MPASExplorer):
            return cls(opts)
        return cls()

    return pipeline_element

def mpas_add_pipeline(datasets, desc, **kwargs):
    sim = kwargs.copy()

    if 'view' not in sim:
        sim['view'] = {}
    if 'filters' not in sim:
        sim['filters'] = []
    if 'writers' not in sim:
        sim['writers'] = []

    # Grid setup
    sim['grid'] = desc['grid']
    sim['image_pattern'] = desc['image_pattern']
    sim['fields'] = desc['fields']

    if 'view_size' in desc:
        size = desc['view_size']
        sim['view']['width'] = size[0]
        sim['view']['height'] = size[1]

    if 'earth_core' in desc:
        sim['filters'].append(desc['earth_core'])

    if desc.get('write_grid', False):
        sim['writers'].append({
            'source': 'simulation',
            'function', XMLUnstructuredGridWriter,
            'pattern': desc.get('grid_pattern', '%(grid)s/%(grid)s_%%t.pvtu' % sim),
            'frequency': desc.get('grid_frequency', 5)
        })

    if 'explorer' in desc:
        explorer = desc['explorer']
        luts = desc['explorer_contours']
        explorer_luts = buildLookupTables(luts)
        sim['writers'].append({
            'source': 'simulation',
            'function': writers(simulation, desc.get(desc.get('explorer_options', {}))),
            'frequency', desc.get('explorer_frequency', 5)
            'properties': {
                'luts', explorer_luts
            }
        })

    datasets.append(sim)
