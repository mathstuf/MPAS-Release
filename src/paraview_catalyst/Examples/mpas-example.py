try: paraview.simple
except: from paraview.simple import *

from mpas_common import *

# A dictionary of datasets to export is used to determine what is wanted when
# running a simulation. Each entry in the dictionary describes a grid to export
# along with information describing how it should be exported.
datasets = {
    # This is the name associated with the output. It is used in the default
    # filename patterns if another name is not given.
    'example': {
        ########################################################################
        # Core information
        #-----------------------------------------------------------------------
        # REQUIRED: This is required and is the name of the grid to export from
        # the simulation. Each entry may only export a single grid, but the
        # same grid may be used multiple times in different outputs.
        'layer': 'LON_LAT_1LAYER-primal',
        # OPTIONAL: The list of fields to use from the generated data. If
        # empty, all fields will be used. (default is empty).
        'fields': ['salinity', 'temperature']

        ########################################################################
        # Image exporting
        #-----------------------------------------------------------------------
        # REQUIRED: How often an image should be written. If not provided, no
        # images will be created.
        'image_frequency': 5,
        # OPTIONAL: The filename pattern to use when writing images. Use '%t'
        # to use the timestep of the image. The file format is automatically
        # determined from the file extension. The default is 'NAME_%t.png'
        # where NAME is the key of the output (the key in 'datasets').
        'image_pattern': 'example_%t.jpg',

        ########################################################################
        # Filters
        #-----------------------------------------------------------------------
        # Filters can be used to modify the incoming data before writing it
        # out. Filters are applied before *any* writers, so if the data type is
        # modified, remember to change the writer used.
        #-----------------------------------------------------------------------
        # OPTIONAL: The list of filters to apply, in order.
        'filters': [
            {
                ################################################################
                # Filters may be used to modify the simulation data before
                # writing it out. In this example, a filter to do the following:
                #     filter = PassArrays(guiName='PassArrays1',
                #                         CellDataArrays=['salinity'])
                #     filter.property = 4
                #---------------------------------------------------------------
                # REQUIRED: A callable value which constructs a filter
                'function': PassArrays,
                # OPTIONAL: Positional arguments to be passed to the function.
                'args': [],
                # OPTIONAL: Keyword arguments to be passed to the function.
                'kwargs': {
                    'guiName': 'PassArrays1',
                    'CellDataArrays': ['salinity']
                },
                # OPTIONAL: Properties to set on the filter.
                'properties': {
                    'property': 4
                }
            }
        ],

        ########################################################################
        # Grid exporting
        #-----------------------------------------------------------------------
        # REQUIRED: How often a grid should be written. If not provided, no
        # grids will be created.
        'grid_frequency': 5,
        # OPTIONAL: The writer class to use when writing out the grid. The
        # default is to use an XMLPUnstructuredGridWriter since all MPAS
        # outputs are unstructured grids.
        'grid_class': XMLPUnstructuredGridWriter,
        # OPTIONAL: The filename pattern to use when writing grid files. Use
        # '%t' to use the timestep of the grid. The default is 'NAME_%t.pvtu'
        # where NAME is the key of the output (the key in 'datasets').
        'grid_pattern': 'example_lonlat1_%t.pvtu',

        ########################################################################
        # In-situ exporting
        #-----------------------------------------------------------------------
        # These options will configure how images for use in the in-situ web
        # viewer will be exported.
        #-----------------------------------------------------------------------
        # REQUIRED: How often in-situ data should be written using the slice
        # exporter. If not provided, no slice data will be exported.
        'in_situ_slice_frequency': 5,
        # REQUIRED: The directory to export images into. It must exist prior to
        # execution. Relative paths are interpreted from the run directory of
        # the simulation.
        'in_situ_slice_dir': 'in_situ-slice',
        # REQUIRED: The file pattern to use when writing images. The string is
        # formatted using Python's str.format() method. The 'sliceColor' and
        # 'slicePosition' keys are required.
        'in_situ_slice_pattern': '{time}_{sliceColor}_{slicePosition}.png',
        # REQUIRED: The colors to use when rendering images.
        # TODO: Have Sebastion document these fields.
        'in_situ_slice_colors': {
            # The name of the lookup table entry.
            'temperature': {
                'lut': GetLookupTableForArray(
                        # The name of the array.
                        'temperature',
                        # The number of components.
                        1,
                        # TODO: Have Sebastion document these fields.
                        RGBPoints=[0, 0.0, 0.0, 1.0, 28.658824920654297, 1.0, 0.0, 0.0],
                        VectorMode='Magnitude',
                        NanColor=[0.498039, 0.498039, 0.498039],
                        ColorSpace='HSV',
                        ScalarRangeInitialized=1.0),
                "type": 'CELL_DATA'
            }
        },
        # OPTIONAL: The number of slices to export (default 10).
        'in_situ_slice_slices': 30,
        # OPTIONAL: The normal vector for the slice plane (default [0, 0, 1]).
        # TODO: Sebastian?
        'in_situ_slice_normal': [0, 0, 1],
        # OPTIONAL: The up vector for the camera (default [0, 1, 0]).
        # TODO: Sebastian?
        'in_situ_slice_viewup': [0, 1, 0],
        # OPTIONAL: The bounds range for the data (default [0, 1]).
        # TODO: Sebastian?
        'in_situ_slice_bound_range': [0, 1],
        # OPTIONAL: Undocumented (default 2).
        # TODO: Sebastian?
        'in_situ_slice_scale_ratio': 2,

        #-----------------------------------------------------------------------
        # XXX: Needs testing.
        # ThreeSixtyImageExporter.
        #-----------------------------------------------------------------------
        # REQUIRED: How often in-situ data should be written using the rotation
        # exporter. If not provided, no rotation data will be exported.
        'in_situ_rotate_frequency': 5,
        # REQUIRED: The directory to export images into. It must exist prior to
        # execution. Relative paths are interpreted from the run directory of
        # the simulation.
        'in_situ_rotate_dir': 'in_situ-rotate',
        # REQUIRED: The file pattern to use when writing images. The string is
        # formatted using Python's str.format() method. The 'phi' and 'theta'
        # keys are required.
        'in_situ_rotate_pattern': '{time}_{phi}_{theta}.png',
        # OPTIONAL: Where the camera should stare (default [0, 0, 0]).
        'in_situ_rotate_focal_point': [0, 0, 0],
        # OPTIONAL: The distance the camera orbits around the focal point
        # (default 100).
        'in_situ_rotate_distance': 100,
        # OPTIONAL: The axis around which the camera rotates.
        'in_situ_rotate_axis': [0, 0, 1],
        # OPTIONAL: The phi and theta steps (in degrees) the camera rotates
        # each step (phi is the angle around the main axis and theta is the
        # angle around the axis between the camera and the focal point)
        # (default [10, 15]).
        'in_situ_rotate_step': [10, 15]
    }
}

# Remove the example from the datasets in case it is wanted for documentation.
if 'example' in datasets:
    del datasets['example']

coprocessor = MPASCreateCoProcessor(datasets)

# These functions is required and is called from Catalyst without arguments.
# Instead, pass the datasets we want to export to MPASCreateCoProcessor.
def RequestDataDescription(datadescription):
    return MPASRequestDataDescription(coprocessor, datadescription)

def DoCoProcessing(datadescription):
    return MPASDoCoProcessing(coprocessor, datadescription)
