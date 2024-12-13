#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some useful function to deal with vtk plotting:
 - proj to perform coordinate conversions (between geographical and vtk)
 - modify_grid to convert structured grid with point values into
               unstructured grid and/or with cell values
 - write_grid to write a grid on disk

A class to gather all the usefull objects (projection, window, renderer...):
 - Usevtk 
"""

from epygram.util import as_numpy_array
from epygram import epygramError

def proj(hCoord, z_factor, offset, reverseZ, geoid):
    """
    Returns a function to transform true coordinates into vtk coordinates and
    the reverse operation 
    :param hCoord: 'll': horizontal coordinates are the lon/lat values
                   'geoid': horizontal coordinates are the lon/lat values on a geoid
                   a basemap/pyproj-like object: horizontal coordinates are set according to this projection
                   a cartopy.crs object
    :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
    :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
    :param reverseZ: True if z coordinate must be reversed
    :param geoid: geoid to use (only used with hCoord='geoid')
    """
    
    def proj3d(X, Y, Z, inverse=False):
        """
        Compute coordinates to and from vtk world
        :param X, Y: longitude, latitude
        :param Z: height, altitude or pressure
        :param inverse: performs the reverse operation if reverse is True
        :return: (x, y, z) in vtk world
        """

        (x_offset, y_offset, z_offset) = offset
        z_offset = z_offset * z_factor

        if not inverse:
            lons = as_numpy_array(X)
            lats = as_numpy_array(Y)
            z = as_numpy_array(Z)
            
            z = z * z_factor
            if hCoord == 'll':
                x, y = lons.flatten(), lats.flatten()
            elif hCoord == 'geoid':
                import pyproj
                ecef = pyproj.Proj(proj='geocent', **geoid)
                lla = pyproj.Proj(proj='latlong', **geoid)
                x, y, z = pyproj.transform(lla, ecef, lons.flatten(), lats.flatten(), z.flatten(), radians=False)
                x_offset, y_offset, z_offset = pyproj.transform(lla, ecef, x_offset, y_offset, z_offset, radians=False)
            elif callable(hCoord):
                #basemap or pyproj-like object
                x, y = hCoord(lons.flatten(), lats.flatten())
                x_offset, y_offset = hCoord(x_offset, y_offset)
            else:
                from cartopy import crs as ccrs
                if isinstance(hCoord, ccrs.CRS):
                    x, y = list(hCoord.transform_points(ccrs.PlateCarree(globe=hCoord.globe), lons.flatten(), lats.flatten()).T)[0:2]
                    x_offset, y_offset = hCoord.transform_point(x_offset, y_offset, ccrs.PlateCarree(globe=hCoord.globe))
                else:
                    raise ValueError("hCoord must be 'geoid', 'll', a basemap/pyproj-like object or a cartopy.crs")
            
            #x, y, z modification to change plot aspect and axis position
            x -= x_offset
            y -= y_offset
            z -= z_offset
            if reverseZ:
                z = -z
            
            return (x, y, z)

        else:
            x = as_numpy_array(X)
            y = as_numpy_array(Y)
            z = as_numpy_array(Z)

            if reverseZ:
                z = -z

            if hCoord == 'll':
                lons, lats = x.flatten(), y.flatten()
            elif hCoord == 'geoid':
                import pyproj
                ecef = pyproj.Proj(proj='geocent', **geoid)
                lla = pyproj.Proj(proj='latlong', **geoid)
                x_offset, y_offset, z_offset = pyproj.transform(lla, ecef, x_offset, y_offset, z_offset, radians=False)
                x += x_offset
                y += y_offset
                z += z_offset
                lons, lats, z = pyproj.transform(ecef, lla, x.flatten(), y.flatten(), z.flatten(), radians=False)
            elif callable(hCoord):
                #basemap or pyproj-like object
                x_offset, y_offset = hCoord(x_offset, y_offset)
                x += x_offset
                y += y_offset
                z += z_offset
                lons, lats = hCoord(x.flatten(), y.flatten(), inverse=True)
            else:
                from cartopy import crs as ccrs
                if isinstance(hCoord, ccrs.CRS):
                    x_offset, y_offset = hCoord.transform_point(x_offset, y_offset, ccrs.PlateCarree(globe=hCoord.globe))
                    x += x_offset
                    y += y_offset
                    z += z_offset
                    lons, lats = list(ccrs.PlateCarree(globe=hCoord.globe).transform_points(hCoord, x.flatten(), y.flatten()).T)[0:2]
                else:
                    raise ValueError("hCoord must be 'geoid', 'll', a basemap/pyproj-like object or a cartopy.crs")

            z = z / z_factor

            return (lons, lats, z)

    return proj3d

def modify_grid(grid, grid_type, datamin=None):
    """
    Modifies the kind of grid
    Input grid must be an sgrid_point
    :param grid_type: can be:
        - sgrid_point: structured grid filled with points
        - sgrid_cell: structured grid filled with hexahedron
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
        - ugrid_point: unstructured grid filled with points
        - ugrid_cell: unstructured grid build filled with cells
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
    :param datamin: for an unknown reason, we need the minimum of the data
                    to transform grid into an unstructured one

    If grid_type is 'sgrid_point', the result is the grid; otherwise
    the result is the function of the last filter used.
    """
    import vtk  # @UnresolvedImport

    if grid_type not in ('sgrid_point', 'sgrid_cell', 'ugrid_point', 'ugrid_cell'):
        raise ValueError("Unknown grid type: " + grid_type)

    if grid_type in ('sgrid_cell', 'ugrid_cell'):
        interp = vtk.vtkPointDataToCellData()
        interp.SetInputData(grid)
        interp.Update()
        grid = interp
        # Values are now associated to cells

    if grid_type in ('ugrid_point', 'ugrid_cell'):
        fil = vtk.vtkThreshold()
        if grid_type in ('sgrid_cell', 'ugrid_cell'):
            fil.SetInputConnection(grid.GetOutputPort())
        else:
            fil.SetInputData(grid)
        # minScalar = grid.GetPointData().GetScalars().GetRange()[0] does not work every time
        if datamin is None:
            raise epygramError("datamin must be provided for unstructured grid types")
        minScalar = datamin
        if hasattr(fil, 'ThresholdByUpper'):
            fil.ThresholdByUpper(minScalar - 1.)
        else:
            fil.SetLowerThreshold(minScalar - 1.)
        grid = fil
        # Grid is now unstructured

    return grid


def write_grid(grid, filename, version='XML', binary=True,
               compression='ZLib', compression_level=5):
    """
    Writes a grid into a file
    :param grid: a vtk grid
    :param filename: filename to save in
    :param version: must be 'legacy' or 'XML'
    :param binary: True (default) for a binary file
    :param compression: must be None, 'LZ4' or 'ZLib'
                        only used for binary XML
    :param compression_level: between 1 and 9 (only for XML binary files Zlib-compressed)
    """
    import vtk  # @UnresolvedImport
    if not isinstance(grid, (vtk.vtkUnstructuredGrid,
                             vtk.vtkStructuredGrid)):
        grid_test = grid.GetOutput()
    else:
        grid_test = grid
    if version == 'legacy':
        if isinstance(grid_test, vtk.vtkUnstructuredGrid):
            writer = vtk.vtkUnstructuredGridWriter()
        elif isinstance(grid_test, vtk.vtkStructuredGrid):
            writer = vtk.vtkStructuredGridWriter()
        else:
            raise epygramError('Unknown grid type: ' + str(grid.__class__))
        if binary:
            writer.SetFileTypeToBinary()
        else:
            writer.SetFileTypeToASCII()
    elif version == 'XML':
        if isinstance(grid_test, vtk.vtkUnstructuredGrid):
            writer = vtk.vtkXMLUnstructuredGridWriter()
        elif isinstance(grid_test, vtk.vtkStructuredGrid):
            writer = vtk.vtkXMLStructuredGridWriter()
        else:
            raise epygramError('Unknown grid type: ' + str(grid.__class__))
        if binary:
            writer.SetDataModeToBinary()
            if compression == 'LZ4':
                writer.SetCompressorTypeToLZ4()
            elif compression == 'ZLib':
                writer.SetCompressorTypeToZLib()
                writer.GetCompressor().SetCompressionLevel(compression_level)
            #elif compression == 'LZMA':
            #    writer.SetCompressorTypeToLZMA()
            elif compression is None:
                writer.SetCompressorTypeToNone()
            else:
                raise epygramError('Unknown compression: ' + str(compression))
        else:
            writer.SetDataModeToAscii()
    else:
        raise epygramError('Unknown vtk file version: ' + str(version))
    writer.SetFileName(filename)
    if not isinstance(grid, (vtk.vtkUnstructuredGrid,
                             vtk.vtkStructuredGrid)):
        writer.SetInputConnection(grid.GetOutputPort())
    else:
        writer.SetInputData(grid)
    writer.Write()


class Usevtk(object):
    def __init__(self, background_color, window_size,
                 hCoord, z_factor, offset,
                 reverseZ, geoid=None,
                 hide_axes=False, offscreen=False,
                 interactor_style=None, title=None,
                 maximum_number_of_peels=None,
                 from_field=None,
                 existing_rendering=None, viewport_pos=None):
        """
        Creates a simple vtk environment
        :param background_color: must be a color name or a 3-tuple with RGB values
        :param window_size: must be a tuple (width, height)
        :param hCoord: 'll': horizontal coordinates are the lon/lat values
                       'geoid': horizontal coordinates are the lon/lat values on a geoid
                       a basemap/pyproj-like object: horizontal coordinates are set according to this projection
                       a cartopy.crs object
        :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
        :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
        :param reverseZ: True if z coordinate must be reversed
        :param geoid: geoid to use (only used with hCoord='geoid')
        :param hide_axes: True to hide the axes representation
        :param offscreen: True to hide window (useful when we only
                          want to produce png file instead of
                          interactively viewing the window)
        :param interactor_style: interactor style class to use (defaults
                                 to vtkInteractorStyleTrackball)
        :param title: title to give to the window
        :param maximum_number_of_peels: if None, renderer will not make use of peeling
                                        if not None, this value represent the maximum
                                        number of peels used by the renderer
                                        (peeling can be used as a workaround when different
                                         layers with opacity are not rendered in the correct
                                         order).
        :param from_field: if not None, must be an epygram field. In this case, if some of the
                           parameters among (hCoord, reverseZ, z_factor, offset, geoid)
                           are set to None, they will be guessed from the field
        :param existing_rendering: Usevtk instance. Usefull if one wants to add a viewport on a
                                   previously created Usevtk instance
        :param viewport_pos: position of the viewport in the window as (xmin, ymin, xmax, ymax)
                             in 0-1 coordinate
        """
        import vtk  # @UnresolvedImport
        
        #window and interactor
        if existing_rendering is None:
            self.window = vtk.vtkRenderWindow()
            if title is not None:
                self.window.SetWindowName(title)
            if offscreen:
                self.window.SetOffScreenRendering(True)
                self.interactor = None
                self.interactorStyle = None
            else:
                self.interactor = vtk.vtkRenderWindowInteractor()
                self.interactor.SetRenderWindow(self.window)
                if interactor_style is None:
                    self.interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
                else:
                    self.interactorStyle = interactor_style()
                self.interactor.SetInteractorStyle(self.interactorStyle)
            
                def exitCheck(obj, event):
                    if obj.GetEventPending() != 0:
                        obj.SetAbortRender(1)
                self.window.AddObserver("AbortCheckEvent", exitCheck)
                self.interactor.Initialize()
            if window_size is not None:
                self.window.SetSize(*window_size)
        else:
            self.window = existing_rendering.window
            self.interactor = existing_rendering.interactor
            self.interactorStyle = existing_rendering.interactorStyle

        #Renderer
        self.renderer = vtk.vtkRenderer()
        self.window.AddRenderer(self.renderer)
        if viewport_pos is not None:
            self.renderer.SetViewport(*viewport_pos)
        if isinstance(background_color, tuple):
            self.renderer.SetBackground(*background_color)
        else:
            self.renderer.SetBackground(vtk.vtkNamedColors().GetColor3d(background_color))
        if maximum_number_of_peels is not None:
            self.renderer.SetUseDepthPeeling(1)
            self.renderer.SetMaximumNumberOfPeels(maximum_number_of_peels)

        #Axes    
        axes = vtk.vtkAxesActor()
        axes.SetTotalLength([30., 30., 30.])
        # axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(colors.GetColor3d("Red"));
        # axes.SetXAxisLabelText("test");
        # axes.GetYAxisShaftProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Yellow"))
        # axes.GetYAxisTipProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Orange"))
        if hide_axes:
            axes.VisibilityOff()
        self.renderer.AddActor(axes)

        #Other attributes
        self.hCoord = hCoord
        self.z_factor = z_factor
        self.offset = offset
        self.reverseZ = reverseZ
        self.geoid = geoid
        
        if from_field is not None:
            #Guess missing values from field
            result = from_field.vtk_guess_param_from_field(reverseZ=self.reverseZ,
                                                           z_factor=self.z_factor, hCoord=self.hCoord,
                                                           offset=self.offset, geoid=self.geoid)
            self.reverseZ = result['reverseZ']
            self.z_factor, self.hCoord = result['z_factor'], result['hCoord']
            self.offset, self.geoid = result['offset'], result['geoid']

    @property
    def proj3d(self):
        return proj(self.hCoord, self.z_factor, self.offset, self.reverseZ, self.geoid)
    
    def print_text_in_renderer(self, text, pos, fontsize=20, color='Black'):
        """
        This function write text in a vtk renderer
        :param text: the text to render
        :param pos: tuple (x, y) of the position of the text (pixel unit)
        :param fontsize: fontsize to use for the text
        :param color: color of the text
        """
        import vtk  # @UnresolvedImport
        textActor = vtk.vtkTextActor()
        textActor.SetInput(text)
        textActor.SetDisplayPosition(*pos)  # SetPosition2(pos)
        textActor.GetTextProperty().SetFontSize(fontsize)
        textActor.GetTextProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
        self.renderer.AddActor2D(textActor)
        return textActor
    
    
    def write_png(self, filename, resolution_increase=1, enable_alpha=True):
        """
        Writes a png file with the present view on vtk window
        :param filename: name of the file to produce
        :param resolution_increase: vtk window resolution is multiplied
                                    by this factor to get the png resolution
        """
        import vtk  # @UnresolvedImport
        windowToImageFilter = vtk.vtkWindowToImageFilter()
        windowToImageFilter.SetInput(self.window)
        try:
            windowToImageFilter.SetScale(resolution_increase, resolution_increase)
        except AttributeError:
            windowToImageFilter.SetMagnification(resolution_increase)
        if enable_alpha:
            windowToImageFilter.SetInputBufferTypeToRGBA()  # also record the alpha (transparency) channel
        # windowToImageFilter.ReadFrontBufferOff() #Do not know what this means but we get bad images when uncomented
        PNGWriter = vtk.vtkPNGWriter()
        PNGWriter.SetFileName(filename)
        windowToImageFilter.Update()
        PNGWriter.SetInputConnection(windowToImageFilter.GetOutputPort())
        PNGWriter.Write()

    def plotBorder(self, color, width=2):
        """
        Plots a border around the viewport
        :param color: color of the border
        :param width: width of the border
        :return: (actor, mapper)
        """
        import vtk  # @UnresolvedImport

        #points start at upper right and proceed anti-clockwise
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(4)
        points.InsertPoint(0, 1, 1, 0)
        points.InsertPoint(1, 0, 1, 0)
        points.InsertPoint(2, 0, 0, 0)
        points.InsertPoint(3, 1, 0, 0)

        #create cells, and lines
        cells = vtk.vtkCellArray()
        cells.Initialize()

        lines = vtk.vtkPolyLine()
        
        lines.GetPointIds().SetNumberOfIds(5)
        for i in range(4):
            lines.GetPointIds().SetId(i, i)
        lines.GetPointIds().SetId(4, 0)
        cells.InsertNextCell(lines);

        #now make the polydata and display it
        poly = vtk.vtkPolyData()
        poly.Initialize()
        poly.SetPoints(points)
        poly.SetLines(cells)

        #use normalized viewport coordinates since
        #they are independent of window size
        coordinate = vtk.vtkCoordinate()
        coordinate.SetCoordinateSystemToNormalizedViewport()

        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(poly)
        mapper.SetTransformCoordinate(coordinate)

        actor = vtk.vtkActor2D()
        actor.SetMapper(mapper)
        
        if isinstance(color, tuple):
            actor.GetProperty().SetColor(color)
        else:
            actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
        actor.GetProperty().SetLineWidth(width)

        self.renderer.AddViewProp(actor)
        
        return (actor, mapper)

    def close(self):
        """
        Close the window, release memory
        
        Caution: it seems that the only way to close the window and to release memory
                 is to garbage collect the variables that hold the window, the interactor
                 and the renderer. If you have set a variable to hold one of this attribute
                 this method will not work
        """
        self.window.Finalize()
        if self.interactor is not None:
            self.interactor.TerminateApp()
        self.window = self.renderer = self.interactor = self.interactorStyle = None
        
        
