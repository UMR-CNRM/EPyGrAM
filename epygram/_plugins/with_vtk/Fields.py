#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes that handle the vtk plotting on fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import numpy

from epygram import epygramError
from epygram.util import as_numpy_array

from usevtk import proj, modify_grid, write_grid

from epygram.fields import _D3CommonField, D3VectorField

def as_vtkGrid(self, rendering, grid_type,
               subzone=None,
               filename=None, name='scalar', grid=None,
               version='XML', binary=True, compression='ZLib',
               compression_level=5):
    """
    Returns a vtkStructuredGrid filled with the field
    :param rendering: a usevtk.Usevtk instance
    :param grid_type: can be:
        - sgrid_point: structured grid filled with points
        - sgrid_cell: structured grid filled with hexahedron
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
        - ugrid_point: unstructured grid filled with points
        - ugrid_cell: unstructured grid build filled with cells
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :param filename: if not None, resulting grid will be written into filename
    :param name: name to give to the scalar array (useful with the grid option)
    :param grid: if grid is not None, the method will add the data to it.
    :param version: must be 'legacy' or 'XML', used with filename
    :param binary: True (default) for a binary file, used with filename
    :param compression: must be None, 'LZ4' or 'ZLib'
                        only used for binary XML
    :param compression_level: between 1 and 9, only used for binary XML Zlib-compressed

    If grid_type is 'sgrid_point', the result is the grid; otherwise
    the result is the function is the last filter used.
    """
    import vtk # @UnresolvedImport
    from vtk.numpy_interface import dataset_adapter as dsa  # @UnresolvedImport

    if len(self.validity) != 1:
        raise NotImplementedError("For now, animation are not possible, only one validity allowed.")
    if self.spectral:
        raise epygramError("Spectral field, please use sp2gp() before.")

    data = self.getdata(d4=True, subzone=subzone).astype(numpy.float32)
    data = data[0, ...].flatten()

    if grid is None:
        grid = self.geometry.make_vtkGrid(rendering, subzone=subzone)
        grid.epygram = dict(geometry=self.geometry, subzone=subzone)
    else:
        if grid.epygram['geometry'] != self.geometry or grid.epygram['subzone'] != subzone:
            raise epygramError("To add a value to an existing grid, geometries must be the same")
        names = [grid.GetPointData().GetArrayName(i) for i in range(grid.GetPointData().GetNumberOfArrays())]
        if name in names:
            raise epygramError("There already is an array with same name: " + name)

    grid.GetPointData().AddArray(dsa.numpyTovtkDataArray(data, name, array_type=vtk.VTK_FLOAT))
    grid.GetPointData().SetActiveScalars(name)
    
    if isinstance(data, numpy.ma.masked_array):
        for index in numpy.nonzero(numpy.ma.getmaskarray(data.flatten()))[0]:
            grid.BlankPoint(index)

    grid = modify_grid(grid, grid_type, datamin=data.min())

    if filename is not None:
        write_grid(grid, filename, version=version, binary=binary,
                   compression=compression, compression_level=compression_level)
    return grid
_D3CommonField.as_vtkGrid = as_vtkGrid

def plot3DOutline(self, rendering,
                  color='Black',
                  subzone=None):
    """
    Plots the outline of the data
    :param rendering:  a usevtk.Usevtk instance
    :param color: color name to use for the outline
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, None
    """
    import vtk # @UnresolvedImport

    grid = self.as_vtkGrid(rendering, 'sgrid_point', subzone)

    # outline = vtk.vtkStructuredGridOutlineFilter()
    outline = vtk.vtkOutlineFilter()
    outline.SetInputData(grid)

    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())

    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))

    rendering.renderer.AddActor(outlineActor)
    return (outlineActor, outlineMapper, None)
_D3CommonField.plot3DOutline = plot3DOutline

def plot3DContour(self, rendering,
                  levels, color='Black', opacity=1.,
                  smoothing=None, colorbar=True,
                  subzone=None):
    """
    This method adds contour lines. If the field is 3D, contours appear as isosurface.
    :param rendering: a usevtk.Usevtk instance
    :param levels: list of values to use to compute the contour lines
    :param color: color name or lookup table or color transfer function
                  for coloring the contour lines
    :param opacity: opacity value for the contour lines
    :param smoothing: None to prevent smoothing
                      otherwise a dictionary with optional keys
                      (cf. vtkWindowedSincPolyDataFilter doc on internet):
                        - nonManifoldSmoothing (boolean)
                        - numberOfIterations (integer)
                        - featureEdgeSmoothing (boolean)
                        - boundarySmoothing (boolean)
                        - featureAngle (float)
                        - edgeAngle (float)
                        - passBand (float)
    :param colorbar: True to plot a colorbar
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, colorbaractor
    """
    import vtk # @UnresolvedImport

    ugrid = self.as_vtkGrid(rendering, 'ugrid_point', subzone)
    iso = vtk.vtkContourFilter()
    # iso.GenerateTrianglesOff()
    iso.SetInputConnection(ugrid.GetOutputPort())
    iso.ComputeScalarsOn()
    for i, value in enumerate(levels):
        iso.SetValue(i, value)

    if smoothing is not None:
        smooth = vtk.vtkWindowedSincPolyDataFilter()
        smooth.SetInputConnection(iso.GetOutputPort())
        smooth.SetNumberOfIterations(min(max(smooth.GetNumberOfIterationsMinValue(),
                                             smoothing.get('numberOfIterations',
                                                           smooth.GetNumberOfIterations())),
                                         smooth.GetNumberOfIterationsMaxValue()))
        smooth.SetFeatureEdgeSmoothing(smoothing.get('featureEdgeSmoothing',
                                                     smooth.GetFeatureEdgeSmoothing()))
        smooth.SetBoundarySmoothing(smoothing.get('boundarySmoothing',
                                                  smooth.GetBoundarySmoothing()))
        smooth.SetNonManifoldSmoothing(smoothing.get('nonManifoldSmoothing',
                                                     smooth.GetNonManifoldSmoothing()))
        smooth.SetFeatureAngle(min(max(smooth.GetFeatureAngleMinValue(),
                                       smoothing.get('featureAngle',
                                                     smooth.GetFeatureAngle())),
                                   smooth.GetFeatureAngleMaxValue()))
        smooth.SetEdgeAngle(min(max(smooth.GetEdgeAngleMinValue(),
                                    smoothing.get('edgeAngle',
                                                  smooth.GetEdgeAngle())),
                                smooth.GetEdgeAngleMaxValue()))
        smooth.SetPassBand(min(max(smooth.GetPassBandMinValue(),
                                   smoothing.get('passBand',
                                                 smooth.GetPassBand())),
                               smooth.GetPassBandMaxValue()))
        iso = smooth

    isoMapper = vtk.vtkPolyDataMapper()
    isoMapper.SetInputConnection(iso.GetOutputPort())
    
    isoActor = vtk.vtkActor()
    isoActor.SetMapper(isoMapper)

    scalarBar = None
    if isinstance(color, vtk.vtkLookupTable) or isinstance(color, vtk.vtkColorTransferFunction):
        isoMapper.ScalarVisibilityOn()
        isoMapper.SetLookupTable(color)
        isoMapper.UseLookupTableScalarRangeOn()
        if colorbar:
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(color)
            rendering.renderer.AddActor2D(scalarBar)
    else:
        isoMapper.ScalarVisibilityOff()
        isoActor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
    isoActor.GetProperty().SetOpacity(opacity)

    rendering.renderer.AddActor(isoActor)
    return (isoActor, isoMapper, scalarBar)
_D3CommonField.plot3DContour = plot3DContour

def plot3DColor(self, rendering,
                color, opacity=1.,
                colorbar=True,
                subzone=None):
    """
    This method color the field.
    :param rendering: a usevtk.Usevtk instance
    :param color: look up table or color transfer function to color the contour lines
    :param opacity: opacity value
    :param colorbar: True to plot a colorbar
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, colorbaractor
    """
    import vtk # @UnresolvedImport

    sgrid = self.as_vtkGrid(rendering, 'sgrid_point', subzone)
    sgridGeom = vtk.vtkStructuredGridGeometryFilter()
    sgridGeom.SetInputData(sgrid)

    sgridGeomMap = vtk.vtkPolyDataMapper()
    sgridGeomMap.SetInputConnection(sgridGeom.GetOutputPort())
    sgridGeomMap.SetLookupTable(color)
    sgridGeomMap.UseLookupTableScalarRangeOn()
    sgridGeomMap.InterpolateScalarsBeforeMappingOn()

    sgridGeomMapActor = vtk.vtkActor()
    sgridGeomMapActor.SetMapper(sgridGeomMap)
    sgridGeomMapActor.GetProperty().SetOpacity(opacity)

    if colorbar:
        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetLookupTable(color)
        rendering.renderer.AddActor2D(scalarBar)
    else:
        scalarBar = None

    rendering.renderer.AddActor(sgridGeomMapActor)
    return (sgridGeomMapActor, sgridGeomMap, scalarBar)
_D3CommonField.plot3DColor = plot3DColor

def plot3DVolume(self, rendering,
                 threshold, color, opacity,
                 algo='OpenGLProjectedTetrahedra',
                 colorbar=True,
                 subzone=None):
    """
    Adds a volume to the vtk rendering system
    :param rendering: a usevtk.Usevtk instance
    :param threshold: a tuple (minval, maxval) used to discard values outside of the interval
                      minval and/or maxval can be replace by None value
    :param color: a vtk.vtkColorTransferFunction object or None to describe the color
    :param opacity: a vtk.vtkPiecewiseFunction object or None to describe the alpha channel
    :param algo: among 'RayCast', 'ZSweep', 'ProjectedTetrahedra', 'OpenGLProjectedTetrahedra'
    :param colorbar: True to plot a colorbar
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, colorbaractor
    """
    import vtk # @UnresolvedImport

    minval, maxval = threshold
    grid = self.as_vtkGrid(rendering, 'ugrid_cell', subzone)

    fil = vtk.vtkThreshold()
    fil.SetInputConnection(grid.GetOutputPort())
    if maxval is None:
        fil.ThresholdByUpper(minval)
    elif minval is None:
        fil.ThresholdByLower(maxval)
    elif minval is None and maxval is None:
        pass
    else:
        fil.ThresholdBetween(minval, maxval)

    tri = vtk.vtkDataSetTriangleFilter()
    tri.SetInputConnection(fil.GetOutputPort())

    volumeMapper = {'RayCast':vtk.vtkUnstructuredGridVolumeRayCastMapper,
                    'ZSweep':vtk.vtkUnstructuredGridVolumeZSweepMapper,
                    'ProjectedTetrahedra':vtk.vtkProjectedTetrahedraMapper,
                    'OpenGLProjectedTetrahedra':vtk.vtkOpenGLProjectedTetrahedraMapper}[algo]()
    volumeMapper.SetInputConnection(tri.GetOutputPort())

    volumeProperty = vtk.vtkVolumeProperty()
    scalarBar = None
    if color is not None:
        volumeProperty.SetColor(color)
        if colorbar:
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(color)
            rendering.renderer.AddActor2D(scalarBar)
    if opacity is not None:
        volumeProperty.SetScalarOpacity(opacity)
    volumeProperty.ShadeOff()
    volumeProperty.SetInterpolationTypeToLinear()

    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    rendering.renderer.AddVolume(volume)
    return (volume, volumeMapper, scalarBar)
_D3CommonField.plot3DVolume = plot3DVolume

def vtk_guess_param_from_field(self,
                               reverseZ=None, z_factor=None,
                               hCoord=None, offset=None, geoid=None):
    """
    Guess suitable values for setting up a usevtk.Usevtk instance:
      - hCoord
      - z_factor
      - offset
      - reverseZ
      - geoid
    :param reverseZ, z_factor, hCoord, offset, geoid: if set, \
          this value are not guessed
    :return: dictionary holding all the values
    
    hCoord can be:
      - 'geoid'
      - 'll'
      - None to use the field geometry
      - a basemap object
      - a name of a proj to use among \
            ('kav7', 'ortho', 'cyl', 'moll',\
             'nsper[,sat_height=3000,lon=15.0,lat=55]')
    """
    
    if reverseZ is None:
        reverseZ = self.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100)
    
    if hCoord in ('geoid', 'll'):
        pass
    elif hCoord is not None and not isinstance(hCoord, six.string_types):
        #Should be a basemap but we do not test it explicitly to not
        #introduce a dependency on a deprecated module
        pass
    else:
        #None or a name of a specific proj
        hCoord = self.geometry.make_basemap(specificproj=hCoord)
    
    #Guess z_factor
    if z_factor is None:
        if self.geometry.isglobal and hCoord == 'geoid':
            z_factor = 500.
        else:
            proj3d = proj(hCoord, 1., (0., 0., 0.), reverseZ, self.geometry.geoid)
            if hasattr(self.geometry, 'gimme_corners_ll'):
                pos = numpy.array([proj3d(*(pos[0], pos[1], 0.)) for pos in self.geometry.gimme_corners_ll().values()])
                xpos, ypos = pos[:, 0], pos[:, 1]
                dx, dy = xpos.max() - xpos.min(), ypos.max() - ypos.min()
                dxy = max([dx, dy])
            else:
                lons, lats = self.geometry.get_lonlat_grid()
                z = numpy.zeros_like(lons)
                x, y, z = proj3d(lons, lats, z)
                dxy = numpy.sqrt((x.max() - x.min())**2 + (y.max() - y.min())**2 + (z.max() - z.min())**2)
                del lons, lats, z, x, y
            dz = numpy.array(self.geometry.vcoordinate.levels).max() - numpy.array(self.geometry.vcoordinate.levels).min()
            z_factor = dxy / dz
    
    #Guess offset
    if offset is None:
        if hasattr(self.geometry, 'gimme_corners_ll'):
            if self.geometry.isglobal:
                if hCoord == 'geoid':
                    offset = (0, 0, -self.geometry.geoid['a'] / z_factor) # divided by z_factor because earth internal is not expanded
                else:
                    offset = (0, 0, 0)
            else:
                offset = self.geometry.gimme_corners_ll()['ll']
                offset = (offset[0], offset[1], 0)
        else:
            offset = (0, 0, -self.geometry.geoid['a'] / z_factor) # divided by z_factor because earth internal is not expanded
    
    if geoid is None and self.geometry.name != 'academic':
        geoid = self.geometry.geoid
    
    return dict(hCoord=hCoord, z_factor=z_factor, offset=offset,
                reverseZ=reverseZ, geoid=geoid)
_D3CommonField.vtk_guess_param_from_field = vtk_guess_param_from_field

def as_vtkGrid(self, rendering, grid_type,
               subzone=None,
               filename=None, module_name='module', vector_name='vector',
               grid=None,
               version='XML', binary=True, compression='ZLib',
               compression_level=5):
    """
    Returns a vtkStructuredGrid filled with the field
    :param rendering: a usevtk.Usevtk instance
    :param grid_type: can be:
        - sgrid_point: structured grid filled with points
        - sgrid_cell: structured grid filled with hexahedron
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
        - ugrid_point: unstructured grid filled with points
        - ugrid_cell: unstructured grid build filled with cells
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :param filename: if not None, resulting grid will be written into filename
    :param module_name: name to give to the scalar field containing the module
                        (useful with the grid option)
    :param vector_name': name of the vector field (useful with the grid option)
    :param grid: if grid is not None, the method will add the data to it.
    :param version: must be 'legacy' or 'XML', used with filename
    :param binary: True (default) for a binary file, used with filename
    :param compression: must be None, 'LZ4' or 'ZLib'
                        only used for binary XML
    :param compression_level: between 1 and 9, only used for binary XML Zlib-compressed

    If grid_type is 'sgrid_point', the result is the grid; otherwise
    the result is the function is the last filter used.
    """
    from vtk.numpy_interface import dataset_adapter as dsa  # @UnresolvedImport

    if len(self.validity) != 1:
        raise NotImplementedError("For now, animation are not possible, only one validity allowed.")
    if self.spectral:
        raise epygramError("Spectral field, please use sp2gp() before.")

    data = self.getdata(d4=True, subzone=subzone)
    data = [d[0, ...].flatten().astype(numpy.float32) for d in data]
    while len(data) < 3:  # We need 3 components to form a vtk vector
        data.append(data[0] * 0.)

    # CAUTION: this part is non trivial and is certainly due to C versus F order
    # see: http://vtk.1045678.n5.nabble.com/Array-order-in-VTK-td5740413.html
    data = [d.flatten() for d in data]
    data = as_numpy_array(data).swapaxes(0, -1).flatten()

    if grid is None:
        grid = self.geometry.make_vtkGrid(rendering, subzone=subzone)
        grid.epygram = dict(geometry=self.geometry, subzone=subzone)
    else:
        if grid.epygram['geometry'] != self.geometry and grid.epygram['subzone'] != subzone:
            raise epygramError("To add a value to an existing grid, geometries must be the same")
        names = [grid.GetPointData().GetArrayName(i) for i in range(grid.GetPointData().GetNumberOfArrays())]
        if module_name in names:
            raise epygramError("There already is an array with same name: " + module_name)
        if vector_name in names:
            raise epygramError("There already is an array with same name: " + vector_name)

    module = self.to_module().getdata().flatten()
    if isinstance(module, numpy.ma.masked_array):
        for index in numpy.nonzero(numpy.ma.getmaskarray(module.flatten()))[0]:
            grid.BlankPoint(index)
    grid.GetPointData().AddArray(dsa.numpyTovtkDataArray(module, module_name))
    grid.GetPointData().SetActiveScalars(module_name)
    vector = dsa.numpyTovtkDataArray(data, vector_name)
    vector.SetNumberOfComponents(3)
    grid.GetPointData().AddArray(vector)
    grid.GetPointData().SetActiveVectors(vector_name)

    grid = modify_grid(grid, grid_type, datamin=data.min())

    if filename is not None:
        write_grid(grid, filename, version=version, binary=binary,
                   compression=compression, compression_level=compression_level)
    return grid
D3VectorField.as_vtkGrid = as_vtkGrid

def vtk_guess_param_from_field(self, *args, **kwargs):
    """Cf. D3Field.vtk_guess_param_from_field()"""
    return self.components[0].vtk_guess_param_from_field(*args, **kwargs)
D3VectorField.vtk_guess_param_from_field = vtk_guess_param_from_field

def plot3DOutline(self, *args, **kwargs):
    """Cf. D3Field.plot3DOutline()"""
    return self.components[0].plot3DOutline(*args, **kwargs)
D3VectorField.plot3DOutline = plot3DOutline

def plot3DVector(self, rendering,
                 samplerate=None, arrowScaleFactor=1.,
                 color='Blue', opacity=1.,
                 colorbar=True,
                 subzone=None):
    """
    This method adds contour lines and/or colorize the field. If
    the field is 3D, contours appear as isosurface.
    :param rendering: a usevtk.Usevtk instance
    :param samplerate: if not None, must be a dictionary. Allowed keys are
                    'x', 'y' and 'z' and values are the sample rate in the given
                    direction. For example {'x':3} means take one over 3 points
                    in the x direction
    :param arrowScaleFactor: scale factor used to plot the vector
    :param color: color name or lookup table or color transfer function
                  to associate colors to the vector norms
    :param opacity: opacity value
    :param colorbar: True to plot a colorbar
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, colorbaractor
    """
    import vtk # @UnresolvedImport

    # generate grid and seed grid
    if samplerate is None:
        samplerate = dict()
    grid = self.as_vtkGrid(rendering, 'sgrid_point', subzone)
    seedGrid = vtk.vtkExtractGrid()
    seedGrid.SetInputData(grid)
    seedGrid.SetSampleRate(samplerate.get('x', 1),
                           samplerate.get('y', 1),
                           samplerate.get('z', 1))

    arrowSource = vtk.vtkArrowSource()
    glyphMapper = vtk.vtkGlyph3DMapper()
    glyphMapper.SetSourceConnection(arrowSource.GetOutputPort())
    glyphMapper.SetInputConnection(seedGrid.GetOutputPort())
    glyphMapper.SetScaleFactor(arrowScaleFactor)
    glyphActor = vtk.vtkActor()
    glyphActor.SetMapper(glyphMapper)
    scalarBar = None
    if isinstance(color, (vtk.vtkLookupTable, vtk.vtkColorTransferFunction)):
        glyphMapper.SetColorModeToMapScalars()
        glyphMapper.UseLookupTableScalarRangeOn()
        glyphMapper.SetLookupTable(color)
        if colorbar:
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(color)
            rendering.renderer.AddActor2D(scalarBar)
    else:
        glyphMapper.ScalarVisibilityOff()
        glyphActor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
    glyphActor.GetProperty().SetOpacity(opacity)

    rendering.renderer.AddActor(glyphActor)
    return (glyphActor, glyphMapper, scalarBar)
D3VectorField.plot3DVector = plot3DVector

def plot3DStream(self, rendering,
                 samplerate=None,
                 maxTime=None, tubesRadius=0.1,
                 color='Blue',
                 opacity=1.,
                 plot_tube=False,
                 colorbar=True,
                 subzone=None):
    """
    This method adds contour lines and/or colorize the field. If
    the field is 3D, contours appear as isosurface.
    :param rendering: a usevtk.Usevtk instance
    :param samplerate: if not None, must be a dictionary. Allowed keys are
                    'x', 'y' and 'z' and values are the sample rate in the given
                    direction. For example {'x':3} means take one over 3 points
                    in the x direction
    :param maxTime: integration time to build the stream lines and tubes
    :param tubesRadius: radius of the tubes
    :param color: a color name, a vtk.vtkColorTransferFunction or a vtk.vtkLookupTable
                  to associate colors to the stream lines or tubes
    :param opacity: opacity value
    :param plot_tube: True to plot the tubes instead of lines
    :param colorbar: True to plot a colorbar
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :return: actor, mapper, colorbaractor
    """
    import vtk # @UnresolvedImport

    # generate grid and seed grid
    if samplerate is None:
        samplerate = dict()
    grid = self.as_vtkGrid(rendering, 'sgrid_point', subzone)
    seedGrid = vtk.vtkExtractGrid()
    seedGrid.SetInputData(grid)
    seedGrid.SetSampleRate(samplerate.get('x', 1),
                           samplerate.get('y', 1),
                           samplerate.get('z', 1))

    # Using directly seedGrid sometimes work
    # But using a vtkStructuredGridGeometryFilter in between
    # helps at suppressing some error messages about extent
    seedGeom = vtk.vtkStructuredGridGeometryFilter()
    seedGeom.SetInputConnection(seedGrid.GetOutputPort())

    scalarRange = list(grid.GetPointData().GetScalars().GetRange())

    if maxTime is None:
        maxVelocity = grid.GetPointData().GetVectors().GetMaxNorm()
        maxTime = 4.0 * grid.GetLength() / maxVelocity
    streamers = vtk.vtkStreamTracer()
    streamers.DebugOn()
    streamers.SetInputData(grid)
    streamers.SetSourceConnection(seedGeom.GetOutputPort())
    streamers.SetMaximumPropagation(maxTime)
    streamers.SetInitialIntegrationStep(.2)
    streamers.SetMinimumIntegrationStep(.01)
    streamers.Update()

    if plot_tube:
        tubes = vtk.vtkTubeFilter()
        tubes.SetInputConnection(streamers.GetOutputPort())
        tubes.SetRadius(tubesRadius)
        tubes.SetNumberOfSides(6)
        tubes.SetVaryRadius(0)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection((tubes if plot_tube else streamers).GetOutputPort())
    mapper.SetScalarRange(scalarRange[0], scalarRange[1])
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    scalarBar = None
    if isinstance(color, vtk.vtkLookupTable) or isinstance(color, vtk.vtkColorTransferFunction):
        mapper.UseLookupTableScalarRangeOn()
        mapper.InterpolateScalarsBeforeMappingOn()
        mapper.SetLookupTable(color)
        if colorbar:
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(color)
            rendering.renderer.AddActor2D(scalarBar)
    else:
        mapper.ScalarVisibilityOff()
        actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
    actor.GetProperty().SetOpacity(opacity)

    rendering.renderer.AddActor(actor)
    return (actor, mapper, scalarBar)
D3VectorField.plot3DStream = plot3DStream
