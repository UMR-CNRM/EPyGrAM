#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend _D3CommonField with plotting methods using 'vtk'.
"""

import numpy

from epygram import epygramError
from epygram.extra.usevtk import proj, modify_grid, write_grid


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([as_vtkGrid, plot3DOutline, plot3DContour,
                                plot3DColor, plot3DVolume,
                                vtk_guess_param_from_field],
                               plugin_name)
    from epygram.fields import _D3CommonField
    _D3CommonField.as_vtkGrid = as_vtkGrid
    _D3CommonField.plot3DOutline = plot3DOutline
    _D3CommonField.plot3DContour = plot3DContour
    _D3CommonField.plot3DColor = plot3DColor
    _D3CommonField.plot3DVolume = plot3DVolume
    _D3CommonField.vtk_guess_param_from_field = vtk_guess_param_from_field


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
    if hasattr(fil, 'ThresholdByUpper'):
        if maxval is None:
            fil.ThresholdByUpper(minval)
        elif minval is None:
            fil.ThresholdByLower(maxval)
        elif minval is None and maxval is None:
            pass
        else:
            fil.ThresholdBetween(minval, maxval)
    else:
        if minval is not None:
            fil.SetLowerThreshold(minval)
        if maxval is not None:
            fil.SetUpperThreshold(maxval)

    tri = vtk.vtkDataSetTriangleFilter()
    tri.SetInputConnection(fil.GetOutputPort())

    volumeMapper = {'RayCast': vtk.vtkUnstructuredGridVolumeRayCastMapper,
                    'ZSweep': vtk.vtkUnstructuredGridVolumeZSweepMapper,
                    'ProjectedTetrahedra': vtk.vtkProjectedTetrahedraMapper,
                    'OpenGLProjectedTetrahedra': vtk.vtkOpenGLProjectedTetrahedraMapper}[algo]()
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
      - None to use the field projection
      - a basemap/pyproj-like object
      - a cratopy.crs object
    """
    
    if reverseZ is None:
        reverseZ = self.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100)
    
    if hCoord in ('geoid', 'll'):
        pass
    elif callable(hCoord):
        #Should be a basemap/pyproj-like object or a cartopy.crs object
        #but we do not test it explicitly to not
        #introduce unnecessary dependency
        pass
    elif hCoord is None:
        #hCoord = self.geometry.default_cartopy_CRS()
        def hCoord(x, y, inverse=False):
            return self.geometry.xy2ll(x, y) if inverse else  self.geometry.ll2xy(x, y)
    else:
        raise ValueError("hCoord must be 'geoid', 'll', None, a basemap/pyproj-like object or a cartopy.crs")
    
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
