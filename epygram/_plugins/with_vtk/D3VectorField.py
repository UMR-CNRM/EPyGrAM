#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend D3VectorField with vtk vector plotting.
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import numpy

from epygram import epygramError
from epygram.util import as_numpy_array
from usevtk import modify_grid, write_grid


def activate():
    """Activate extension of _D3CommonField."""
    from epygram.fields import D3VectorField
    D3VectorField.as_vtkGrid = as_vtkGrid
    D3VectorField.vtk_guess_param_from_field = vtk_guess_param_from_field
    D3VectorField.plot3DOutline = plot3DOutline
    D3VectorField.plot3DVector = plot3DVector
    D3VectorField.plot3DStream = plot3DStream


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


def vtk_guess_param_from_field(self, *args, **kwargs):
    """Cf. D3Field.vtk_guess_param_from_field()"""
    return self.components[0].vtk_guess_param_from_field(*args, **kwargs)


def plot3DOutline(self, *args, **kwargs):
    """Cf. D3Field.plot3DOutline()"""
    return self.components[0].plot3DOutline(*args, **kwargs)


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
