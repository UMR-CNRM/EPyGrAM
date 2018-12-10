#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for a 3D field.

Plus a function to create a Vector field from 2 scalar fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy
import sys

from footprints import proxy as fpx, FPList

from epygram import epygramError
from epygram.base import Field, FieldValidityList
from epygram.util import vtk_modify_grid, vtk_check_transform, vtk_write_grid
from . import D3Field


def make_vector_field(fX, fY):
    """
    Creates a new :class:`epygram.D3VectorField` or subclass from
    two :class:`epygram.D3Field` or subclass *fX, fY* representing resp.
    the X and Y components of the vector in the field geometry.
    """
    if not isinstance(fX, D3Field) or not isinstance(fY, D3Field):
        raise epygramError("'fX', 'fY' must be (subclass of) D3Field.")
    if fX.geometry.dimensions != fY.geometry.dimensions:
        raise epygramError("'fX', 'fY' must be share their gridpoint" +
                           " dimensions.")
    if fX.spectral_geometry != fY.spectral_geometry:
        raise epygramError("'fX', 'fY' must be share their spectral" +
                           " geometry.")
    if fX.structure != fY.structure:
        raise epygramError("'fX', 'fY' must share their structure.")

    f = fpx.field(fid={'op':'make_vector()'},
                  structure=fX.structure,
                  validity=fX.validity.copy(),
                  processtype=fX.processtype,
                  vector=True,
                  components=[fX, fY])
    return f


def psikhi2uv(psi, khi):
    """
    Compute wind (on the grid) as a D3VectorField (or subclass)
    from streamfunction **psi** and velocity potential **khi**.
    """
    (dpsidx, dpsidy) = psi.compute_xy_spderivatives()
    (dkhidx, dkhidy) = khi.compute_xy_spderivatives()
    u = dkhidx - dpsidy
    v = dkhidy + dpsidx
    u.fid = {'derivative':'u-wind'}
    v.fid = {'derivative':'v-wind'}
    u.processtype = psi.processtype
    u.processtype = psi.processtype
    u.validity = psi.validity
    v.validity = psi.validity
    return make_vector_field(u, v)


class D3VectorField(Field):
    """
    3-Dimensions Vector field class.

    This is a wrapper to a list of D3Field(s), representing the components
    of a vector projected on its geometry (the grid axes).
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                info="Type of Field geometry.",
                values=set(['3D'])),
            vector=dict(
                info="Intrinsic vectorial nature of the field.",
                type=bool,
                values=set([True])),
            validity=dict(
                info="Validity of the field.",
                type=FieldValidityList,
                optional=True,
                access='rwx',
                default=FieldValidityList()),
            components=dict(
                info="List of Fields that each compose a component of the vector.",
                type=FPList,
                optional=True,
                default=FPList([])),
            processtype=dict(
                optional=True,
                info="Generating process.")
        )
    )

##############
# ABOUT DATA #
##############

    @property
    def spectral_geometry(self):
        return self.components[0].spectral_geometry

    @property
    def spectral(self):
        """Returns True if the field is spectral."""
        return self.spectral_geometry is not None

    @property
    def geometry(self):
        return self.components[0].geometry

    def attach_components(self, *components):
        """
        Attach components of the vector to the VectorField.
        *components* must be a series of D3Field.
        """
        for f in components:
            if not isinstance(f, D3Field):
                raise epygramError("*components* must inherit from D3Field(s).")
            if f.structure != self.structure:
                raise epygramError("*components* must share the same strucuture.")
        for f in components:
            self.components.append(f)

    def sp2gp(self):
        """
        Transforms the spectral field into gridpoint, according to its spectral
        geometry. Replaces data in place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *sp2gp()* method.
        """
        for f in self.components:
            f.sp2gp()

    def gp2sp(self, spectral_geometry=None):
        """
        Transforms the gridpoint field into spectral space, according to the
        *spectral_geometry* mandatorily passed as argument. Replaces data in
        place.

        :param spectral_geometry: instance of SpectralGeometry, actually
                                  containing spectral transform subroutine (in
                                  in its own *gp2sp()* method).
        """
        for f in self.components:
            f.gp2sp(spectral_geometry=spectral_geometry)

    def getdata(self, subzone=None, **kwargs):
        """
        Returns the field data, with 2D shape if the field is not spectral,
        1D if spectral, as a tuple with data for each component.

        :param subzone: optional, among ('C', 'CI'), for LAM fields only, returns
          the data resp. on the C or C+I zone.
          Default is no subzone, i.e. the whole field.

        Shape of 2D data: (x (0) being the X component, y (1) the Y one) \n
        - Rectangular grids:\n
          grid[0,0,x] is SW, grid[-1,-1,x] is NE \n
          grid[0,-1,x] is SE, grid[-1,0,x] is NW
        - Gauss grids:\n
          grid[0,:Nj,x] is first (Northern) band of latitude, masked
          after Nj = number of longitudes for latitude j \n
          grid[-1,:Nj,x] is last (Southern) band of latitude (idem).
        """
        return [f.getdata(subzone=subzone, **kwargs) for f in self.components]

    def setdata(self, data):
        """
        Sets data to its components.

        :param data: [data_i for i components]
        """
        if len(data) != len(self.components):
            raise epygramError("data must have as many components as VectorField.")
        for i in range(len(self.components)):
            self.components[i].setdata(data[i])

    def deldata(self):
        """Empties the data."""
        for i in range(len(self.components)):
            self.components[i].deldata()

    data = property(getdata, setdata, deldata, "Accessor to the field data.")

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.

        :param level: the requested level expressed in coordinate value (Pa, m...)
        :param k: the index of the requested level
        """
        components = [comp.getlevel(level=level, k=k) for comp in self.components]
        return fpx.field(fid={'op':'getlevel()'},
                         structure=components[0].structure,
                         validity=components[0].validity.copy(),
                         processtype=components[0].processtype,
                         vector=True,
                         components=components)

    def to_module(self):
        """
        Returns a :class:`epygram.D3Field` (or subclass) whose data
        is the module of the Vector field.
        """
        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata(d4=True)
        else:
            datagp = self.getdata(d4=True)
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
        else:
            loc_sqrt = numpy.sqrt
        module = 0.
        for i in range(len(self.components)):
            module += datagp[i] ** 2
        module = loc_sqrt(module)
        f = fpx.field(geometry=self.geometry.copy(),
                      structure=self.structure,
                      fid={'op':'VectorField.to_module()'},
                      validity=self.validity.copy(),
                      processtype=self.processtype)
        f.setdata(module)
        if self.spectral:
            f.gp2sp(self.spectral_geometry)

        return f

    def compute_direction(self):
        """
        Returns a :class:`epygram.D3Field` or subclass whose data
        is the direction of the horizontal part of the Vector field
        (the two firsts components), in degrees.
        """
        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata()
        else:
            datagp = self.getdata()
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
            loc_arccos = numpy.ma.arccos
        else:
            loc_sqrt = numpy.sqrt
            loc_arccos = numpy.arccos
        module = loc_sqrt(datagp[0] ** 2 + datagp[1] ** 2)
        module_cal = numpy.where(module < 1.E-15, 1.E-15, module)
        u_norm = -datagp[0] / module_cal
        v_norm = -datagp[1] / module_cal
        numpy.clip(v_norm, -1, 1, out=v_norm)
        dd1 = loc_arccos(v_norm)
        dd2 = 2. * numpy.pi - dd1
        direction = numpy.degrees(numpy.where(u_norm >= 0., dd1, dd2))
        direction = numpy.where(module < 1.E-15, 0., direction)
        f = fpx.field(geometry=self.geometry.copy(),
                      structure=self.structure,
                      fid={'op':'VectorField.compute_direction()'},
                      validity=self.validity.copy(),
                      processtype=self.processtype)
        f.setdata(direction)
        if self.spectral:
            f.gp2sp(self.spectral_geometry)

        return f

    def reproject_wind_on_lonlat(self,
                                 map_factor_correction=True,
                                 reverse=False):
        """
        Reprojects a wind vector (u, v) from the grid axes onto real
        sphere, i.e. with components on true zonal/meridian axes.
        Other components are kept untouched.

        :param map_factor_correction: if True, apply a correction of magnitude
                                      due to map factor.
        :param reverse: if True, apply the reverse reprojection.
        """
        (lon, lat) = self.geometry.get_lonlat_grid()
        assert not self.spectral
        u = self.components[0].getdata()
        v = self.components[1].getdata()
        if self.geometry.name == 'rotated_reduced_gauss':
            (u, v) = self.geometry.reproject_wind_on_lonlat(u.compressed(),
                                                            v.compressed(),
                                                            lon.compressed(),
                                                            lat.compressed(),
                                                            map_factor_correction=map_factor_correction,
                                                            reverse=reverse)
            u = self.geometry.reshape_data(u, first_dimension='T')
            v = self.geometry.reshape_data(v, first_dimension='T')
        else:
            (u, v) = self.geometry.reproject_wind_on_lonlat(u, v, lon, lat,
                                                            map_factor_correction=map_factor_correction,
                                                            reverse=reverse)
        self.setdata([u, v] + self.components[2:])

    def map_factorize(self, reverse=False):
        """
        Multiply the field by its map factor.
        Only the first two components are affected.

        :param reverse: if True, divide.
        """
        if self.spectral:
            spgeom = self.spectral_geometry
            self.sp2gp()
            was_spectral = True
        else:
            was_spectral = False
        m = self.geometry.map_factor_field()
        if reverse:
            op = '/'
        else:
            op = '*'
        self.components[0].operation_with_other(op, m)
        self.components[1].operation_with_other(op, m)
        if was_spectral:
            self.gp2sp(spgeom)

    def compute_vordiv(self, divide_by_m=False):
        """
        Compute vorticity and divergence fields from the vector field.

        :param divide_by_m: if True, apply f = f/m beforehand, where m is the
                            map factor.
        """
        if divide_by_m:
            field = self.deepcopy()
            field.map_factorize(reverse=True)
        else:
            field = self
        (dudx, dudy) = field.components[0].compute_xy_spderivatives()
        (dvdx, dvdy) = field.components[1].compute_xy_spderivatives()
        vor = dvdx - dudy
        div = dudx + dvdy
        vor.fid = {'derivative':'vorticity'}
        div.fid = {'derivative':'divergence'}
        vor.validity = dudx.validity
        div.validity = dudx.validity

        return (vor, div)

    def remove_level(self, *args, **kwargs):
        """Cf. D3Field.remove_level()"""
        for component in self.components:
            component.remove_level(*args, **kwargs)

    def extract_subdomain(self, *args, **kwargs):
        """Cf. D3Field.extract_subdomain()"""
        result = make_vector_field(self.components[0].extract_subdomain(*args, **kwargs),
                                   self.components[1].extract_subdomain(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subdomain(*args, **kwargs))
        return result

    def extract_zoom(self, *args, **kwargs):
        """Cf. D3Field.extract_zoom()"""
        result = make_vector_field(self.components[0].extract_zoom(*args, **kwargs),
                                   self.components[1].extract_zoom(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_zoom(*args, **kwargs))
        return result

    def extract_subarray(self, *args, **kwargs):
        """Cf. D3Field.extract_subarray()"""
        result = make_vector_field(self.components[0].extract_subarray(*args, **kwargs),
                                   self.components[1].extract_subarray(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subarray(*args, **kwargs))
        return result
    
    def extract_subsample(self, *args, **kwargs):
        """Cf. D3Field.extract_subsample"""
        result = make_vector_field(self.components[0].extract_subsample(*args, **kwargs),
                                   self.components[1].extract_subsample(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subsample(*args, **kwargs))
        return result
        
    def resample(self, *args, **kwargs):
        """Cf. D3Field.resample()"""
        result = make_vector_field(self.components[0].resample(*args, **kwargs),
                                   self.components[1].resample(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.resample(*args, **kwargs))
        return result

    def resample_on_regularll(self, *args, **kwargs):
        """Cf. D3Field.resample_on_regularll()"""
        result = make_vector_field(self.components[0].resample_on_regularll(*args, **kwargs),
                                   self.components[1].resample_on_regularll(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.resample_on_regularll(*args, **kwargs))
        return result
    
    def center(self, *args, **kwargs):
        """Cf. D3Field.center()"""
        for component in self.components:
            component.center(*args, **kwargs)
    
    def select_subzone(self, *args, **kwargs):
        """Cf. D3Field.select_subzone()"""
        for component in self.components:
            component.select_subzone(*args, **kwargs)
    
    def use_field_as_vcoord(self, *args, **kwargs):
        """Cf. D3Field.use_field_as_vcoord()"""
        for component in self.components:
            component.use_field_as_vcoord(*args, **kwargs)

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def getvalue_ij(self, *args, **kwargs):
        """
        Returns the value of the different components of the field from indexes.
        """
        return [f.getvalue_ij(*args, **kwargs) for f in self.components]

    def getvalue_ll(self, *args, **kwargs):
        """
        Returns the value of the different components of the field from coordinates.
        """
        return [f.getvalue_ll(*args, **kwargs) for f in self.components]

    def min(self, subzone=None):
        """Returns the minimum value of data."""
        return [f.min(subzone=subzone) for f in self.components]

    def max(self, subzone=None):
        """Returns the maximum value of data."""
        return [f.max(subzone=subzone) for f in self.components]

    def mean(self, subzone=None):
        """Returns the mean value of data."""
        return [f.mean(subzone=subzone) for f in self.components]

    def std(self, subzone=None):
        """Returns the standard deviation of data."""
        return [f.std(subzone=subzone) for f in self.components]

    def quadmean(self, subzone=None):
        """Returns the quadratic mean of data."""
        return [f.quadmean(subzone=subzone) for f in self.components]

    def nonzero(self, subzone=None):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        return [f.nonzero(subzone=subzone) for f in self.components]

    def global_shift_center(self, longitude_shift):
        """
        Shifts the center of the geometry (and the data accordingly) by
        *longitude_shift* (in degrees). *longitude_shift* has to be a multiple
        of the grid's resolution in longitude.

        For global RegLLGeometry grids only.
        """
        if self.geometry.name != 'regular_lonlat':
            raise epygramError("only for regular lonlat geometries.")
        for f in self.components:
            f.global_shift_center(longitude_shift)

    def what(self, out=sys.stdout,
             vertical_geometry=True,
             cumulativeduration=True):
        """
        Writes in file a summary of the field.

        :param out: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        :param vertical_geometry: if True, writes the vertical geometry of the
          field.
        """
        for f in self.components:
            f.what(out,
                   vertical_geometry=vertical_geometry,
                   cumulativeduration=cumulativeduration)

    def as_vtkGrid(self, hCoord, grid_type, z_factor, offset,
                   filename=None, module_name='module', vector_name='vector',
                   grid=None):
        """
        Returns a vtkStructuredGrid filled with the field
        :param hCoord: 'll': horizontal coordinates are the lon/lat values
                       a basemap: horizontal coordinates are set according to this basemap
        :param grid_type: can be:
            - sgrid_point: structured grid filled with points
            - sgrid_cell: structured grid filled with hexahedron
                          If the field is 2D, a zero thickness is used.
                          If the field is 3D, thickness are approximately computed
            - ugrid_point: unstructured grid filled with points
            - ugrid_cell: unstructured grid build filled with cells
                          If the field is 2D, a zero thickness is used.
                          If the field is 3D, thickness are approximately computed
        :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
        :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
        :param filename: if not None, resulting grid will be written into filename
        :param module_name: name to give to the scalar field containing the module
                            (useful with the grid option)
        :param vector_name': name of the vector field (useful with the grid option)
        :param grid: if grid is not None, the method will add the data to it.
        
        If grid_type is 'sgrid_point', the output is directly the grid. Otherwise,
        the output is the last used filter.
        """
        import vtk
        from vtk.numpy_interface import dataset_adapter as dsa

        if len(self.validity) != 1:
            raise NotImplementedError("For now, animation are not possible, only one validity allowed.")
        if self.spectral:
            raise epygramError("Spectral field, please use sp2gp() before.")
        
        #data = numpy.array(self.getdata(d4=True)).astype(numpy.float64)
        data = self.getdata(d4=True)
        if not all([d.shape[0] == 1 for d in data]):
            raise NotImplementedError("For now, animation are not possible.")
        data = [d[0, ...].astype(numpy.float32) for d in data]
        while len(data) < 3: #We need 3 components to form a vtk vector
            data.append(data[0] * 0.)
            
        #CAUTION: this part is non trivial and is certainly due to C versus F order
        #see: http://vtk.1045678.n5.nabble.com/Array-order-in-VTK-td5740413.html
        data = [d.flatten() for d in data]
        data = numpy.array(data).swapaxes(0, -1).flatten()

        if grid is None:
            grid = self.geometry.make_vtkGrid(hCoord, z_factor, offset)
            grid.epygram = dict(geometry=self.geometry)
        else:
            if grid.epygram['geometry'] != self.geometry:
                raise epygramError("To add a value to an existing grid, geometries must be the same")
            names = [grid.GetPointData().GetArrayName(i) for i in range(grid.GetPointData().GetNumberOfArrays())]
            if module_name in names:
                raise epygramError("There already is an array with same name: " + module_name)
            if vector_name in names:
                raise epygramError("There already is an array with same name: " + vector_name)
        
        grid.GetPointData().AddArray(dsa.numpyTovtkDataArray(self.to_module().getdata().flatten(), module_name))
        grid.GetPointData().SetActiveScalars(module_name)
        vector = dsa.numpyTovtkDataArray(data, vector_name)
        vector.SetNumberOfComponents(3)
        grid.GetPointData().AddArray(vector)
        grid.GetPointData().SetActiveVectors(vector_name)
        
        grid = vtk_modify_grid(grid, grid_type, datamin=data.min())

        if filename is not None:
            vtk_write_grid(grid, filename)
        return grid

    def plot3DOutline(self, *args, **kwargs):
        """Cf. D3Field.plot3DOutline()"""
        return self.components[0].plot3DOutline(*args, **kwargs)

    def plot3DVector(self, rendering,
                     samplerate=None, arrowScaleFactor=1., color=None,
                     hCoord=None, z_factor=None, offset=None):
        """
        This method adds contour lines and/or colorize the field. If
        the field is 3D, contours appear as isosurface.
        :param rendering: a dictionary containing, at least, the renderer key
        :param samplerate: if not None, must be a dictionary. Allowed keys are
                        'x', 'y' and 'z' and values are the sample rate in the given
                        direction. For example {'x':3} means take one over 3 points
                        in the x direction
        :param arrowScaleFactor: scale factor used to plot the vector
        :param color: lookup table to associate colors to the vectors
        :param hCoord: 'll': horizontal coordinates are the lon/lat values
                       a basemap: horizontal coordinates are set according to this basemap
        :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
        :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
        """
        import vtk
        
        hCoord, z_factor, offset = vtk_check_transform(rendering,
                                                       self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       hCoord, z_factor, offset)

        #generate grid and seed grid
        if samplerate is None:
            samplerate = dict()
        grid = self.as_vtkGrid(hCoord, 'sgrid_point', z_factor, offset)
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
        if color is not None:
            glyphMapper.SetColorModeToMapScalars()
            glyphMapper.UseLookupTableScalarRangeOn()
            glyphMapper.SetLookupTable(color)
        glyphActor = vtk.vtkActor()
        glyphActor.SetMapper(glyphMapper)
        rendering['renderer'].AddActor(glyphActor)
        return (glyphActor, glyphMapper)
        
    def plot3DStream(self, rendering,
                     samplerate=None,
                     maxTime=None, tubesRadius=0.1,
                     color=None,
                     plot_tube=False,
                     hCoord=None, z_factor=None, offset=None):
        """
        This method adds contour lines and/or colorize the field. If
        the field is 3D, contours appear as isosurface.
        :param rendering: a dictionary containing, at least, the renderer key
        :param samplerate: if not None, must be a dictionary. Allowed keys are
                        'x', 'y' and 'z' and values are the sample rate in the given
                        direction. For example {'x':3} means take one over 3 points
                        in the x direction
        :param maxTime: integration time to build the stream lines and tubes
        :param tubesRadius: radius of the tubes
        :param color: a vtk.vtkColorTransferFunction or a vtk.vtkLookupTable
                      to associate colors to the stream lines or tubes
        :param alpha: a vtk.vtkPiecewiseFunction object or None to describe the alpha channel
        :param plot_tube: True to plot the tubes instead of lines
        :param hCoord: 'll': horizontal coordinates are the lon/lat values
                       a basemap: horizontal coordinates are set according to this basemap
        :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
        :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
        """
        import vtk
        
        hCoord, z_factor, offset = vtk_check_transform(rendering,
                                                       self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       hCoord, z_factor, offset)

        #generate grid and seed grid
        if samplerate is None:
            samplerate = dict()
        grid = self.as_vtkGrid(hCoord, 'sgrid_point', z_factor, offset)
        seedGrid = vtk.vtkExtractGrid()
        seedGrid.SetInputData(grid)
        seedGrid.SetSampleRate(samplerate.get('x', 1),
                               samplerate.get('y', 1),
                               samplerate.get('z', 1))
        
        #Using directly seedGrid sometimes work
        #But using a vtkStructuredGridGeometryFilter in between
        #helps at suppressing some error messages about extent
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
        if color is not None:
            mapper.UseLookupTableScalarRangeOn()
            mapper.InterpolateScalarsBeforeMappingOn()
            mapper.SetLookupTable(color)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        rendering['renderer'].AddActor(actor)
        return (actor, mapper)
        
#############
# OPERATORS #
#############

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """
        if 'vector' not in other._attributes:
            raise epygramError("cannot operate a Vector field with a" +
                               " non-Vector one.")
        else:
            if isinstance(other, self.__class__):
                if len(self.components) != len(other.components):
                    raise epygramError("vector fields must have the same" +
                                       " number of components.")
            super(D3VectorField, self)._check_operands(other)

    def __add__(self, other):
        """
        Definition of addition, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'+'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] + other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] + other for i in range(len(self.components))]

        newid = {'op':'+'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)

        return newfield

    def __mul__(self, other):
        """
        Definition of multiplication, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'*'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] * other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] * other for i in range(len(self.components))]
        newid = {'op':'*'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __sub__(self, other):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] - other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] - other for i in range(len(self.components))]
        newid = {'op':'-'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __div__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] / other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] / other for i in range(len(self.components))]
        newid = {'op':'/'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, other):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [other.components[i] - self.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [other - self.components[i] for i in range(len(self.components))]
        newid = {'op':'-'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __rdiv__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [other.components[i] / self.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [other / self.components[i] for i in range(len(self.components))]
        newid = {'op':'/'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield
