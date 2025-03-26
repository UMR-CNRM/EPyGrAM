#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for 3D geometries of fields.
"""

from functools import lru_cache
import numpy
import math
import copy
import sys

import footprints
from footprints import proxy as fpx

from epygram import epygramError, config
from epygram.config import rounding_decimal as _rd
from epygram.util import (degrees_nearest_mod, Angle,
                          write_formatted,
                          as_numpy_array,
                          is_scalar, Comparator)

from .VGeometry import VGeometry
from .AbstractGeometry import RectangularGridGeometry

epylog = footprints.loggers.getLogger(__name__)



class ProjectedGeometry(RectangularGridGeometry):
    """
    Handles the geometry for a Projected 3-Dimensions Field.
    """

    _ghost_attributes = RectangularGridGeometry._ghost_attributes + ['_proj']

    def __init__(self, name, grid, dimensions, vcoordinate, projection,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be among ['lambert', 'mercator', 'polar_stereographic',
                                         'space_view']
        :param grid: Handles description of the horizontal grid.
        :param dimensions: Handles grid dimensions.
        :param vcoordinate: Handles vertical geometry parameters.
        :param position_on_horizontal_grid: Position of points w/r to the horizontal.
                                            among: ['upper-right', 'upper-left',
                                                    'lower-left', 'lower-right',
                                                    'center-left', 'center-right',
                                                    'lower-center', 'upper-center',
                                                    'center', '__unknown__']
        :param geoid: To specify geoid shape.
        """
        self.add_attr_inlist('name', ['lambert', 'mercator', 'polar_stereographic',
                                      'space_view'])
        self.add_attr_dict('projection')

        self.name = name
        self.projection = projection

        super(ProjectedGeometry, self).__init__(grid, dimensions, vcoordinate,
                                                      position_on_horizontal_grid, geoid)

    @property
    @lru_cache
    def _get_proj(self):
        """
        Returns (proj, K, center_lon, center_lat)
        """
        import pyproj

        def compute_center_proj(p, center):
            if center == self.grid['input_position']:
                _center_lon = self.grid['input_lon']
                _center_lat = self.grid['input_lat']
            else:
                # x1, y1: coordinates in non rotated proj of input point
                x1, y1 = p(float(self.grid['input_lon'].get('degrees')),
                           float(self.grid['input_lat'].get('degrees')))
                # offset between center and input points is known in rotated proj
                # dx, dy is the offset in non rotated proj
                (dx, dy) = self._rotate_axis(
                    round((center[0] - self.grid['input_position'][0]) *
                          self.grid['X_resolution'],
                          _rd),
                    round((center[1] - self.grid['input_position'][1]) *
                          self.grid['Y_resolution'],
                          _rd),
                    'xy2ll')
                # xc, yc: coordinates of center point in non rotated proj
                xc = x1 + dx
                yc = y1 + dy
                center_lon, center_lat = p(xc, yc, inverse=True)
                _center_lon = Angle(center_lon, 'degrees')
                _center_lat = Angle(center_lat, 'degrees')
            return _center_lon, _center_lat

        projdict = {'lambert':'lcc',
                    'mercator':'merc',
                    'polar_stereographic':'stere',
                    'space_view':'geos'}
        proj = projdict[self.name]
        # build proj
        if self.grid['LAMzone'] is not None:
            nx = self.dimensions['X_CIzone']
            ny = self.dimensions['Y_CIzone']
            io = self.dimensions.get('X_CIoffset', 0)
            jo = self.dimensions.get('Y_CIoffset', 0)
            centerPoint = (io + (float(nx) - 1) / 2.,
                           jo + (float(ny) - 1) / 2.)  # Coordinates of center point
        else:
            nx = self.dimensions['X']
            ny = self.dimensions['Y']
            centerPoint = ((float(nx) - 1) / 2.,
                           (float(ny) - 1) / 2.)  # Coordinates of center point
        # CAUTION: x0, y0 are deeply linked with ij2xy and xy2ij methods:
        # the origin of the grid is the center of the domain !
        if self.name == 'lambert':
            if self.secant_projection:
                lat_1 = self.projection['secant_lat1'].get('degrees')
                lat_2 = self.projection['secant_lat2'].get('degrees')
                m1 = math.cos(math.radians(lat_1))
                m2 = math.cos(math.radians(lat_2))
                t1 = math.tan(math.pi / 4. - math.radians(lat_1) / 2.)
                t2 = math.tan(math.pi / 4. - math.radians(lat_2) / 2.)
                _K = (math.log(m1) - math.log(m2)) / \
                     (math.log(t1) - math.log(t2))
            else:
                lat_1 = self.projection['reference_lat'].get('degrees')
                lat_2 = self.projection['reference_lat'].get('degrees')
                _K = abs(self.projection['reference_lat'].get('cos_sin')[1])
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_1=lat_1, lat_2=lat_2,
                            **self.geoid)
            _center_lon, _center_lat = compute_center_proj(p, centerPoint)
            x0, y0 = p(_center_lon.get('degrees'),
                       _center_lat.get('degrees'))
            _proj = pyproj.Proj(proj=proj,
                                lon_0=self.projection['reference_lon'].get('degrees'),
                                lat_1=lat_1, lat_2=lat_2,
                                x_0=-x0, y_0=-y0,
                                **self.geoid)
        elif self.name == 'mercator':
            _K = None
            if self.secant_projection:
                lat_ts = self.projection['secant_lat'].get('degrees')
            else:
                lat_ts = 0.
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_ts=lat_ts,
                            **self.geoid)
            _center_lon, _center_lat = compute_center_proj(p, centerPoint)
            x0, y0 = p(_center_lon.get('degrees'),
                       _center_lat.get('degrees'))
            _proj = pyproj.Proj(proj=proj,
                                lon_0=self.projection['reference_lon'].get('degrees'),
                                lat_ts=lat_ts,
                                x_0=-x0, y_0=-y0,
                                **self.geoid)
        elif self.name == 'polar_stereographic':
            _K = None
            lat_0 = self.projection['reference_lat'].get('degrees')
            if self.secant_projection:
                lat_ts = self.projection['secant_lat'].get('degrees')
            else:
                lat_ts = self.projection['reference_lat'].get('degrees')
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_0=lat_0, lat_ts=lat_ts,
                            **self.geoid)
            _center_lon, _center_lat = compute_center_proj(p, centerPoint)
            x0, y0 = p(_center_lon.get('degrees'),
                       _center_lat.get('degrees'))
            _proj = pyproj.Proj(proj=proj,
                                lon_0=self.projection['reference_lon'].get('degrees'),
                                lat_0=lat_0, lat_ts=lat_ts,
                                x_0=-x0, y_0=-y0,
                                **self.geoid)
        elif self.name == 'space_view':
            _K = None
            latSat = self.projection['satellite_lat'].get('degrees')
            lonSat = self.projection['satellite_lon'].get('degrees')
            height = self.projection['satellite_height']  # Height above ellipsoid
            if latSat != 0:
                raise epygramError("Only space views with satellite_lat=0 are allowed")
            p = pyproj.Proj(proj=proj,
                            h=height,
                            lon_0=lonSat,
                            **self.geoid)
            _center_lon, _center_lat = compute_center_proj(p, centerPoint)
            x0, y0 = p(_center_lon.get('degrees'),
                       _center_lat.get('degrees'))
            _proj = pyproj.Proj(proj=proj,
                                h=height,
                                lon_0=lonSat,
                                x_0=-x0, y_0=-y0,
                                **self.geoid)
        else:
            raise NotImplementedError("projection: " + self.name)

        return (_proj, _K, _center_lon, _center_lat)

    @property
    def _center_lon(self):
        """Get the center's longitude"""
        proj, K, center_lon, center_lat = self._get_proj
        return center_lon

    @property
    def _center_lat(self):
        """Get the center's latitude"""
        proj, K, center_lon, center_lat = self._get_proj
        return center_lat

    @property
    def _K(self):
        """Get the K parameter"""
        proj, K, center_lon, center_lat = self._get_proj
        return K

    @property
    def _proj(self):
        """Get the projection"""
        proj, K, center_lon, center_lat = self._get_proj
        return proj

    @property
    def secant_projection(self):
        """ Is the projection secant to the sphere ? (or tangent)"""
        return ('secant_lat' in self.projection or
                'secant_lat1' in self.projection)

    def tolerant_equal(self, other, tolerance=config.epsilon):
        if self.__class__ == other.__class__:
            # create copies of inner objects to filter some ghost attributes
            almost_self = {k:copy.deepcopy(self.__dict__[k])
                           for k in self.__dict__.keys()
                           if k not in self._ghost_attributes}
            almost_other = {k: copy.deepcopy(other.__dict__[k])
                            for k in other.__dict__.keys()
                            if k not in other._ghost_attributes}
            # (the same grid could be defined from different inputs)
            for almost in (almost_self, almost_other):
                for k in ('input_lon', 'input_lat', 'input_position'):
                    almost['grid'].pop(k)
            almost_self['grid']['center'] = self.getcenter()
            almost_other['grid']['center'] = other.getcenter()
            return Comparator.are_equal(almost_self, almost_other, tolerance)
        else:
            return False

    def _rotate_axis(self, x, y, direction):
        """
        Internal method used to rotate x/y coordinates to handle rotated geometries.

        :param direction:, if evals to 'xy2ll', direction used when converting (x, y) to (lat, lon).
                           if evals to 'll2xy', direction used when converting (lat, lon) to (x, y).
        """
        if self.projection['rotation'].get('degrees') == 0:
            return (x, y)
        else:
            beta = self.projection['rotation'].get('radians')
            if direction == 'xy2ll':
                return (numpy.array(x) * numpy.cos(-beta) + numpy.array(y) * numpy.sin(-beta),
                        - numpy.array(x) * numpy.sin(-beta) + numpy.array(y) * numpy.cos(-beta))
            elif direction == 'll2xy':
                return (numpy.array(x) * numpy.cos(beta) + numpy.array(y) * numpy.sin(beta),
                        - numpy.array(x) * numpy.sin(beta) + numpy.array(y) * numpy.cos(beta))
            else:
                raise epygramError('Wrong direction of rotation.')

    def _consistency_check(self):
        """Check that the geometry is consistent."""
        if self.name == 'lambert':
            sets_of_keys = (['reference_lon', 'reference_lat',
                             'rotation'],
                            ['reference_lon', 'secant_lat1', 'secant_lat2',
                             'rotation'])
        elif self.name in ('polar_stereographic', 'mercator'):
            sets_of_keys = (['reference_lon', 'reference_lat',
                             'rotation'],
                            ['reference_lon', 'reference_lat', 'secant_lat',
                             'rotation'])
        elif self.name == 'space_view':
            sets_of_keys = ['satellite_lat', 'satellite_lon',
                            'satellite_height',
                            'rotation',
                            'reference_lon']
            sets_of_keys = (sets_of_keys, sets_of_keys)
        else:
            raise NotImplementedError("projection: " + self.name)
        if set(self.projection.keys()) != set(sets_of_keys[0]) and \
           set(self.projection.keys()) != set(sets_of_keys[1]):
            # print self.projection.keys()
            raise epygramError("attributes for projection " + self.name +
                               " must consist in keys: " +
                               str(sets_of_keys[0]) + " or " +
                               str(sets_of_keys[1]))
        for a in ['reference_lon', 'reference_lat', 'rotation',
                  'secant_lat1', 'secant_lat2', 'secant_lat',
                  'satellite_lat', 'satellite_lon']:
            assert isinstance(self.projection.get(a, Angle(0., 'degrees')),
                              Angle)
        grid_keys = ['LAMzone', 'X_resolution', 'Y_resolution',
                     'input_lon', 'input_lat', 'input_position']
        if set(self.grid.keys()) != set(grid_keys):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys))
        assert isinstance(self.grid['input_lon'], Angle)
        assert isinstance(self.grid['input_lat'], Angle)
        LAMzone_values = [None, 'CI', 'CIE']
        if self.grid['LAMzone'] not in LAMzone_values:
            raise epygramError("grid['LAMzone'] must be among " +
                               str(LAMzone_values))
        dimensions_keys = ['X', 'Y']
        if self.grid['LAMzone'] in ('CI', 'CIE'):
            dimensions_keys.extend(['X_Iwidth', 'Y_Iwidth',
                                    'X_Czone', 'Y_Czone',
                                    'X_CIzone', 'Y_CIzone'])
        if self.grid['LAMzone'] == 'CIE':
            dimensions_keys.extend(['X_CIoffset', 'Y_CIoffset'])
        if set(self.dimensions.keys()) != set(dimensions_keys):
            raise epygramError("dimensions attribute must consist in keys: " +
                               str(dimensions_keys))
        # if self.projection['rotation'].get('degrees') != 0.0:
        #    epylog.warning('*rotation* != 0. may not have been thoroughly tested...')  # TOBECHECKED: here and there, ...

    def select_subzone(self, subzone):
        """
        If a LAMzone defines the geometry, select only the *subzone* from it
        and return a new geometry object.

        :param subzone: among ('C', 'CI').
        """
        assert subzone in ('C', 'CI'), \
               'unknown subzone : ' + subzone
        if self.grid.get('LAMzone') is not None:
            newgeom = self.deepcopy()
            if subzone == 'CI':
                io = self.dimensions.get('X_CIoffset', 0)
                jo = self.dimensions.get('Y_CIoffset', 0)
                centerPoint = (io + (float(self.dimensions['X_CIzone']) - 1) / 2.,
                               jo + (float(self.dimensions['Y_CIzone']) - 1) / 2.)  # Coordinates of center point
                newgeom.grid['LAMzone'] = subzone
                newgeom.dimensions = {'X':self.dimensions['X_CIzone'],
                                      'Y':self.dimensions['Y_CIzone'],
                                      'X_CIzone':self.dimensions['X_CIzone'],
                                      'Y_CIzone':self.dimensions['Y_CIzone'],
                                      'X_Iwidth':self.dimensions['X_Iwidth'],
                                      'Y_Iwidth':self.dimensions['Y_Iwidth'],
                                      'X_Czone':self.dimensions['X_Czone'],
                                      'Y_Czone':self.dimensions['Y_Czone']}
            elif subzone == 'C':
                centerPoint = ((float(self.dimensions['X_Czone']) - 1) / 2.,
                               (float(self.dimensions['Y_Czone']) - 1) / 2.)  # Coordinates of center point
                newgeom.grid['LAMzone'] = None
                newgeom.dimensions = {'X':self.dimensions['X_Czone'],
                                      'Y':self.dimensions['Y_Czone']}
            newgeom.grid['input_lon'] = self._center_lon
            newgeom.grid['input_lat'] = self._center_lat
            newgeom.grid['input_position'] = centerPoint
        else:
            newgeom = self

        return newgeom

    def getcenter(self):
        """
        Returns the coordinate of the grid center as a tuple of Angles
        (center_lon, center_lat).
        """
        return (self._center_lon, self._center_lat)

    def ij2xy(self, i, j, position=None):
        """
        Return the (x, y) coordinates of point *(i,j)*, in the projection.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """
        if isinstance(i, list) or isinstance(i, tuple):
            i = numpy.array(i)
        if isinstance(j, list) or isinstance(j, tuple):
            j = numpy.array(j)
        Xresolution = self.grid['X_resolution']
        Yresolution = self.grid['Y_resolution']
        if self.grid.get('LAMzone', None) is None:
            Xpoints = self.dimensions['X']
            Ypoints = self.dimensions['Y']
        else:
            Xpoints = self.dimensions['X_CIzone']
            Ypoints = self.dimensions['Y_CIzone']
        Xorigin = 0.
        Yorigin = 0.
        # origin of coordinates is the center of CI domain
        i0 = self.dimensions.get('X_CIoffset', 0) + float(Xpoints - 1) / 2.
        j0 = self.dimensions.get('Y_CIoffset', 0) + float(Ypoints - 1) / 2.
        (oi, oj) = self._getoffset(position)
        x = Xorigin + (i - i0 + oi) * Xresolution
        y = Yorigin + (j - j0 + oj) * Yresolution
        return (x, y)

    def xy2ij(self, x, y, position=None):
        """
        Return the (i, j) indexes of point *(x, y)*,
        in the 2D matrix of gridpoints.

        :param x: X coordinate of point in the projection
        :param y: Y coordinate of point in the projection
        :param position: position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float (the nearest grid point is the nearest
        integer).

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)
        Xresolution = self.grid['X_resolution']
        Yresolution = self.grid['Y_resolution']
        if self.grid.get('LAMzone', None) is None:
            Xpoints = self.dimensions['X']
            Ypoints = self.dimensions['Y']
        else:
            Xpoints = self.dimensions['X_CIzone']
            Ypoints = self.dimensions['Y_CIzone']
        Xorigin = 0.
        Yorigin = 0.
        # origin of coordinates is the center of CI domain
        i0 = self.dimensions.get('X_CIoffset', 0) + float(Xpoints - 1) / 2.
        j0 = self.dimensions.get('Y_CIoffset', 0) + float(Ypoints - 1) / 2.
        (oi, oj) = self._getoffset(position)
        i = i0 + (x - Xorigin) / Xresolution - oi
        j = j0 + (y - Yorigin) / Yresolution - oj
        return (i, j)

    def ll2xy(self, lon, lat):
        """
        Return the (x, y) coordinates of point *(lon, lat)* in degrees.

        :param lon: longitude of point in degrees
        :param lat: latitude of point in degrees

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """
        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)
        lon = degrees_nearest_mod(lon, self.projection['reference_lon'].get('degrees'))
        xy = self._proj(lon, lat)
        return self._rotate_axis(*xy, direction='ll2xy')

    def xy2ll(self, x, y):
        """
        Return the (lon, lat) coordinates of point *(x, y)* in the 2D matrix
        of gridpoints*(i,j)*, in degrees.

        :param x: X coordinate of point in the projection
        :param y: Y coordinate of point in the projection

        Note that origin of coordinates in projection is the center of the
        C+I domain.
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)
        (a, b) = self._rotate_axis(x, y, direction='xy2ll')
        ll = self._proj(a, b, inverse=True)
        # mask invalid values
        lons, lats = ll
        lons_is_scalar, lats_is_scalar = is_scalar(lons), is_scalar(lats)
        lons, lats = as_numpy_array(lons), as_numpy_array(lats)
        mask = numpy.logical_or(lons == config.pyproj_invalid_values,
                                lats == config.pyproj_invalid_values)
        mask = numpy.logical_or(mask,
                                numpy.logical_or(lons == numpy.inf,
                                                 lats == numpy.inf))
        lons = numpy.ma.masked_where(mask, lons)
        lats = numpy.ma.masked_where(mask, lats)
        if lons_is_scalar:
            lons = lons[0]
        if lats_is_scalar:
            lats = lats[0]
        ll = (lons, lats)
        return ll

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        return self.xy2ll(*self.ij2xy(i, j, position))

    def ll2ij(self, lon, lat, position=None):
        """
        Return the (i, j) indexes of point *(lon, lat)* in degrees,
        in the 2D matrix of gridpoints.

        :param lon: longitude of point in degrees
        :param lat: latitude of point in degrees
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float.
        """
        return self.xy2ij(*self.ll2xy(lon, lat), position=position)

    def map_factor(self, lat):
        """Returns the map factor at the given latitude(s) *lat* in degrees."""
        lat = numpy.radians(lat)
        if self.name == 'mercator':
            if self.secant_projection:
                lat_0 = self.projection['secant_lat'].get('radians')
            else:
                lat_0 = 0.
            m = numpy.cos(lat_0) / numpy.cos(lat)
        elif self.name == 'polar_stereographic':
            if self.secant_projection:
                lat_0 = self.projection['secant_lat'].get('radians')
            else:
                lat_0 = self.projection['reference_lat'].get('radians')
            m = (1. + numpy.copysign(1., lat_0) * numpy.sin(lat_0)) / \
                (1. + numpy.copysign(1., lat_0) * numpy.sin(lat))
        elif self.name == 'lambert':
            if self.secant_projection:
                lat_0 = self.projection['secant_lat2'].get('radians')
            else:
                lat_0 = self.projection['reference_lat'].get('radians')
            m = (numpy.cos(lat_0) / numpy.cos(lat)) ** (1. - self._K) * \
                ((1. + numpy.copysign(1., lat_0) * numpy.sin(lat_0)) /
                 (1. + numpy.copysign(1., lat_0) * numpy.sin(lat))) ** self._K
        else:
            raise epygramError("projection " + self.name + " not implemented.")
        return m

    def map_factor_field(self, position=None):
        """
        Returns a new field whose data is the map factor over the field.

        :param position: grid position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        kwargs_vcoord = {'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [0]}
        vcoordinate = VGeometry(**kwargs_vcoord)
        geometry = self.deepcopy()
        geometry.vcoordinate=vcoordinate
        f = fpx.field(structure='H2D',
                      geometry=geometry,
                      fid={'geometry':'Map Factor'},
                      units='-')
        lats = self.get_lonlat_grid(position=position)[1]
        data = self.map_factor(lats)
        f.setdata(data)
        return f

    def mesh_area(self, lat):
        """
        Compute the area of a mesh/gridpoint, given its latitude.
        """
        return self.grid['X_resolution'] * self.grid['Y_resolution'] / (self.map_factor(lat)**2)

    def mesh_area_field(self, position=None):
        """
        Returns a new field whose data is the mesh area of gridpoints,
        i.e. X_resolution x Y_resolution / m^2, where m is the local map factor.

        :param position: grid position with respect to the model cell.
                         Defaults to self.position_on_horizontal_grid.
        """
        m = self.map_factor_field(position=position)
        dxdy = self.grid['X_resolution'] * self.grid['Y_resolution']
        m.setdata(dxdy / m.data**2)
        m.fid['geometry'] = 'Mesh Area'
        m.units = 'm^2'
        return m

    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up in the geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        :param num: number of points, including point1 and point2.
        """
        if num < 2:
            raise epygramError("'num' must be at least 2.")
        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        xy_linspace = list(zip(numpy.linspace(x1, x2, num=num),
                               numpy.linspace(y1, y2, num=num)))
        return [tuple(numpy.around(self.xy2ll(*xy), 8)) for xy in xy_linspace]

    def distance(self, end1, end2):
        """
        Computes the distance between two points along a straight line in the
        geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        """
        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        # distance on the projected surface
        distance = numpy.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        # map factor computed as the mean over 100 points along the line
        # between point1 and point2
        mean_map_factor = numpy.array([self.map_factor(lat) for (_, lat) in self.linspace(end1, end2, 10)]).mean()
        return distance / mean_map_factor

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.

        :param lon: longitude of the point in degrees
        :param lat: latitude of the point in degrees
        """
        return self.resolution_ij(*self.ll2ij(lon, lat))

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        """
        (iint, jint) = (numpy.rint(i).astype('int'),
                        numpy.rint(j).astype('int'))
        points_list = [(iint + oi, jint + oj)
                       for oi in [-1, 0, 1]
                       for oj in [-1, 0, 1]
                       if (oi, oj) != (0, 0)]
        return numpy.array([self.distance(self.ij2ll(iint, jint),
                                          self.ij2ll(*p))
                            for p in points_list]).min()

    def plane_azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points in plane local referential
        geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        """
        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        beta = self.projection['rotation'].get('degrees')
        return (numpy.degrees(numpy.arctan2(x2 - x1, y2 - y1)) - beta + 180.) % 360. - 180.

    def make_section_geometry(self, end1, end2,
                              points_number=None,
                              resolution=None,
                              position=None):
        """
        Returns a projected Geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param position: defines the position of data in the grid (defaults to 'center')
        """
        if resolution is not None and points_number is not None:
            raise epygramError("only one of resolution and " +
                               " points_number can be given.")

        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        distance = numpy.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        if resolution is None and points_number is None:
            resolution = 0.5 * (self.resolution_ll(*end1) +
                                self.resolution_ll(*end2))
            points_number = int(numpy.rint(distance / resolution)) + 1
            if resolution > distance:
                raise epygramError("'ends' are too near: pure" +
                                   " interpolation between two gridpoints.")
        if points_number is not None:
            if points_number < 2:
                raise epygramError("'points_number' must be at least 2.")
            resolution = distance / (points_number - 1)
        else:
            points_number = int(numpy.rint(distance / resolution)) + 1
        rotation = numpy.arctan2(y2 - y1, x2 - x1) + self.projection['rotation'].get('radians')
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[])
        grid = {'LAMzone':None,
                'X_resolution':resolution,
                'Y_resolution':resolution,
                'input_lon':Angle(end1[0], 'degrees'),
                'input_lat':Angle(end1[1], 'degrees'),
                'input_position':(0, 0)}
        projection = dict(self.projection)
        projection['rotation'] = Angle(rotation, 'radians')
        dimensions = {'X':points_number, 'Y':1}
        kwargs_geom = dict(name=self.name,
                           grid=grid,
                           dimensions=dimensions,
                           projection=projection,
                           geoid=self.geoid,
                           position_on_horizontal_grid='center' if position is None else position,
                           vcoordinate=vcoordinate)
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return ProjectedGeometry(**kwargs_geom)

    def compass_grid(self, subzone=None, position=None):
        """
        Get the compass grid, i.e. the angle between Y-axis and North for each
        gridpoint.

        :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        :param position: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        (lons, _) = self.get_lonlat_grid(subzone=subzone, position=position)
        if self.name == 'lambert':
            if self.secant_projection:
                k = numpy.copysign(self._K,
                                   self.projection['secant_lat1'].get('degrees') +
                                   self.projection['secant_lat1'].get('degrees'))
            else:
                k = numpy.copysign(self._K, self.projection['reference_lat'].get('degrees'))
        elif self.name == 'mercator':
            k = 0.
        elif self.name == 'polar_stereographic':
            k = self.projection['reference_lat'].get('cos_sin')[1]
        deviation = lons - self.projection['reference_lon'].get('degrees')
        deviation = ((deviation + 180) % 360) - 180  # Between -180 and 180
        theta = k * deviation - self.projection.get('rotation', 0.).get('degrees')
        return theta

    def reproject_wind_on_lonlat(self, u, v,
                                 lon=None, lat=None,
                                 map_factor_correction=True,
                                 reverse=False):
        """
        Reprojects a wind vector (u, v) on the grid onto real
        axes, i.e. with components on true zonal/meridian axes.

        :param u: the u == zonal-on-the-grid component of wind
        :param v: the v == meridian-on-the-grid component of wind
        :param lon: longitudes of points in degrees, if u/v are not vectors
          on the whole grid
        :param lat: latitudes of points in degrees, if u/v are not vectors
          on the whole grid
        :param map_factor_correction:, applies a correction of magnitude due
          to map factor.
        :param reverse: if True, apply the reverse reprojection.
        """
        theta = (-1) ** int(reverse) * self.compass_grid()
        costheta = numpy.cos(-theta * numpy.pi / 180.)
        sintheta = numpy.sin(-theta * numpy.pi / 180.)
        if numpy.shape(u) == costheta.shape:
            # whole grid
            if map_factor_correction:
                m = self.map_factor_field().getdata()
            else:
                m = numpy.ones(costheta.shape)
            if reverse:
                m = 1. / m
            ru = m * (u * costheta - v * sintheta)
            rv = m * (u * sintheta + v * costheta)
        else:
            # some points
            ru = numpy.ndarray(numpy.shape(u))
            rv = numpy.ndarray(numpy.shape(v))
            for k in range(len(u)):
                if map_factor_correction:
                    m = self.map_factor(lat[k])
                else:
                    m = 1.
                if reverse:
                    m = 1. / m
                (i, j) = self.ll2ij(lon[k], lat[k])
                i = int(i)
                j = int(j)
                c = costheta[j, i]
                s = sintheta[j, i]
                ru[k] = m * (u[k] * c - v[k] * s)
                rv[k] = m * (u[k] * s + v[k] * c)

        return (ru, rv)

    def _what_projection(self, out=sys.stdout, arpifs_var_names=False):
        """
        Writes in file a summary of the projection of the field's grid.

        :param out: the output open file-like object
        :param arpifs_var_names: if True, prints the equivalent 'arpifs' variable
          names.
        """
        projection = self.projection
        grid = self.grid
        dimensions = self.dimensions
        varname = ''
        projmap = {'regular_lonlat':'Regular Lon/Lat',
                   'lambert':'Lambert (conformal conic)',
                   'mercator':'Mercator',
                   'polar_stereographic':'Polar Stereographic',
                   'space_view':'Space View'}
        write_formatted(out, "Kind of projection",
                        projmap[self.name])
        if self.name == 'space_view':
            write_formatted(out, "Satellite Longitude in deg",
                            projection['satellite_lon'].get('degrees'))
            write_formatted(out, "Satellite Latitude in deg",
                            projection['satellite_lat'].get('degrees'))
            write_formatted(out, "Reference Longitude in deg",
                            projection['reference_lon'].get('degrees'))
        else:
            if self.secant_projection:
                if self.name == 'lambert':
                    write_formatted(out, "Secant Latitude 1 in deg",
                                    projection['secant_lat1'].get('degrees'))
                    write_formatted(out, "Secant Latitude 2 in deg",
                                    projection['secant_lat2'].get('degrees'))
                else:
                    write_formatted(out, "Secant Latitude in deg",
                                    projection['secant_lat'].get('degrees'))
            else:
                write_formatted(out, "Sinus of Reference Latitude",
                                projection['reference_lat'].get('cos_sin')[1])
                if arpifs_var_names:
                    varname = ' (ELAT0)'
                write_formatted(out, "Reference Latitude in deg" + varname,
                                projection['reference_lat'].get('degrees'))
            if arpifs_var_names:
                varname = ' (ELON0)'
            write_formatted(out, "Reference Longitude in deg" + varname,
                            projection['reference_lon'].get('degrees'))
        write_formatted(out, "Angle of rotation in deg", projection['rotation'].get('degrees'))
        if self.grid.get('LAMzone', False) in ('C', False):
            (lons, lats) = self.get_lonlat_grid()
            corners = self.gimme_corners_ll()
        else:
            (lons, lats) = self.get_lonlat_grid(subzone='CI')
            corners = self.gimme_corners_ll(subzone='CI')
        if arpifs_var_names:
            varname = ' (ELONC)'
        write_formatted(out, "Center Longitude (of C+I) in deg" + varname,
                        self._center_lon.get('degrees'))
        if arpifs_var_names:
            varname = ' (ELATC)'
        write_formatted(out, "Center Latitude (of C+I) in deg" + varname,
                        self._center_lat.get('degrees'))
        if arpifs_var_names:
            varname = ' (EDELX)'
        write_formatted(out, "Resolution in X, in metres" + varname,
                        grid['X_resolution'])
        if arpifs_var_names:
            varname = ' (EDELY)'
        write_formatted(out, "Resolution in Y, in metres" + varname,
                        grid['Y_resolution'])
        if self.grid['LAMzone'] is None:
            dimX = dimensions['X']
            dimY = dimensions['Y']
        else:
            dimX = dimensions['X_CIzone']
            dimY = dimensions['Y_CIzone']
        if arpifs_var_names:
            varname = ' (ELX)'
        write_formatted(out,
                        "Domain width (of C+I) in X, in metres" + varname,
                        grid['X_resolution'] * (dimX - 1))
        if arpifs_var_names:
            varname = ' (ELY)'
        write_formatted(out,
                        "Domain width (of C+I) in Y, in metres" + varname,
                        grid['Y_resolution'] * (dimY - 1))
        write_formatted(out, "Max Longitude (of C+I) in deg", lons.max())
        write_formatted(out, "Min Longitude (of C+I) in deg", lons.min())
        write_formatted(out, "Max Latitude (of C+I) in deg", lats.max())
        write_formatted(out, "Min Latitude (of C+I) in deg", lats.min())
        write_formatted(out,
                        "Low-Left corner (of C+I) Longitude in deg",
                        corners['ll'][0])
        write_formatted(out,
                        "Low-Left corner (of C+I) Latitude in deg",
                        corners['ll'][1])
        write_formatted(out,
                        "Low-Right corner (of C+I) Longitude in deg",
                        corners['lr'][0])
        write_formatted(out,
                        "Low-Right corner (of C+I) Latitude in deg",
                        corners['lr'][1])
        write_formatted(out,
                        "Upper-Left corner (of C+I) Longitude in deg",
                        corners['ul'][0])
        write_formatted(out,
                        "Upper-Left corner (of C+I) Latitude in deg",
                        corners['ul'][1])
        write_formatted(out,
                        "Upper-Right corner (of C+I) Longitude in deg",
                        corners['ur'][0])
        write_formatted(out,
                        "Upper-Right corner (of C+I) Latitude in deg",
                        corners['ur'][1])

    # FIXME: cleanme def __eq__(self, other):
    #    """Test of equality by recursion on the object's attributes."""
    """    if self.__class__ == other.__class__ and \
           set(self._attributes.keys()) == set(other._attributes.keys()):
            for attr in self._attributes.keys():
                if attr == 'grid':
                    # the same grid could be defined from different inputs
                    selfgrid = {k:v for k, v in self.grid.items() if
                                k not in ['input_lon',
                                          'input_lat',
                                          'input_position']}
                    othergrid = {k:v for k, v in self.grid.items() if
                                 k not in ['input_lon',
                                           'input_lat',
                                           'input_position']}
                    ok = (selfgrid == othergrid and
                          self.getcenter() == other.getcenter())
                else:
                    ok = self._attributes[attr] == other._attributes[attr]
                if not ok:
                    break
        else:
            ok = False
        return ok

    def __hash__(self):
        # known issue __eq__/__hash__ must be defined both or none, else inheritance is broken
        return super(ProjectedGeometry, self).__hash__()"""
