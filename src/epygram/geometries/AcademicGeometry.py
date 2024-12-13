#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for 3D geometries of fields.
"""

import numpy
import copy
import sys

import footprints

from epygram import epygramError, config
from epygram.util import Angle, write_formatted, Comparator

from .VGeometry import VGeometry
from .AbstractGeometry import RectangularGridGeometry

epylog = footprints.loggers.getLogger(__name__)


class AcademicGeometry(RectangularGridGeometry):
    """Handles the geometry for an academic 3-Dimensions Field."""

    def __init__(self, name, grid, dimensions, vcoordinate, projection,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be 'academic'
        :param grid: Handles description of the horizontal grid.
        :param dimensions: Handles grid dimensions.
        :param vcoordinate: Handles vertical geometry parameters.
        :param position_on_horizontal_grid: Position of points w/r to the horizontal.
                                            among: ['upper-right', 'upper-left',
                                                    'lower-left', 'lower-right',
                                                    'center-left', 'center-right',
                                                    'lower-center', 'upper-center',
                                                    'center', '__unknown__']
        :param geoid: To specify geoid shape (of no meaning in this geometry).
        """
        self.add_attr_inlist('name', ['academic'])
        self.add_attr_dict('projection')

        self.name = name
        self.projection = projection

        super(AcademicGeometry, self).__init__(grid, dimensions, vcoordinate,
                                                      position_on_horizontal_grid, geoid)

        if self.grid['input_position'] != (0, 0):
            raise NotImplementedError("For now, only input_position = (0, 0) is allowed for academic geometries.")

    @property
    def _center_lon(self):
        """Get the center's longitude"""
        return (self.dimensions['X'] - 1) / 2.

    @property
    def _center_lat(self):
        """Get the center's latitude"""
        return (self.dimensions['Y'] - 1) / 2.

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
                    almost['_attributes']['grid'].pop(k)
            almost_self['_attributes']['grid']['center'] = self.getcenter()
            almost_other['_attributes']['grid']['center'] = other.getcenter()
            return Comparator.are_equal(almost_self, almost_other, tolerance)
        else:
            return False

    def _rotate_axis(self, x, y, direction):
        """
        Internal method used to rotate x/y coordinates to handle rotated geometries.
        *direction*, if evals to 'xy2ll', direction used when converting (x, y) to (lat, lon).
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
        grid_keys = ['LAMzone', 'X_resolution', 'Y_resolution',
                     'input_lat', 'input_lon', 'input_position']
        if set(self.grid.keys()) != set(grid_keys) and \
           set(self.grid.keys()) != set(grid_keys + ['longitude', 'latitude']):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys) + " or " +
                               str(grid_keys + ['longitude', 'latitude']))
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
        projection_keys = ['reference_dX', 'reference_dY', 'rotation']
        if set(self.projection.keys()) != set(projection_keys):
            raise epygramError("projection attribute must consist in keys: " +
                               str(projection_keys))

    def _getoffset(self, position=None):
        """
        Returns the offset to use for this position.
        Replaces the method defined in RectangularGridGeometry to deal with
        1D or 2D simulations.
        """
        offset = super(AcademicGeometry, self)._getoffset(position)
        if self.dimensions['X'] == 1:
            offset = (0, offset[1])
        if self.dimensions['Y'] == 1:
            offset = (offset[0], 0)
        return offset

    def getcenter(self):
        """
        Returns the coordinate of the grid center as a tuple
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

        Note that origin of coordinates in projection is the center of the C+I domain.
        """
        if isinstance(i, list) or isinstance(i, tuple):
            i = numpy.array(i)
        if isinstance(j, list) or isinstance(j, tuple):
            j = numpy.array(j)
        (oi, oj) = self._getoffset(position)
        return ((i + oi) * self.grid['X_resolution'],
                (j + oj) * self.grid['Y_resolution'])

    def xy2ij(self, x, y, position=None):
        """
        Return the (i, j) indexes of point *(x, y)*,
        in the 2D matrix of gridpoints.

        :param x: X coordinate of point in the academic projection
        :param y: Y coordinate of point in the academic projection
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
        (oi, oj) = self._getoffset(position)
        return (x / self.grid['X_resolution'] - oi,
                y / self.grid['Y_resolution'] - oj)

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        This routine has a special meaning for this geometry. It returns the (i, j)
        position in the original grid (especially after a section extraction) with
        an offset of one (for old tools compatibility).
        """
        if isinstance(i, list) or isinstance(i, tuple):
            i = numpy.array(i)
        if isinstance(j, list) or isinstance(j, tuple):
            j = numpy.array(j)
        return self.xy2ll(*self.ij2xy(i, j, position))

    def ll2ij(self, lon, lat, position=None):
        """
        Return the (i, j) indexes of point *(lon, lat)* in degrees,
        in the 2D matrix of gridpoints.

        :param lon: longitude of point in degrees
        :param lat: latitude of point in degrees
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: the returned (*i,j*) are float.

        This routine has a special meaning for this geometry. It returns the (lon, lat)
        position in the original grid (especially after a section extraction) with
        an offset of one (for old tools compatibility).
        """
        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)
        return self.xy2ij(*self.ll2xy(lon, lat), position=position)

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
        # TOBECHECKED: position of self.ij2xy(0, 0) inside the formula
        xy = ((lon - self.grid['input_lon']) * self.projection['reference_dX'] + self.ij2xy(0, 0, 'center')[0],
              (lat - self.grid['input_lat']) * self.projection['reference_dY'] + self.ij2xy(0, 0, 'center')[1])
        return self._rotate_axis(*xy, direction='ll2xy')

    def xy2ll(self, x, y):
        """
        Return the (lon, lat) coordinates of point *(x, y)* in the 2D matrix
        of gridpoints*(i,j)*, in degrees.

        :param x: X coordinate of point in the academic projection
        :param y: Y coordinate of point in the academic projection

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)

        (a, b) = self._rotate_axis(x, y, direction='xy2ll')
        # TOBECHECKED: position of self.ij2xy(0, 0) inside the formula
        return ((a - self.ij2xy(0, 0, 'center')[0]) / self.projection['reference_dX'] + self.grid['input_lon'],
                (b - self.ij2xy(0, 0, 'center')[1]) / self.projection['reference_dY'] + self.grid['input_lat'])

    def distance(self, end1, end2):
        """
        Computes the distance between two points along a straight line in the
        geometry.
        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        """
        return numpy.sqrt(((end1[0] - end2[0]) * self.grid['X_resolution']) ** 2 +
                          ((end1[1] - end2[1]) * self.grid['Y_resolution']) ** 2)

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
        return list(zip(numpy.linspace(end1[0], end2[0], num=num),
                        numpy.linspace(end1[1], end2[1], num=num)))

    def resolution_ll(self, *_, **__):
        """Returns the minimum of X and Y resolution."""
        return min(self.grid['X_resolution'], self.grid['Y_resolution'])

    def resolution_ij(self, *_, **__):
        """Returns the minimum of X and Y resolution."""
        return min(self.grid['X_resolution'], self.grid['Y_resolution'])

    def plane_azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points in plane local referential
        geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        """
        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        return (numpy.degrees(numpy.arctan2(x2 - x1, y2 - y1)) + 180.) % 360. - 180.

    def azimuth(self, end1, end2):
        """Same as plane_azimuth in this geometry."""
        return self.plane_azimuth(end1, end2)

    def _what_position(self, out=sys.stdout):
        """
        Writes in file a summary of the position of the field.

        :param out: the output open file-like object
        """
        write_formatted(out, "Kind of Geometry", 'Academic')
        if 'latitude' in self.grid:
            write_formatted(out, "Latitude", self.grid['latitude'])
        if 'longitude' in self.grid:
            write_formatted(out, "Longitude", self.grid['longitude'])

    def _what_projection(self, out=sys.stdout, **_):
        """
        Writes in file a summary of the projection of the field's grid.

        :param out: the output open file-like object
        """
        projection = self.projection
        write_formatted(out, "X resolution of the reference grid",
                        projection['reference_dX'])
        write_formatted(out, "Y resolution of the reference grid",
                        projection['reference_dY'])

    def make_section_geometry(self, end1, end2,
                              points_number=None,
                              resolution=None,
                              position=None):
        """
        Returns a academic Geometry.

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
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[])
        grid = {'LAMzone':'CIE',
                'X_resolution':resolution,
                'Y_resolution':resolution,
                'input_lat':end1[1],
                'input_lon':end1[0],
                'input_position':(0, 0)}
        if 'longitude' in self.grid:
            grid['longitude'] = self.grid['longitude']
        if 'latitude' in self.grid:
            grid['latitude'] = self.grid['latitude']

        dimensions = {'X':points_number, 'Y':1,
                      'X_Iwidth':0, 'Y_Iwidth':0,
                      'X_Czone':0, 'Y_Czone':0,
                      'X_CIzone':points_number, 'Y_CIzone':1,
                      'X_CIoffset':0, 'Y_CIoffset':0}

        rotation = numpy.arctan2(y2 - y1, x2 - x1)
        projection = {'rotation':Angle(rotation, 'radians'),
                      'reference_dX':self.projection['reference_dX'],
                      'reference_dY':self.projection['reference_dY']}

        kwargs_geom = dict(name=self.name,
                           grid=grid,
                           projection=projection,
                           dimensions=dimensions,
                           position_on_horizontal_grid='center' if position is None else position,
                           vcoordinate=vcoordinate
                           )
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return AcademicGeometry(**kwargs_geom)
