#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for 3D geometries of fields.
"""

import numpy
import sys

import footprints

from epygram import epygramError
from epygram.config import rounding_decimal as _rd
from epygram.util import degrees_nearest_mod, Angle, write_formatted

from .AbstractGeometry import LLGeometry

epylog = footprints.loggers.getLogger(__name__)



class RotLLGeometry(LLGeometry):
    """
    Handles the geometry for a Rotated Lon/Lat 3-Dimensions Field.

    TODO: Is this class really necessary. Maybe we could make only one class
          with LLGeometry, RegLLGeometry and RotLLGEometry?
          Is a RegLLGeometry equivalent to a RotLLGeometry with a null rotation angle?
    """

    def __init__(self, name, grid, dimensions, vcoordinate,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be 'rotated_lonlat'
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
        self.add_attr_inlist('name', ['rotated_lonlat'])

        self.name = name

        super(RotLLGeometry, self).__init__(grid, dimensions, vcoordinate,
                                            position_on_horizontal_grid, geoid)

        if self.grid.get('rotation', Angle(0., 'degrees')).get('degrees') != 0.:
            raise NotImplementedError("rotation != Angle(0.)")

    @property
    def _center_rlon(self):
        """Get the rotation center's longitude"""
        if self.grid['input_position'] == ((float(self.dimensions['X']) - 1) / 2.,
                                           (float(self.dimensions['Y']) - 1) / 2.):
            return self.grid['input_lon']
        elif self.grid['input_position'][0] == 0:
            return Angle(round(self.grid['input_lon'].get('degrees') +
                               self.grid['X_resolution'].get('degrees') *
                               (self.dimensions['X'] - 1) / 2, _rd),
                         'degrees')
        elif self.grid['input_position'][0] == self.dimensions['X'] - 1:
            return Angle(round(self.grid['input_lon'].get('degrees') -
                               self.grid['X_resolution'].get('degrees') *
                               (self.dimensions['X'] - 1) / 2, _rd),
                         'degrees')
        else:
            raise NotImplementedError("this 'input_position': " +
                                      str(self.grid['input_position']))

    @property
    def _center_rlat(self):
        """Get the rotation center's latitude"""
        if self.grid['input_position'] == ((float(self.dimensions['X']) - 1) / 2.,
                                           (float(self.dimensions['Y']) - 1) / 2.):
            return self.grid['input_lat']
        elif self.grid['input_position'][1] == 0:
            return Angle(round(self.grid['input_lat'].get('degrees') +
                               self.grid['Y_resolution'].get('degrees') *
                               (self.dimensions['Y'] - 1) / 2, _rd),
                         'degrees')
        elif self.grid['input_position'][1] == self.dimensions['Y'] - 1:
            return Angle(round(self.grid['input_lat'].get('degrees') -
                               self.grid['Y_resolution'].get('degrees') *
                               (self.dimensions['Y'] - 1) / 2, _rd),
                         'degrees')
        else:
            raise NotImplementedError("this 'input_position': " +
                                      str(self.grid['input_position']))

    @property
    def _center_lon(self):
        """Get the center's longitude"""
        lon, lat = self.xy2ll(self._center_rlon.get('degrees'),
                              self._center_rlat.get('degrees'))
        return Angle(lon, 'degrees')

    @property
    def _center_lat(self):
        """Get the center's latitude"""
        lon, lat = self.xy2ll(self._center_rlon.get('degrees'),
                              self._center_rlat.get('degrees'))
        return Angle(lat, 'degrees')

    def getcenter(self):
        """
        Returns the coordinate of the grid center as a tuple of Angles
        (center_lon, center_lat) in true lon/lat coordinates.
        """
        return (self._center_lon, self._center_lat)

    def _consistency_check(self):
        """
        Check that the geometry is consistent.

        Note:
        **input_lon** and **input_lat** parameters are supposed to be
        longitude/latitude of input point in the Rotated Lon/lat referential.
        """
        grid_keys = ['input_lon', 'input_lat', 'input_position',
                     'X_resolution', 'Y_resolution',
                     'southern_pole_lon', 'southern_pole_lat',
                     'rotation']
        if set(self.grid.keys()) != set(grid_keys):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys))
        assert isinstance(self.grid['input_lon'], Angle)
        assert isinstance(self.grid['input_lat'], Angle)
        assert isinstance(self.grid['X_resolution'], Angle)
        assert isinstance(self.grid['Y_resolution'], Angle)
        assert isinstance(self.grid['southern_pole_lon'], Angle)
        assert isinstance(self.grid['southern_pole_lat'], Angle)
        assert isinstance(self.grid['rotation'], Angle)
        dimensions_keys = ['X', 'Y']
        if set(self.dimensions.keys()) != set(dimensions_keys):
            raise epygramError("dimensions attribute must consist in keys: " + str(dimensions_keys))

    def ij2xy(self, i, j, position=None):
        """
        Return the (x, y) == (rot'lon, rot'lat) coordinates of point *(i,j)*,
        in the rotated coordinates system.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        if isinstance(i, list) or isinstance(i, tuple):
            i = numpy.array(i)
        if isinstance(j, list) or isinstance(j, tuple):
            j = numpy.array(j)
        (oi, oj) = self._getoffset(position)
        Xpoints = self.dimensions['X']
        Ypoints = self.dimensions['Y']
        Xresolution = self.grid['X_resolution'].get('degrees')
        Yresolution = self.grid['Y_resolution'].get('degrees')
        Xorigin = self._center_rlon.get('degrees')
        Yorigin = self._center_rlat.get('degrees')
        # origin of coordinates is the center of domain
        i0 = float(Xpoints - 1) / 2.
        j0 = float(Ypoints - 1) / 2.
        x = Xorigin + (i - i0 + oi) * Xresolution
        y = Yorigin + (j - j0 + oj) * Yresolution
        return (x, y)

    def xy2ij(self, x, y, position=None):
        """
        Return the (i, j) indexes of point *(x, y)* == (rot'lon, rot'lat),
        in the 2D matrix of gridpoints.

        :param x: rotated lon coordinate of point
        :param y: rotated lat coordinate of point
        :param position: position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float (the nearest grid point is the nearest
        integer).
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)

        (oi, oj) = self._getoffset(position)
        Xpoints = self.dimensions['X']
        Ypoints = self.dimensions['Y']
        Xresolution = self.grid['X_resolution'].get('degrees')
        Yresolution = self.grid['Y_resolution'].get('degrees')
        Xorigin = self._center_rlon.get('degrees')
        Yorigin = self._center_rlat.get('degrees')
        # origin of coordinates is the center of domain
        i0 = float(Xpoints - 1) / 2.
        j0 = float(Ypoints - 1) / 2.
        i = i0 + (x - Xorigin) / Xresolution - oi
        j = j0 + (y - Yorigin) / Yresolution - oj
        return (i, j)

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) true coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        return self.xy2ll(*self.ij2xy(i, j, position=position))

    def ll2ij(self, lon, lat, position=None):
        """
        Return the (i, j) indexes of point *(lon, lat)* in degrees,
        in the 2D matrix of gridpoints.

        :param lon: true longitude of point in degrees
        :param lat: true latitude of point in degrees
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float.
        """
        lon = degrees_nearest_mod(lon, self._center_lon.get('degrees'))
        return self.xy2ij(*self.ll2xy(lon, lat), position=position)

    def ll2xy(self, lon, lat):
        """
        Return the (x, y) == (rot'lon, rot'lat) coordinates of
        point *(lon, lat)* in degrees, in the rotated system.

        :param lon: true longitude of point in degrees
        :param lat: true latitude of point in degrees
        """
        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)
        return self._lonlat_to_rotatedlonlat(lon, lat, reverse=False)

    def xy2ll(self, x, y):
        """
        Return the (lon, lat) coordinates of point *(x, y)* == (rot'lon, rot'lat)
        in the rotated system, in degrees.

        :param x: rotated lon coordinate of point in the projection
        :param y: rotated lat coordinate of point in the projection
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)
        return self._lonlat_to_rotatedlonlat(x, y, reverse=True)

    def _lonlat_to_rotatedlonlat(self, lon, lat, reverse=False):
        """
        Conversion formula from true lon/lat to rotated lon/lat, or **reverse**.
        Inspired from https://gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
        """
        lon = numpy.radians(lon)
        lat = numpy.radians(lat)

        theta = numpy.pi / 2. + self.grid['southern_pole_lat'].get('radians')
        phi = self.grid['southern_pole_lon'].get('radians')
        # spherical to cartesian
        x = numpy.cos(lon) * numpy.cos(lat)
        y = numpy.sin(lon) * numpy.cos(lat)
        z = numpy.sin(lat)
        # conversion
        if not reverse:
            x_new = (numpy.cos(theta) * numpy.cos(phi) * x +
                     numpy.cos(theta) * numpy.sin(phi) * y +
                     numpy.sin(theta) * z)
            y_new = -numpy.sin(phi) * x + numpy.cos(phi) * y
            z_new = (-numpy.sin(theta) * numpy.cos(phi) * x -
                     numpy.sin(theta) * numpy.sin(phi) * y +
                     numpy.cos(theta) * z)
        else:
            phi = -phi
            theta = -theta
            x_new = (numpy.cos(theta) * numpy.cos(phi) * x +
                     numpy.sin(phi) * y +
                     numpy.sin(theta) * numpy.cos(phi) * z)
            y_new = (-numpy.cos(theta) * numpy.sin(phi) * x +
                     numpy.cos(phi) * y -
                     numpy.sin(theta) * numpy.sin(phi) * z)
            z_new = -numpy.sin(theta) * x + numpy.cos(theta) * z
        # cartesian back to spherical coordinates
        lon_new = numpy.degrees(numpy.arctan2(y_new, x_new))
        lat_new = numpy.degrees(numpy.arcsin(z_new))
        return (lon_new, lat_new)

    # TODO: def default_cartopy_CRS(self):

    def _what_grid(self, out=sys.stdout):
        """
        Writes in file a summary of the grid of the field.

        :param out: the output open file-like object
        """
        grid = self.grid
        dimensions = self.dimensions
        (lons, lats) = self.get_lonlat_grid()
        corners = self.gimme_corners_ll()
        write_formatted(out, "Kind of Geometry", 'Rotated Lon/Lat')
        write_formatted(out, "Southern Pole Longitude in deg",
                        grid['southern_pole_lon'].get('degrees'))
        write_formatted(out, "Southern Pole Latitude in deg",
                        grid['southern_pole_lat'].get('degrees'))
        write_formatted(out, "Center Longitude (in rotated referential) in deg",
                        self._center_rlon.get('degrees'))
        write_formatted(out, "Center Latitude (in rotated referential) in deg",
                        self._center_rlat.get('degrees'))
        write_formatted(out, "Center Longitude (true) in deg",
                        self._center_lon.get('degrees'))
        write_formatted(out, "Center Latitude (true) in deg",
                        self._center_lat.get('degrees'))
        write_formatted(out, "Resolution in X, in deg",
                        grid['X_resolution'].get('degrees'))
        write_formatted(out, "Resolution in Y, in deg",
                        grid['Y_resolution'].get('degrees'))
        write_formatted(out, "Domain width in X, in deg",
                        grid['X_resolution'].get('degrees') *
                        (dimensions['X'] - 1))
        write_formatted(out, "Domain width in Y, in deg",
                        grid['Y_resolution'].get('degrees') *
                        (dimensions['Y'] - 1))
        write_formatted(out, "Max Longitude in deg", lons.max())
        write_formatted(out, "Min Longitude in deg", lons.min())
        write_formatted(out, "Max Latitude in deg", lats.max())
        write_formatted(out, "Min Latitude in deg", lats.min())
        write_formatted(out, "Low-Left corner Longitude in deg",
                        corners['ll'][0])
        write_formatted(out, "Low-Left corner Latitude in deg",
                        corners['ll'][1])
        write_formatted(out, "Low-Right corner Longitude in deg",
                        corners['lr'][0])
        write_formatted(out, "Low-Right corner Latitude in deg",
                        corners['lr'][1])
        write_formatted(out, "Upper-Left corner Longitude in deg",
                        corners['ul'][0])
        write_formatted(out, "Upper-Left corner Latitude in deg",
                        corners['ul'][1])
        write_formatted(out, "Upper-Right corner Longitude in deg",
                        corners['ur'][0])
        write_formatted(out, "Upper-Right corner Latitude in deg",
                        corners['ur'][1])
