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

from epygram import epygramError, config
from epygram.config import rounding_decimal as _rd
from epygram.util import degrees_nearest_mod, Angle, write_formatted

from .AbstractGeometry import LLGeometry

epylog = footprints.loggers.getLogger(__name__)



class RegLLGeometry(LLGeometry):
    """
    Handles the geometry for a Regular Lon/Lat 3-Dimensions Field.

    TODO: Is this class really necessary. Maybe we could make only one class
          with LLGeometry, RegLLGeometry and RotLLGEometry?
          Is a RegLLGeometry equivalent to a RotLLGeometry with a null rotation angle?
    """

    def __init__(self, name, grid, dimensions, vcoordinate,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be 'regular_lonlat'
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
        self.add_attr_inlist('name', ['regular_lonlat'])

        self.name = name

        super(RegLLGeometry, self).__init__(grid, dimensions, vcoordinate,
                                            position_on_horizontal_grid, geoid)

        # earth-round grids: wrap TODO:
        corners = self.gimme_corners_ll()
        if abs(abs(degrees_nearest_mod(corners['ul'][0], corners['ur'][0]) -
                   corners['ur'][0]) -
               self.grid['X_resolution'].get('degrees')) <= config.epsilon:
            self._earthround = True
        else:
            self._earthround = False

    @property
    def _center_lon(self):
        """Get the center's longitude"""
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
    def _center_lat(self):
        """Get the center's latitude"""
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

    def getcenter(self):
        """
        Returns the coordinate of the grid center as a tuple of Angles
        (center_lon, center_lat).
        """
        return (self._center_lon, self._center_lat)

    def _consistency_check(self):
        """Check that the geometry is consistent."""
        if self.name == 'regular_lonlat':
            grid_keys = ['input_lon', 'input_lat', 'input_position',
                         'X_resolution', 'Y_resolution']
            if set(self.grid.keys()) != set(grid_keys):
                raise epygramError("grid attribute must consist in keys: " +
                                   str(grid_keys))
            assert isinstance(self.grid['input_lon'], Angle)
            assert isinstance(self.grid['input_lat'], Angle)
            assert isinstance(self.grid['X_resolution'], Angle)
            assert isinstance(self.grid['Y_resolution'], Angle)
            dimensions_keys = ['X', 'Y']
            if set(self.dimensions.keys()) != set(dimensions_keys):
                raise epygramError("dimensions attribute must consist in keys: " + str(dimensions_keys))

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

        (oi, oj) = self._getoffset(position)
        Xpoints = self.dimensions['X']
        Ypoints = self.dimensions['Y']
        Xresolution = self.grid['X_resolution'].get('degrees')
        Yresolution = self.grid['Y_resolution'].get('degrees')
        Xorigin = self._center_lon.get('degrees')
        Yorigin = self._center_lat.get('degrees')
        # origin of coordinates is the center of domain
        i0 = float(Xpoints - 1) / 2.
        j0 = float(Ypoints - 1) / 2.
        x = numpy.round(Xorigin + (i - i0 + oi) * Xresolution, _rd)
        y = numpy.round(Yorigin + (j - j0 + oj) * Yresolution, _rd)
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

        (oi, oj) = self._getoffset(position)
        Xpoints = self.dimensions['X']
        Ypoints = self.dimensions['Y']
        Xresolution = self.grid['X_resolution'].get('degrees')
        Yresolution = self.grid['Y_resolution'].get('degrees')
        Xorigin = self._center_lon.get('degrees')
        Yorigin = self._center_lat.get('degrees')
        # origin of coordinates is the center of domain
        i0 = float(Xpoints - 1) / 2.
        j0 = float(Ypoints - 1) / 2.
        i = i0 + (x - Xorigin) / Xresolution - oi
        j = j0 + (y - Yorigin) / Yresolution - oj
        if self._earthround and isinstance(i, numpy.ndarray):
            for idx, ii in enumerate(i):
                if ii < -0.5:
                    ii = self.dimensions['X'] - 1 + (1 - ii)
                elif ii > (self.dimensions['X'] - 1) + 0.5:
                    ii = ii - (self.dimensions['X'] - 1) - 1
                i[idx] = ii
        return (i, j)

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        return self.ij2xy(i, j, position)

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
        lon = degrees_nearest_mod(lon, self._center_lon.get('degrees'))
        return self.xy2ij(lon, lat, position)

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
        return self.ij2xy(*self.ll2ij(lon, lat))

    def xy2ll(self, x, y):
        """
        Return the (lon, lat) coordinates of point *(x, y)* in the 2D matrix
        of gridpoints*(i,j)*, in degrees.

        :param x: X coordinate of point in the projection
        :param y: Y coordinate of point in the projection

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """
        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)
        return self.ij2ll(*self.xy2ij(x, y))

    def _what_grid(self, out=sys.stdout, arpifs_var_names=False):
        """
        Writes in file a summary of the grid of the field.

        :param out: the output open file-like object
        :param arpifs_var_names: if True, prints the equivalent 'arpifs' variable
          names.
        """
        grid = self.grid
        dimensions = self.dimensions
        varname = ''
        (lons, lats) = self.get_lonlat_grid()
        corners = self.gimme_corners_ll()
        write_formatted(out, "Kind of Geometry", 'Regular Lon/Lat')
        if arpifs_var_names:
            varname = ' (ELONC)'
        write_formatted(out, "Center Longitude in deg" + varname,
                        self._center_lon.get('degrees'))
        if arpifs_var_names:
            varname = ' (ELATC)'
        write_formatted(out, "Center Latitude in deg" + varname,
                        self._center_lat.get('degrees'))
        if arpifs_var_names:
            varname = ' (EDELX)'
        write_formatted(out, "Resolution in X, in deg" + varname,
                        grid['X_resolution'].get('degrees'))
        if arpifs_var_names:
            varname = ' (EDELY)'
        write_formatted(out, "Resolution in Y, in deg" + varname,
                        grid['Y_resolution'].get('degrees'))
        if arpifs_var_names:
            varname = ' (ELX)'
        write_formatted(out, "Domain width in X, in deg" + varname,
                        grid['X_resolution'].get('degrees') *
                        (dimensions['X'] - 1))
        if arpifs_var_names:
            varname = ' (ELY)'
        write_formatted(out, "Domain width in Y, in deg" + varname,
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
