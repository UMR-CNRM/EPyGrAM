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
from bronx.syntax.arrays import stretch_array

from epygram import epygramError
from epygram.util import (Angle,
                          positive_longitudes, longitudes_between_minus180_180,
                          write_formatted,
                          as_numpy_array, moveaxis)

from .AbstractGeometry import RectangularGridGeometry

epylog = footprints.loggers.getLogger(__name__)


class UnstructuredGeometry(RectangularGridGeometry):
    """Handles the geometry for an unstructured 3-Dimensions Field."""

    def __init__(self, name, grid, dimensions, vcoordinate,
                 position_on_horizontal_grid='center', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be 'unstructured'
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
        self.add_attr_inlist('name', ['unstructured',
                                      'DDH:point', 'DDH:ij_point', 'DDH:quadrilateral',
                                      'DDH:rectangle', 'DDH:globe', 'DDH:zonal_bands'])

        self.name = name

        super(UnstructuredGeometry, self).__init__(grid, dimensions, vcoordinate,
                                                   position_on_horizontal_grid, geoid)

    def _consistency_check(self):
        """Check that the geometry is consistent."""
        grid_keys = ['latitudes', 'longitudes']
        if set(self.grid.keys()) != set(grid_keys) and \
           set(self.grid.keys()) != set(['DDH_domain']):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys) + " or " +
                               str(['DDH_domain']))

    def get_lonlat_grid(self,
                        subzone=None,
                        position=None,
                        d4=False,
                        nb_validities=0,
                        force_longitudes=None):
        """
        Returns a tuple of two tables containing one the longitude of each
        point, the other the latitude, with 2D shape.

        :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        :param position: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with respect to geometry.
                       d4=True requires nb_validities > 0
        :param nb_validities: number of validities represented in data values
        :param force_longitudes: if 'positive', the longitudes will be forced positive
                                 if ']-180,180]', the longitudes will be in the ]-180, 180] interval

        Shape of 2D data on Rectangular grids: \n
          - grid[0,0] is SW, grid[-1,-1] is NE \n
          - grid[0,-1] is SE, grid[-1,0] is NW
        """
        if self._getoffset(position) != (0., 0.):
            raise epygramError('We can only retrieve latitude and longitude of mass point on an unstructured grid')
        lons = self.grid['longitudes']
        if not isinstance(lons, numpy.ndarray):
            lons = numpy.array(lons)
        lats = self.grid['latitudes']
        if not isinstance(self.grid['latitudes'], numpy.ndarray):
            lats = numpy.array(lats)
        if len(lons.shape) == 1:
            lons = self.reshape_data(lons)
            lats = self.reshape_data(lats)
        if subzone and self.grid.get('LAMzone', None):
            lons = self.extract_subzone(lons, subzone)
            lats = self.extract_subzone(lats, subzone)

        if d4:
            lons, lats = self._reshape_lonlat_4d(lons, lats, nb_validities)
        elif not d4 and nb_validities != 0:
            raise ValueError("*nb_validities* must be 0 when d4==False")
        else:
            lons = lons.squeeze()
            lats = lats.squeeze()
        if force_longitudes == 'positive':
            lons = positive_longitudes(lons)
        elif force_longitudes == ']-180,180]':
            lons = longitudes_between_minus180_180(lons)
        return (lons, lats)

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        (lons, lats) = self.get_lonlat_grid(position=position)
        if lons.ndim == 1:
            assert numpy.all(j == 0), "j must be 0 when geometry does not contain j dimension"
            return (lons[i], lats[i])
        else:
            return (lons[j, i], lats[j, i])

    def ll2ij(self, lon, lat, position=None):
        """
        Return the (i, j) indexes of point *(lon, lat)* in degrees,
        in the 2D matrix of gridpoints.

        :param lon: longitude of point in degrees
        :param lat: latitude of point in degrees
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: the returned (*i,j*) are float.
        """
        lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()

        i = numpy.zeros(len(lon))
        j = numpy.zeros(len(lon))
        (lons, lats) = self.get_lonlat_grid(position=position)

        for n in range(len(lon)):
            where = (numpy.abs(lons - lon[n]) + numpy.abs(lats - lat[n]) == 0).nonzero()
            if len(where[0]) == 0:
                raise epygramError("No point found with these coordinates.")
            elif len(where[0]) > 1:
                raise epygramError("Several points have the same coordinates.")
            if self.datashape['j']:
                i[n], j[n] = where[::-1]
            else:
                i[n] = where[0]
        return (i.squeeze(), j.squeeze())

    def nearest_points(self, lon, lat, request,
                       position=None,
                       external_distance=None,
                       squeeze=True):
        """
        Returns the (i, j) position of the points needed to perform an
        interpolation. This is a list of (lon, lat) tuples.

        :param lon: longitude of point in degrees.
        :param lat: latitude of point in degrees.
        :param request: criteria for selecting the points, among:
               * {'n':'1'} - the nearest point
        :param position: position in the model cell of the lat lon position.
          Defaults to self.position_on_horizontal_grid.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':a_3DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = |target_value - external_field.data|
        :param squeeze: True to suppress useless dimensions

        :rtype: general output form is [list, list, ..., list]
                with as many list items as the length of lon/lat.
                Each list item is of the form [tuple, tuple, ..., tuple]
                with as many tuples as the request implies. A tuple
                represents one of the nearest points associated with one
                value taken from lon/lat. Each tuple as the form
                (i, j).

                Dimensions with a length of one are removed except if
                squeeze is False. If squeeze is True and if request
                implies only one nearest point, the list item of the general
                output form is replaced by the tuple item; if length of
                lon/lat is one, the output is directly the list item of
                the general output form. Hence, if length of lon/lat is
                one and the request implies only one point, the output is
                a tuple.

                In case of a simple square request, output is actually
                an array. Otherwise, the output is as described (it cannot
                be an array because the number of nearest points can vary
                with the entry point).
        """
        if request != {'n':'1'}:
            raise NotImplementedError("*request* != '{'n':'1'}' for UnstructuredGeometry.")
        if external_distance is not None:
            raise NotImplementedError("*external_distance* is not None for UnstructuredGeometry.")

        lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()

        (lons, lats) = self.get_lonlat_grid(position=position)

        result = []
        for (one_lon, one_lat) in zip(lon, lat):
            dist = (one_lon - lons) ** 2 + (one_lat - lats) ** 2
            i = dist.argmin()
            if lons.ndim == 2:
                result.append(numpy.unravel_index(i, lons.shape)[::-1])
            else:
                result.append((i, 0))
        result = numpy.array(result)
        if squeeze:
            result = result.squeeze()
        return result

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.

        :param lon: longitude of point in degrees.
        :param lat: latitude of point in degrees.
        """
        near_point = moveaxis(self.nearest_points(lon, lat, request={'n':'1'}), 0, -1)
        return self.resolution_ij(*near_point)

    def resolution_ij(self, i, j, position=None):
        """
        Returns the distance to the nearest point of (i,j) point.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        """
        (lons, lats) = self.get_lonlat_grid(position=position)
        i, j = as_numpy_array(i), as_numpy_array(j)
        assert len(i) == len(j), "Both coordinates must have the same length"
        result = numpy.zeros(len(i))
        for k in range(len(i)):
            dist = self.distance((lons[j[k], i[k]], lats[j[k], i[k]]),
                                 (stretch_array(lons), stretch_array(lats)))
            result[k] = dist[dist != 0].min()
        return result.squeeze()

    def getcenter(self, position=None):
        """
        Returns the coordinate of the grid center as a tuple of Angles
        (center_lon, center_lat).

        Caution: this is computed as the raw average of all grid points.
        A barycentric computation would be more adequate.
        """
        (lons, lats) = self.get_lonlat_grid(position=position)
        return (Angle(lons.mean(), 'degrees'),
                Angle(lats.mean(), 'degrees'))

    def _what_grid(self, out=sys.stdout):
        """
        Writes in file a summary of the grid of the field.

        :param out: the output open file-like object
        """
        (lons, lats) = self.get_lonlat_grid()
        corners = self.gimme_corners_ll()
        write_formatted(out, "Kind of Geometry", 'Unstructured')
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
