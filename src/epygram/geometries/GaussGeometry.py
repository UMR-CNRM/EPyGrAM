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
import re

import footprints
from footprints import proxy as fpx
from bronx.syntax.arrays import stretch_array

from epygram import epygramError, config
from epygram.util import (degrees_nearest_mod, Angle,
                          positive_longitudes, longitudes_between_minus180_180,
                          write_formatted,
                          nearlyEqual,
                          as_numpy_array, moveaxis)

from .VGeometry import VGeometry
from .AbstractGeometry import Geometry
from . import _need_pyproj_geod

epylog = footprints.loggers.getLogger(__name__)
_re_nearest_sq = re.compile(r'(?P<n>\d+)\*(?P<m>\d+)')



class GaussGeometry(Geometry):
    """
    Handles the geometry for a Global Gauss grid 3-Dimensions Field.
    """
    _ghost_attributes = Geometry._ghost_attributes + ['_buffered_gauss_grid']

    def __init__(self, name, grid, dimensions, vcoordinate,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
        :param name: Name of geometrical type of representation of points on the Globe.
                     Name must be among ['rotated_reduced_gauss', 'reduced_gauss', 'regular_gauss'
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
        self.add_attr_inlist('name', ['rotated_reduced_gauss', 'reduced_gauss', 'regular_gauss'])

        self.name = name

        super(GaussGeometry, self).__init__(grid, dimensions, vcoordinate,
                                            position_on_horizontal_grid, geoid)

    @property
    def isglobal(self):
        """
        :return: True if geometry is global
        """
        return True

    def _consistency_check(self):
        """Check that the geometry is consistent."""
        grid_keys = ['dilatation_coef', 'latitudes']
        if self.name == 'rotated_reduced_gauss':
            grid_keys.extend(['pole_lon', 'pole_lat'])
        if set(self.grid.keys()) != set(grid_keys):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys))
        assert isinstance(self.grid['latitudes'][0], Angle), \
            "'latitudes' attribute of grid must be a list of Angle objects."
        if self.name == 'rotated_reduced_gauss':
            assert isinstance(self.grid['pole_lon'], Angle)
            assert isinstance(self.grid['pole_lat'], Angle)
        dimensions_keys = ['max_lon_number', 'lat_number', 'lon_number_by_lat']
        if set(self.dimensions.keys()) != set(dimensions_keys):
            raise epygramError("dimensions attribute must consist in keys: " +
                               str(dimensions_keys))
        if self.name == 'rotated_reduced_gauss' \
           and nearlyEqual(self.grid['pole_lon'].get('degrees'), 0.) \
           and nearlyEqual(self.grid['pole_lat'].get('degrees'), 90.):
            self._attributes['name'] = 'reduced_gauss'
            self.grid.pop('pole_lon')
            self.grid.pop('pole_lat')

    def suggested_GRIB2_sample(self, spectral=False):
        if self.structure == 'H2D':
            prefix = {'rotated_reduced_gauss':'reduced_rotated_gg',
                      'reduced_gauss':'reduced_gg',
                      'regular_gauss':'regular_gg'}[self.name]
            return self._GRIB2_sample('sh' if spectral else prefix)
        else:
            raise NotImplementedError()

    def ij2ll(self, i, j, position=None):
        """
        Return the (lon, lat) coordinates of point *(i,j)*, in degrees.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        if self._getoffset(position) != (0., 0.):
            raise NotImplementedError("horizontal staggered grids for reduced\
                                       gauss grid are not implemented.")

        if isinstance(i, list) or isinstance(i, tuple) or\
           isinstance(i, numpy.ndarray):
            latitudes = as_numpy_array([self.grid['latitudes'][n].get('degrees') for n in range(len(self.grid['latitudes']))])
            dimensions = as_numpy_array(self.dimensions['lon_number_by_lat'])
            lat = latitudes[j]
            lon = (numpy.pi * 2 * as_numpy_array(i)) / dimensions[j]
        else:
            lat = self.grid['latitudes'][j].get('degrees')
            lon = (numpy.pi * 2. * i) / self.dimensions['lon_number_by_lat'][j]

        lon = numpy.degrees(lon)
        (lon, lat) = self._rotate_stretch(lon, lat, reverse=True)
        return (lon, lat)

    def ll2ij(self, lon, lat, position=None):
        """
        Return the (*i, j*) coordinates in the 2D matrix of gridpoints,
        of the gridpoint nearest to (*lon*, *lat*).

        :param lon: longitude of point in degrees
        :param lat: latitude of point in degrees
        :param position: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        lon, lat = as_numpy_array(lon), as_numpy_array(lat)
        ij = moveaxis(self.nearest_points(lon, lat, {'n': '1'}, position), 0, -1).squeeze()
        return (ij[0], ij[1])

    def _allocate_colocation_grid(self, compressed=False, as_float=False):
        """
        Creates the array for lonlat grid.
        Just a trick to avoid recomputing the array for several fields that
        share their geometry.

        :param compressed: if True, return 1D arrays, else 2D masked arrays.
        :param as_float: if True, return arrays with dtype float64, else int64.
        """
        Jmax = self.dimensions['lat_number']
        Imax = self.dimensions['lon_number_by_lat']
        igrid = []
        jgrid = []
        for j in range(Jmax):
            for i in range(Imax[j]):
                igrid.append(i)
                jgrid.append(j)
        if as_float:
            igrid = numpy.array(igrid, dtype=numpy.float64)
            jgrid = numpy.array(jgrid, dtype=numpy.float64)
        else:
            igrid = numpy.array(igrid)
            jgrid = numpy.array(jgrid)
        if not compressed:
            igrid = self.reshape_data(igrid)
            jgrid = self.reshape_data(jgrid)
        return (igrid, jgrid)

    def _clear_buffered_gauss_grid(self):
        """Deletes the buffered lonlat grid if any."""
        if hasattr(self, '_buffered_gauss_grid'):
            del self._buffered_gauss_grid

    @property
    def gridpoints_number(self, **_):
        """Returns the number of gridpoints of the grid."""
        return sum(self.dimensions['lon_number_by_lat'])

    def get_lonlat_grid(self,
                        position=None,
                        d4=False,
                        nb_validities=0,
                        force_longitudes=None,
                        **_):
        """
        Returns a tuple of two tables containing one the longitude of each
        point, the other the latitude, with 2D shape.

        Shape of 2D data in Gauss grids: \n
          - grid[0, 0:Nj] is first (Northern) band of latitude, masked after
            Nj = number of longitudes for latitude j \n
          - grid[-1, 0:Nj] is last (Southern) band of latitude (idem).

        :param position: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with respect to geometry.
                      d4=True requires nb_validities > 0
        :param nb_validities: number of validities represented in data values
        :param force_longitudes: if 'positive', the longitudes will be forced positive
                                 if ']-180,180]', the longitudes will be in the ]-180, 180] interval
        """
        # !!! **_ enables the method to receive arguments specific to
        #     other geometries but useless here ! Do not remove.
        if hasattr(self, '_buffered_gauss_grid') and \
           self._buffered_gauss_grid.get('filled'):
            lons = self._buffered_gauss_grid['lons']
            lats = self._buffered_gauss_grid['lats']
        else:
            (igrid, jgrid) = self._allocate_colocation_grid(compressed=True,
                                                            as_float=False)
            (lons, lats) = self.ij2ll(igrid, jgrid, position)
            lons = self.reshape_data(lons)
            lats = self.reshape_data(lats)
            if config.FA_buffered_gauss_grid:
                if not hasattr(self, '_buffered_gauss_grid'):
                    self._buffered_gauss_grid = {'lons':lons, 'lats':lats}
                else:
                    # trick: the arrays remain pointers to where they were
                    # created, so that they can be shared by several geometry
                    # objects or fields !
                    self._buffered_gauss_grid['lons'].data[...] = lons[...].data
                    self._buffered_gauss_grid['lons'].mask[...] = lons[...].mask
                    self._buffered_gauss_grid['lats'].data[...] = lats[...].data
                    self._buffered_gauss_grid['lats'].mask[...] = lons[...].mask
                self._buffered_gauss_grid['filled'] = True

        if d4:
            lons, lats = self._reshape_lonlat_4d(lons, lats, nb_validities)
        elif not d4 and nb_validities != 0:
            raise ValueError("*nb_validities* must be 0 when d4==False")
        if force_longitudes == 'positive':
            lons = positive_longitudes(lons)
        elif force_longitudes == ']-180,180]':
            lons = longitudes_between_minus180_180(lons)
        return (lons, lats)

    def get_datashape(self, force_dimZ=None, dimT=None, d4=False, **_):
        """
        Returns the data shape according to the geometry.

        :param force_dimZ: if supplied, force the Z dimension instead of that
          of the vertical geometry
        :param dimT: if supplied, is the time dimension to be added to the
          data shape
        :param d4: - if True,  shape is 4D (need to specify *dimT*)
                   - if False, shape is 3D if dimZ > 1 else 2D
        """
        dimY = self.dimensions['lat_number']
        dimX = self.dimensions['max_lon_number']
        if force_dimZ is not None:
            dimZ = force_dimZ
        else:
            dimZ = len(self.vcoordinate.levels)
        if d4:
            assert dimT is not None, \
                "*dimT* must be supplied with *d4*=True"
            shape = [dimT, dimZ, dimY, dimX]
        else:
            shape = []
            if self.datashape['k'] or dimZ > 1:
                shape.append(dimZ)
            shape.append(dimY)
            shape.append(dimX)
        return tuple(shape)

    def reshape_data(self, data, first_dimension=None, d4=False):
        """
        Returns a 2D data (horizontal dimensions) reshaped from 1D,
        according to geometry.

        :param data: the 1D data (or 3D with a T and Z dimensions,
          or 2D with either a T/Z dimension, to be specified),
          of dimension concording with geometry. In case data is 3D, T must be
          first dimension and Z the second.
        :param first_dimension: in case data is 2D, specify what is the first
          dimension of data among ('T', 'Z')
        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with respect to geometry
        """
        assert 1 <= len(data.shape) <= 3
        shp_in = data.shape
        nb_levels = 1
        nb_validities = 1
        if len(shp_in) == 2:
            assert first_dimension in ('T', 'Z'), \
                "*first_dimension* must be among ('T', 'Z') if *data*.shape == 2"
            if first_dimension == 'T':
                nb_validities = shp_in[0]
            elif first_dimension == 'Z':
                nb_levels = shp_in[0]
        elif len(shp_in) == 3:
            nb_validities = shp_in[0]
            nb_levels = shp_in[1]
        assert nb_levels in (1, len(self.vcoordinate.levels)), \
               "vertical dimension of data must be 1 or self.vcoordinate.levels=" + \
               str(self.vcoordinate.levels)
        shp4D = self.get_datashape(dimT=nb_validities, force_dimZ=nb_levels, d4=True)
        data4D = numpy.ma.masked_all(shp4D)
        ind_end = 0
        for j in range(self.dimensions['lat_number']):
            ind_begin = ind_end
            ind_end = ind_begin + self.dimensions['lon_number_by_lat'][j]
            if len(shp_in) == 1:
                buff = data[slice(ind_begin, ind_end)]
                data4D[0, 0, j, slice(0, self.dimensions['lon_number_by_lat'][j])] = buff
            elif len(shp_in) == 2:
                buff = data[:, slice(ind_begin, ind_end)]
                if nb_levels > 1:
                    data4D[0, :, j, slice(0, self.dimensions['lon_number_by_lat'][j])] = buff
                else:
                    data4D[:, 0, j, slice(0, self.dimensions['lon_number_by_lat'][j])] = buff
            elif len(shp_in) == 3:
                buff = data[:, :, slice(ind_begin, ind_end)]
                data4D[:, :, j, slice(0, self.dimensions['lon_number_by_lat'][j])] = buff
        if ind_end != data.shape[-1]:
            raise epygramError("data have a wrong length")
        if d4 or len(shp_in) == 3:
            data_out = data4D
        else:
            if len(shp_in) == 1:
                data_out = data4D[0, 0, :, :]
            elif len(shp_in) == 2:
                if first_dimension == 'T':
                    data_out = data4D[:, 0, :, :]
                elif first_dimension == 'Z':
                    data_out = data4D[0, :, :, :]
        return data_out

    def fill_maskedvalues(self, data, fill_value=None):
        """
        Returns a copy of *data* with 'real' masked values (i.e. not those
        linked to reduced Gauss) filled with *fill_value*.
        *data* must be already 4D for simplicity reasons.
        """
        assert isinstance(data, numpy.ma.masked_array)
        assert len(data.shape) == 4
        data_filled = numpy.ma.masked_array(data.filled(fill_value))
        mask = numpy.zeros(data_filled.data.shape, dtype=bool)
        for j in range(self.dimensions['lat_number']):
            i0 = self.dimensions['lon_number_by_lat'][j]
            mask[:, :, j, i0:] = True
        data_filled.mask = mask
        return data_filled

    def horizontally_flattened(self, data):
        """
        Returns a copy of *data* with horizontal dimensions flattened and
        compressed (cf. numpy.ma.masked_array.compressed).
        *data* must be 4D for simplicity reasons.
        """
        assert len(data.shape) == 4
        assert isinstance(data, numpy.ma.masked_array)
        data3D = numpy.empty(tuple(list(data.shape[:2]) + [self.gridpoints_number]))
        for t in range(data.shape[0]):
            for k in range(data.shape[1]):
                data3D[t, k, :] = data[t, k, :, :].compressed()
        return data3D

    def resolution_ll(self, lon, lat):
        """
        Returns the average meridian resolution (worst directional resolution)
        at point position.

        :param lon: longitude of the point in degrees
        :param lat: latitude of the point in degrees
        """
        return self.resolution_j(self.ll2ij(lon, lat)[1])

    def meridian_resolution_j(self, j):
        """
        Returns the average meridian resolution at longitude circle number *j*.
        """
        jint = numpy.rint(j).astype('int')
        jm1 = jint - 1
        jp1 = jint + 1
        if jm1 == -1:
            dist = self.distance(self.ij2ll(0, jint),
                                 self.ij2ll(0, jp1))
        elif jp1 == self.dimensions['lat_number']:
            dist = self.distance(self.ij2ll(0, jint),
                                 self.ij2ll(0, jm1))
        else:
            dist = (self.distance(self.ij2ll(0, jint),
                                  self.ij2ll(0, jp1)) +
                    self.distance(self.ij2ll(0, jint),
                                  self.ij2ll(0, jm1))) / 2.
        return dist

    def zonal_resolution_j(self, j):
        """
        Returns the average zonal resolution at longitude circle number j.
        """
        jint = numpy.rint(j).astype('int')
        return self.distance(self.ij2ll(0, jint),
                             self.ij2ll(1, jint))

    @_need_pyproj_geod
    def resolution_field_from_stretching(self):
        """
        Returns a field which values are the local resolution computed as the
        nominal resolution stretched locally by the map factor.
        """
        assert self._pyproj_geod.sphere, "Method is not available with a non-spheroid geoid."
        zonal_equatorial_resolution = 2. * numpy.pi * self._pyproj_geod.a / self.dimensions['max_lon_number']
        mf = self.map_factor_field()
        mf.fid['geometry'] = 'resolution_from_stretching'
        mf.setdata(zonal_equatorial_resolution / mf.data)
        return mf

    def resolution_j(self, j):
        """
        Returns the average meridian resolution at longitude circle number j.
        """
        return self.meridian_resolution_j(j)

    def resolution_field(self, direction='meridian'):
        """
        Returns a field whose values are the local resolution in m.

        :param direction: among ('zonal', 'meridian'), direction in which
                          the resolution is computed.
        """
        assert direction in ('zonal', 'meridian')
        resolutions = [getattr(self, direction + '_resolution_j')(j)
                       for j in range(self.dimensions['lat_number'])]
        resol_2d = (numpy.ma.ones(self.get_lonlat_grid()[0].data.shape).transpose() *
                    numpy.array(resolutions)).transpose()
        resol_2d.mask = self.get_lonlat_grid()[0].mask
        f = fpx.field(structure='H2D',
                      geometry=self,
                      fid={'geometry':direction + ' resolution'},
                      units='m')
        f.setdata(resol_2d)
        return f

    def distance_to_nearest_neighbour_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.

        :param lon: longitude of the point in degrees
        :param lat: latitude of the point in degrees
        """
        return self.distance_to_nearest_neighbour_ij(*self.ll2ij(lon, lat))

    def distance_to_nearest_neighbour_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        """
        # FIXME: not sure this is exactly computed
        (iint, jint) = (numpy.rint(i).astype('int'),
                        numpy.rint(j).astype('int'))
        points_list = []
        for oi in [-1, 0, 1]:
            for oj in [-1, 0, 1]:
                oj = 0
                if (oi, oj) != (0, 0):
                    pi = iint + oi
                    pj = jint + oj
                    if pj < self.dimensions['lat_number'] and pj >= 0:
                        pi = pi % self.dimensions['lon_number_by_lat'][pj]
                        points_list.append((pi, pj))

        return numpy.array([self.distance(self.ij2ll(iint, jint),
                                          self.ij2ll(*p))
                            for p in points_list]).min()

    def point_is_inside_domain_ll(self, lon, *_, **__):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.
        This is always the case in Gauss grids, no real meaning.

        :param lon: longitude of the point in degrees
        :param lat: latitude of the point in degrees
        :param margin: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        :param position: position of the grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        return numpy.ones_like(lon, dtype=bool)

    def point_is_inside_domain_ij(self, i, j, margin=-0.1):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.

        :param i: X index of point in the 2D matrix of gridpoints
        :param j: Y index of point in the 2D matrix of gridpoints
        :param margin: DEPRECATED
        """
        i = as_numpy_array(i) if isinstance(i, (list, tuple)) else i
        j = as_numpy_array(j) if isinstance(j, (list, tuple)) else j
        dimensions = as_numpy_array(self.dimensions['lon_number_by_lat'])
        # Firstly we test the validity of j.
        # In case j is invalid result will be False
        # but we need a valid value for j to test the i validity
        j2 = numpy.maximum(0, numpy.minimum(j, self.dimensions['lat_number'] - 1))
        inside = numpy.logical_and(numpy.logical_and(j >= 0,
                                                     j < self.dimensions['lat_number']),
                                   numpy.logical_and(i >= 0,
                                                     i < dimensions[j2]))
        return inside

    def _rotate_stretch(self, lon, lat, reverse=False):
        """
        Internal method used to transform true lon/lat into
        rotated and stretched lon/lat.

        :param lon: longitude in degrees
        :param lat: latitude in degrees
        :param reverse: if True, do the reverse transform.

        Computation adapted from arpifs/transform/trareca.F90 and tracare.F90.
        """
        if self.name == 'rotated_reduced_gauss':
            KTYP = 2
        else:
            KTYP = 1
        PFACDI = self.grid['dilatation_coef']
        ZDBLC = 2.0 * PFACDI
        ZC2P1 = PFACDI * PFACDI + 1.
        ZC2M1 = PFACDI * PFACDI - 1.
        ZCRAP = -ZC2M1 / ZC2P1
        ZEPS = config.epsilon

        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)

        if not reverse:
            if 'rotated' in self.name:
                lon0 = self.grid['pole_lon'].get('degrees')
            else:
                lon0 = 0.0
            lon = degrees_nearest_mod(lon, lon0)
            PSLAR = numpy.sin(numpy.radians(lat))
            PSLOR = numpy.sin(numpy.radians(lon))
            PCLOR = numpy.cos(numpy.radians(lon))

            if KTYP == 2:
                # Complete ARPEGE
                (PCLAP, PSLAP) = self.grid['pole_lat'].get('cos_sin')
                (PCLOP, PSLOP) = self.grid['pole_lon'].get('cos_sin')
                ZCLAR = numpy.sqrt(1. - PSLAR * PSLAR)
                ZA = PSLAP * PSLAR + PCLAP * ZCLAR * \
                     (PCLOP * PCLOR + PSLOP * PSLOR)
                PSLAC = (ZCRAP + ZA) / (1. + ZA * ZCRAP)
                ZB = 1. / numpy.sqrt(numpy.maximum(ZEPS, 1. - ZA * ZA))
                PCLOC = (PCLAP * PSLAR - PSLAP * ZCLAR *
                         (PCLOP * PCLOR + PSLOP * PSLOR)) * ZB
                PSLOC = ZCLAR * (PSLOP * PCLOR - PCLOP * PSLOR) * ZB
                PSLAC = numpy.maximum(-1., numpy.minimum(1., PSLAC))
                PCLOC = numpy.maximum(-1., numpy.minimum(1., PCLOC))
                PSLOC = numpy.maximum(-1., numpy.minimum(1., PSLOC))
                if isinstance(PSLOC, numpy.ndarray):
                    for k in range(0, len(PSLOC)):
                        if PCLOC[k] * PCLOC[k] + PSLOC[k] * PSLOC[k] < 0.5:
                            PSLOC[k] = 1.
                            PCLOC[k] = 0.
                else:
                    if PCLOC * PCLOC + PSLOC * PSLOC < 0.5:
                        PSLOC = 1.
                        PCLOC = 0.
            else:
                # Schmidt
                PSLOC = PSLOR
                PCLOC = PCLOR
                PSLAC = (ZC2P1 * PSLAR - ZC2M1) / (ZC2P1 - ZC2M1 * PSLAR)
                PSLAC = numpy.maximum(-1., numpy.minimum(1., PSLAC))

            lat = numpy.arcsin(PSLAC)
            lon = numpy.arctan2(PSLOC, PCLOC) % (numpy.pi * 2)
        elif reverse:
            PSLAC = numpy.sin(numpy.radians(lat))
            PCLOC = numpy.cos(numpy.radians(lon))
            PSLOC = numpy.sin(numpy.radians(lon))
            if KTYP == 2:
                # Complete ARPEGE
                (PCLAP, PSLAP) = self.grid['pole_lat'].get('cos_sin')
                (PCLOP, PSLOP) = self.grid['pole_lon'].get('cos_sin')
                ZCLAC = numpy.sqrt(1. - PSLAC * PSLAC)
                ZA = 1. / (ZC2P1 + ZC2M1 * PSLAC)
                ZB = ZC2P1 * PSLAC + ZC2M1
                ZC = ZDBLC * PCLAP * ZCLAC * PCLOC + ZB * PSLAP
                PSLAR = ZC * ZA
                ZD = ZA / numpy.sqrt(numpy.maximum(ZEPS, 1. - PSLAR * PSLAR))
                ZE = ZB * PCLAP * PCLOP - ZDBLC * ZCLAC * \
                     (PSLAP * PCLOC * PCLOP - PSLOP * PSLOC)
                ZF = ZB * PCLAP * PSLOP - ZDBLC * ZCLAC * \
                     (PSLAP * PCLOC * PSLOP + PCLOP * PSLOC)
                PCLOR = ZE * ZD
                PSLOR = ZF * ZD
                PSLAR = numpy.maximum(-1., numpy.minimum(1., PSLAR))
                PCLOR = numpy.maximum(-1., numpy.minimum(1., PCLOR))
                PSLOR = numpy.maximum(-1., numpy.minimum(1., PSLOR))
            else:
                # Schmidt
                PSLOR = PSLOC
                PCLOR = PCLOC
                PSLAR = (ZC2P1 * PSLAC + ZC2M1) / (ZC2P1 + ZC2M1 * PSLAC)
                PSLAR = numpy.maximum(-1., numpy.minimum(1., PSLAR))

            lon = numpy.arctan2(PSLOR, PCLOR)
            lat = numpy.arcsin(PSLAR)
        lon = numpy.degrees(lon)
        lat = numpy.degrees(lat)

        return (lon, lat)

    def map_factor(self, lon, lat):
        """
        Returns the map factor at the given longitude/latitude(s)

        :param lon: longitude in degrees
        :param lat: latitude in degrees
        """
        pc = self.grid['dilatation_coef']
        (_, plab) = self._rotate_stretch(lon, lat)
        zlat1 = numpy.radians(plab)
        # From rotated/streched sphere to rotated
        zinterm = 1. / pc * numpy.cos(zlat1) / (1. + numpy.sin(zlat1))
        zlat2 = 2. * numpy.arctan((1. - zinterm) / (1. + zinterm))
        zm = numpy.cos(zlat1) / numpy.cos(zlat2)
        return zm

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
        (lons, lats) = self.get_lonlat_grid(position=position)
        data = self.map_factor(stretch_array(lons), stretch_array(lats))
        data = self.reshape_data(data)
        f.setdata(data)
        return f

    def reproject_wind_on_lonlat(self, u, v, lon, lat,
                                 map_factor_correction=True,
                                 reverse=False):
        """
        Reprojects a wind vector (u, v) on rotated/stretched sphere onto real
        sphere, i.e. with components on true zonal/meridian axes.

        :param u: the u == zonal-on-the-grid component of wind
        :param v: the v == meridian-on-the-grid component of wind
        :param lon: longitudes of points in degrees
        :param lat: latitudes of points in degrees
        :param map_factor_correction: applies a correction of magnitude due
                                      to map factor.
        :param reverse: if True, apply the reverse reprojection.

        lon/lat are coordinates on real sphere.
        """
        pc = self.grid['dilatation_coef']
        if 'rotated' in self.name:
            plac = self.grid['pole_lat'].get('degrees')
        else:
            plac = 90.
        (plob, plab) = self._rotate_stretch(lon, lat)
        # the below formulas seem to be working with
        # lon_on_rotated_sphere begining below pole_of_rotation,
        # hence a 180° rotation.
        plob += 180.
        pust = u
        pvst = v

        # Adapted from J.M.Piriou's pastrv.F90
        zlon1 = numpy.radians(plob)
        zlat1 = numpy.radians(plab)
        zlatp = numpy.radians(plac)

        # From rotated/streched sphere to rotated
        zinterm = 1. / pc * numpy.cos(zlat1) / (1. + numpy.sin(zlat1))
        zlat2 = 2. * numpy.arctan((1. - zinterm) / (1. + zinterm))
        zlon2 = zlon1

        # Map factor
        if map_factor_correction:
            zm = numpy.cos(zlat1) / numpy.cos(zlat2)
            # zm = self.map_factor(lon, lat) # but redundant computations
        else:
            if self.grid['dilatation_coef'] != 1.:
                epylog.warning('check carefully *map_factor_correction* w.r.t. dilatation_coef')
            zm = numpy.ones(zlat1.shape)

        # From rotated sphere to real sphere
        # Compute latitude on real sphere
        zsla3 = -numpy.cos(zlat2) * numpy.cos(zlon2) * numpy.cos(zlatp) + numpy.sin(zlat2) * numpy.sin(zlatp)
        if (zsla3 > 1.).any() or (zsla3 < -1.).any():
            epylog.warning('reproject_wind_on_lonlat: zsla3=' + str(zsla3))
            zsla3 = min(1., max(-1., zsla3))
        zlat3 = numpy.arcsin(zsla3)

        # Real North components
        zca = (numpy.cos(zlat2) * numpy.sin(zlatp) +
               numpy.sin(zlat2) * numpy.cos(zlatp) * numpy.cos(zlon2)) / \
               numpy.cos(zlat3)
        zsa = numpy.cos(zlatp) * numpy.sin(zlon2) / numpy.cos(zlat3)

        # Wind transformation
        if reverse:
            zm = 1. / zm
            zsa = -zsa
        pusr = zm * (zca * pust - zsa * pvst)
        pvsr = zm * (zsa * pust + zca * pvst)

        return (pusr, pvsr)

    def nearest_points(self, lon, lat, request,
                       position=None,
                       external_distance=None,
                       squeeze=True):
        """
        Returns a list of the (i, j) position of the points needed to perform
        an interpolation.

        :param lon: longitude of point in degrees.
        :param lat: latitude of point in degrees.
        :param request: criteria for selecting the points, among:
               * {'n':'1'} - the nearest point
               * {'n':'2*2'} - the 2*2 square points around the position
               * {'n':'4*4'} - the 4*4 square points around the position
               * {'n':'N*N'} - the N*N square points around the position: N must be even
               * {'radius':xxxx, 'shape':'square'} - the points which are xxxx metres
                 around the position in X or Y direction
               * {'radius':xxxx, 'shape':'circle'} - the points within xxxx metres
                 around the position. (default shape == circle)
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
        if self._getoffset(position) != (0., 0.):
            raise NotImplementedError("horizontal staggered grids for " +
                                      "reduced gauss grid are not implemented.")

        lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()

        # ## internal functions
        def nearest_lats(latrs, num):
            """
            Internal method used to find the nearest latitude circles.

            *latrs* is the rotated and streched latitude.
            *num* is:
            - 2 to find the two surrounding latitude cicles,
            - 3 to find the nearest latitude circles plus the one above
                and the one under.
            - 4 to find the two surrounding latitude circles plus the
                preceding one and the following one.
            - and so on...

            returns an array of shape (len(latrs), num)
            """
            if not numpy.all(numpy.array(latitudes[1:]) <= numpy.array(latitudes[:-1])):
                raise ValueError('latitudes must be in descending order')

            latrs = as_numpy_array(latrs)
            nearest = numpy.ndarray((len(latrs), ), dtype=int)
            distmin = numpy.ndarray((len(latrs), ), dtype=float)
            # Slicing is needed to prevent memory error
            batch_size = 100  # on a particular example, increasing this value implies a longuer execution time
            for imin in range(0, len(latrs), batch_size):
                imax = min(imin + batch_size, len(latrs))
                dist = latrs[imin:imax, numpy.newaxis] - as_numpy_array(latitudes)[numpy.newaxis, :]
                nearest[imin:imax] = numpy.argmin(numpy.abs(dist), axis=1)
                distmin[imin:imax] = dist[numpy.arange(dist.shape[0]), nearest[imin:imax]]
            result = numpy.zeros((len(latrs), num), dtype=int)
            if num == 2:
                result[:, 1] = - numpy.copysign(1, distmin).astype('int')
            elif num == 3:
                result[:, 0] = -1
                result[:, 2] = 1
            elif num % 2:  # odd
                for k in range(1, num // 2 + 1):
                    result[:, num // 2 - k] = - k
                    result[:, num // 2 + k] = k
            elif not num % 2:  # even
                for k in range(1, num // 2):
                    result[:, num // 2 - k - 1] = - k
                    result[:, num // 2 + k - 1] = k
                result[:, -1] = num // 2
            result = result + nearest[:, numpy.newaxis]
            return result

        def nearest_lons(lonrs, latnum, num):
            """
            Internal method used to find the nearest points on a latitude
            circle.

            Args:\n
            - *lonrs* is the rotated and streched longitude.
            - *num* is:
              - 1 to find the nearest point,
              - 2 to find the two surrounding points,
              - 4 to find the two surrounding points plus the preceding one and
                  the following one.
              - and so on
            - *latnum*: if -1 (resp. -2), we search for the opposite longitude
              on the first (resp. second) latitude circle.
              The same is true for the other pole.

            latnum must be of shape (len(lonrs), x) where x is the number of latitude
            to search for each longitude

            result shape is (len(lonrs), x, num, 2)
            """
            lonrs = as_numpy_array(lonrs)
            latnum = as_numpy_array(latnum)
            assert len(lonrs.shape) == 1 and lonrs.shape == latnum.shape[:1]
            result = numpy.ndarray(tuple(list(latnum.shape) + [num, 2]), dtype=int)
            lonrs = numpy.radians(lonrs)
            lonnummax = numpy.ndarray(latnum.shape, dtype=int)
            j = numpy.ndarray(latnum.shape, dtype=int)
            i = numpy.ndarray(latnum.shape, dtype=float)
            lonrs2d = numpy.repeat(lonrs[:, numpy.newaxis], latnum.shape[1], axis=1)
            # Near the pole, have to look-up on the first latitudes
            # circles, symmetrically with regards to the pole
            mask1 = latnum < 0
            j[mask1] = abs(latnum[mask1]) - 1
            lonnummax[mask1] = as_numpy_array(self.dimensions['lon_number_by_lat'])[j[mask1]]
            i[mask1] = ((lonrs2d[mask1] - numpy.pi) % (numpy.pi * 2)) * \
                       lonnummax[mask1] / (numpy.pi * 2)

            mask2 = latnum >= len(latitudes)
            j[mask2] = len(latitudes) - (latnum[mask2] - len(latitudes)) - 1  # TOBECHECKED: next circle past the pole
            lonnummax[mask2] = as_numpy_array(self.dimensions['lon_number_by_lat'])[j[mask2]]
            i[mask2] = ((lonrs2d[mask2] - numpy.pi) % (numpy.pi * 2)) * \
                            lonnummax[mask2] / (numpy.pi * 2)

            mask = numpy.logical_not(numpy.logical_or(mask1, mask2))
            j[mask] = latnum[mask]
            lonnummax[mask] = as_numpy_array(self.dimensions['lon_number_by_lat'])[latnum[mask]]
            i[mask] = lonrs2d[mask] * lonnummax[mask] / (numpy.pi * 2)

            if num == 1:
                result[:, :, 0, 0] = numpy.rint(i).astype('int') % lonnummax
                result[:, :, 0, 1] = j
            elif num == 2:
                result[:, :, 0, 0] = numpy.floor(i).astype('int') % lonnummax
                result[:, :, 0, 1] = j
                result[:, :, 1, 0] = (numpy.floor(i).astype('int') + 1) % lonnummax
                result[:, :, 1, 1] = j
            elif num % 2:  # odd
                raise NotImplementedError("but is it necessary ?")
            elif not num % 2:  # even
                ii = numpy.floor(i).astype('int')
                for k in range(1, num // 2):
                    result[:, :, num // 2 - k - 1, 0] = (ii - k) % lonnummax
                    result[:, :, num // 2 - k - 1, 1] = j
                result[:, :, num // 2 - 1, 0] = ii % lonnummax
                result[:, :, num // 2 - 1, 1] = j
                for k in range(1, num // 2 + 1):
                    result[:, :, num // 2 + k - 1, 0] = (ii + num // 2 + 1 - k) % lonnummax  # reverse order? why?
                    result[:, :, num // 2 + k - 1, 1] = j

            return result

        def nearest(lon, lat, lonrs, latrs):
            """
            Internal method used to find the nearest point.
            lon/lat are the true coordinate, lonrs/latrs are rotated and
            streched coordinates.

            Returns an array of points
            """
            lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()
            lonrs, latrs = as_numpy_array(lonrs).flatten(), as_numpy_array(latrs).flatten()
            assert len(lon) == len(lat) == len(lonrs) == len(latrs)

            all_nearest_lats = nearest_lats(latrs, 3)
            all_nearest_lons = nearest_lons(lonrs, all_nearest_lats, 1)
            lon2d = numpy.repeat(lon[:, numpy.newaxis], 3, axis=1)
            lat2d = numpy.repeat(lat[:, numpy.newaxis], 3, axis=1)
            ll = self.ij2ll(*[all_nearest_lons[:, :, 0, 0].flatten(),
                              all_nearest_lons[:, :, 0, 1].flatten()])
            all_dist = self.distance((lon2d.flatten(), lat2d.flatten()), ll)
            all_dist = all_dist.reshape((len(lonrs), 3))
            result = all_nearest_lons[numpy.arange(all_dist.shape[0]), numpy.argmin(all_dist, axis=1), 0, :]
            return result

        def nearests(lonrs, latrs, n):
            """
            Internal methods used to find the n*n points surrunding the
            point lonrs/latrs, lonrs/latrs are rotated and stretched
            coordinates.

            lonrs and latrs can be arrays
            """
            lonrs = as_numpy_array(lonrs)
            latrs = as_numpy_array(latrs)
            assert lonrs.shape == latrs.shape and len(lonrs.shape) == 1, "only scalar or 1d arrays with same length"

            all_nearest_lats = nearest_lats(latrs, n)
            all_nearest_lons = nearest_lons(lonrs, all_nearest_lats, n)
            result = all_nearest_lons.reshape((len(lonrs), n**2, 2))

            return result

        # ## actual algorithm
        # initializations
        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)
        latitudes = [self.grid['latitudes'][n].get('degrees')
                     for n in range(0, self.dimensions['lat_number'])]
        # compute rotated/stretched lon/lat
        (lonrs, latrs) = self._rotate_stretch(lon, lat)
        # 1.A: nearest point is needed only
        nsquare_match = _re_nearest_sq.match(request.get('n', ''))
        if nsquare_match:
            assert nsquare_match.group('n') == nsquare_match.group('m'), \
                "anisotropic request {'n':'N*M'} is not supported."
        if external_distance is not None:
            assert request == {'n':'1'}

        def _increments(n):
            def _rng(n):
                return numpy.arange(-n // 2 + 1, n // 2 + 1)
            if self.name == 'academic' and self.dimensions['X'] == 1:
                i_incr = _rng(1)
            else:
                i_incr = _rng(n)
            if self.name == 'academic' and self.dimensions['Y'] == 1:
                j_incr = _rng(1)
            else:
                j_incr = _rng(n)
            return (i_incr, j_incr)
        if request == {'n':'1'} and not external_distance:
            result = nearest(lon, lat, lonrs, latrs)
        # 2.: several points are needed
        else:
            # 2.1: how many ?
            if external_distance:
                n = 1
            elif request.get('radius'):
                if isinstance(lat, numpy.ndarray) and len(lat) > 1:
                    raise NotImplementedError("request:{'radius':..., 'shape':'radius'} with several points.")
                resolution = self.resolution_ll(lon, lat)
                n = max(numpy.around(float(request['radius']) / float(resolution)).astype('int') * 2,
                        1)
            elif nsquare_match:
                n = int(nsquare_match.group('n'))
            else:
                raise epygramError("unrecognized **request**: " + str(request))
            # 2.2: get points
            result = nearests(lonrs, latrs, n)

            # 2.3: only select the nearest with regards to external_distance
            if external_distance:
                result = list(result)  # We transform into list to be able to modify the length
                for ipt, points in enumerate(result):
                    mindistance = None
                    for p in points:
                        dist = abs(external_distance['external_field'].getvalue_ij(*p, one=True) - external_distance['target_value'])
                        if mindistance is None or dist < mindistance:
                            result[ipt] = [p]
                            mindistance = dist
                result = as_numpy_array(result)  # all item must now have the same length

        if request.get('radius'):
            result = list(result)  # We transform into list to be able to modify the length
            for ipt, points in enumerate(result):
                if request.get('shape', 'circle') == 'circle':
                    #for i,j in points:
                    #    print(i,j,self.distance((lon[ipt], lat[ipt]), self.ij2ll(i, j)))
                    result[ipt] = [(i, j)
                                   for (i, j) in points
                                   if self.distance((lon[ipt], lat[ipt]), self.ij2ll(i, j)) <= request['radius']]
                elif request.get('shape') == 'square':
                    result[ipt] = [(i, j)
                                   for (i, j) in points
                                   if all(abs(numpy.array(self.ll2xy(lon[ipt], lat[ipt])) - numpy.array(self.ij2xy(i, j))) <= request['radius'])]
                assert len(result[ipt]) > 0, "no points found: radius may be too small."
            if squeeze and len(result) == 1:
                result = result[0]
        else:
            if squeeze:
                result = result.squeeze()

        return result

    def _what_grid_dimensions(self, out=sys.stdout, spectral_geometry=None):
        """
        Writes in file a summary of the grid & dimensions of the field.

        :param out: the output open file-like object
        :param spectral_geometry: an optional dict containing the spectral
          truncature {'max':}.
        """

        grid = self.grid
        dimensions = self.dimensions
        gridmap = {'reduced_gauss':'Reduced Gauss',
                   'rotated_reduced_gauss':'Rotated Reduced Gauss',
                   'regular_gauss':'Regular Gauss'}
        write_formatted(out, "Grid", gridmap[self.name])
        if self.name == 'rotated_reduced_gauss':
            write_formatted(out, "Pole Longitude",
                            grid['pole_lon'].get('degrees'))
            write_formatted(out, "Pole Latitude",
                            grid['pole_lat'].get('degrees'))
        write_formatted(out, "Dilatation coefficient",
                        grid['dilatation_coef'])
        write_formatted(out, "Number of latitudes",
                        dimensions['lat_number'])
        write_formatted(out, "Maximum number of longitudes on a parallel",
                        dimensions['max_lon_number'])
        if spectral_geometry is not None:
            write_formatted(out, "Truncation",
                            spectral_geometry['max'])
