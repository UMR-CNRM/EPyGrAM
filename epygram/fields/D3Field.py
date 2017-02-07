#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import copy
import numpy
import sys

import footprints
from footprints import FPDict, FPList, proxy as fpx
from epygram import epygramError
from epygram.util import write_formatted, stretch_array, Angle
from epygram.base import Field, FieldSet, FieldValidity, FieldValidityList, Resource
from epygram.geometries import D3Geometry, SpectralGeometry


class D3CommonField(Field):
    """
    3-Dimensions common field class.
    """

    _collector = ('field',)
    _abstract = True
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['3D']),
                info="Type of Field geometry.")
        )
    )

    @property
    def spectral(self):
        """Returns True if the field is spectral."""
        return self.spectral_geometry is not None and self.spectral_geometry != '__unknown__'


##############
# ABOUT DATA #
##############

    def getvalue_ll(self, lon=None, lat=None, level=None, validity=None,
                    interpolation='nearest',
                    neighborinfo=False,
                    one=True,
                    external_distance=None):
        """
        Returns the value of the field on point of coordinates (*lon, lat, level*): \n
        - if *interpolation == 'nearest'* (default), returns the value of the
          nearest neighboring gridpoint;
        - if *interpolation == 'linear'*, computes and returns the field value
          with linear spline interpolation;
        - if *interpolation == 'cubic'*, computes and returns the field value
          with cubic spline interpolation.
        *level* is the True level not the index of the level. Depending on the
        vertical coordinate, it could be expressed in Pa, m.
        *validity* is a FieldValidity or a FieldValidityList instance

        If *neighborinfo* is set to **True**, returns a tuple
        *(value, (lon, lat))*, with *(lon, lat)* being the actual coordinates
        of the neighboring gridpoint (only for *interpolation == 'nearest'*).

        *lon* and *lat* may be longer than 1.
        If *one* is False and len(lon) is 1, returns [value] instead of value.

        *external_distance* can be a dict containing the target point value
        and an external field on the same grid as self, to which the distance
        is computed within the 4 horizontally nearest points; e.g.
        {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
        If so, the nearest point is selected with
        distance = |target_value - external_field.data|

        Warning: for interpolation on Gauss geometries, requires the
        :mod:`pyproj` module.
        """

        if isinstance(validity, FieldValidity):
            myvalidity = FieldValidityList(validity)
        else:
            myvalidity = validity
        if self.spectral:
            raise epygramError("field must be gridpoint to get value of a" +
                               " lon/lat point.")
        if len(self.validity) > 1 and myvalidity is None:
            raise epygramError("*validity* is mandatory when there are several validities")
        if self.geometry.datashape['k'] and level is None:
            raise epygramError("*level* is mandatory when field has a vertical coordinate")
        if (self.geometry.datashape['j'] or self.geometry.datashape['i']) and \
           (lon is None or lat is None):
            raise epygramError("*lon* and *lat* are mandatory when field has an horizontal extension")

        maxsize = numpy.array([numpy.array(dim).size for dim in [lon, lat, level] if dim is not None]).max()
        if myvalidity is not None:
            maxsize = max(maxsize, len(myvalidity))

        # We look for indexes for vertical and time coordinates (no interpolation)
        if myvalidity is None:
            my_t = numpy.zeros(maxsize, dtype=int)
        else:
            my_t = []
            for v in myvalidity:
                for t in range(len(self.validity)):
                    if v == self.validity[t]:
                        my_t.append(t)
                        break
            my_t = numpy.array(my_t)
            if my_t.size != maxsize:
                if my_t.size != 1:
                    raise epygramError("validity must be uniq or must have the same length as other indexes")
                my_t = numpy.array([my_t.item()] * maxsize)
        if len(self.geometry.vcoordinate.levels) != len(set(self.geometry.vcoordinate.levels)):
            raise epygramError('Some levels are represented twice in levels list.')
        if level is None:
            my_k = numpy.zeros(maxsize, dtype=int)
        else:
            my_level = numpy.array(level)
            if my_level.size == 1:
                my_level = numpy.array([my_level.item()])
            my_k = []
            for l in my_level:
                my_k.append(self.geometry.vcoordinate.levels.index(l))
            my_k = numpy.array(my_k)
            if my_k.size != maxsize:
                if my_k.size != 1:
                    raise epygramError("k must be scalar or must have the same length as other indexes")
                my_k = numpy.array([my_k.item()] * maxsize)
        my_lon = numpy.array(lon)
        if my_lon.size == 1:
            my_lon = numpy.array([my_lon.item()])
        my_lat = numpy.array(lat)
        if my_lat.size == 1:
            my_lat = numpy.array([my_lat.item()])
        if my_lon.size != maxsize or my_lat.size != maxsize:
            raise epygramError("lon and lat must have the same length and the same length as level and validity")

        if interpolation == 'nearest':
            (ri, rj) = self.geometry.nearest_points(lon, lat, {'n':'1'},
                                                    external_distance=external_distance)
            value = self.getvalue_ij(ri, rj, my_k, my_t, one=one)
            if neighborinfo:
                (lon, lat) = self.geometry.ij2ll(ri, rj)
                if numpy.shape(lon) in ((1,), ()):
                    lon = float(lon)
                    lat = float(lat)
                value = (value, (lon, lat))
        elif interpolation in ('linear', 'cubic'):
            from scipy.interpolate import interp1d, interp2d
            nvalue = numpy.zeros(maxsize)
            interp_points = []
            for n in range(maxsize):
                if maxsize > 1:
                    lonn = lon[n]
                    latn = lat[n]
                    my_kn = my_k[n]
                    my_tn = my_t[n]
                else:
                    lonn = my_lon.item()
                    latn = my_lat.item()
                    my_kn = my_k.item()
                    my_tn = my_t.item()
                interp_points.append(self.geometry.nearest_points(lonn, latn,
                                                                  {'linear':{'n':'2*2'},
                                                                   'cubic':{'n':'4*4'}}[interpolation]))
            # depack
            all_i = []
            all_j = []
            for points in interp_points:
                for p in points:
                    all_i.append(p[0])
                    all_j.append(p[1])
            # get values and lons/lats
            flat_values_at_interp_points = list(self.getvalue_ij(all_i, all_j, my_kn, my_tn))
            all_lonslats = self.geometry.ij2ll(all_i, all_j)
            all_lonslats = (list(all_lonslats[0]), list(all_lonslats[1]))
            # repack and interpolate
            for n in range(maxsize):
                loc_values = [flat_values_at_interp_points.pop(0) for _ in range(len(interp_points[n]))]
                loc_lons = [all_lonslats[0].pop(0) for _ in range(len(interp_points[n]))]
                loc_lats = [all_lonslats[1].pop(0) for _ in range(len(interp_points[n]))]
                if self.geometry.name == 'academic' and \
                   1 in (self.geometry.dimensions['X'], self.geometry.dimensions['Y'] == 1):
                    if self.geometry.dimensions['X'] == 1:
                        f = interp1d(loc_lats, loc_values, kind=interpolation)
                        value = f(latn)
                    else:
                        f = interp1d(loc_lons, loc_values, kind=interpolation)
                        value = f(lonn)
                else:
                    f = interp2d(loc_lons, loc_lats, loc_values, kind=interpolation)
                    value = f(lonn, latn)
                nvalue[n] = value
            value = nvalue
        if one:
            try:
                value = float(value)
            except (ValueError, TypeError):
                pass

        return copy.copy(value)

    def as_lists(self, order='C', subzone=None):
        """
        Export values as a dict of lists (in fact numpy arrays).
        - *order*: whether to flatten arrays in 'C' (row-major) or
                   'F' (Fortran, column-major) order.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
                     among: 'C', 'CI'.
        """

        if self.spectral:
            raise epygramError("as_lists method needs a grid-point field, not a spectral one.")

        lons4d, lats4d = self.geometry.get_lonlat_grid(d4=True, nb_validities=len(self.validity), subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)
        data4d = self.getdata(d4=True, subzone=subzone)

        dates4d = []
        times4d = []
        for i, t in enumerate(self.validity):
            dates4d += [t.get().year * 10000 + t.get().month * 100 + t.get().day] * levels4d[i].size
            times4d += [t.get().hour * 100 + t.get().minute] * levels4d[i].size
        dates4d = numpy.array(dates4d).reshape(data4d.shape).flatten(order=order)
        times4d = numpy.array(times4d).reshape(data4d.shape).flatten(order=order)
        result = dict(values=data4d.flatten(order=order),
                      latitudes=lats4d.flatten(order=order),
                      longitudes=lons4d.flatten(order=order),
                      levels=levels4d.flatten(order=order),
                      dates=dates4d.flatten(order=order),
                      times=times4d.flatten(order=order))
        return result

    def as_dicts(self, subzone=None):
        """
        Export values as a list of dicts.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
                     among: 'C', 'CI'.
        """

        if self.spectral:
            raise epygramError("as_dicts method needs a grid-point field, not a spectral one.")

        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone)
        data4d = self.getdata(d4=True, subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)

        result = []
        for t in range(data4d.shape[0]):
            validity = self.validity[t]
            date = validity.get().year * 10000 + validity.get().month * 100 + validity.get().day
            time = validity.get().hour * 100 + validity.get().minute
            for k in range(data4d.shape[1]):
                for j in range(data4d.shape[2]):
                    for i in range(data4d.shape[3]):
                        result.append(dict(value=data4d[t, k, j, i],
                                           date=date,
                                           time=time,
                                           latitude=lats[j, i],
                                           longitude=lons[j, i],
                                           level=levels4d[t, k, j, i]))
        return result

    def as_points(self, subzone=None):
        """
        Export values as a fieldset of points.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
                     among: 'C', 'CI'.
        """

        if self.spectral:
            raise epygramError("as_points method needs a grid-point field, not a spectral one.")

        field_builder = fpx.field
        geom_builder = fpx.geometry
        vcoord_builer = fpx.geometry

        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone)
        data4d = self.getdata(d4=True, subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)

        result = FieldSet()
        kwargs_vcoord = copy.deepcopy(self.geometry.vcoordinate.footprint_as_dict())
        for t in range(data4d.shape[0]):
            validity = self.validity[t]
            for k in range(data4d.shape[1]):
                for j in range(data4d.shape[2]):
                    for i in range(data4d.shape[3]):
                        kwargs_vcoord['levels'] = [levels4d[t, k, j, i]]
                        vcoordinate = vcoord_builer(**copy.deepcopy(kwargs_vcoord))
                        geometry = geom_builder(structure='Point',
                                                dimensions={'X':1, 'Y':1},
                                                vcoordinate=vcoordinate,
                                                grid={'longitudes':[lons[j, i]],
                                                      'latitudes':[lats[j, i]]},
                                                position_on_horizontal_grid='center'
                                                )
                        pointfield = field_builder(structure='Point',
                                                   fid=dict(copy.deepcopy(self.fid)),
                                                   geometry=geometry,
                                                   validity=validity.copy())
                        pointfield.setdata(data4d[t, k, j, i])
                        result.append(pointfield)
        return result

    def as_profiles(self, subzone=None):
        """
        Export values as a fieldset of profiles.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
                     among: 'C', 'CI'.
        """

        if self.spectral:
            raise epygramError("as_profiles method needs a grid-point field, not a spectral one.")

        field_builder = fpx.field
        geom_builder = fpx.geometry
        vcoord_builer = fpx.geometry

        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone)
        data4d = self.getdata(d4=True, subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)

        result = FieldSet()
        kwargs_vcoord = copy.deepcopy(self.geometry.vcoordinate.footprint_as_dict())
        for t in range(data4d.shape[0]):
            validity = self.validity[t]
            for j in range(data4d.shape[2]):
                for i in range(data4d.shape[3]):
                    kwargs_vcoord['levels'] = [levels4d[t, :, j, i]]
                    vcoordinate = vcoord_builer(**copy.deepcopy(kwargs_vcoord))
                    geometry = geom_builder(structure='V1D',
                                            dimensions={'X':1, 'Y':1},
                                            vcoordinate=vcoordinate,
                                            grid={'longitudes':[lons[j, i]],
                                                  'latitudes':[lats[j, i]]},
                                            position_on_horizontal_grid='center'
                                            )
                    profilefield = field_builder(structure='V1D',
                                                 fid=dict(copy.deepcopy(self.fid)),
                                                 geometry=geometry,
                                                 validity=validity.copy())
                    profilefield.setdata(data4d[t, :, j, i])
                    result.append(profilefield)
        return result

    def extract_subdomain(self, geometry, interpolation='nearest',
                          external_distance=None,
                          exclude_extralevels=True):
        """
        Extracts a subdomain from a field, given a new geometry.

        Args: \n
        - *geometry* defines the geometry on which extract data
        - *interpolation* defines the interpolation function used to compute
          the profile at requested lon/lat from the fields grid:
          - if 'nearest' (default), extracts profile at the horizontal nearest
            neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear interpolation;
          - if 'cubic', computes profile with horizontal cubic interpolation.
        - *external_distance* can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = |target_value - external_field.data|
        - *exclude_extralevels* if True levels with no physical meaning are
          suppressed.
        """

        # build subdomain fid
        subdomainfid = {key:(FPDict(value)
                             if isinstance(value, dict)
                             else value)
                        for (key, value) in self.fid.items()}

        # build vertical geometry
        kwargs_vcoord = {'structure':'V',
                         'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid}
        if self.geometry.vcoordinate.typeoffirstfixedsurface == 119:
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
            kwargs_vcoord['levels'] = copy.copy(self.geometry.vcoordinate.levels)
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 118:
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
            kwargs_vcoord['levels'] = copy.copy(self.geometry.vcoordinate.levels)
            # Suppression of levels above or under physical domain
            if exclude_extralevels:
                for level in kwargs_vcoord['levels']:
                    if level < 1 or level > len(self.geometry.vcoordinate.grid['gridlevels']) - 1:
                        kwargs_vcoord['levels'].remove(level)
        elif self.geometry.vcoordinate.typeoffirstfixedsurface in [100, 103, 109, 1, 106, 255, 160, 200]:
            kwargs_vcoord['levels'] = copy.copy(self.geometry.vcoordinate.levels)
        else:
            raise NotImplementedError("type of first surface level: " + str(self.geometry.vcoordinate.typeoffirstfixedsurface))
        if geometry.vcoordinate.typeoffirstfixedsurface not in [255, kwargs_vcoord['typeoffirstfixedsurface']]:
            raise epygramError("extract_subdomain cannot change vertical coordinate.")
        if geometry.vcoordinate.position_on_grid not in [None, '__unknown__', kwargs_vcoord['position_on_grid']]:
            raise epygramError("extract_subdomain cannot change position on vertical grid.")
        if geometry.vcoordinate.grid != {} and geometry.vcoordinate.grid != kwargs_vcoord['grid']:
            # One could check if requested grid is a subsample of field grid
            raise epygramError("extract_subdomain cannot change vertical grid")
        if geometry.vcoordinate.levels != []:
            for level in geometry.vcoordinate.levels:
                if level not in kwargs_vcoord['levels']:
                    raise epygramError("extract_subdomain cannot do vertical interpolations.")
            kwargs_vcoord['levels'] = geometry.vcoordinate.levels
        vcoordinate = fpx.geometry(**kwargs_vcoord)
        # build geometry
        kwargs_geom = {'structure': geometry.structure,
                       'name': geometry.name,
                       'grid': dict(geometry.grid),  # do not remove dict(), it is usefull for unstructured grid
                       'dimensions': copy.copy(geometry.dimensions),
                       'vcoordinate': vcoordinate,
                       'position_on_horizontal_grid': 'center'}
        if geometry.projected_geometry:
            kwargs_geom['projection'] = copy.copy(geometry.projection)
            kwargs_geom['projtool'] = geometry.projtool
            kwargs_geom['geoid'] = geometry.geoid
        if geometry.position_on_horizontal_grid not in [None, '__unknown__', kwargs_geom['position_on_horizontal_grid']]:
            raise epygramError("extract_subdomain cannot deal with position_on_horizontal_grid other than 'center'")
        newgeometry = fpx.geometry(**kwargs_geom)

        # location & interpolation
        lons, lats = newgeometry.get_lonlat_grid()
        if len(lons.shape) > 1:
            lons = stretch_array(lons)
            lats = stretch_array(lats)
        for (lon, lat) in numpy.nditer([lons, lats]):
            if not self.geometry.point_is_inside_domain_ll(lon, lat):
                raise ValueError("point (" + str(lon) + ", " +
                                 str(lat) + ") is out of field domain.")
        comment = None
        if interpolation == 'nearest':
            if lons.size == 1:
                true_loc = self.geometry.ij2ll(*self.geometry.nearest_points(lons, lats,
                                                                             {'n':'1'},
                                                                             external_distance=external_distance))
                distance = self.geometry.distance((float(lons), float(lats)),
                                                  (float(true_loc[0]),
                                                   float(true_loc[1])))
                az = self.geometry.azimuth((float(lons), float(lats)),
                                           (float(true_loc[0]),
                                            float(true_loc[1])))
                if -22.5 < az <= 22.5:
                    direction = 'N'
                elif -77.5 < az <= -22.5:
                    direction = 'NW'
                elif -112.5 < az <= -77.5:
                    direction = 'W'
                elif -157.5 < az <= -112.5:
                    direction = 'SW'
                elif -180.5 <= az <= -157.5 or 157.5 < az <= 180.:
                    direction = 'S'
                elif 22.5 < az <= 77.5:
                    direction = 'NE'
                elif 77.5 < az <= 112.5:
                    direction = 'E'
                elif 112.5 < az <= 157.5:
                    direction = 'SE'
                gridpointstr = "(" + \
                               '{:.{precision}{type}}'.format(float(true_loc[0]),
                                                              type='F',
                                                              precision=4) + \
                               ", " + \
                               '{:.{precision}{type}}'.format(float(true_loc[1]),
                                                              type='F',
                                                              precision=4) + \
                               ")"
                comment = "Profile @ " + str(int(distance)) + "m " + \
                          direction + " from " + str((float(lons), float(lats))) + \
                          "\n" + "( = nearest gridpoint: " + gridpointstr + ")"
        elif interpolation in ('linear', 'cubic'):
            if interpolation == 'linear':
                interpstr = 'linearly'
            elif interpolation == 'cubic':
                interpstr = 'cubically'
            if lons.size == 1:
                comment = "Profile " + interpstr + " interpolated @ " + \
                          str((float(lons), float(lats)))
        else:
            raise NotImplementedError(interpolation + "interpolation.")

        # Values
        shp = newgeometry.get_datashape(dimT=len(self.validity), d4=True)
        data = numpy.ndarray(shp)
        for t in range(len(self.validity)):
            for k in range(len(newgeometry.vcoordinate.levels)):
                level = newgeometry.vcoordinate.levels[k]
                extracted = self.getvalue_ll(lons, lats, level, self.validity[t],
                                             interpolation=interpolation,
                                             external_distance=external_distance,
                                             one=False)
                data[t, k, :, :] = newgeometry.reshape_data(extracted)

        # Field
        newfield = fpx.field(fid=FPDict(subdomainfid),
                             structure=newgeometry.structure,
                             geometry=newgeometry,
                             validity=self.validity,
                             processtype=self.processtype,
                             comment=comment)
        newfield.setdata(data)

        return newfield

    def extract_zoom(self, zoom):
        """
        Extract an unstructured field with the gridpoints contained in *zoom*,
        *zoom* being a dict(lonmin=, lonmax=, latmin=, latmax=).
        """

        assert not self.spectral, \
               "spectral field: convert to gridpoint beforehand"

        (lons, lats) = self.geometry.get_lonlat_grid()
        kwargs_zoomgeom = {'structure':self.geometry.structure,
                           'vcoordinate':self.geometry.vcoordinate,
                           'position_on_horizontal_grid':self.geometry.position_on_horizontal_grid,
                           'geoid':self.geometry.geoid}
        if self.geometry.name == 'regular_lonlat':
            imin, jmin = self.geometry.ll2ij(zoom['lonmin'], zoom['latmin'])
            imax, jmax = self.geometry.ll2ij(zoom['lonmax'], zoom['latmax'])
            imin = int(numpy.ceil(imin))
            imax = int(numpy.floor(imax))
            jmin = int(numpy.ceil(jmin))
            jmax = int(numpy.floor(jmax))
            kwargs_zoomgeom['dimensions'] = {'X':imax - imin + 1,
                                             'Y':jmax - jmin + 1}
            lonmin, latmin = self.geometry.ij2ll(imin, jmin)
            kwargs_zoomgeom['name'] = self.geometry.name
            kwargs_zoomgeom['grid'] = {'input_position':(0, 0),
                                       'input_lon':Angle(lonmin, 'degrees'),
                                       'input_lat':Angle(latmin, 'degrees'),
                                       'X_resolution':self.geometry.grid['X_resolution'],
                                       'Y_resolution':self.geometry.grid['Y_resolution']}
        else:
            lons = lons.flatten()
            lats = lats.flatten()
            zoomlons = FPList([])
            zoomlats = FPList([])
            flat_indexes = []
            for i in range(len(lons)):
                if zoom['lonmin'] <= lons[i] <= zoom['lonmax'] and \
                   zoom['latmin'] <= lats[i] <= zoom['latmax']:
                    zoomlons.append(lons[i])
                    zoomlats.append(lats[i])
                    flat_indexes.append(i)
            assert len(zoomlons) > 0, "zoom not in domain."
            kwargs_zoomgeom['dimensions'] = {'X':len(zoomlons),
                                             'Y':1}
            kwargs_zoomgeom['name'] = 'unstructured'
            kwargs_zoomgeom['grid'] = {'longitudes':zoomlons,
                                       'latitudes':zoomlats}
        zoom_geom = fpx.geometry(**kwargs_zoomgeom)

        # Serait plus élégant mais pb d'efficacité (2x plus lent)
        # car extract_subdomain fait une recherche des plus proches points:
        # zoom_field = self.extract_subdomain(zoom_geom)
        # zoom_field.fid = fid
        shp = zoom_geom.get_datashape(dimT=len(self.validity), d4=True)
        data = numpy.empty(shp)
        values = self.getdata(d4=True)
        for t in range(len(self.validity)):
            for k in range(len(self.geometry.vcoordinate.levels)):
                if self.geometry.name == 'regular_lonlat':
                    data[t, k, :, :] = values[t, k, jmin:jmax + 1, imin:imax + 1]
                else:
                    vals = values[t, k, :, :].flatten()
                    zoomvals = []
                    for i in flat_indexes:
                        zoomvals.append(vals[i])
                    data[t, k, :, :] = numpy.array(zoomvals).reshape(shp[2:])

        fid = {k:v for k, v in self.fid.items()}
        for k, v in fid.items():
            if isinstance(v, dict):
                fid[k] = FPDict(v)
        zoom_field = fpx.field(fid=fid,
                                            structure=self.structure,
                                            geometry=zoom_geom,
                                            validity=self.validity,
                                            spectral_geometry=None,
                                            processtype=self.processtype)
        zoom_field.setdata(data)

        return zoom_field

    def extract_subarray(self,
                         first_i, last_i,
                         first_j, last_j):
        """
        Extract a rectangular sub-array from the field, given the i,j index limits
        of the sub-array, and return the extracted field.
        """

        newgeom = self.geometry.make_subarray_geometry(first_i, last_i,
                                                       first_j, last_j)
        # select data
        subdata = self.getdata(d4=True)[:, :, first_j:last_j, first_i:last_i]
        # copy the field object and set the new geometry, then data
        field_kwargs = copy.deepcopy(self._attributes)
        field_kwargs['geometry'] = newgeom
        newfield = fpx.field(**field_kwargs)
        newfield.setdata(subdata)

        return newfield

    def extend(self, another_field_with_time_dimension):
        """
        Extend the field with regard to time dimension with the field given as
        argument.

        Be careful no check is done for consistency between the two fields
        geometry (except that dimensions match) nor their validities.
        """

        another = another_field_with_time_dimension
        d1 = self.getdata(d4=True)
        d2 = another.getdata(d4=True)
        d = numpy.concatenate([d1, d2], axis=0)
        self.validity.extend(another.validity)
        self.setdata(d)

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def plotfield(self):
        raise NotImplementedError("plot of 3D field is not implemented")

    def stats(self, subzone=None):
        """
        Computes some basic statistics on the field, as a dict containing:
        {'min', 'max', 'mean', 'std', 'quadmean', 'nonzero'}.

        See each of these methods for details.

        - *subzone*: optional, among ('C', 'CI'), for LAM fields only, plots
          the data resp. on the C or C+I zone. \n
          Default is no subzone, i.e. the whole field.
        """

        return {'min':self.min(subzone=subzone),
                'max':self.max(subzone=subzone),
                'mean':self.mean(subzone=subzone),
                'std':self.std(subzone=subzone),
                'quadmean':self.quadmean(subzone=subzone),
                'nonzero':self.nonzero(subzone=subzone)}

    def min(self, subzone=None):
        """Returns the minimum value of data."""
        return super(D3CommonField, self).min(subzone=subzone)

    def max(self, subzone=None):
        """Returns the maximum value of data."""
        return super(D3CommonField, self).max(subzone=subzone)

    def mean(self, subzone=None):
        """Returns the mean value of data."""
        return super(D3CommonField, self).mean(subzone=subzone)

    def std(self, subzone=None):
        """Returns the standard deviation of data."""
        return super(D3CommonField, self).std(subzone=subzone)

    def quadmean(self, subzone=None):
        """Returns the quadratic mean of data."""
        return super(D3CommonField, self).quadmean(subzone=subzone)

    def nonzero(self, subzone=None):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        return super(D3CommonField, self).nonzero(subzone=subzone)

    def dctspectrum(self, level_index=None, validity_index=None, subzone=None):
        """
        Returns the DCT spectrum of the field, as a
        :class:`epygram.spectra.Spectrum` instance.
        *k* is the level index to use to compute de DCT
        """
        import epygram.spectra as esp

        if level_index is None:
            if self.geometry.datashape['k']:
                raise epygramError("must provide *level_index* for a 3D field.")
            else:
                level_index = 0
        elif level_index > 0 and not self.geometry.datashape['k']:
            raise epygramError("invalid *level_index* with regard to vertical levels.")
        if validity_index is None:
            if len(self.validity) > 1:
                raise epygramError("must provide *validity_index* for a time-dimensioned field.")
            else:
                validity_index = 0
        elif validity_index > 0 and len(self.validity) == 1:
            raise epygramError("invalid *validity_index* with regard to time dimension.")

        if self.geometry.datashape['k']:
            field2d = self.getlevel(k=level_index)
        else:
            field2d = self
        if len(self.validity) > 1:
            field2dt0 = field2d.getvalidity(validity_index)
        else:
            field2dt0 = field2d
        variances = esp.dctspectrum(field2dt0.getdata(subzone=subzone))
        spectrum = esp.Spectrum(variances[1:],
                                name=str(self.fid),
                                resolution=self.geometry.grid['X_resolution'] / 1000.,
                                mean2=variances[0])

        return spectrum

    def global_shift_center(self, longitude_shift):
        """
        For global RegLLGeometry grids only !
        Shifts the center of the geometry (and the data accordingly) by
        *longitude_shift* (in degrees). *longitude_shift* has to be a multiple
        of the grid's resolution in longitude.
        """

        if self.geometry.name != 'regular_lonlat':
            raise epygramError("only for regular lonlat geometries.")
        self.geometry.global_shift_center(longitude_shift)
        n = int(longitude_shift / self.geometry.grid['X_resolution'].get('degrees'))
        data = self.getdata(d4=True)
        data[:, :, :, :] = numpy.concatenate((data[:, :, :, n:], data[:, :, :, 0:n]),
                                             axis=3)
        self.setdata(data)

    def what(self, out=sys.stdout,
             validity=True,
             vertical_geometry=True,
             cumulativeduration=True,
             arpifs_var_names=False,
             fid=True):
        """
        Writes in file a summary of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *vertical_geometry*: if True, writes the validity of the
          field.
        - *vertical_geometry*: if True, writes the vertical geometry of the
          field.
        - *cumulativeduration*: if False, not written.
        - *arpifs_var_names*: if True, prints the equivalent 'arpifs' variable
          names.
        - *fid*: if True, prints the fid.
        """

        if self.spectral:
            spectral_geometry = self.spectral_geometry.truncation
        else:
            spectral_geometry = None
        write_formatted(out, "Kind of producting process", self.processtype)
        if validity:
            self.validity.what(out, cumulativeduration=cumulativeduration)
        self.geometry.what(out,
                           vertical_geometry=vertical_geometry,
                           arpifs_var_names=arpifs_var_names,
                           spectral_geometry=spectral_geometry)
        if fid:
            for key in self.fid:
                write_formatted(out, "fid " + key, self.fid[key])

    def dump_to_nc(self, filename, variablename=None, fidkey=None):
        """
        Dumps the field in a netCDF file named *filename*, with variable
        name being either *variablename* or self.fid[*fidkey*] or
        self.fid['netCDF'] if existing.
        """
        from epygram.formats import resource

        _fid = self.fid.get('netCDF')
        assert None in (variablename, fidkey), \
               "only one of *variablename*, *fidkey* can be provided."
        if variablename is not None:
            self.fid['netCDF'] = variablename
        elif fidkey is not None:
            self.fid['netCDF'] = self.fid[fidkey]
        else:
            assert 'netCDF' in self.fid, \
                   ' '.join(["must provide *variablename* or *fidkey* for",
                             "determining variable name in netCDF."])
        with resource(filename, 'w', fmt='netCDF') as r:
            r.writefield(self)
        if _fid is not None:
            self.fid['netCDF'] = _fid

#############
# OPERATORS #
#############

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """

        if isinstance(other, self.__class__):
            assert self.spectral == other.spectral, \
                   "cannot operate a spectral field with a non-spectral field."
            assert self.geometry.dimensions == other.geometry.dimensions, \
                   ' '.join(["operations on fields cannot be done if fields do",
                             "not share their gridpoint dimensions."])
            assert self.spectral_geometry == other.spectral_geometry, \
                   ' '.join(["operations on fields cannot be done if fields do",
                             "not share their spectral geometry."])
            assert len(self.validity) == len(other.validity), \
                   ' '.join(["operations on fields cannot be done if fields do",
                             "not share their time dimension."])
        else:
            super(D3CommonField, self)._check_operands(other)

    def __add__(self, other):
        """
        Definition of addition, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'+'} and null validity.
        """

        newfield = self._add(other,
                             structure=self.structure,
                             geometry=self.geometry,
                             spectral_geometry=self.spectral_geometry,
                             validity=FieldValidityList(length=len(self.validity)))
        return newfield

    def __mul__(self, other):
        """
        Definition of multiplication, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'*'} and null validity.
        """

        newfield = self._mul(other,
                             structure=self.structure,
                             geometry=self.geometry,
                             spectral_geometry=self.spectral_geometry,
                             validity=FieldValidityList(length=len(self.validity)))
        return newfield

    def __sub__(self, other):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """

        newfield = self._sub(other,
                             structure=self.structure,
                             geometry=self.geometry,
                             spectral_geometry=self.spectral_geometry,
                             validity=FieldValidityList(length=len(self.validity)))
        return newfield

    def __div__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """

        newfield = self._div(other,
                             structure=self.structure,
                             geometry=self.geometry,
                             spectral_geometry=self.spectral_geometry,
                             validity=FieldValidityList(length=len(self.validity)))
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

        newfield = self._rsub(other,
                              structure=self.structure,
                              geometry=self.geometry,
                              spectral_geometry=self.spectral_geometry,
                              validity=FieldValidityList(length=len(self.validity)))
        return newfield

    def __rdiv__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """

        newfield = self._rdiv(other,
                              structure=self.structure,
                              geometry=self.geometry,
                              spectral_geometry=self.spectral_geometry,
                              validity=FieldValidityList(length=len(self.validity)))
        return newfield


class D3Field(D3CommonField):
    """
    3-Dimensions field class.
    A field is defined by its identifier 'fid',
    its data, its geometry (gridpoint and optionally spectral),
    and its validity.

    The natural being of a field is gridpoint, so that:
    a field always has a gridpoint geometry, but it has a spectral geometry only
    in case it is spectral.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['3D'])),
            geometry=dict(
                type=D3Geometry,
                info="Geometry defining the position of the field gridpoints."),
            validity=dict(
                type=FieldValidityList,
                info="Validity of the field.",
                optional=True,
                access='rwx',
                default=FieldValidityList()),
            spectral_geometry=dict(
                info="For a spectral field, its spectral geometry handles \
                      spectral transforms and dimensions.",
                type=SpectralGeometry,
                optional=True),
            processtype=dict(
                optional=True,
                info="Generating process.")
        )
    )


##############
# ABOUT DATA #
##############

    def sp2gp(self):
        """
        Transforms the spectral field into gridpoint, according to its spectral
        geometry. Replaces data in place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *sp2gp()* method.
        """

        if self.spectral:
            gpdims = self._get_gpdims_for_spectral_transforms()
            if self.geometry.rectangular_grid:
                # LAM
                gpdata = numpy.empty(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
            else:
                # global
                gpdata = numpy.ma.zeros(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
            for t in range(len(self.validity)):
                for k in range(len(self.geometry.vcoordinate.levels)):
                    spdata_i = self.getdata(d4=True)[t, k, :]
                    gpdata_i = self.spectral_geometry.sp2gp(spdata_i, gpdims)
                    gpdata_i = self.geometry.reshape_data(gpdata_i)
                    gpdata[t, k, :, :] = gpdata_i[:, :]

            self._attributes['spectral_geometry'] = None
            self.setdata(gpdata)

    def gp2sp(self, spectral_geometry):
        """
        Transforms the gridpoint field into spectral space, according to the
        *spectral_geometry* mandatorily passed as argument. Replaces data in
        place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *gp2sp()* method.
        """

        assert isinstance(spectral_geometry, SpectralGeometry)

        if not self.spectral:
            gpdims = self._get_gpdims_for_spectral_transforms()
            spdata = None
            for t in range(len(self.validity)):
                for k in range(len(self.geometry.vcoordinate.levels)):
                    gpdata_i = stretch_array(self.getdata(d4=True)[t, k, :, :])
                    spdata_i = spectral_geometry.gp2sp(gpdata_i, gpdims)
                    n = len(spdata_i)
                    if spdata is None:
                        spdata = numpy.empty((len(self.validity),
                                              len(self.geometry.vcoordinate.levels),
                                              n))
                    spdata[t, k, :] = spdata_i[:]

            self._attributes['spectral_geometry'] = spectral_geometry
            self.setdata(spdata)

    def _get_gpdims_for_spectral_transforms(self):
        """
        Build a dictionary containing gridpoint dimensions for the call to
        spectral transforms.
        """

        if self.geometry.rectangular_grid:
            # LAM
            gpdims = {}
            for dim in ['X', 'Y', 'X_CIzone', 'Y_CIzone']:
                gpdims[dim] = self.geometry.dimensions[dim]
            for item in ['X_resolution', 'Y_resolution']:
                gpdims[item] = self.geometry.grid[item]
        else:
            # global
            gpdims = {}
            for dim in ['lat_number', 'lon_number_by_lat']:
                gpdims[dim] = self.geometry.dimensions[dim]

        return gpdims

    def compute_xy_spderivatives(self):
        """
        Compute the derivatives of field in spectral space, then come back in
        gridpoint space.
        Returns the two derivative fields.

        The spectral transform and derivatives subroutines are actually included
        in the spectral geometry's *compute_xy_spderivatives()* method.
        """

        if self.spectral:
            gpdims = self._get_gpdims_for_spectral_transforms()
            if self.geometry.rectangular_grid:
                # LAM
                gpderivX = numpy.empty(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
                gpderivY = numpy.empty(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
            else:
                # global
                gpderivX = numpy.ma.zeros(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
                gpderivY = numpy.ma.zeros(self.geometry.get_datashape(dimT=len(self.validity), d4=True))
            for t in range(len(self.validity)):
                for k in range(len(self.geometry.vcoordinate.levels)):
                    spdata_i = self.getdata(d4=True)[t, k, :]
                    (dx, dy) = self.spectral_geometry.compute_xy_spderivatives(spdata_i, gpdims)
                    dx = self.geometry.reshape_data(dx)
                    dy = self.geometry.reshape_data(dy)
                    gpderivX[t, k, :, :] = dx[:, :]
                    gpderivY[t, k, :, :] = dy[:, :]

            field_dX = copy.deepcopy(self)
            field_dY = copy.deepcopy(self)
            for field in (field_dX, field_dY):
                field._attributes['spectral_geometry'] = None
                field.fid = {'derivative':'x'}
            field_dX.setdata(gpderivX)
            field_dY.setdata(gpderivY)
        else:
            raise epygramError('field must be spectral to compute its spectral derivatives.')

        return (field_dX, field_dY)

    def getdata(self, subzone=None, d4=False):
        """
        Returns the field data, with 3D shape if the field is not spectral,
        2D if spectral.

        - *subzone*: optional, among ('C', 'CI'), for LAM fields only, returns
          the data resp. on the C or C+I zone.
          Default is no subzone, i.e. the whole field.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry

        Shape of 4D data: \n
        - Rectangular grids:\n
          grid[t,k,0,0] is SW, grid[t,k,-1,-1] is NE \n
          grid[t,k,0,-1] is SE, grid[t,k,-1,0] is NW \n
          with k the level, t the temporal dimension
        - Gauss grids:\n
          grid[t,k,0,:Nj] is first (Northern) band of latitude, masked after
          Nj = number of longitudes for latitude j \n
          grid[t,k,-1,:Nj] is last (Southern) band of latitude (idem). \n
          with k the level, t the temporal dimension
        """

        data = self._data
        if not self.spectral and subzone is not None:
            if self.geometry.grid.get('LAMzone') is not None:
                data = self.geometry.extract_subzone(data, subzone)
            else:
                raise epygramError("*subzone* cannot be provided for this field.")
        if not d4:
            data = data.squeeze()

        return data

    def setdata(self, data):
        """
        Sets field data, checking *data* to have the good shape according to
        geometry.
        *data* dimensions should in any case be ordered in a subset of
        (t,z,y,x), or (t,z,n) if spectral (2D spectral coefficients must be 1D
        with ad hoc ordering, n being the total number of spectral
        coefficients).

        *data* may be 4D (3D if spectral) even if the field is not, as long as
        the above dimensions ordering is respected.
        """

        if not isinstance(data, numpy.ndarray):
            data = numpy.array(data)

        if self.spectral:
            shp = (len(self.validity),
                   len(self.geometry.vcoordinate.levels),
                   data.shape[-1])
        elif 'gauss' in self.geometry.name:
            shp = (len(self.validity),
                   len(self.geometry.vcoordinate.levels),
                   self.geometry.dimensions['lat_number'],
                   self.geometry.dimensions['max_lon_number'])
        else:
            shp = (len(self.validity),
                   len(self.geometry.vcoordinate.levels),
                   self.geometry.dimensions['Y'],
                   self.geometry.dimensions['X'])

        if self.spectral:
            d4 = len(data.shape) == 3
        else:
            d4 = len(data.shape) == 4
        if d4:
            assert data.shape == shp, \
                   ' '.join(['data', str(data.shape),
                             'should have shape', str(shp)])
        else:
            # find indexes corresponding to dimensions
            dimensions = 0
            indexes = {'t':0, 'z':1, 'y':2, 'x':3}
            # t, z
            if len(self.validity) > 1:
                dimensions += 1
            else:
                indexes['t'] = None
                for i in ('z', 'y', 'x'):
                    indexes[i] = indexes[i] - 1
            if self.geometry.datashape['k']:
                dimensions += 1
            else:
                indexes['z'] = None
                for i in ('y', 'x'):
                    indexes[i] = indexes[i] - 1
            # y, x or spectral ordering
            if self.spectral:
                dimensions += 1
                dataType = "spectral"
            else:
                if self.geometry.datashape['j']:
                    dimensions += 1
                else:
                    indexes['y'] = None
                    for i in ('x',):
                        indexes[i] = indexes[i] - 1
                if self.geometry.datashape['i']:
                    dimensions += 1
                else:
                    indexes['x'] = None
                dataType = "gridpoint"
            # check dimensions
            assert len(numpy.shape(data)) == dimensions \
                   or numpy.shape(data) == (1,), \
                   dataType + " data should be " + str(dimensions) + "D array."
            if indexes['t'] is not None:
                assert data.shape[0] == len(self.validity), \
                       ' == '.join(['data.shape[0] should be len(self.validity)',
                                    str(len(self.validity))])
            if self.geometry.datashape['k']:
                assert data.shape[indexes['z']] == len(self.geometry.vcoordinate.levels), \
                       ' == '.join(['data.shape[' + str(indexes['z']) +
                                    '] should be len(self.geometry.vcoordinate.levels)',
                                    str(len(self.geometry.vcoordinate.levels))])
            if not self.spectral:
                if 'gauss' in self.geometry.name:
                    if self.geometry.datashape['j']:
                        assert data.shape[indexes['y']] == self.geometry.dimensions['lat_number'], \
                               ' == '.join(['data.shape[' + str(indexes['y']) +
                                            "] should be self.geometry.dimensions['lat_number']",
                                            str(self.geometry.dimensions['lat_number'])])
                    if self.geometry.datashape['i']:
                        assert data.shape[indexes['x']] == self.geometry.dimensions['max_lon_number'], \
                               ' == '.join(['data.shape[' + str(indexes['x']) +
                                            "] should be self.geometry.dimensions['max_lon_number']",
                                            str(self.geometry.dimensions['max_lon_number'])])
                else:
                    if self.geometry.datashape['j']:
                        assert data.shape[indexes['y']] == self.geometry.dimensions['Y'], \
                               ' == '.join(['data.shape[' + str(indexes['y']) +
                                            "] should be self.geometry.dimensions['Y']",
                                            str(self.geometry.dimensions['Y'])])
                    if self.geometry.datashape['i']:
                        assert data.shape[indexes['x']] == self.geometry.dimensions['X'], \
                               ' == '.join(['data.shape[' + str(indexes['x']) +
                                            "] should be self.geometry.dimensions['X']",
                                            str(self.geometry.dimensions['X'])])
            # reshape to 4D
            data = data.reshape(shp)
        super(D3Field, self).setdata(data)

    data = property(getdata, setdata, Field.deldata, "Accessor to the field data.")

    def select_subzone(self, subzone):
        """
        If a LAMzone defines the field, select only the *subzone* from it.
        *subzone* among ('C', 'CI').
        Warning: modifies the field and its geometry in place !
        """

        if self.geometry.grid.get('LAMzone') is not None:
            data = self.getdata(subzone=subzone)
            self._attributes['geometry'] = self.geometry.select_subzone(subzone)
            self.setdata(data)

    def getvalue_ij(self, i=None, j=None, k=None, t=None,
                    one=True):
        """
        Returns the value of the field on point of indices (*i, j, k, t*).
        Take care (*i, j, k, t*) is python-indexing, ranging from 0 to dimension - 1.
        *k* is the index of the level (not a value in Pa or m...)
        *t* is the index of the temporal dimension (not a validity object)
        *k* and *t* can be scalar even if *i* and *j* are arrays.

        If *one* is False, returns [value] instead of value.
        """

        if len(self.validity) > 1 and t is None:
            raise epygramError("*t* is mandatory when there are several validities")
        if self.geometry.datashape['k'] and k is None:
            raise epygramError("*k* is mandatory when field has a vertical coordinate")
        if self.geometry.datashape['j'] and j is None:
            raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
        if self.geometry.datashape['i'] and j is None:
            raise epygramError("*i* is mandatory when field has one horizontal dimension")

        if not self.geometry.point_is_inside_domain_ij(i, j):
            raise ValueError("point is out of field domain.")

        maxsize = numpy.array([numpy.array(dim).size for dim in [i, j, k, t] if dim is not None]).max()
        if t is None:
            my_t = numpy.zeros(maxsize, dtype=int)
        else:
            my_t = numpy.array(t, dtype=int)
            if my_t.size != maxsize:
                if my_t.size != 1:
                    raise epygramError("t must be scalar or must have the same length as other indexes")
                my_t = numpy.array([my_t.item()] * maxsize, dtype=int)
        if k is None:
            my_k = numpy.zeros(maxsize, dtype=int)
        else:
            my_k = numpy.array(k, dtype=int)
            if my_k.size != maxsize:
                if my_k.size != 1:
                    raise epygramError("k must be scalar or must have the same length as other indexes")
                my_k = numpy.array([my_k.item()] * maxsize, dtype=int)
        if j is None:
            my_j = numpy.zeros(maxsize, dtype=int)
        else:
            my_j = numpy.array(j, dtype=int)
            if my_j.size != maxsize:
                raise epygramError("j must have the same length as other indexes")
        if i is None:
            my_i = numpy.zeros(maxsize, dtype=int)
        else:
            my_i = numpy.array(i, dtype=int)
            if my_i.size != maxsize:
                raise epygramError("i must have the same length as other indexes")

        value = numpy.copy(self.getdata(d4=True)[my_t, my_k, my_j, my_i])
        if value.size == 1 and one:
            value = value.item()
        return value

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.
        *level* is the requested level expressed in coordinate value (Pa, m...)
        *k* is the index of the requested level
        """

        if k is None and level is None:
            raise epygramError("You must give k or level.")
        if k is not None and level is not None:
            raise epygramError("You cannot give, at the same time, k and level")
        if level is not None:
            if level not in self.geometry.vcoordinate.levels:
                raise epygramError("The requested level does not exist.")
            my_k = self.geometry.vcoordinate.levels.index(level)
        else:
            my_k = k

        if self.structure == '3D':
            newstructure = 'H2D'
        elif self.structure == 'V2D':
            newstructure = 'H2D'
        elif self.structure == 'V1D':
            newstructure = 'Point'
        else:
            raise epygramError("It's not possible to extract a level from a " + self.structure + " field.")

        kwargs_vcoord = {'structure': 'V',
                         'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid,
                         'levels':[level]}
        if self.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
        newvcoordinate = fpx.geometry(**kwargs_vcoord)
        kwargs_geom = {'structure':newstructure,
                       'name': self.geometry.name,
                       'grid': dict(self.geometry.grid),
                       'dimensions': copy.copy(self.geometry.dimensions),
                       'vcoordinate': newvcoordinate,
                       'position_on_horizontal_grid': self.geometry.position_on_horizontal_grid}
        if self.geometry.projected_geometry:
            kwargs_geom['projection'] = copy.copy(self.geometry.projection)
            kwargs_geom['projtool'] = self.geometry.projtool
            kwargs_geom['geoid'] = self.geometry.geoid
        newgeometry = fpx.geometry(**kwargs_geom)
        generic_fid = self.fid.get('generic', {})
        generic_fid['level'] = level
        kwargs_field = {'structure':newstructure,
                        'validity':self.validity.copy(),
                        'processtype':self.processtype,
                        'geometry':newgeometry,
                        'fid':{'generic':generic_fid}}
        if self.spectral_geometry is not None:
            kwargs_field['spectral_geometry'] = self.spectral_geometry.copy()
        newfield = fpx.field(**kwargs_field)
        newfield.setdata(self.getdata(d4=True)[:, my_k:my_k + 1, :, :])

        return newfield

    def getvalidity(self, index_or_validity):
        """
        Returns the field restrained to one of its temporal validity as a new
        field.

        *index_or_validity* can be either a :class:`epygram.base.FieldValidity`
        instance or the index of the requested validity in the field's
        FieldValidityList.
        """

        assert isinstance(index_or_validity, FieldValidity) or isinstance(index_or_validity, int), \
               "index_or_validity* should be either a FieldValidity instance or an int."
        if isinstance(index_or_validity, FieldValidity):
            validity = index_or_validity
            try:
                index = self.validity.index(index_or_validity)
            except ValueError:
                raise ValueError('*validity* not in field validity.')
        else:
            index = index_or_validity
            validity = self.validity[index]
        newfield = self.deepcopy()
        newfield.validity = validity
        newfield.setdata(self.getdata(d4=True)[index:index + 1, :, :, ])

        return newfield


class D3VirtualField(D3CommonField):
    """
    3-Dimensions Virtual field class.

    Data is taken from other fields, either:
    - a given *fieldset*
    - a *resource* in which are stored fields defined by *resource_fids*.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['3D'])),
            fieldset=dict(
                info="Set of real fields that can compose the Virtual Field.",
                type=FieldSet,
                optional=True,
                default=FieldSet()),
            resource=dict(
                info="Resource in which is stored the fields defined by \
                      resource_fids.",
                type=Resource,
                optional=True),
            resource_fids=dict(
                info="Definition of the fields in resource that compose the \
                      virtual field.",
                type=FPList,
                optional=True,
                default=FPList()),
        )
    )

    def __init__(self, *args, **kwargs):
        """
        Constructor. See its footprint for arguments.
        """

        super(D3VirtualField, self).__init__(*args, **kwargs)

        if self.fieldset != FieldSet():
            if self.resource is not None or self.resource_fids != FPList():
                raise epygramError("You cannot set fieldset and (resource or resource_fids) at the same time.")
            raise NotImplementedError("D3VirtualField from a fieldset is not yet implemented")
        else:
            if self.resource is None or self.resource_fids == FPList():
                raise epygramError("If you do not set fieldset, you need to provide resource and resource_fids.")
            fidlist = self.resource.find_fields_in_resource(seed=self.resource_fids,
                                                            fieldtype=['H2D', '3D'])
            if len(fidlist) == 0:
                raise epygramError("There is no field in resource matching with resource_fids")
            first = True
            self._fidList = []
            levelList = []
            for fid in fidlist:
                field = self.resource.readfield(fid, getdata=False)
                if field.structure != 'H2D':
                    raise epygramError("3D virtual fields must be build from H2D fields only")
                if first:
                    self._geometry = field.geometry.copy()
                    self._validity = field.validity.copy()
                    if field.spectral_geometry is not None:
                        self._spectral_geometry = field.spectral_geometry.copy()
                    else:
                        self._spectral_geometry = None
                    self._processtype = field.processtype
                else:
                    if self._geometry.structure != field.geometry.structure or \
                       self._geometry.name != field.geometry.name or \
                       self._geometry.grid != field.geometry.grid or \
                       self._geometry.dimensions != field.geometry.dimensions or \
                       self._geometry.position_on_horizontal_grid != field.geometry.position_on_horizontal_grid:
                        raise epygramError("All H2D fields must share the horizontal geometry")
                    if self._geometry.projected_geometry or field.geometry.projected_geometry:
                        if self._geometry.projection != field.geometry.projection or \
                           self._geometry.geoid != field.geometry.geoid:
                            raise epygramError("All H2D fields must share the geometry projection")
                    if self._geometry.vcoordinate.typeoffirstfixedsurface != field.geometry.vcoordinate.typeoffirstfixedsurface or \
                       self._geometry.vcoordinate.position_on_grid != field.geometry.vcoordinate.position_on_grid:
                        raise epygramError("All H2D fields must share the vertical geometry")
                    if self._geometry.vcoordinate.grid is not None or field.geometry.vcoordinate.grid is not None:
                        if self._geometry.vcoordinate.grid != field.geometry.vcoordinate.grid:
                            raise epygramError("All H2D fields must share the vertical grid")
                    if self._validity != field.validity:
                        raise epygramError("All H2D fields must share the validity")
                    if self._spectral_geometry != field.spectral_geometry:
                        raise epygramError("All H2D fields must share the spectral geometry")
                    if self._processtype != field.processtype:
                        raise epygramError("All H2D fields must share the sprocesstype")
                if len(field.geometry.vcoordinate.levels) != 1:
                    raise epygramError("H2D fields must have only one level")
                if field.geometry.vcoordinate.levels[0] in levelList:
                    raise epygramError("This level have already been found")
                levelList.append(field.geometry.vcoordinate.levels[0])
                self._fidList.append(fid)
            kwargs_vcoord = dict(structure='V',
                                 typeoffirstfixedsurface=self._geometry.vcoordinate.typeoffirstfixedsurface,
                                 position_on_grid=self._geometry.vcoordinate.position_on_grid)
            if self._geometry.vcoordinate.grid is not None:
                kwargs_vcoord['grid'] = self._geometry.vcoordinate.grid
            kwargs_vcoord['levels'], self._fidList = (list(t) for t in zip(*sorted(zip(levelList, self._fidList))))  # TOBECHECKED
            self._geometry.vcoordinate = fpx.geometry(**kwargs_vcoord)
            self._spgpOpList = []

    @property
    def geometry(self):
        return self._geometry

    @property
    def validity(self):
        return self._validity

    @property
    def spectral_geometry(self):
        return self._spectral_geometry

    @property
    def processtype(self):
        return self._processtype


##############
# ABOUT DATA #
##############
    def sp2gp(self):
        """
        Transforms the spectral field into gridpoint, according to its spectral
        geometry. Replaces data in place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *sp2gp()* method.
        """
        self._spgpOpList.append(('sp2gp', {}))
        self._spectral_geometry = None

    def gp2sp(self, spectral_geometry):
        """
        Transforms the gridpoint field into spectral space, according to the
        *spectral geometry* mandatorily passed as argument. Replaces data in
        place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *gp2sp()* method.
        """
        self._spgpOpList.append(('gp2sp', {'spectral_geometry':spectral_geometry}))

    def getdata(self, subzone=None, d4=False):
        """
        Returns the field data, with 3D shape if the field is not spectral,
        2D if spectral.

        - *subzone*: optional, among ('C', 'CI'), for LAM fields only, returns
          the data resp. on the C or C+I zone.
          Default is no subzone, i.e. the whole field.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry

        Shape of 3D data: \n
        - Rectangular grids:\n
          grid[k,0,0] is SW, grid[k,-1,-1] is NE \n
          grid[k,0,-1] is SE, grid[k,-1,0] is NW \n
          with k the level
        - Gauss grids:\n
          grid[k,0,:Nj] is first (Northern) band of latitude, masked after
          Nj = number of longitudes for latitude j \n
          grid[k,-1,:Nj] is last (Southern) band of latitude (idem). \n
          with k the level
        """

        dataList = []
        for k in range(len(self.geometry.vcoordinate.levels)):
            dataList.append(self.getlevel(k).getdata(subzone=subzone, d4=d4))

        return numpy.array(dataList)

    def setdata(self, data):
        """setdata() not implemented on virtual fields."""
        raise epygramError("setdata cannot be implemented on virtual fields")

    def deldata(self):
        """deldata() not implemented on virtual fields."""
        raise epygramError("deldata cannot be implemented on virtual fields")

    data = property(getdata)

    def getvalue_ij(self, i=None, j=None, k=None, t=None,
                    one=True):
        """
        Returns the value of the field on point of indices (*i, j, k, t*).
        Take care (*i, j, k, t*) is python-indexing, ranging from 0 to dimension - 1.
        *k* is the index of the level (not a value in Pa or m...)
        *t* is the index of the temporal dimension (not a validity object)
        *k* and *t* can be scalar even if *i* and *j* are arrays.

        If *one* is False, returns [value] instead of value.
        """

        if len(self.validity) > 1 and t is None:
            raise epygramError("*t* is mandatory when there are several validities")
        if self.geometry.datashape['k'] and k is None:
            raise epygramError("*k* is mandatory when field has a vertical coordinate")
        if self.geometry.datashape['j'] and j is None:
            raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
        if self.geometry.datashape['i'] and j is None:
            raise epygramError("*i* is mandatory when field has one horizontal dimension")

        if not self.geometry.point_is_inside_domain_ij(i, j):
            raise ValueError("point is out of field domain.")

        maxsize = numpy.array([numpy.array(dim).size for dim in [i, j, k, t] if dim is not None]).max()
        if t is None:
            my_t = numpy.zeros(maxsize, dtype=int)
        else:
            my_t = numpy.array(t)
            if my_t.size != maxsize:
                if my_t.size != 1:
                    raise epygramError("t must be scalar or must have the same length as other indexes")
                my_t = numpy.array([my_t.item()] * maxsize)
        if k is None:
            my_k = numpy.zeros(maxsize, dtype=int)
        else:
            my_k = numpy.array(k)
            if my_k.size != maxsize:
                if my_k.size != 1:
                    raise epygramError("k must be scalar or must have the same length as other indexes")
                my_k = numpy.array([my_k.item()] * maxsize)
        if j is None:
            my_j = numpy.zeros(maxsize, dtype=int)
        else:
            my_j = numpy.array(j)
            if my_j.size != maxsize:
                raise epygramError("j must have the same length as other indexes")
        if i is None:
            my_i = numpy.zeros(maxsize, dtype=int)
        else:
            my_i = numpy.array(i)
            if my_i.size != maxsize:
                raise epygramError("i must have the same length as other indexes")

        value = []
        oldk = None
        for x in range(my_k.size):
            thisk = my_k[x] if my_k.size > 1 else my_k.item()
            if thisk != oldk:
                field2d = self.getlevel(k=thisk)
                oldk = thisk
            if my_t.size == 1:
                pos = (my_t.item(), 0, my_j.item(), my_i.item())
            else:
                pos = (my_t[x], 0, my_j[x], my_i[x])
            value.append(field2d.getdata(d4=True)[pos])
        value = numpy.array(value)
        if value.size == 1 and one:
            value = value.item()
        return value

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.
        *level* is the requested level expressed in coordinate value (Pa, m...)
        *k* is the index of the requested level
        """

        if k is None and level is None:
            raise epygramError("You must give k or level.")
        if k is not None and level is not None:
            raise epygramError("You cannot give, at the same time, k and level")
        if level is not None:
            if level not in self.geometry.vcoordinate.levels:
                raise epygramError("The requested level does not exist.")
            my_k = self.geometry.vcoordinate.levels.index(level)
        else:
            my_k = k

        result = self.resource.readfield(self._fidList[my_k])

        for op, kwargs in self._spgpOpList:
            if op == 'sp2gp':
                result.sp2gp(**kwargs)
            elif op == 'gp2sp':
                result.gp2sp(**kwargs)
            else:
                raise epygramError("operation not known")

        return result


footprints.collectors.get(tag='fields').fasttrack = ('structure',)
