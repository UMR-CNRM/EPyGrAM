#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

import copy
import numpy
import sys

import footprints
from footprints import FPDict, FPList, proxy as fpx
from bronx.graphics.axes import set_figax
from bronx.syntax.arrays import stretch_array

from epygram import epygramError, config
from epygram.util import (write_formatted, Angle,
                          degrees_nearest_mod,
                          as_numpy_array,
                          moveaxis)
from epygram.base import Field, FieldSet, FieldValidity, FieldValidityList, Resource
from epygram.geometries import (Geometry, SpectralGeometry, VGeometry, ProjectedGeometry, 
                                RectangularGridGeometry, RegLLGeometry, UnstructuredGeometry)


class _D3CommonField(Field):
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

###########################
# ABOUT VERTICAL GEOMETRY #
###########################
    def use_field_as_vcoord(self, field, force_kind=None, validity_check=True):
        """
        Use values from the field as vertical coordinates.
        One example of use: field1 is a temperature field
        extracted from a resource on hybrid levels whereas
        field2 is the pressure field expressed on the same
        grid. field1.use_field_as_vcoord(field2) transforms
        the vertical coordinate of field1 to obtain a field
        on pressure levels.

        :param field: must be a field on the same geometry and
                      validity as self that we will use to replace
                      the vertical coordinate. 'generic' fid must
                      allow to recognize the vertical coordinate
                      type (except if force_kind is set)
        :param force_kind: if not None, is used as the vertical
                           coordinate kind instead of trying to
                           guess it from the fid
        :param validity_check: if True (default) checks if validity
                               of current object and of field
                               parameter are the same
        """
        if self.geometry != field.geometry:
            raise epygramError("geometries must be identical")
        if validity_check and self.validity != field.validity:
            raise epygramError("validities must be identical, " + \
                               "or validity_check must be False")
        if field.spectral:
            raise epygramError("field must not be spectral")

        self.geometry.vcoordinate = field.as_vcoordinate(force_kind)

    def as_vcoordinate(self, force_kind=None):
        """
        Returns a VGeometry build from the values of field.
        Field must have a 'generic' fid allowing to
        recognize the vertical coordinate type (except if
        force_kind is set)

        :param force_kind: if not None, is used as the vertical
                           coordinate kind instead of trying to
                           guess it from the fid
        """
        if force_kind is None:
            if 'generic' not in self.fid:
                raise epygramError("fid must contain a 'generic' key or force_kind must be set")
            for k in ['discipline', 'parameterCategory', 'parameterNumber']:
                if k not in self.fid['generic']:
                    raise epygramError("fid['generic'] must contain the " + k + " key.")
            fid = (self.fid['generic']['discipline'],
                   self.fid['generic']['parameterCategory'],
                   self.fid['generic']['parameterNumber'])
            if fid == (0, 3, 0):
                vcoord = 100
            elif fid == (0, 3, 6):
                vcoord = 102
            else:
                raise NotImplementedError("Do not know which kind of vertical coordinate fid corresponds to")
        else:
            vcoord = force_kind

        levels = self.getdata().copy()
        if vcoord == 100:
            levels = levels / 100.
        if len(self.validity) == 1:
            # data is an array z, y, x with z being optional
            if len(self.geometry.vcoordinate.levels) == 1:
                # in case z dimension does not exist in data
                levels = [levels]
        else:
            # data is an array t, z, y, x with z being optional
            if len(self.geometry.vcoordinate.levels) == 1:
                # in case z dimension does not exist in data
                levels = [levels]
            else:
                # if z dimension exist, z and t dimensions must be exchanged
                levels = levels.swapaxes(0, 1)

        kwargs_vcoord = {'typeoffirstfixedsurface': vcoord,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid,
                         'levels': list(levels)
                         }
        return VGeometry(**kwargs_vcoord)

##############
# ABOUT DATA #
##############

    def getvalue_ll(self, lon=None, lat=None,
                    level=None,
                    k=None,
                    validity=None,
                    interpolation='nearest',
                    neighborinfo=False,
                    one=True,
                    external_distance=None):
        """
        Returns the value of the field on point of coordinates
        (*lon, lat, level*).
        *lon* and *lat* may be list/numpy.array.

        :param interpolation: if:\n
                              - 'nearest' (default), returns the value of the
                                nearest neighboring gridpoint;
                              - 'linear', computes and returns the field value
                                with bilinear or linear spline interpolation;
                              - 'cubic', computes and returns the field value
                                with cubic spline interpolation;
        :param level: is the True level not the index of the level. Depending on
                      the vertical coordinate, it could be expressed in Pa, m.
        :param k: is the index of the level.
        :param validity: is a FieldValidity or a FieldValidityList instance
        :param neighborinfo: if set to **True**, returns a tuple
                             *(value, (lon, lat))*, with *(lon, lat)* being
                             the actual coordinates of the neighboring gridpoint
                             (only for *interpolation == 'nearest'*).
        :param one: if False and len(lon) is 1, returns [value] instead of value.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = abs(target_value - external_field.data)

        When linear interpolation method is choosen, if the grid is rectangular,
        the interpolation is a bilinear one on the model grid, otherwise it
        is a spline interpolation on the lat/lon domain.

        Warning: for interpolation on Gauss geometries, requires the
        :mod:`pyproj` module.
        """
        if self.spectral:
            raise epygramError("field must be gridpoint to get value of a" +
                               " lon/lat point.")
        if len(numpy.array(self.geometry.vcoordinate.levels).shape) == 1 and \
           len(self.geometry.vcoordinate.levels) != len(set(self.geometry.vcoordinate.levels)):
            raise epygramError('Some levels are represented twice in levels list.')

        # We look for indexes for time coordinate (no interpolation)
        if validity is None:
            if len(self.validity) > 1:
                raise epygramError("*validity* is mandatory when there are several validities")
            t = 0
        else:
            if isinstance(validity, FieldValidity):
                validity = FieldValidityList(validity)
            t = [[ind for ind , val in enumerate(self.validity) if val == v][0] for v in validity]
        t = as_numpy_array(t).flatten()

        # We look for indexes for vertical coordinate (no interpolation)
        if k is not None and level is not None:
            raise epygramError("*level* and *k* cannot be different from None together")
        elif k is level is None:
            if self.geometry.datashape['k']:
                raise epygramError("*level* or *k* is mandatory when field has a vertical coordinate")
            k = 0
        elif k is None:
            # level is not None
            k = [self.geometry.vcoordinate.levels.index(l) for l in as_numpy_array(level)]
        k = as_numpy_array(k)

        if lon is None or lat is None:
            raise epygramError("*lon* and *lat* are mandatory")
        lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()
        if len(lon) != len(lat):
            raise epygramError("*lon* and *lat* must have the same length")

        sizes = set([len(x) for x in [lon, lat, k, t]])
        if len(sizes) > 2 or (len(sizes) == 2 and 1 not in sizes):
            raise epygramError("each of lon, lat, k/level and validity must be scalar or have the same length as the others")

        if interpolation == 'linear':
            if isinstance(self.geometry, RectangularGridGeometry):
                method = 'bilinear'
            else:
                method = 'linear_spline'
        else:
            method = interpolation

        if method == 'nearest':
            (ri, rj) = moveaxis(numpy.array(self.geometry.nearest_points(lon, lat, {'n':'1'},
                                            external_distance=external_distance)), 0, -1)
            value = self.getvalue_ij(ri, rj, k, t, one=one)
            if neighborinfo:
                (lon, lat) = self.geometry.ij2ll(ri, rj)
                if numpy.shape(lon) in ((1,), ()):
                    lon = float(lon)
                    lat = float(lat)
                value = (value, (lon, lat))

        elif method in ('linear_spline', 'cubic', 'bilinear'):
            lon = lon if len(lon) == max(sizes) else lon.repeat(max(sizes))
            lat = lat if len(lat) == max(sizes) else lat.repeat(max(sizes))
            k = k if len(k) == max(sizes) else k.repeat(max(sizes))
            t = t if len(t) == max(sizes) else t.repeat(max(sizes))
            value = numpy.zeros(max(sizes))
            interp_points = self.geometry.nearest_points(lon, lat,
                                                         {'linear':{'n':'2*2'},
                                                          'linear_spline':{'n':'2*2'},
                                                          'bilinear':{'n':'2*2'},
                                                          'cubic':{'n':'4*4'}}[interpolation],
                                                         squeeze=False)
            # depack
            all_i = interp_points[:, :, 0]
            all_j = interp_points[:, :, 1]
            all_k = numpy.repeat(k[:, numpy.newaxis], interp_points.shape[1], axis=1)
            all_t = numpy.repeat(t[:, numpy.newaxis], interp_points.shape[1], axis=1)

            # get values and lons/lats
            values_at_interp_points = self.getvalue_ij(all_i.flatten(),
                                                       all_j.flatten(),
                                                       all_k.flatten(),
                                                       all_t.flatten()).reshape(all_i.shape)

            if method in ('linear_spline', 'cubic'):
                from scipy.interpolate import interp1d, interp2d
                all_lons, all_lats = self.geometry.ij2ll(all_i.flatten(), all_j.flatten())
                all_lons = all_lons.reshape(all_i.shape)
                all_lats = all_lats.reshape(all_i.shape)
                for n in range(max(sizes)):
                    if self.geometry.name == 'academic' and \
                       1 in (self.geometry.dimensions['X'], self.geometry.dimensions['Y'] == 1):
                        if self.geometry.dimensions['X'] == 1:
                            f = interp1d(all_lats[n], values_at_interp_points[n], kind=interpolation)
                            value[n] = f(lat[n])
                        else:
                            f = interp1d(all_lons[n], values_at_interp_points[n], kind=interpolation)
                            value[n] = f(lon[n])
                    else:
                        f = interp2d(all_lons[n], all_lats[n], values_at_interp_points[n], kind=interpolation)
                        value[n] = f(as_numpy_array(lon)[n], as_numpy_array(lat)[n])

            elif method == 'bilinear':
                def simple_inter(x1, q1, x2, q2, x):
                    """Simple linear interpolant"""
                    return (x2 - x) / (x2 - x1) * q1 + (x - x1) / (x2 - x1) * q2

                if 'gauss' in self.geometry.name:
                    test = numpy.all(all_j[:, 0] == all_j[:, 1]) and numpy.all(all_j[:, 2] == all_j[:, 3])
                    if not test:
                        raise RuntimeError("Points are not in the awaited order...")

                    lonrs, latrs = self.geometry._rotate_stretch(lon, lat) #waited point
                    lonrs_near, latrs_near = self.geometry._rotate_stretch(*self.geometry.ij2ll(all_i.flatten(), all_j.flatten()))
                    lonrs_near, latrs_near = lonrs_near.reshape(all_i.shape), latrs_near.reshape(all_i.shape)
                    value = simple_inter(latrs_near[:, 0], simple_inter(lonrs_near[:, 0], values_at_interp_points[:, 0],
                                                                        lonrs_near[:, 1], values_at_interp_points[:, 1], lonrs),
                                         latrs_near[:, 2], simple_inter(lonrs_near[:, 2], values_at_interp_points[:, 2],
                                                                        lonrs_near[:, 3], values_at_interp_points[:, 3], lonrs),
                                         latrs)
                else:
                    square = numpy.all(all_i[:, 0] == all_i[:, 1]) and numpy.all(all_i[:, 2] == all_i[:, 3]) and \
                             numpy.all(all_j[:, 0] == all_j[:, 2]) and numpy.all(all_j[:, 1] == all_j[:, 3])
                    if not square:
                        raise NotImplementedError("Bilinear interpolation on non gaussian, irregular grid")
                    
                    # Position of wanted point
                    wanted_i, wanted_j = self.geometry.ll2ij(lon, lat)
                    
                    value = simple_inter(all_j[:, 0], simple_inter(all_i[:, 0], values_at_interp_points[:, 0],
                                                                    all_i[:, 2], values_at_interp_points[:, 2], wanted_i),
                                         all_j[:, 1], simple_inter(all_i[:, 1], values_at_interp_points[:, 1],
                                                                   all_i[:, 3], values_at_interp_points[:, 3], wanted_i),
                                         wanted_j)
        else:
            raise NotImplementedError('*interpolation*=' + interpolation)
        if one:
            try:
                value = float(value)
            except (ValueError, TypeError):
                pass

        return copy.copy(value)


    def as_lists(self, order='C', subzone=None):
        """
        Export values as a dict of lists (in fact numpy arrays).

        :param order: whether to flatten arrays in 'C' (row-major) or
                      'F' (Fortran, column-major) order.
        :param subzone: defines the LAM subzone to be included, in LAM case,
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

        :param subzone: defines the LAM subzone to be included, in LAM case,
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

        :param subzone: defines the LAM subzone to be included, in LAM case,
                     among: 'C', 'CI'.
        """
        if self.spectral:
            raise epygramError("as_points method needs a grid-point field, not a spectral one.")

        field_builder = fpx.field

        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone)
        data4d = self.getdata(d4=True, subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)

        result = FieldSet()
        vcoordinate = self.geometry.vcoordinate.deepcopy()
        for t in range(data4d.shape[0]):
            validity = self.validity[t]
            for k in range(data4d.shape[1]):
                for j in range(data4d.shape[2]):
                    for i in range(data4d.shape[3]):
                        vcoordinate.levels = [levels4d[t, k, j, i]]
                        geometry = UnstructuredGeometry(name='unstructured',
                                                        dimensions={'X':1, 'Y':1},
                                                        vcoordinate=vcoordinate.deepcopy(),
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

        :param subzone: defines the LAM subzone to be included, in LAM case,
                        among: 'C', 'CI'.
        """
        if self.spectral:
            raise epygramError("as_profiles method needs a grid-point field, not a spectral one.")

        field_builder = fpx.field

        lons4d, lats4d = self.geometry.get_lonlat_grid(subzone=subzone, d4=True, nb_validities=1)
        data4d = self.getdata(d4=True, subzone=subzone)
        levels4d = self.geometry.get_levels(d4=True, nb_validities=len(self.validity), subzone=subzone)

        result = FieldSet()
        vcoordinate = self.geometry.vcoordinate.deepcopy()
        for j in range(data4d.shape[2]):
            for i in range(data4d.shape[3]):
                if all([numpy.all(levels4d[0, :, j, i] == levels4d[t, :, j, i]) for t in range(len(self.validity))]):
                    vcoordinate.levels = list(levels4d[0, :, j, i])
                else:
                    vcoordinate.levels = list(levels4d[:, :, j, i].swapaxes(0, 1))
                geometry = UnstructuredGeometry(name='unstructured',
                                                dimensions={'X':1, 'Y':1},
                                                vcoordinate=vcoordinate.deepcopy(),
                                                grid={'longitudes':[lons4d[0, 0, j, i].tolist()],
                                                      'latitudes':[lats4d[0, 0, j, i].tolist()],
                                                      },
                                                position_on_horizontal_grid='center',
                                                )
                profilefield = field_builder(structure='V1D',
                                             fid=dict(copy.deepcopy(self.fid)),
                                             geometry=geometry,
                                             validity=self.validity.copy())
                profilefield.setdata(data4d[:, :, j, i])
                result.append(profilefield)
        return result

    def extract_physicallevels(self, getdata=True):
        """
        Returns a new field excluding levels with no physical meaning.

        :param getdata: if False returns a field without data
        """
        if self.geometry.vcoordinate.typeoffirstfixedsurface == 118:
            # build fid
            fid = {key:(FPDict(value)
                        if isinstance(value, dict)
                        else value)
                   for (key, value) in self.fid.items()}
            #build geometry
            geometry = self.geometry.make_physicallevels_geometry()
            # Field
            newfield = fpx.field(fid=FPDict(fid),
                                 structure=geometry.structure,
                                 geometry=geometry,
                                 validity=self.validity.deepcopy(),
                                 processtype=self.processtype,
                                 comment=self.comment)
            if getdata:
                shp = geometry.get_datashape(dimT=len(self.validity), d4=True)
                data = numpy.ndarray(shp)
                for t in range(len(self.validity)):
                    for k in range(len(geometry.vcoordinate.levels)):
                        k_index = self.geometry.vcoordinate.levels.index(geometry.vcoordinate.levels[k])
                        data[t, k, :, :] = self.getdata(d4=True)[t, k_index, :, :]
                newfield.setdata(data)
    
            return newfield
        else:
            return self


    def extractprofile(self, lon, lat,
                       interpolation='nearest',
                       external_distance=None,
                       exclude_extralevels=True,
                       getdata=True):
        """
        Extracts a vertical profile from the field, given
        the geographic location (*lon*/*lat*) of the profile.

        :param lon: the longitude of the desired point.
        :param lat: the latitude of the desired point.
          If both None, extract a horizontally-averaged profile.
        :param interpolation: defines the interpolation function used to compute
          the profile at requested lon/lat from the fields grid:\n
          - if 'nearest' (default), extracts profile at the horizontal nearest neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear spline interpolation;
          - if 'cubic', computes profile with horizontal cubic spline interpolation.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = abs(target_value - external_field.data)
        :param exclude_extralevels: if True levels with no physical meaning are
                                    suppressed.
        :param getdata: if False returns a field without data
        """
        pointG = self.geometry.make_profile_geometry(lon, lat)
        profile = self.extract_subdomain(pointG,
                                         interpolation=interpolation,
                                         external_distance=external_distance,
                                         getdata=getdata)
        if exclude_extralevels:
            profile = profile.extract_physicallevels()
        return profile

    def extractsection(self, end1, end2,
                       points_number=None,
                       resolution=None,
                       interpolation='linear',
                       exclude_extralevels=True,
                       getdata=True):
        """
        Extracts a vertical section from the field, given
        the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param exclude_extralevels: if True levels with no physical meaning are
                                    suppressed.
        :param getdata: if False returns a field without data
        """
        sectionG = self.geometry.make_section_geometry(end1, end2,
                                                       points_number=points_number,
                                                       resolution=resolution)
        section = self.extract_subdomain(sectionG,
                                         interpolation=interpolation,
                                         getdata=getdata)
        if exclude_extralevels:
            section = section.extract_physicallevels()
        return section

    def extract_subdomain(self, geometry,
                          interpolation='nearest',
                          external_distance=None,
                          getdata=True):
        """
        Extracts a subdomain from a field, given a new geometry.

        :param geometry: defines the geometry on which extract data
        :param interpolation: defines the interpolation function used to compute
                              the profile at requested lon/lat from the fields
                              grid:\n
          - if 'nearest' (default), extracts profile at the horizontal nearest
            neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear interpolation;
          - if 'cubic', computes profile with horizontal cubic interpolation.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = abs(target_value - external_field.data)
        :param getdata: if False returns a field without data
        """
        # build subdomain fid
        subdomainfid = {key:(FPDict(value)
                             if isinstance(value, dict)
                             else value)
                        for (key, value) in self.fid.items()}

        # build vertical geometry
        # k_index: ordered list of the level indexes of self used to build the new field
        kwargs_vcoord = {'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid}
        if self.geometry.vcoordinate.typeoffirstfixedsurface == 119:
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
            kwargs_vcoord['levels'] = copy.copy(self.geometry.vcoordinate.levels)
            k_index = range(len(kwargs_vcoord['levels']))
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 118:
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
            kwargs_vcoord['levels'] = copy.copy(self.geometry.vcoordinate.levels)
            k_index = list(range(len(kwargs_vcoord['levels'])))
        elif self.geometry.vcoordinate.typeoffirstfixedsurface in [100, 102, 103, 109, 1, 106, 255, 160, 200]:
            kwargs_vcoord['levels'] = list(copy.deepcopy(self.geometry.vcoordinate.levels))
            k_index = range(len(kwargs_vcoord['levels']))
        else:
            raise NotImplementedError("type of first surface level: " + str(self.geometry.vcoordinate.typeoffirstfixedsurface))
        if geometry.vcoordinate.typeoffirstfixedsurface not in [255, kwargs_vcoord['typeoffirstfixedsurface']]:
            raise epygramError("extract_subdomain cannot change vertical coordinate.")
        if geometry.vcoordinate.position_on_grid not in [None, '__unknown__', kwargs_vcoord['position_on_grid']]:
            raise epygramError("extract_subdomain cannot change position on vertical grid.")
        if geometry.vcoordinate.grid != {} and geometry.vcoordinate.grid != kwargs_vcoord['grid']:
            raise epygramError("extract_subdomain cannot change vertical grid")
        if geometry.vcoordinate.levels != []:
            k_index = []
            for level in geometry.vcoordinate.levels:
                if level not in kwargs_vcoord['levels']:
                    raise epygramError("extract_subdomain cannot do vertical interpolations.")
                k_index.append(self.geometry.vcoordinate.levels.index(level))
            kwargs_vcoord['levels'] = geometry.vcoordinate.levels
        vcoordinate = VGeometry(**kwargs_vcoord)
        # build geometry
        newgeometry = geometry.deepcopy()
        newgeometry.vcoordinate = vcoordinate
        newgeometry.position_on_horizontal_grid = 'center'
        if geometry.position_on_horizontal_grid not in [None, '__unknown__', newgeometry.position_on_horizontal_grid]:
            raise epygramError("extract_subdomain cannot deal with position_on_horizontal_grid other than 'center'")

        if any([isinstance(level, numpy.ndarray) for level in newgeometry.vcoordinate.levels]):
            # level value is not constant on the domain, we must extract a subdomain.
            # We build a new field with same structure/geometry/validity as self
            # and we fill its data with the level values. We replace levels by index
            # in the geometry to stop recursion and use the extract_subdomain method
            # to get the level values on the new geometry levels.
            temp_field = self.geometry.vcoord_as_field(self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       self.validity)
            temp_geom = geometry.deepcopy()
            del temp_geom.vcoordinate.levels[:]
            temp_geom.vcoordinate.levels.extend(k_index)  # We extract only the requested levels
            extracted = temp_field.extract_subdomain(temp_geom,
                                                     interpolation=interpolation,
                                                     external_distance=external_distance,
                                                     getdata=True).getdata(d4=True)
            extracted = extracted.swapaxes(0, 1)  # z first, then validity for levels (opposite of fields)
            if len(self.validity) > 1:
                # We keep all dims when validity is not unique
                extracted = list(extracted)
            else:
                # We suppress the validity dimension
                extracted = list(extracted[:, 0, ...])
            del newgeometry.vcoordinate.levels[:]
            newgeometry.vcoordinate.levels.extend(extracted)
            del temp_field, temp_geom

        # location & interpolation
        lons, lats = newgeometry.get_lonlat_grid()
        if len(lons.shape) > 1:
            lons = stretch_array(lons)
            lats = stretch_array(lats)
        if not numpy.all(self.geometry.point_is_inside_domain_ll(lons, lats)):
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
        elif interpolation in ('linear', 'cubic', 'bilinear'):
            if interpolation == 'linear':
                interpstr = 'linearly'
            elif interpolation == 'cubic':
                interpstr = 'cubically'
            elif interpolation == 'bilinear':
                interpstr = 'bilinearly'
            if lons.size == 1:
                comment = "Profile " + interpstr + " interpolated @ " + \
                          str((float(lons), float(lats)))
        else:
            raise NotImplementedError(interpolation + " interpolation.")

        # Values
        if getdata:
            shp = newgeometry.get_datashape(dimT=len(self.validity), d4=True)
            data = numpy.ndarray(shp)
            for t in range(len(self.validity)):
                #This is the old version, to delete later
                #for k in range(len(newgeometry.vcoordinate.levels)):
                #    # level = newgeometry.vcoordinate.levels[k]
                #    extracted = self.getvalue_ll(lons, lats,
                #                                 # level=level,  # equivalent to k=k_index[k] except that the k option
                #                                 k=k_index[k],  # still works when level is an array
                #                                 validity=self.validity[t],
                #                                 interpolation=interpolation,
                #                                 external_distance=external_distance,
                #                                 one=False)
                #    data[t, k, :, :] = newgeometry.reshape_data(extracted)
                data[t, ...] = self.getvalue_ll(as_numpy_array(lons)[numpy.newaxis, :].repeat(len(k_index), axis=0).flatten(),
                                                as_numpy_array(lats)[numpy.newaxis, :].repeat(len(k_index), axis=0).flatten(),
                                                k=as_numpy_array(k_index).repeat(lons.size),
                                                validity=self.validity[t],
                                                interpolation=interpolation,
                                                external_distance=external_distance,
                                                one=False).reshape(shp[1:])

        # Field
        newfield = fpx.field(fid=FPDict(subdomainfid),
                             structure=newgeometry.structure,
                             geometry=newgeometry,
                             validity=self.validity.deepcopy(),
                             processtype=self.processtype,
                             comment=comment)
        if getdata:
            newfield.setdata(data)

        return newfield


    def extract_zoom(self, zoom, extra_10th=False):
        """
        Extract an unstructured field with the gridpoints contained in *zoom*.

        :param zoom: a dict(lonmin=, lonmax=, latmin=, latmax=).
        :param extra_10th: if True, add 1/10th of the X/Y extension of the zoom
                           (regular_lonlat grid case only).
        """
        assert not self.spectral, \
               "spectral field: convert to gridpoint beforehand"

        kwargs_zoomgeom = {'vcoordinate':self.geometry.vcoordinate.deepcopy(),
                           'position_on_horizontal_grid':self.geometry.position_on_horizontal_grid,
                           'geoid':self.geometry.geoid}
        if self.geometry.name == 'regular_lonlat':
            if extra_10th:
                dx = (degrees_nearest_mod(zoom['lonmax'], 0.) -
                      degrees_nearest_mod(zoom['lonmin'], 0.)) / 10.
                dy = (zoom['latmax'] - zoom['latmin']) / 10.
                zoom = {'lonmin':zoom['lonmin'] - dx,
                        'lonmax':zoom['lonmax'] + dx,
                        'latmin':zoom['latmin'] - dy,
                        'latmax':zoom['latmax'] + dy}
            imin, jmin = self.geometry.ll2ij(zoom['lonmin'], zoom['latmin'])
            imax, jmax = self.geometry.ll2ij(zoom['lonmax'], zoom['latmax'])
            if imin > imax:
                gridmin = self.geometry.gimme_corners_ll()['ll'][0]
                diff_lonmin = (gridmin - degrees_nearest_mod(zoom['lonmin'],
                                                             gridmin))
                Xres = self.geometry.grid['X_resolution'].get('degrees')
                shift = (diff_lonmin // Xres + 1) * Xres
                shifted_self = self.deepcopy()
                shifted_self.global_shift_center(-shift)
                return shifted_self.extract_zoom(zoom)
            elif imin == imax:  # means 360deg wide
                imin = 0
                imax = self.geometry.dimensions['X'] - 1
            imin = max(int(numpy.ceil(imin)),
                       0)
            imax = min(int(numpy.floor(imax)),
                       self.geometry.dimensions['X'] - 1)
            jmin = max(int(numpy.ceil(jmin)),
                       0)
            jmax = min(int(numpy.floor(jmax)),
                       self.geometry.dimensions['Y'] - 1)
            kwargs_zoomgeom['dimensions'] = {'X':imax - imin + 1,
                                             'Y':jmax - jmin + 1}
            lonmin, latmin = self.geometry.ij2ll(imin, jmin)
            geometryclass = RegLLGeometry
            kwargs_zoomgeom['name'] = 'regular_lonlat'
            kwargs_zoomgeom['grid'] = {'input_position':(0, 0),
                                       'input_lon':Angle(lonmin, 'degrees'),
                                       'input_lat':Angle(latmin, 'degrees'),
                                       'X_resolution':self.geometry.grid['X_resolution'],
                                       'Y_resolution':self.geometry.grid['Y_resolution']}
        else:
            (lons, lats) = self.geometry.get_lonlat_grid()
            lons = lons.flatten()
            lats = lats.flatten()
            zoomlons = []
            zoomlats = []
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
            geometryclass = UnstructuredGeometry
            kwargs_zoomgeom['name'] = 'unstructured'
            kwargs_zoomgeom['grid'] = {'longitudes':zoomlons,
                                       'latitudes':zoomlats}
        zoom_geom = geometryclass(**kwargs_zoomgeom)
        # Il serait mieux d'utiliser
        # zoom_geom = self.geometry.make_zoom_geometry(zoom, extra_10th)
        # mais il faudrait s'assurer d'appliquer le global_shift_center de
        # la même manière sur la géométrie et sur le champ.
        # Le problème ne se poserait plus en utilisant extract_subdomain

        if any([isinstance(level, numpy.ndarray) for level in zoom_geom.vcoordinate.levels]):
            # level value is not constant on the domain, we must extract a subdomain.
            # We build a new field with same structure/geometry/validity as self
            # and we fill its data with the level values. We replace levels by index
            # in the geometry to stop recursion and use the extract_subdomain method
            # to get the level values on the new geometry levels.
            temp_field = self.geometry.vcoord_as_field(self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       self.validity)
            extracted = temp_field.extract_zoom(zoom, extra_10th).getdata(d4=True)
            extracted = extracted.swapaxes(0, 1)  # z first, then validity for levels (opposite of fields)
            if len(self.validity) > 1:
                # We keep all dims when validity is not unique
                extracted = list(extracted)
            else:
                # We suppress the validity dimension
                extracted = list(extracted[:, 0, ...])
            del zoom_geom.vcoordinate.levels[:]
            zoom_geom.vcoordinate.levels.extend(extracted)
            del temp_field

        # Serait plus élégant mais pb d'efficacité (2x plus lent)
        # car extract_subdomain fait une recherche des plus proches points:
        # zoom_field = self.extract_subdomain(zoom_geom)
        # zoom_field.fid = fid
        
        shp = zoom_geom.get_datashape(dimT=len(self.validity), d4=True)
        data = numpy.ma.empty(shp)
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
                    data[t, k, :, :] = numpy.ma.array(zoomvals).reshape(shp[2:])

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


    #This method is intended to replace the extract_zoom one
    #but is still longer at execution
    #def _extract_zoom_with_subdo(self, zoom, extra_10th=False):
    #    """
    #    Extract an unstructured field with the gridpoints contained in *zoom*.
    #
    #    :param zoom: a dict(lonmin=, lonmax=, latmin=, latmax=).
    #    :param extra_10th: if True, add 1/10th of the X/Y extension of the zoom
    #                       (regular_lonlat grid case only).
    #    """
    #    assert not self.spectral, \
    #           "spectral field: convert to gridpoint beforehand"
    #
    #    zoom_geom = self.geometry.make_zoom_geometry(zoom, extra_10th)
    #    if any([isinstance(level, numpy.ndarray) for level in self.geometry.vcoordinate.levels]):
    #        del zoom_geom.vcoordinate.levels[:]
    #    zoom_field = self.extract_subdomain(zoom_geom, interpolation='nearest')
    #
    #    return zoom_field

    def extract_subarray(self,
                         first_i, last_i,
                         first_j, last_j,
                         getdata=True,
                         deepcopy=True):
        """
        Extract a rectangular sub-array from the field, given the i,j index limits
        of the sub-array, and return the extracted field.
        If deepcopy is False, current field and returned field will share
        some attributes (like validity).
        """
        newgeom = self.geometry.make_subarray_geometry(first_i, last_i,
                                                       first_j, last_j)
        if any([isinstance(level, numpy.ndarray) for level in newgeom.vcoordinate.levels]):
            # level value is not constant on the domain, we must extract a subdomain.
            # We build a new field with same structure/geometry/validity as self
            # and we fill its data with the level values. We replace levels by index
            # in the geometry to stop recursion and use the extract_subdomain method
            # to get the level values on the new geometry levels.
            temp_field = self.geometry.vcoord_as_field(self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       self.validity)
            extracted = temp_field.extract_subarray(first_i, last_i,
                                                    first_j, last_j,
                                                    getdata=True,
                                                    deepcopy=True).getdata(d4=True)
            extracted = extracted.swapaxes(0, 1)  # z first, then validity for levels (opposite of fields)
            if len(self.validity) > 1:
                # We keep all dims when validity is not unique
                extracted = list(extracted)
            else:
                # We suppress null dims
                extracted = list(extracted.squeeze())
            del newgeom.vcoordinate.levels[:]
            newgeom.vcoordinate.levels.extend(extracted)
            del temp_field

        # copy the field object and set the new geometry, then data
        if hasattr(self, 'as_real_field'):
            field_kwargs = self.as_real_field()._attributes
        else:
            field_kwargs = self._attributes
        if deepcopy:
            field_kwargs = copy.deepcopy(field_kwargs)
        field_kwargs['geometry'] = newgeom
        newfield = fpx.field(**field_kwargs)
        if getdata:
            assert not self.spectral
            # select data
            subdata = self.getdata(d4=True)[:, :, first_j:last_j + 1, first_i:last_i + 1]
            newfield.setdata(subdata)

        return newfield

    def extract_subsample(self,
                          sample_x, sample_y, sample_z,
                          getdata=True,
                          deepcopy=True):
        """
        Extract a subsample from field by decreasing resolution.
        
        :param sample_x: take one over <sample_x> points in the x direction
        :param sample_y: same for the y direction
        :param sample_z: same for the z direction
        If deepcopy is False, current field and returned field will share
        some attributes (like validity).

        The extension zone of the original field is kept in the new field.
        Hence, it's recommended to call this method on field without extension zone.

        Because the resolution change during this operation, field must be
        localize on the grid center (use center method before if not).
        """
        assert self.geometry.position_on_horizontal_grid == 'center', \
               "Field must be centered on the grid"

        newgeom = self.geometry.make_subsample_geometry(sample_x,
                                                        sample_y,
                                                        sample_z)
        
        if any([isinstance(level, numpy.ndarray) for level in newgeom.vcoordinate.levels]):
            # level value is not constant on the domain, we must extract a subdomain.
            # We build a new field with same structure/geometry/validity as self
            # and we fill its data with the level values. We replace levels by index
            # in the geometry to stop recursion and use the extract_subdomain method
            # to get the level values on the new geometry levels.
            temp_field = self.geometry.vcoord_as_field(self.geometry.vcoordinate.typeoffirstfixedsurface,
                                                       self.validity)
            extracted = temp_field.extract_subsample(sample_x, sample_y, sample_z,
                                                     getdata=True,
                                                     deepcopy=True).getdata(d4=True)
            extracted = extracted.swapaxes(0, 1)  # z first, then validity for levels (opposite of fields)
            if len(self.validity) > 1:
                # We keep all dims when validity is not unique
                extracted = list(extracted)
            else:
                # We suppress null dims
                extracted = list(extracted.squeeze())
            del newgeom.vcoordinate.levels[:]
            newgeom.vcoordinate.levels.extend(extracted)
            del temp_field
        
        # copy the field object and set the new geometry, then data
        field_kwargs = self._attributes
        if deepcopy:
            field_kwargs = copy.deepcopy(field_kwargs)
        field_kwargs['geometry'] = newgeom
        newfield = fpx.field(**field_kwargs)
        if getdata:
            assert not self.spectral
            # select data
            subdata = self.getdata(d4=True)[:, ::sample_z,
                                               ::sample_y,
                                               ::sample_x]
            newfield.setdata(subdata)

        return newfield

    def shave(self, minval=None, maxval=None):
        """
        Shave data: cut values greater (resp. lower) than **maxval** to
        **maxval** (resp. **minval**).

        :param maxval: upper threshold
        :param minval: lower threshold
        """
        data = self.getdata()
        if isinstance(data, numpy.ma.masked_array):
            if maxval is not None:
                data.data[data.data > maxval] = maxval
            if minval is not None:
                data.data[data.data < minval] = minval
        else:
            if maxval is not None:
                data[data > maxval] = maxval
            if minval is not None:
                data[data < minval] = minval
        self.data = data

    def resample(self, target_geometry,
                 weighting='nearest',
                 subzone=None,
                 neighbour_info=None,
                 # neighbouring options
                 radius_of_influence=None,
                 neighbours=8,
                 epsilon=0,
                 reduce_data=True,
                 nprocs=1,
                 segments=None,
                 # resampling options
                 fill_value=None,
                 sigma=None,
                 with_uncert=False):
        """
        Resample data (interpolate) on a **target geometry**, using pyresample.

        :param weighting: among ['nearest',
                                 'gauss' (in which case sigma can be specified,
                                          cf. sigma option below),
                                 lambda r:1/r**2 (or other function, r in m)]
        :param subzone: restrain source field to LAMzone in LAM geometry case
        :param neighbour_info: - if True: skips the resampling, only compute
                                          neighbours and return this
                                          information as a tuple;
                               - if given as a tuple, skips the computation of
                                 neighbours, taking those information
                                 given in argument (from a former call with
                                 neighbour_info=True)

        Options from pyresample, neighbouring:

        :param radius_of_influence: Cut off distance in meters, default value
                                    computed from geometry
        :param neighbours: The number of neigbours to consider for each grid
                           point
        :param epsilon: Allowed uncertainty in meters. Increasing uncertainty
                        reduces execution time
        :param reduce_data: Perform initial coarse reduction of source dataset
                            in order to reduce execution time
        :param nprocs: Number of processor cores to be used
        :param segments: Number of segments to use when resampling.
                         If set to None an estimate will be calculated

        Options from pyresample, resampling:

        :param fill_value: Set undetermined pixels to this value.
                           If fill_value is None a masked array is returned
                           with undetermined pixels masked
        :param sigma: sigma (in m) to use for the gauss weighting
                      w = exp(-dist^2/sigma^2)
                      Default is twice the resolution.
        :param with_uncert: Calculate uncertainty estimates (returned as side
                            fields)
        """
        assert not self.spectral, "field must be gridpoint"
        try:
            from pyresample.geometry import GridDefinition, SwathDefinition
            from pyresample.kd_tree import (get_neighbour_info,
                                            get_sample_from_neighbour_info)
        except ImportError:
            raise ImportError("""pyresample is not installed.
                                 You can install it locally using:\n
                                 'pip install --user pyresample'""")
        # PART 0: computation of parameters for PARTS 1 & 2
        lons, lats = target_geometry.get_lonlat_grid(subzone=subzone, d4=True, nb_validities=1,
                                                     force_longitudes=']-180,180]')
        lons = lons[0, 0, :, :]
        lats = lats[0, 0, :, :]
        if isinstance(lons, numpy.ma.MaskedArray) and lons.mask.any():
            reduce_data = False  # functionality does not work properly
        target_geo = GridDefinition(lons, lats)

        def _resolution():
            if 'gauss' in self.geometry.name:
                # gauss grid case
                resolution = self.geometry.resolution_j(self.geometry.dimensions['lat_number'] - 1)
            elif self.geometry.name in ('regular_lonlat', 'unstructured'):
                resolution = self.geometry.resolution_ij(self.geometry.dimensions['X'] / 2,
                                                         self.geometry.dimensions['Y'] / 2)
            else:
                resolution = (self.geometry.grid.get('X_resolution') +
                              self.geometry.grid.get('Y_resolution')) / 2.
            return resolution

        # PART 1: computation of neighbours
        if neighbour_info in (True, None, False):
            lons, lats = self.geometry.get_lonlat_grid(subzone=subzone, d4=True, nb_validities=1,
                                                       force_longitudes=']-180,180]')
            lons = lons[0,0,:,:]
            lats = lats[0,0,:,:]
            if isinstance(lons, numpy.ma.MaskedArray) and lons.mask.any():
                reduce_data = False  # functionality does not work properly
            source_geo = GridDefinition(lons, lats)
            if radius_of_influence is None:
                radius_of_influence = 4. * _resolution()
            if weighting == 'nearest':
                neighbours = 1
            (valid_input_index,
             valid_output_index,
             index_array,
             distance_array) = get_neighbour_info(source_geo,
                                                  target_geo,
                                                  radius_of_influence=radius_of_influence,
                                                  neighbours=neighbours,
                                                  epsilon=epsilon,
                                                  reduce_data=reduce_data,
                                                  nprocs=nprocs,
                                                  segments=segments)
            if neighbour_info is True:
                return (valid_input_index,
                        valid_output_index,
                        index_array,
                        distance_array)
        else:  # neighbour_info is given: tuple of information about neighbours
            assert len(neighbour_info) == 4
            (valid_input_index,
             valid_output_index,
             index_array,
             distance_array) = neighbour_info
        # PART 2: computation of resampling
        source_data = self.getdata(d4=True, subzone=subzone)
        resampled_data = numpy.ma.zeros((source_data.shape[0],  # t
                                         source_data.shape[1],  # z
                                         target_geo.shape[0],  # y
                                         target_geo.shape[1]))  # x
        if with_uncert:
            stddev = numpy.ma.zeros((source_data.shape[0],  # t
                                     source_data.shape[1],  # z
                                     target_geo.shape[0],  # y
                                     target_geo.shape[1]))  # x
            counts = numpy.ma.zeros((source_data.shape[0],  # t
                                     source_data.shape[1],  # z
                                     target_geo.shape[0],  # y
                                     target_geo.shape[1]))  # x
        if callable(weighting):  # custom: weighting function
            weight_funcs = weighting
            weighting = 'custom'
        elif weighting == 'nearest':
            weighting = 'nn'
            weight_funcs = None
        elif weighting == 'gauss':

            def gauss(s):
                return lambda r: numpy.exp(-r ** 2 / s ** 2)

            weighting = 'custom'
            if sigma is None:
                sigma = 2. * _resolution()
            else:
                sigma = float(sigma)
            weight_funcs = gauss(sigma)
        else:
            raise epygramError("unknown weighting='" + str(weighting))
        for t in range(source_data.shape[0]):
            for z in range(source_data.shape[1]):
                rdata = get_sample_from_neighbour_info(weighting,
                                                       target_geo.shape,
                                                       source_data[t, z, :, :],
                                                       valid_input_index,
                                                       valid_output_index,
                                                       index_array,
                                                       distance_array,
                                                       weight_funcs=weight_funcs,
                                                       fill_value=fill_value,
                                                       with_uncert=with_uncert)
                if with_uncert:
                    resampled_data[t, z, :, :] = rdata[0][:, :]
                    stddev[t, z, :, :] = rdata[1][:, :]
                    counts[t, z, :, :] = rdata[2][:, :]
                else:
                    resampled_data[t, z, :, :] = rdata[:, :]

        # build final field
        field_kwargs = copy.deepcopy(self._attributes)
        field_kwargs['geometry'] = target_geometry
        newfield = fpx.field(**field_kwargs)
        if not resampled_data.mask.any():
            resampled_data = resampled_data.data
        newfield.setdata(resampled_data)
        if with_uncert:
            stddev_field = newfield.deepcopy()
            stddev_field.setdata(stddev)
            stddev_field.comment += ':stddev_from_resampling'
            counts_field = newfield.deepcopy()
            counts_field.setdata(counts)
            stddev_field.comment += ':number_of_points_from_resampling'
            return newfield, stddev_field, counts_field
        else:
            return newfield

    def resample_on_regularll(self, borders, resolution_in_degrees, **kwargs):
        """
        Resample data (interpolate) on a target geometry defined as a
        regular lonlat grid of given resolution.
        Wrapper to resample() method, building target geometry.

        :param dict borders: dict(lonmin=, latmin=, lonmax=, latmax=)
        :param resolution_in_degrees: resolution of the lonlat grid in degrees.
                                      Supposed to be a divisor of lonmax-lonmin
                                      and latmax-latmin.

        Other kwargs to be passed to resample() method.
        """
        DX = borders['lonmax'] - borders['lonmin']
        if DX < 0.:
            DX += 360. 
        DY = borders['latmax'] - borders['latmin']
        X = int(DX / resolution_in_degrees) + 1
        Y = int(DY / resolution_in_degrees) + 1
        grid = {'input_position':(0, 0),
                'input_lon':Angle(borders['lonmin'], 'degrees'),
                'input_lat':Angle(borders['latmin'], 'degrees'),
                'X_resolution':Angle(resolution_in_degrees, 'degrees'),
                'Y_resolution':Angle(resolution_in_degrees, 'degrees')}
        target_geometry = RegLLGeometry(name='regular_lonlat',
                                        vcoordinate=self.geometry.vcoordinate.deepcopy(),
                                        dimensions={'X':X, 'Y':Y},
                                        grid=grid,
                                        position_on_horizontal_grid='center',
                                        geoid=self.geometry.geoid)
        return self.resample(target_geometry, **kwargs)

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

    def remove_level(self, k):
        """
        Remove one level from the current field
        :param k: index of the level to remove
        """
        del self.geometry.vcoordinate.levels[k]
        self.setdata(numpy.delete(self.getdata(d4=True), k, 1))

    def decumulate(self, center=False):
        """
        Decumulate cumulated fields (for a field with temporal dimension !).

        :param center: - if False, values at t are mean values between t-1 and t,
                         except at t=0 where it remains untouched.
                       - if True, values at t are then (v[t-1, t] + v[t,t+1])/2.,
                         except for t=0 where it becomes v[0, 1]
                         and t=last where it remains the same.
        """
        if len(self.validity) > 1:
            data = self.getdata(d4=True)
            data[1:, ...] = data[1:, ...] - data[0:-1, ...]
            if center:
                data[0, ...] = data[1, ...]  # the actual mean between 0 and 1
                if len(self.validity) > 2:
                    data[1:-1, ...] = (data[1:-1, ...] + data[2:, ...]) / 2.
            self.setdata(data)

    def time_smooth(self, length, window='center'):
        """
        Smooth data in time.

        :param length: time span to average on
        :param window: replace data[t] by:\n
          - data[t-length/2:t+length/2].mean() if 'center';
          - data[t-length:t].mean() if 'left';
          - data[t:t+length/2].mean() if 'right'.
        """
        assert len(self.validity) > 1
        data = copy.deepcopy(self.getdata(d4=True))
        if window == 'center':
            tinf = length // 2
            tsup = length // 2
        elif window == 'left':
            tinf = length
            tsup = 0
        elif window == 'right':
            tinf = 0
            tsup = length
        for t in range(len(self.validity)):
            t9 = max(0, t - tinf)
            t1 = min(len(self.validity), t + tsup)
            data[t, ...] = self.getdata(d4=True)[t9:t1, ...].mean(axis=0)
        self.setdata(data)

    def time_reduce(self, reduce_function='mean'):
        """
        Make a time reduction of field,
        i.e. apply the requested **reduce_function** on the time dimension of
        the field. Replace data (and validity) in place.

        :param reduce_function: must be a numpy.ndarray method,
                                e.g. 'mean', 'sum', 'min', 'max', ...
        """
        data = self.getdata(d4=True)
        assert hasattr(data, reduce_function), 'unknown method:{} for numpy arrays'.format(reduce_function)
        reduced = getattr(data, reduce_function)(axis=0)
        validity = self.validity[-1].deepcopy()
        validity.set(cumulativeduration=self.validity[-1].get() -
                     self.validity[0].get(),
                     statistical_process_on_duration=reduce_function)
        self.validity = validity
        self.setdata(reduced.squeeze())

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def stats(self, subzone=None):
        """
        Computes some basic statistics on the field, as a dict containing:
        {'min', 'max', 'mean', 'std', 'quadmean', 'nonzero'}.
        See each of these methods for details.

        :param subzone: optional, among ('C', 'CI'), for LAM fields only, plots
                        the data resp. on the C or C+I zone.
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
        return super(_D3CommonField, self).min(subzone=subzone)

    def max(self, subzone=None):
        """Returns the maximum value of data."""
        return super(_D3CommonField, self).max(subzone=subzone)

    def mean(self, subzone=None):
        """Returns the mean value of data."""
        return super(_D3CommonField, self).mean(subzone=subzone)

    def std(self, subzone=None):
        """Returns the standard deviation of data."""
        return super(_D3CommonField, self).std(subzone=subzone)

    def quadmean(self, subzone=None):
        """Returns the quadratic mean of data."""
        return super(_D3CommonField, self).quadmean(subzone=subzone)

    def nonzero(self, subzone=None):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        return super(_D3CommonField, self).nonzero(subzone=subzone)

    def sgmean(self):
        """Returns an average of the field, weighted by the map factor at each point."""
        map_factor = self.geometry.map_factor_field()  # Extraction of the map_factor field from geometry object
        map_factor_data = map_factor.getdata()
        nonmasked_value = map_factor_data.count(1)  # Number of longitude points for each (pseudo)latitude circle

        myfield_data = self.getdata()
        myfield_data_norm = numpy.ma.divide(myfield_data, map_factor_data)  # Global map factor normalisation
        myfield_data_norm_sumlat = numpy.sum(myfield_data_norm, 1)  # Sum over longitude for each latitude circle
        myfield_data_norm_sumlat = numpy.ma.divide(myfield_data_norm_sumlat, nonmasked_value)  # Normalisation by the number of longitude points for each latitude circle
        Nb_lat = len(myfield_data_norm_sumlat)  # Nb of latitude circles
        myfield_data_sgmean = numpy.sum(myfield_data_norm_sumlat) / Nb_lat  # Sum over latitude normalised by the number of latitude circles
        return myfield_data_sgmean

    def dctspectrum(self, subzone=None, level_index=None, validity_index=None):
        """
        Returns the DCT spectrum of the field, as a
        :class:`epygram.spectra.Spectrum` instance.

        :param subzone: LAM zone on which to compute the DCT,
                        among (None, 'C', 'CI', 'CIE'). Defaults to 'C'.
        :param level_index: the level index on which to compute the DCT
        :param validity_index: the validity index on which to compute the DCT
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
        if subzone is None and self.geometry.grid.get('LAMzone', None) is not None:
            subzone = 'C'
        variances = esp.dctspectrum(field2dt0.getdata(subzone=subzone))
        spectrum = esp.Spectrum(variances[1:],
                                name=str(self.fid),
                                resolution=self.geometry.grid['X_resolution'] / 1000.,
                                mean2=variances[0])

        return spectrum

    def histogram(self,
                  subzone=None,
                  over=(None, None),
                  title=None,
                  get_n_bins_patches=False,
                  together_with=[],
                  mask_threshold=None,
                  center_hist_on_0=False,
                  minmax_in_title=True,
                  figsize=None,
                  **hist_kwargs):
        """
        Build an histogram of the field data.

        :param subzone: LAM zone on which to compute the DCT,
                        among (None, 'C', 'CI', 'CIE').
        :param over: any existing figure and/or ax to be used for the
                     plot, given as a tuple (fig, ax), with None for
                     missing objects. *fig* is the frame of the
                     matplotlib figure, containing eventually several
                     subplots (axes); *ax* is the matplotlib axes on
                     which the drawing is done. When given (is not None),
                     these objects must be coherent, i.e. ax being one of
                     the fig axes.
        :param title: title for the plot
        :param get_n_bins_patches: if True, returns n, bins, patches
                                   (cf. matplotlib hist()) in addition
        :param together_with: another field or list of fields to plot on the
                              same histogram
        :param mask_threshold: dict with min and/or max value(s) to mask outside.
        :param center_hist_on_0: to center the histogram on 0.
        :param minmax_in_title: if True and **range** is not None,
                                adds min and max values in title.
        :param hist_kwargs: any keyword argument to be passed to
                            matplotlib's hist()
        :param figsize: figure sizes in inches, e.g. (5, 8.5).
                        Default figsize is config.plotsizes.
        """
        import matplotlib.pyplot as plt
        # initializations
        plt.rc('font', family='serif')
        if figsize is None:
            figsize = config.plotsizes

        if self.spectral:
            raise epygramError("please convert to gridpoint with sp2gp()" +
                               " method before plotting.")
        mask_outside = {'min':-config.mask_outside,
                        'max':config.mask_outside}
        if mask_threshold is not None:
            mask_outside.update(mask_threshold)

        # get data
        data1d = numpy.ma.masked_outside(self.getdata(subzone=subzone),
                                         mask_outside['min'],
                                         mask_outside['max'])
        data1d = [stretch_array(data1d)]
        fig, ax = set_figax(*over, figsize=figsize)
        if together_with != []:
            if isinstance(together_with, D3Field):
                together_with = [together_with]
            data1d = data1d + [stretch_array(numpy.ma.masked_outside(f.getdata(subzone=subzone),
                                                                     mask_outside['min'],
                                                                     mask_outside['max']))
                               for f in together_with]
        # plot params
        if title is None:
            if len(data1d) == 1:
                ax.set_title("\n".join([str(self.fid[sorted(self.fid.keys())[0]]),
                                        str(self.validity.get())]))
        else:
            ax.set_title(title)
        if 'rwidth' not in hist_kwargs.keys():
            hist_kwargs['rwidth'] = 0.9
        if 'label' not in hist_kwargs.keys() and len(data1d) > 1:
            hist_kwargs['label'] = [self.fid[sorted(self.fid.keys())[0]]] + \
                                   [f.fid[sorted(f.fid.keys())[0]]
                                    for f in together_with]
        if hist_kwargs.get('range') is not None:
            try:
                m = float(hist_kwargs['range'][0])
            except ValueError:
                m = min([d.min() for d in data1d])
            try:
                M = float(hist_kwargs['range'][1])
            except ValueError:
                M = max([d.max() for d in data1d])
            hist_kwargs['range'] = (m, M)
            if minmax_in_title:
                minmax_in_title = '(min: ' + \
                    '{: .{precision}{type}}'.format(min([d.min() for d in data1d]),
                                                    type='E', precision=3) + \
                    ' // max: ' + \
                    '{: .{precision}{type}}'.format(max([d.max() for d in data1d]),
                                                    type='E', precision=3) + ')'
        else:
            minmax_in_title = ''
        if hist_kwargs.get('range') is None:
            if hist_kwargs.get('bins') is None or \
               isinstance(hist_kwargs.get('bins'), int):
                hist_kwargs['range'] = (min([d.min() for d in data1d]),
                                        max([d.max() for d in data1d]))
        if isinstance(hist_kwargs.get('bins'), list):
            if hist_kwargs['bins'][0] == 'min':
                hist_kwargs['bins'][0] = min([d.min() for d in data1d])
            if hist_kwargs['bins'][-1] == 'min':
                hist_kwargs['bins'][-1] = max([d.max() for d in data1d])
        # build histogram
        n, bins, patches = ax.hist(data1d, **hist_kwargs)
        # finalize graphical options
        if len(bins) > 8:
            step = len(bins) // 7
            ticks = [b for b in bins[::step]]
            if bins[-1] != ticks[-1]:
                ticks.append(bins[-1])
        else:
            ticks = bins
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.set_xticks(bins, minor=True)
        ax.set_xticks(ticks)
        ax.set_xticklabels(['{:.3E}'.format(t) for t in ticks],
                           rotation=27.5,
                           horizontalalignment='right')
        if center_hist_on_0:
            M = max(abs(ticks[0]), abs(ticks[-1]))
            ax.set_xlim(-M, M)
            ax.axvline(0., linestyle='-', color='k')  # vertical lines
        else:
            ax.set_xlim(ticks[0], ticks[-1])
        ax.set_ylabel('Number of gridpoints')
        xlabel = 'Gridpoint values'
        if minmax_in_title:
            xlabel += '\n' + minmax_in_title
        ax.set_xlabel(xlabel)
        if 'label' in hist_kwargs.keys():
            ax.legend()
        ax.grid()
        if get_n_bins_patches:
            return fig, ax, n, bins, patches
        else:
            return fig, ax

    def scatter_with(self, other,
                     fig=None,
                     ax=None,
                     figsize=None,
                     mask_outside=config.mask_outside,
                     subzone=None,
                     fidkey=None,
                     xlabel='__auto__',
                     ylabel='__auto__',
                     bisect=True,
                     **scatter_kwargs):
        """
        Make a scatter plot of self data against other data.

        :param fig: a preexisting figure on which to plot
        :param ax: a preexisting axis on which to plot
        :param figsize: figure sizes in inches, e.g. (5, 8.5).
                        Default figsize is config.plotsizes.
        :param mask_outside: threshold to mask values outside (-mask_outside, +mask_outside)
        :param subzone: LAM subzone
        :param fidkey: fid key to be written on axes labels
        :param xlabel: custom x label
        :param ylabel: custom y label
        :param bisect: plot x=y line

        Other kwargs are passed to scatter() method.
        """
        import matplotlib.pyplot  # Workaround for https://github.com/UMR-CNRM/bronx/issues/9
        for f in (self, other):
            if f.spectral:
                f.sp2gp()
        data, otherdata = self._masked_any(other, mask_outside, subzone=subzone)
        otherdata = otherdata.compressed()
        data = data.compressed()
        fig, ax = set_figax(fig, ax, figsize=figsize)
        ax.scatter(data, otherdata, **scatter_kwargs)
        if xlabel is not None:
            if xlabel == '__auto__':
                if fidkey is None:
                    xlabel = str(self.fid)
                else:
                    xlabel = str(self.fid[fidkey])
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            if ylabel == '__auto__':
                if fidkey is None:
                    ylabel = str(other.fid)
                else:
                    ylabel = str(other.fid[fidkey])
            ax.set_ylabel(ylabel)
        if bisect:
            xmin = min(data.min(), otherdata.min())
            xmax = max(data.max(), otherdata.max())
            ax.plot((xmin, xmax), (xmin, xmax), color='red')
        return fig, ax

    def global_shift_center(self, longitude_shift):
        """
        Shifts the center of the geometry (and the data accordingly) by
        *longitude_shift* (in degrees).

        :param longitude_shift: has to be a multiple of the grid's resolution
                                in longitude.

        For global RegLLGeometry grids only.
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

        :param out: the output open file-like object.
        :param vertical_geometry: if True, writes the validity of the field.
        :param vertical_geometry: if True, writes the vertical geometry of the
                                  field.
        :param cumulativeduration: if False, not written.
        :param arpifs_var_names: if True, prints the equivalent 'arpifs'
                                 variable names.
        :param fid: if True, prints the fid.
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
        Dumps the field in a netCDF file.

        :param filename: filename to dump in
        :param variablename: variable name in netCDF
        :param fidkey: forces *variablename* to self.fid[*fidkey*];
                       defaults to 'netCDF'
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
            super(_D3CommonField, self)._check_operands(other)

    def __add__(self, other):
        """
        Definition of addition, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'+'} and null validity in the general case,
        or conserved fid if shared by both fields, and extended validity if
        both consistent; to be checked anyway.
        """
        new_attributes = dict(structure=self.structure,
                              geometry=self.geometry,
                              spectral_geometry=self.spectral_geometry,
                              validity=FieldValidityList(length=len(self.validity)))
        if isinstance(other, self.__class__) and (self.validity.is_valid() and other.validity.is_valid()):
            if len(self.validity) == 1:
                if self.validity[0] == other.validity[0]:
                    new_attributes.update(validity=self.validity.deepcopy())
                else:
                    if None not in (self.validity.cumulativeduration(), other.validity.cumulativeduration()):
                        if self.validity.get() != other.validity.get():
                            if self.validity.get() < other.validity.get():
                                last = other
                                first = self
                            elif other.validity.get() < self.validity.get():
                                last = self
                                first = other
                            if last.validity.get() - last.validity.cumulativeduration() == first.validity.get():
                                # the beginning of other cumulation is the end of self cumulation
                                validity = last.validity.deepcopy()
                                validity.set(cumulativeduration=self.validity.cumulativeduration() + other.validity.cumulativeduration())
                                new_attributes.update(validity=validity)
            if self.fid == other.fid:
                new_attributes.update(fid=self.fid)
        newfield = self._add(other, **new_attributes)
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
        with 'fid' = {'op':'-'} and null validity in the general case,
        or conserved fid if shared by both fields, and extended validity if
        both consistent; to be checked anyway.
        """
        new_attributes = dict(structure=self.structure,
                              geometry=self.geometry,
                              spectral_geometry=self.spectral_geometry,
                              validity=FieldValidityList(length=len(self.validity)))
        if isinstance(other, self.__class__) and (self.validity.is_valid() and other.validity.is_valid()):
            if len(self.validity) == 1:
                if self.validity[0] == other.validity[0]:
                    new_attributes.update(validity=self.validity.deepcopy())
                else:
                    if None not in (self.validity.cumulativeduration(), other.validity.cumulativeduration()):
                        if self.validity.get() != other.validity.get():
                            if other.validity.get() < self.validity.get():
                                if self.validity.get() - self.validity.cumulativeduration() == \
                                   other.validity.get() - other.validity.cumulativeduration():
                                    # same start of cumulation
                                    validity = self.validity.deepcopy()
                                    validity.set(cumulativeduration=self.validity.get() - other.validity.get())
                                    new_attributes.update(validity=validity)
                    else:
                        if self.validity.get() != other.validity.get():
                            if None in (self.validity.get(), other.validity.get()):
                                validity = FieldValidityList(length=len(self.validity))
                            elif self.validity.get() < other.validity.get():
                                validity = other.validity.deepcopy()
                                validity.set(cumulativeduration=other.validity.get() - self.validity.get(),
                                             statistical_process_on_duration=8)
                            elif other.validity.get() < self.validity.get():
                                validity = self.validity.deepcopy()
                                validity.set(cumulativeduration=self.validity.get() - other.validity.get(),
                                             statistical_process_on_duration=4)
                            new_attributes.update(validity=validity)
            if self.fid == other.fid:
                new_attributes.update(fid=self.fid)
        newfield = self._sub(other, **new_attributes)
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


class D3Field(_D3CommonField):
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
                type=Geometry,
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
                access='rwx',
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

        :param spectral_geometry: instance of SpectralGeometry, actually
                                  containing spectral transform subroutine (in
                                  in its own *gp2sp()* method).
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
            field_dX.fid = {'derivative':'x'}
            field_dY.fid = {'derivative':'y'}
            field_dX.setdata(gpderivX)
            field_dY.setdata(gpderivY)
        else:
            raise epygramError('field must be spectral to compute its spectral derivatives.')

        return (field_dX, field_dY)

    def getdata(self, subzone=None, d4=False):
        """
        Returns the field data, with 3D shape if the field is not spectral,
        2D if spectral.

        :param subzone: optional, among ('C', 'CI'), for LAM fields only,
                        returns the data resp. on the C or C+I zone.
                        Default is no subzone, i.e. the whole field.
        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with
                     respect to geometry.

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

        :param data: dimensions should in any case be ordered in a subset of
        (t,z,y,x), or (t,z,n) if spectral (2D spectral coefficients must be 1D
        with ad hoc ordering, n being the total number of spectral
        coefficients).

        *data* may be 4D (3D if spectral) even if the field is not, as long as
        the above dimensions ordering is respected.
        """
        if not isinstance(data, numpy.ndarray):
            data = numpy.array(data)

        # determination of field shape
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

        if (len(data.shape) == 3 and self.spectral) or len(data.shape) == 4:
            # data is given as 4D: check compatibility
            assert data.shape == shp, \
                ' '.join(['data', str(data.shape),
                          'should have shape', str(shp)])
        else:
            # data is squeezed: reshape
            try:
                data = data.reshape(shp)
            except ValueError:
                print("'data' shape has to be compatible with field shape: {}".format(shp))
                raise
            """# find indexes corresponding to dimensions
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
            assert (len(numpy.shape(data)) == dimensions
                        or numpy.shape(data) == (1,)), \
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
            data = data.reshape(shp)"""
        super(D3Field, self).setdata(data)

    data = property(getdata, setdata, Field.deldata, "Accessor to the field data.")

    def select_subzone(self, subzone):
        """
        If a LAMzone defines the field, select only the *subzone* from it.

        :param subzone: among ('C', 'CI').

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
        Take care (*i, j, k, t*) is python-indexing, ranging from 0 to
        dimension-1.

        :param i: the index in X direction
        :param j: the index in Y direction
        :param k: the index of the level (not a value in Pa or m...)
        :param t: the index of the temporal dimension (not a validity object)

        *k* and *t* can be scalar even if *i* and *j* are arrays.

        If *one* is False, returns [value] instead of value.
        """
        if t is None:
            if len(self.validity) > 1:
                raise epygramError("*t* is mandatory when there are several validities")
            t = 0
        if k is None:
            if self.geometry.datashape['k']:
                raise epygramError("*k* is mandatory when field has a vertical coordinate")
            k = 0
        if j is None:
            if self.geometry.datashape['j']:
                raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
            j = 0
        if i is None:
            if self.geometry.datashape['i']:
                raise epygramError("*i* is mandatory when field has one horizontal dimension")
            i = 0

        i, j = as_numpy_array(i).flatten(), as_numpy_array(j).flatten()
        k, t = as_numpy_array(k).flatten(), as_numpy_array(t).flatten()

        if not numpy.all(self.geometry.point_is_inside_domain_ij(i, j)):
            raise ValueError("point is out of field domain.")

        sizes = set([len(x) for x in [i, j, k, t]])
        if len(sizes) > 2 or (len(sizes) == 2 and 1 not in sizes):
            raise epygramError("each of i, j, k and t must be scalar or have the same length as the others")

        value = (self.getdata(d4=True)[t, k, j, i]).copy()

        if value.size == 1 and one:
            value = value.item()
        return value

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.

        :param level: the requested level expressed in coordinate value (Pa, m...)
        :param k: the index of the requested level
        """
        if k is None and level is None:
            raise epygramError("You must give k or level.")
        if k is not None and level is not None:
            raise epygramError("You cannot give, at the same time, k and level")
        if level is not None:
            if level not in self.geometry.vcoordinate.levels:
                raise epygramError("The requested level does not exist.")
            my_level = level
            my_k = self.geometry.vcoordinate.levels.index(level)
        else:
            my_k = k
            my_level = self.geometry.vcoordinate.levels[k]

        if self.structure not in ['3D', 'V2D', 'V1D']:
            raise epygramError("It's not possible to extract a level from a " + self.structure + " field.")

        kwargs_vcoord = {'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid,
                         'levels':[my_level]}
        if self.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
            kwargs_vcoord['grid'] = copy.copy(self.geometry.vcoordinate.grid)
        newvcoordinate = VGeometry(**kwargs_vcoord)
        newgeometry = self.geometry.deepcopy()
        newgeometry.vcoordinate = newvcoordinate
        generic_fid = self.fid.get('generic', {})
        generic_fid['level'] = my_level if not isinstance(my_level, numpy.ndarray) else my_k  # to avoid arrays in fid
        kwargs_field = {'structure':newgeometry.structure,
                        'validity':self.validity.copy(),
                        'processtype':self.processtype,
                        'geometry':newgeometry,
                        'fid':{'generic':generic_fid}}
        if self.spectral_geometry is not None:
            kwargs_field['spectral_geometry'] = self.spectral_geometry.copy()
        newfield = fpx.field(**kwargs_field)
        newfield.setdata(self.getdata(d4=True)[:, my_k:my_k + 1, ...])

        return newfield

    def getvalidity(self, index_or_validity):
        """
        Returns the field restrained to one of its temporal validity as a new
        field.

        :param index_or_validity: can be either a
                                  :class:`epygram.base.FieldValidity` instance
                                  or the index of the requested validity in the
                                  field's FieldValidityList.
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
        if len(numpy.array(newfield.geometry.vcoordinate.levels).shape) > 1:
            # levels are, at least, dependent on position
            levels4d = newfield.geometry.get_levels(d4=True, nb_validities=len(self.validity))
            newfield.geometry.vcoordinate.levels = list(levels4d[index, ...].squeeze())

        return newfield

    def center(self):
        """
        Performs an averaging on data values to center the field on the grid (aka shuman)
        and modify geometry.
        """
        posSpl = self.geometry.position_on_horizontal_grid.split('-')
        if len(posSpl) == 2:  # not center
            posY, posX = posSpl
            if (posY == 'lower' or posX == 'left'):
                if not isinstance(self.geometry, ProjectedGeometry):
                    raise NotImplementedError("Centering is not available with non-projected geometries")
                data = self.getdata(d4=True)
                if data is not None:  # This test is useful if field has been retrieved with getdata=False
                    if posY == 'lower':
                        data[:, :, :-1, :] = 0.5 * (data[:, :, :-1, :] + data[:, :, 1:, :])
                        data[:, :, -1, :] = data[:, :, -2, :]
                    if posX == 'left':
                        data[:, :, :, :-1] = 0.5 * (data[:, :, :, :-1] + data[:, :, :, 1:])
                        data[:, :, :, -1] = data[:, :, :, -2]
                    self.setdata(data)
        self.geometry.position_on_horizontal_grid = 'center'


class D3VirtualField(_D3CommonField):
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
        super(D3VirtualField, self).__init__(*args, **kwargs)
        if len(self.fieldset) != 0:
            if self.resource is not None or self.resource_fids != FPList():
                raise epygramError("You cannot set fieldset and (resource or resource_fids) at the same time.")

            def fieldGenerator(getdata):
                for field in self.fieldset:
                    yield field.fid, field

            self._fieldGenerator = fieldGenerator

            def getFieldByFid(fid, getdata):
                result = [field for field in self.fieldset if field.fid == fid]
                if len(result) != 1:
                    raise epygramError("There is no, or more than one, field(s) matching with this fid.")
                return result[0]

            self._getFieldByFid = getFieldByFid
            self._mode = 'fieldset'
        else:
            if self.resource is None or self.resource_fids == FPList():
                raise epygramError("If you do not set fieldset, you need to provide resource and resource_fids.")
            if all([k in self.resource.listfields() for k in self.resource_fids]):
                fidlist = self.resource_fids
            else:
                #There could be pattern to search for instead of true fid 
                fidlist = self.resource.find_fields_in_resource(seed=self.resource_fids,
                                                                fieldtype=['Point', 'V1D', 'V2D', 'H2D', '3D'])
            if len(fidlist) == 0:
                raise epygramError("There is no field in resource matching with resource_fids")

            def fieldGenerator(getdata):
                for fid in fidlist:
                    field = self.resource.readfield(fid, getdata=getdata)
                    yield fid, field

            self._fieldGenerator = fieldGenerator

            def getFieldByFid(fid, getdata):
                return self.resource.readfield(fid, getdata=getdata)

            self._getFieldByFid = getFieldByFid
            self._mode = 'resource'

        first = True
        self._fidList = []
        levelList = []
        levelIsArray = False

        for fid, field in self._fieldGenerator(getdata=False):
            if field.structure not in ['Point', 'V1D', 'V2D', 'H1D', 'H2D', '3D']:
                raise epygramError("3D virtual fields must be build from 'physical' fields only")
            if first:
                self._structure = field.structure
                self._geometry = field.geometry.deepcopy()  # We need a deepcopy as we modify it
                self._validity = field.validity.copy()
                if field.spectral_geometry is not None:
                    self._spectral_geometry = field.spectral_geometry.copy()
                else:
                    self._spectral_geometry = None
                self._processtype = field.processtype
            else:
                if self._structure != field.structure:
                    raise epygramError("All fields must share the structure")
                if self._geometry.structure != field.geometry.structure or \
                   self._geometry.name != field.geometry.name or \
                   self._geometry.grid != field.geometry.grid or \
                   self._geometry.dimensions != field.geometry.dimensions or \
                   self._geometry.position_on_horizontal_grid != field.geometry.position_on_horizontal_grid:
                    raise epygramError("All fields must share the horizontal geometry")
                if self._geometry.projected_geometry or field.geometry.projected_geometry or \
                   self._geometry.name == 'academic' or field.geometry.name == 'academic':
                    if self._geometry.projection != field.geometry.projection:
                        raise epygramError("All fields must share the geometry projection")
                if self._geometry.projected_geometry or field.geometry.projected_geometry:
                    if self._geometry.projection != field.geometry.projection or \
                       self._geometry.geoid != field.geometry.geoid:
                        raise epygramError("All fields must share the geometry geoid")
                if self._geometry.vcoordinate.typeoffirstfixedsurface != field.geometry.vcoordinate.typeoffirstfixedsurface or \
                   self._geometry.vcoordinate.position_on_grid != field.geometry.vcoordinate.position_on_grid:
                    raise epygramError("All fields must share the vertical geometry")
                if self._geometry.vcoordinate.grid is not None or field.geometry.vcoordinate.grid is not None:
                    if self._geometry.vcoordinate.grid != field.geometry.vcoordinate.grid:
                        raise epygramError("All fields must share the vertical grid")
                if self._validity != field.validity:
                    raise epygramError("All fields must share the validity")
                if self._spectral_geometry != field.spectral_geometry:
                    raise epygramError("All fields must share the spectral geometry")
                if self._processtype != field.processtype:
                    raise epygramError("All fields must share the processtype")
            if len(field.geometry.vcoordinate.levels) != 1:
                raise epygramError("fields must have only one level")
            levelIsArray = levelIsArray or isinstance(field.geometry.vcoordinate.levels[0], numpy.ndarray)
            if (not isinstance(field.geometry.vcoordinate.levels[0], numpy.ndarray)) and field.geometry.vcoordinate.levels[0] in levelList:
                raise epygramError("This level have already been found")
            levelList.append(field.geometry.vcoordinate.levels[0])
            if fid in self._fidList:
                raise epygramError("fields must have different fids")
            self._fidList.append(fid)

        kwargs_vcoord = dict(typeoffirstfixedsurface=self._geometry.vcoordinate.typeoffirstfixedsurface,
                             position_on_grid=self._geometry.vcoordinate.position_on_grid)
        if self._geometry.vcoordinate.grid is not None:
            kwargs_vcoord['grid'] = self._geometry.vcoordinate.grid
        if levelIsArray:
            kwargs_vcoord['levels'] = levelList
        else:
            kwargs_vcoord['levels'], self._fidList = (list(t) for t in zip(*sorted(zip(levelList, self._fidList))))  # TOBECHECKED
        newvcoordinate = VGeometry(**kwargs_vcoord)

        newstructure = {'3D': '3D',
                        'H2D': '3D',
                        'V1D': 'V1D',
                        'Point': 'V1D',
                        'V2D': 'V2D',
                        'H1D': 'V2D'}[field.structure]
        assert newstructure == self.structure, \
               "Individual fields structure do not match the field strucuture"

        self._geometry = self.geometry.deepcopy()
        self._geometry.vcoordinate = newvcoordinate
        self._spgpOpList = []

    def as_real_field(self, getdata=True):
        field3d = fpx.field(fid=dict(self.fid),
                            structure={'H2D':'3D',
                                       'Point':'V1D',
                                       'H1D':'V2D',
                                       'V1D':'V1D',
                                       'V2D':'V2D',
                                       '3D':'3D'}[self._structure],
                            geometry=self.geometry,
                            validity=self.validity,
                            spectral_geometry=self.spectral_geometry,
                            processtype=self.processtype)
        if getdata:
            field3d.setdata(self.getdata(d4=True))
        return field3d

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
        self._spectral_geometry = None
        if self._mode == 'resource':
            self._spgpOpList.append(('sp2gp', {}))
        else:
            for field in self.fieldset:
                field.sp2gp()

    def gp2sp(self, spectral_geometry):
        """
        Transforms the gridpoint field into spectral space, according to the
        *spectral_geometry* mandatorily passed as argument. Replaces data in
        place.

        :param spectral_geometry: instance of SpectralGeometry, actually
                                  containing spectral transform subroutine (in
                                  in its own *gp2sp()* method).
        """
        self._spectral_geometry = spectral_geometry
        if self._mode == 'resource':
            self._spgpOpList.append(('gp2sp', {'spectral_geometry':spectral_geometry}))
        else:
            for field in self.fieldset:
                field.gp2sp(spectral_geometry)

    def getdata(self, subzone=None, d4=False):
        """
        Returns the field data, with 3D shape if the field is not spectral,
        2D if spectral.

        :param subzone: optional, among ('C', 'CI'), for LAM fields only,
                        returns the data resp. on the C or C+I zone.
                        Default is no subzone, i.e. the whole field.
        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with
                     respect to geometry

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
        concat = numpy.concatenate
        arr = numpy.array
        missing = set()
        for k in range(len(self.geometry.vcoordinate.levels)):
            data = self.getlevel(k=k).getdata(subzone=subzone, d4=d4)
            if isinstance(data, numpy.ma.masked_array):
                concat = numpy.ma.concatenate
                arr = numpy.ma.array
                missing.add(data.fill_value)
            dataList.append(data)
        if d4:
            # vertical dimension already exists, and is the second one
            result = concat(dataList, axis=1)
        else:
            if len(self.validity) > 1:
                # vertical dimension does not exist and
                # must be the second one of the resulting array
                result = numpy.stack(dataList, axis=1)
            else:
                # vertical dimension does not exist and
                # must be the first one of the resulting array
                result = arr(dataList)
        if len(missing) == 1:
            result.set_fill_value(missing.pop())
        return result

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
        Take care (*i, j, k, t*) is python-indexing, ranging from 0 to
        dimension-1.

        :param i: the index in X direction
        :param j: the index in Y direction
        :param k: the index of the level (not a value in Pa or m...)
        :param t: the index of the temporal dimension (not a validity object)

        *k* and *t* can be scalar even if *i* and *j* are arrays.

        If *one* is False, returns [value] instead of value.
        """
        if t is None:
            if len(self.validity) > 1:
                raise epygramError("*t* is mandatory when there are several validities")
            t = 0
        if k is None:
            if self.geometry.datashape['k']:
                raise epygramError("*k* is mandatory when field has a vertical coordinate")
            k = 0
        if j is None:
            if self.geometry.datashape['j']:
                raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
            j = 0
        if i is None:
            if self.geometry.datashape['i']:
                raise epygramError("*i* is mandatory when field has one horizontal dimension")
            i = 0

        i, j = as_numpy_array(i).flatten(), as_numpy_array(j).flatten()
        k, t = as_numpy_array(k).flatten(), as_numpy_array(t).flatten()

        if not numpy.all(self.geometry.point_is_inside_domain_ij(i, j)):
            raise ValueError("point is out of field domain.")

        sizes = set([len(x) for x in [i, j, k, t]])
        if len(sizes) > 2 or (len(sizes) == 2 and 1 not in sizes):
            raise epygramError("each of i, j, k and t must be scalar or have the same length as the others")

        if max(sizes) > 1 and len(k) == 1:
            k = k.repeat(max(sizes))

        oldk = k[0]
        data = self.getlevel(k=oldk).getdata(d4=True)
        value = numpy.ndarray((max(sizes),), dtype=data.dtype)

        for thisk in numpy.unique(k):
            if thisk != oldk:
                data = self.getlevel(k=thisk).getdata(d4=True)
                oldk = thisk
            mask = k == thisk
            value[mask] = data[t[mask] if len(t) > 1 else t,
                               0,
                               j[mask] if len(j) > 1 else j,
                               i[mask] if len(i) > 1 else i]

        if value.size == 1 and one:
            value = value.item()
        return value

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.

        :param level: the requested level expressed in coordinate value (Pa, m...)
        :param k: the index of the requested level
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

        result = self._getFieldByFid(self._fidList[my_k], True)
        
        #update vccordinate in case of this attribute was changed in the geometry
        vcoord = self.geometry.vcoordinate.deepcopy()
        levels = vcoord.levels[my_k]
        for _ in range(len(vcoord.levels)):
            vcoord.levels.pop()
        vcoord.levels.append(levels)
        result.geometry.vcoordinate = vcoord

        for op, kwargs in self._spgpOpList:
            if op == 'sp2gp':
                result.sp2gp(**kwargs)
            elif op == 'gp2sp':
                result.gp2sp(**kwargs)
            else:
                raise epygramError("operation not known")

        return result


footprints.collectors.get(tag='fields').fasttrack = ('structure',)
