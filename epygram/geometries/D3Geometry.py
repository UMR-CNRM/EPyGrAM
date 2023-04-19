#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for 3D geometries of fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy
import math
import copy
import sys
import re

import footprints
from footprints import FootprintBase, FPDict, FPList, proxy as fpx
from bronx.graphics.axes import set_figax
from bronx.syntax.arrays import stretch_array
from bronx.syntax.decorators import nicedeco

from epygram import epygramError, config
from epygram.config import rounding_decimal as _rd
from epygram.util import (RecursiveObject, degrees_nearest_mod, Angle,
                          positive_longitudes, longitudes_between_minus180_180,
                          separation_line, write_formatted,
                          nearlyEqual,
                          as_numpy_array, moveaxis,
                          is_scalar, Comparator)

from .VGeometry import VGeometry

epylog = footprints.loggers.getLogger(__name__)
_re_nearest_sq = re.compile('(?P<n>\d+)\*(?P<m>\d+)')


@nicedeco
def _need_pyproj_geod(mtd):
    """
    Decorator for Geometry object: if the method needs a pyproj.Geod
    object to be set.
    """
    def with_geod(self, *args, **kwargs):
        if not hasattr(self, '_pyproj_geod'):
            self._set_geoid()
        return mtd(self, *args, **kwargs)
    return with_geod


class D3Geometry(RecursiveObject, FootprintBase):
    """
    Handles the geometry for a 3-Dimensions Field.
    Abstract mother class.
    """

    _abstract = True
    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                info="Type of geometry.",
                values=set(['3D'])),
            name=dict(
                info="Name of geometrical type of representation of points on \
                      the Globe.",
                values=set(['lambert', 'mercator', 'polar_stereographic',
                            'regular_lonlat', 'rotated_lonlat',
                            'rotated_reduced_gauss', 'reduced_gauss', 'regular_gauss',
                            'unstructured'])),
            grid=dict(
                type=FPDict,
                info="Handles description of the horizontal grid."),
            dimensions=dict(
                type=FPDict,
                info="Handles grid dimensions."),
            vcoordinate=dict(
                access='rwx',
                type=VGeometry,
                info="Handles vertical geometry parameters."),
            position_on_horizontal_grid=dict(
                optional=True,
                access='rwx',
                info="Position of points w/r to the horizontal.",
                default='__unknown__',
                values=set(['upper-right', 'upper-left',
                            'lower-left', 'lower-right',
                            'center-left', 'center-right',
                            'lower-center', 'upper-center',
                            'center', '__unknown__'])),
            geoid=dict(
                type=FPDict,
                optional=True,
                default=FPDict(config.default_geoid),
                info="To specify geoid shape.")
        )
    )

    # ghost attributes are ignored when comparing 2 objects between them
    _ghost_attributes = RecursiveObject._ghost_attributes + ['_puredict', '_observer']  # footprints special attributes

    def __init__(self, *args, **kwargs):
        super(D3Geometry, self).__init__(*args, **kwargs)
        # Checks !
        self._consistency_check()

    def _consistency_check(self):
        # implemented in child classes
        pass

    @property
    def rectangular_grid(self):
        """ Is the grid rectangular ? """
        return isinstance(self, D3RectangularGridGeometry)

    @property
    def projected_geometry(self):
        """ Is the geometry a projection ? """
        return 'projection' in self._attributes

    @property
    def datashape(self):
        """Returns the data shape requested by this geometry."""
        if self.structure == "3D":
            return {'k': True, 'j':True, 'i':True}
        elif self.structure == "V2D":
            return {'k': True, 'j':False, 'i':True}
        elif self.structure == "H2D":
            return {'k': False, 'j':True, 'i':True}
        elif self.structure == "H1D":
            return {'k':False, 'j':False, 'i':True}
        elif self.structure == "V1D":
            return {'k': True, 'j':False, 'i':False}
        elif self.structure == "Point":
            return {'k': False, 'j':False, 'i':False}
        else:
            raise epygramError('This structure is unknown')

    def get_levels(self, d4=False, nb_validities=0, subzone=None):
        """
        Returns an array containing the level for each data point.

        :param d4: - if True,  returned values are shaped in a 4 dimensions array
                   - if False, shape of returned values is determined with respect to geometry
                     d4=True requires nb_validities > 0
        :param nb_validities: the number of validities represented in data values
        :param subzone: optional, among ('C', 'CI'), for LAM grids only, extracts
          the levels resp. from the C or C+I zone off the C+I(+E) zone.

        Levels are internally stored with the vertical dimension first whereas this
        method puts the time in first dimension.
        """
        if d4:
            if nb_validities < 1:
                raise ValueError("nb_validities must be >=1 when d4==True")

        levels = as_numpy_array(self.vcoordinate.levels)

        # We add the horizontal axis
        h_shape2D = self.get_datashape(d4=True,
                                       dimT=nb_validities,
                                       force_dimZ=1,
                                       subzone=subzone)[-2:]
        if len(levels.shape) == 1:
            # level values constant over the horizontal domain and time
            # We add the horizontal dimension
            original_has_time = False
            if self.datashape['j'] or d4:
                shape = tuple(list(levels.shape) + [h_shape2D[0]])
                levels = levels.repeat(h_shape2D[0]).reshape(shape)
            if self.datashape['i'] or d4:
                shape = tuple(list(levels.shape) + [h_shape2D[1]])
                levels = levels.repeat(h_shape2D[1]).reshape(shape)
        else:
            # level values with horizontal variations
            h_shape = self.get_datashape(force_dimZ=2)[1:]  # workaround for inconsistency between rectangular and gauss, self.get_datashape(force_dimZ=1)
            if len(levels.shape) == 1 + len(h_shape):
                original_has_time = False
            elif len(levels.shape) == 2 + len(h_shape):
                original_has_time = True
            else:
                raise epygramError("Wrong number of dimensions")
            if levels.shape[len(levels.shape) - len(h_shape):] != h_shape:  # OK for h_shape=tuple()
                raise epygramError("Shape of self.vcoordinate.levels does not agree with horizontal dimensions")
            if subzone is not None:
                levels = self.extract_subzone(levels, subzone)
            if d4 and ((not self.datashape['i']) or (not self.datashape['j'])):
                shape = levels.shape[:len(levels.shape) - len(h_shape)]  # shape without the horizontal dimensions, OK for h_shape=tuple()
                shape = tuple(list(shape) + list(h_shape2D))  # shape with the new horizontal dimensions
                levels = levels.reshape(shape)
        # We suppress the vertical dimension if we do not need it
        if len(self.vcoordinate.levels) == 1 and not d4 and nb_validities <= 1:
            levels = levels[0]

        # We add the time axis
        if original_has_time:
            if levels.shape[1] != nb_validities:
                raise epygramError("Shape of self.vcoordinate.levels does not agree with nb_validities")
            shape_range = range(len(levels.shape))
            levels = levels.transpose(tuple([1, 0] + list(shape_range[2:])))
        elif d4 or nb_validities >= 2:
            shape = tuple(list(levels.shape) + [nb_validities])
            levels = levels.repeat(nb_validities).reshape(shape)
            shape_range = list(range(len(shape)))
            levels = levels.transpose(tuple([shape_range[-1]] + list(shape_range[0:-1])))  # last axis in first

        return levels

    def vcoord_as_field(self, surface_type, validity=None, levels=None):
        """
        Returns a field filled with the level values associated to a fake geometry
        :param validity: validities to associate with the returned field
                         if None, we try without setting validity
        :param surface_type: typeOfFirstFixedSurface to associate to the fake geometry
        :param levels: list of values to use as the levels of the fake geometry
                       if None, levels will be replaced by a range
        """

        geometry = self.deepcopy()
        vgeom_kwargs = copy.deepcopy(self.vcoordinate._attributes)
        vgeom_kwargs['typeOfFirstFixedSurface'] = surface_type
        vgeom_kwargs['levels'] = range(len(self.vcoordinate.levels)) if levels is None else levels
        geometry.vcoordinate = fpx.geometry(**vgeom_kwargs)
        field_kwargs = dict(fid=dict(temp='temp'),
                            structure=geometry.structure,
                            geometry=geometry)
        if validity is not None:
            field_kwargs['validity'] = validity
        field = fpx.field(**field_kwargs)
        field.setdata(self.get_levels(d4=True, nb_validities=1 if validity is None else len(validity)))
        return field

    def _getoffset(self, position=None):
        """Returns the offset to use for this **position**."""
        if position is not None:
            pos = position
        else:
            pos = self.position_on_horizontal_grid
        if pos == '__unknown__' or pos is None:
            raise epygramError("position_on_horizontal_grid must be" +
                               " defined.")

        return {'upper-right' : (.5, .5),
                'upper-left'  : (-.5, .5),
                'lower-left'  : (-.5, -.5),
                'lower-right' : (.5, -.5),
                'center-left' : (-.5, 0.),
                'center-right': (.5, 0.),
                'lower-center': (0., -.5),
                'upper-center': (0., .5),
                'center'      : (0., 0.)}[pos]

    def _set_geoid(self):
        import pyproj
        self._pyproj_geod = pyproj.Geod(**self.geoid)

    @_need_pyproj_geod
    def distance(self, end1, end2):
        """
        Computes the distance between two points along a Great Circle.

        :param end1: first point, must be a tuple (lon, lat) in degrees.
        :param end2: second point, must be a tuple (lon, lat) in degrees.

        If one of (end1, end2) is a tuple of arrays (like (<array of lon>, <array of lat>))
        the other one must be a scalar tuple or a tuple of arrays with same dimensions.

        Warning: requires the :mod:`pyproj` module.
        """
        scalar = all([is_scalar(val) for val in [end1[0], end1[1], end2[0], end2[1]]])
        end1_0 = as_numpy_array(end1[0]).flatten()
        end1_1 = as_numpy_array(end1[1]).flatten()
        end2_0 = as_numpy_array(end2[0]).flatten()
        end2_1 = as_numpy_array(end2[1]).flatten()
        assert len(end1_0) == len(end1_1) and len(end2_0) == len(end2_1), 'pb with dims'
        assert len(end1_0.shape) == len(end1_1.shape) == len(end2_0.shape) == len(end2_1.shape) == 1, 'only scalars or 1D arrays'

        if len(end1_0) == len(end2_0):
            pass
        else:
            if len(end1_0) == 1:
                end1_0 = numpy.repeat(end1_0, len(end2_0), axis=0)
                end1_1 = numpy.repeat(end1_1, len(end2_0), axis=0)
            elif len(end2_0) == 1:
                end2_0 = numpy.repeat(end2_0, len(end1_0), axis=0)
                end2_1 = numpy.repeat(end2_1, len(end1_0), axis=0)
            else:
                raise epygramError('at least one point must be fixed or both arrays must have the same length')
        if True in [isinstance(a, numpy.ma.MaskedArray) for a in [end1_0, end1_1, end2_0, end2_1]]:
            # inv method does not like masked arrays
            # and can raise a ValueError: undefined inverse geodesic (may be an antipodal point)
            # on masked points
            mask = numpy.logical_not(numpy.ma.getmaskarray(end1_0 + end1_1 + end2_0 + end2_1))
            distance = numpy.empty_like(end1_0, dtype=float)
            distance[mask] = self._pyproj_geod.inv(end1_0[mask], end1_1[mask], end2_0[mask], end2_1[mask])[2]
        else:
            distance = self._pyproj_geod.inv(end1_0, end1_1, end2_0, end2_1)[2]
        return distance[0] if scalar else distance

    @_need_pyproj_geod
    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up along a Great Circle.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        :param num: the number of points, including point1 and point2.

        Warning: requires the :mod:`pyproj` module.
        """
        if num < 2:
            raise epygramError("'num' must be at least 2.")
        transect = self._pyproj_geod.npts(end1[0], end1[1],
                                          end2[0], end2[1],
                                          num - 2)
        transect.insert(0, end1)
        transect.append(end2)
        return transect

    @_need_pyproj_geod
    def azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points following a Great Circle.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.

        Warning: requires the :mod:`pyproj` module.
        """
        return self._pyproj_geod.inv(end1[0], end1[1], end2[0], end2[1])[0]

    def make_point_geometry(self, lon, lat):
        """Returns a PointGeometry at coordinates *(lon,lat)* in degrees."""
        vcoordinate = fpx.geometry(structure='V',
                                   typeoffirstfixedsurface=255,
                                   levels=[0])
        return fpx.geometry(structure='Point',
                            name='unstructured',
                            vcoordinate=vcoordinate,
                            dimensions={'X':1, 'Y':1},
                            grid={'longitudes':[lon],
                                  'latitudes':[lat]},
                            position_on_horizontal_grid='center'
                            )

    def make_profile_geometry(self, lon, lat):
        """Returns a V1DGeometry at coordinates *(lon,lat)* in degrees."""
        vcoordinate = fpx.geometry(structure='V',
                                   typeoffirstfixedsurface=255,
                                   levels=[])
        return fpx.geometry(structure='V1D',
                            name='unstructured',
                            vcoordinate=vcoordinate,
                            dimensions={'X':1, 'Y':1},
                            grid={'longitudes':[lon],
                                  'latitudes':[lat]},
                            position_on_horizontal_grid='center'
                            )

    def make_section_geometry(self, end1, end2,
                              points_number=None,
                              resolution=None,
                              position=None):
        """
        Returns a V2DGeometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param position: defines the position of data in the grid (for projected grids only)
        """
        assert isinstance(end1, tuple) or isinstance(end1, list)
        assert isinstance(end2, tuple) or isinstance(end2, list)
        if position not in [None, 'center']:
            raise epygramError("position can only be None or 'center' for non-projected geometries.")
        if resolution is not None and points_number is not None:
            raise epygramError("only one of resolution and " +
                               " points_number can be given.")

        distance = self.distance(end1, end2)
        if resolution is None and points_number is None:
            resolution = 0.5 * (self.resolution_ll(*end1) +
                                self.resolution_ll(*end2))
            if resolution > distance:
                raise epygramError("'ends' are too near: pure" +
                                   " interpolation between two gridpoints.")
        elif points_number is not None and points_number < 2:
            raise epygramError("'points_number' must be at least 2.")
        if resolution is not None:
            points_number = int(numpy.rint(distance / resolution)) + 1
        if points_number >= 3:
            transect = self.linspace(end1, end2, points_number)
        elif points_number == 2:
            transect = [end1, end2]
        else:
            raise epygramError("cannot make a section with less than" +
                               " 2 points.")
        vcoordinate = fpx.geometry(structure='V',
                                   typeoffirstfixedsurface=255,
                                   levels=[])
        kwargs_geom = {'structure':'V2D',
                       'name':'unstructured',
                       'vcoordinate':vcoordinate,
                       'dimensions':{'X':len(transect), 'Y':1},
                       'grid':{'longitudes':[p[0] for p in transect],
                               'latitudes':[p[1] for p in transect]},
                       'position_on_horizontal_grid':'center' if position is None else position}
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return fpx.geometry(**kwargs_geom)

    def make_physicallevels_geometry(self):
        """
        Returns a new geometry excluding levels with no physical meaning.

        :param getdata: if False returns a field without data
        """
        geometry = self.deepcopy()
        if self.vcoordinate.typeoffirstfixedsurface == 118:
            # build vertical geometry
            kwargs_vcoord = {'structure':'V',
                             'typeoffirstfixedsurface': self.vcoordinate.typeoffirstfixedsurface,
                             'position_on_grid': self.vcoordinate.position_on_grid,
                             'grid': copy.copy(self.vcoordinate.grid),
                             'levels': copy.copy(self.vcoordinate.levels)}
            # Suppression of levels above or under physical domain
            for level in list(kwargs_vcoord['levels']):
                if level < 1 or level > len(self.vcoordinate.grid['gridlevels']) - 1:
                    kwargs_vcoord['levels'].remove(level)
            # build geometry
            geometry.vcoordinate = fpx.geometry(**kwargs_vcoord)
        return geometry

    def make_zoom_geometry(self, zoom, extra_10th=False):
        """
        Returns an unstructured geometry with the points contained in *zoom*.

        :param zoom: a dict(lonmin=, lonmax=, latmin=, latmax=).
        :param extra_10th: if True, add 1/10th of the X/Y extension of the zoom
                           (only usefull for the regular_lonlat grid implementation).
        """
        (lons, lats) = self.get_lonlat_grid()
        kwargs_zoomgeom = {'structure':self.structure,
                           'vcoordinate':self.vcoordinate.deepcopy(),
                           'position_on_horizontal_grid':self.position_on_horizontal_grid,
                           'geoid':self.geoid}
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
        return fpx.geometry(**kwargs_zoomgeom)

    def _reshape_lonlat_4d(self, lons, lats, nb_validities):
        """Make lons, lats grids 4D."""
        if nb_validities < 1:
            raise ValueError("nb_validities must be >=1 when d4==True")
        # We add vertical dimension, and missing horizontal dimension
        shape = list(lons.shape)
        if len(shape) == 1:
            shape = [1] + shape
        shape = tuple(shape + [len(self.vcoordinate.levels)])
        lons = lons.repeat(len(self.vcoordinate.levels)).reshape(shape)
        lats = lats.repeat(len(self.vcoordinate.levels)).reshape(shape)
        shape_range = list(range(len(shape)))
        lons = lons.transpose(tuple([shape_range[-1]] + list(shape_range[0:-1])))  # last axis in first
        lats = lats.transpose(tuple([shape_range[-1]] + list(shape_range[0:-1])))  # last axis in first
        # We add validities
        shape = tuple(list(lons.shape) + [nb_validities])
        lons = lons.repeat(nb_validities).reshape(shape)
        lats = lats.repeat(nb_validities).reshape(shape)
        shape_range = list(range(len(shape)))
        lons = lons.transpose(tuple([shape_range[-1]] + list(shape_range[0:-1])))  # last axis in first
        lats = lats.transpose(tuple([shape_range[-1]] + list(shape_range[0:-1])))  # last axis in first
        return lons, lats

    def eq_Hgeom(self, other):
        """
        Tests if the horizontal part of the geometry is equal to
        the horizontal part of another geometry.
        :param: other: other geometry to use in the comparison
        """
        assert isinstance(other, D3Geometry), "Other must be a geometry object"
        vc_self = self.vcoordinate
        self.vcoordinate = other.vcoordinate
        result = self == other
        self.vcoordinate = vc_self
        return result

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def make_field(self, fid=None):
        """Make a field out of the geometry."""
        if fid is None:
            fid = self.name
        ll_field = fpx.field(structure='H2D',
                             geometry=self,
                             fid=FPDict({'geometry': fid}))
        data = numpy.zeros((self.dimensions['Y'], self.dimensions['X']))
        data[1:-1, 1:-1, ...] = 1.
        ll_field.setdata(data)
        return ll_field

    def plotgeometry(self,
                     plotlib='cartopy',
                     **kwargs):
        """
        Makes a simple plot of the geometry, with a number of options.

        :param plotlib: library to be used for plotting: 'basemap' is DEPRECATED;
            'cartopy' (default) is recommended !
        """
        if plotlib == 'cartopy':
            return self.cartoplot_geometry(**kwargs)
        else:
            raise NotImplementedError("'basemap' plotlib has been removed, only remains 'cartopy'")

    def cartoplot_geometry(self, **kwargs):
        """
        Makes a simple plot of the geometry, using cartopy.
        For kwargs please refer to epygram.geometries.domain_making.output.cartoplot_rect_geometry
        """
        from epygram.geometries.domain_making.output import cartoplot_rect_geometry
        if 'color' in kwargs:  # compatibility dirty-fix
            kwargs['contourcolor'] = kwargs.pop('color')
        fig, ax = cartoplot_rect_geometry(self, **kwargs)
        return fig, ax

    def what(self, out=sys.stdout,
             vertical_geometry=True,
             arpifs_var_names=False,
             spectral_geometry=None):
        """
        Writes in file a summary of the geometry.

        :param out: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        :param vertical_geometry: if True, writes the vertical geometry of the
          field.
        :param arpifs_var_names: if True, prints the equivalent 'arpifs' variable
          names.
        :param spectral_geometry: an optional dict containing the spectral
          truncatures {'in_X':, 'in_Y':}  (LAM) or {'max':} (global).
        """
        out.write("###########################\n")
        out.write("### HORIZONTAL GEOMETRY ###\n")
        out.write("###########################\n")
        write_formatted(out, "Rectangular grid ( = LAM or reg. Lon/Lat)",
                        self.rectangular_grid)

        self._what_grid_dimensions(out, spectral_geometry=spectral_geometry)
        if self.rectangular_grid:
            if self.projected_geometry:
                self._what_projection(out, arpifs_var_names=arpifs_var_names)
            elif self.name == 'regular_lonlat':
                self._what_grid(out, arpifs_var_names=arpifs_var_names)
            elif self.name == 'rotated_lonlat':
                self._what_grid(out)
            elif self.name == 'academic':
                self._what_position(out)
            elif self.name == 'unstructured':
                self._what_grid(out)
        out.write(separation_line)
        out.write("\n")

        if vertical_geometry:
            self.vcoordinate.what(out)


class D3RectangularGridGeometry(D3Geometry):
    """
    Handles the geometry for a rectangular 3-Dimensions Field.
    Abstract.
    """

    _abstract = True
    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic',
                            'regular_lonlat', 'rotated_lonlat',
                            'academic', 'unstructured']))
        )
    )

    @property
    def isglobal(self):
        """
        :return: True if geometry is global
        """
        return False  # Apart for global lon-lat

    def _get_grid(self, indextype, subzone=None, position=None):
        """
        Returns a tuple of two tables containing the two indexes of each
        point, with 2D shape.

        :param indextype: either 'ij', 'xy' or 'll' to get
          i,j indexes, x,y coordinates or lon,lat coordinates
        :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        :param position: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        igrid, jgrid = numpy.meshgrid(range(self.dimensions['X']),
                                      range(self.dimensions['Y']))
        igrid = igrid.flatten()
        jgrid = jgrid.flatten()
        if indextype == 'xy':
            (a, b) = self.ij2xy(igrid, jgrid, position)
        elif indextype == 'll':
            (a, b) = self.ij2ll(igrid, jgrid, position)
        elif indextype == 'ij':
            pass
        else:
            raise ValueError('*indextype*== ' + indextype)
        a = self.reshape_data(a)
        b = self.reshape_data(b)
        if subzone and self.grid.get('LAMzone') is not None:
            a = self.extract_subzone(a, subzone)
            b = self.extract_subzone(b, subzone)

        return (a, b)

    @property
    def gridpoints_number(self, subzone=None):
        """
        Returns the number of gridpoints of the grid.

        :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        """
        shp = self.get_datashape(dimT=1, force_dimZ=1, subzone=subzone)
        return shp[0] * shp[1]

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
                   - if False, shape of returned values is determined with
                     respect to geometry. d4=True requires nb_validities > 0
        :param nb_validities: number of validities represented in data values
        :param force_longitudes: if 'positive', the longitudes will be forced positive
                                 if ']-180,180]', the longitudes will be in the ]-180, 180] interval

        Shape of 2D data on Rectangular grids: \n
          - grid[0,0] is SW, grid[-1,-1] is NE \n
          - grid[0,-1] is SE, grid[-1,0] is NW
        """
        (lons, lats) = self._get_grid('ll', subzone=subzone, position=position)
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

    def extract_subzone(self, data, subzone):
        """
        Extracts the subzone C or CI from a LAM field.

        :param data: the data values with shape concording with geometry.
        :param subzone: optional, among ('C', 'CI'), for LAM grids only,
          extracts the data resp. from the C or C+I zone off the C+I(+E) zone.
        """
        if subzone not in ('C', 'CI'):
            raise epygramError("only possible values for 'subzone' are 'C'" +
                               " or 'CI'.")
        if 'LAMzone' not in self.grid:
            raise epygramError("method for LAM grids only.")

        selectionE = []  # To remove E-zone
        selectionI = []  # To eventually remove I-zone
        if len(data.shape) == 4:
            selectionE.extend([slice(None)] * 2)
            selectionI.extend([slice(None)] * 2)
        elif len(data.shape) == 3:
            selectionE.append(slice(None))
            selectionI.append(slice(None))
        if self.datashape['j']:
            y1 = self.dimensions.get('Y_CIoffset', 0)
            y2 = self.dimensions.get('Y_CIoffset', 0) + self.dimensions['Y_CIzone']
            selectionE.append(slice(y1, y2))
            y1 = self.dimensions['Y_Iwidth']
            y2 = -self.dimensions['Y_Iwidth']
            selectionI.append(slice(y1, y2))
        if self.datashape['i']:
            x1 = self.dimensions.get('X_CIoffset', 0)
            x2 = self.dimensions.get('X_CIoffset', 0) + self.dimensions['X_CIzone']
            selectionE.append(slice(x1, x2))
            x1 = self.dimensions['X_Iwidth']
            x2 = -self.dimensions['X_Iwidth']
            selectionI.append(slice(x1, x2))

        if self.grid['LAMzone'] == 'CIE':
            # remove E-zone
            edata = data[tuple(selectionE)]
        else:
            edata = data
        if subzone == 'C':
            # remove I-zone
            edata = edata[tuple(selectionI)]
        return edata

    def make_subarray_geometry(self,
                               first_i, last_i,
                               first_j, last_j):
        """
        Make a modified geometry consisting in a subarray of the grid, defined
        by the indexes given as argument.
        """
        if first_i < 0 or first_i >= self.dimensions['X'] or \
           last_i < 0 or last_i >= self.dimensions['X'] or \
           first_j < 0 or first_j >= self.dimensions['Y'] or \
           last_j < 0 or last_j >= self.dimensions['Y']:
            raise epygramError("first_i, last_i, first_j and last_j must be inside the geometry")

        geom_kwargs = copy.deepcopy(self._attributes)
        geom_kwargs.pop('dimensions')
        if 'LAMzone' in geom_kwargs['grid']:
            geom_kwargs['grid']['LAMzone'] = None
        if 'input_position' in geom_kwargs['grid']:
            coords_00 = self.ij2ll(first_i, first_j, position='center')
            geom_kwargs['grid']['input_position'] = (0, 0)
            input_lon = Angle(coords_00[0], 'degrees') if isinstance(self.grid['input_lon'], Angle) else coords_00[0]
            input_lat = Angle(coords_00[1], 'degrees') if isinstance(self.grid['input_lat'], Angle) else coords_00[1]
            geom_kwargs['grid']['input_lon'] = input_lon
            geom_kwargs['grid']['input_lat'] = input_lat
        geom_kwargs['dimensions'] = {'X':last_i - first_i + 1, 'Y':last_j - first_j + 1}
        newgeom = fpx.geometry(**geom_kwargs)  # create new geometry object
        return newgeom

    def make_subsample_geometry(self, sample_x, sample_y, sample_z):
        """
        Make a sample geometry by decreasing resolution.
        :param sample_x: take one over <sample_x> points in the x direction
        :param sample_y: same for the y direction
        :param sample_z: same for the z direction

        CAUTION: if your grid contains non physical point, these points can be
        retain in subsample. Use select_zone beforehand to suppress these points.
        """
        assert isinstance(sample_x, int) and \
               isinstance(sample_y, int) and \
               isinstance(sample_z, int), "sample_x, y and z must be integers"
        geom_kwargs = copy.deepcopy(self._attributes)
        if 'LAMzone' in geom_kwargs['grid']:
            geom_kwargs['grid']['LAMzone'] = None
        if 'input_position' in geom_kwargs['grid']:
            coords_00 = self.ij2ll(0, 0, position='center')
            geom_kwargs['grid']['input_position'] = (0, 0)
            input_lon = Angle(coords_00[0], 'degrees') if isinstance(self.grid['input_lon'], Angle) else coords_00[0]
            input_lat = Angle(coords_00[1], 'degrees') if isinstance(self.grid['input_lat'], Angle) else coords_00[1]
            geom_kwargs['grid']['input_lon'] = input_lon
            geom_kwargs['grid']['input_lat'] = input_lat
        geom_kwargs['dimensions'] = {'X':int(math.ceil(geom_kwargs['dimensions']['X'] / float(sample_x))),
                                     'Y':int(math.ceil(geom_kwargs['dimensions']['Y'] / float(sample_y)))}
        geom_kwargs['grid']['X_resolution'] = geom_kwargs['grid']['X_resolution'] * sample_x
        geom_kwargs['grid']['Y_resolution'] = geom_kwargs['grid']['Y_resolution'] * sample_y
        levels = geom_kwargs['vcoordinate'].levels[::sample_z]
        del geom_kwargs['vcoordinate'].levels[:]
        geom_kwargs['vcoordinate'].levels.extend(levels)
        newgeom = fpx.geometry(**geom_kwargs)  # create new geometry object
        return newgeom

    def get_datashape(self,
                      dimT=1,
                      force_dimZ=None,
                      d4=False,
                      subzone=None):
        """
        Returns the data shape according to the geometry.

        :param force_dimZ: if supplied, force the Z dimension instead of that
          of the vertical geometry
        :param dimT: if supplied, is the time dimension to be added to the
          data shape
        :param d4: - if True,  shape is 4D
                   - if False, shape has only those > 1
        :param subzone: optional, among ('C', 'CI'), for LAM grids only, informes that
          data is resp. on the C or C+I zone off the C+I(+E) zone.
        """
        if subzone is None:
            dimX = self.dimensions['X']
            dimY = self.dimensions['Y']
        else:
            assert self.grid.get('LAMzone', False), \
                "*subzone* cannot be requested for this geometry"
            assert subzone in ('C', 'CI')
            if self.grid['LAMzone'] == 'CIE':
                if subzone == 'CI':
                    dimX = self.dimensions['X_CIzone']
                    dimY = self.dimensions['Y_CIzone']
                elif subzone == 'C':
                    dimX = self.dimensions['X_CIzone'] - 2 * self.dimensions['X_Iwidth']
                    dimY = self.dimensions['Y_CIzone'] - 2 * self.dimensions['Y_Iwidth']
            elif self.grid['LAMzone'] == 'CI':
                dimX = self.dimensions['X'] - 2 * self.dimensions['X_Iwidth']
                dimY = self.dimensions['Y'] - 2 * self.dimensions['Y_Iwidth']

        if force_dimZ is not None:
            dimZ = force_dimZ
        else:
            dimZ = len(self.vcoordinate.levels)
        if d4:
            shape = [dimT, dimZ, dimY, dimX]
        else:
            shape = []
            if dimT > 1:
                shape.append(dimT)
            if dimZ > 1:
                shape.append(dimZ)
            if self.datashape['j']:
                shape.append(dimY)
            if self.datashape['i']:
                shape.append(dimX)
        return tuple(shape)

    def reshape_data(self, data,
                     first_dimension=None,
                     d4=False,
                     subzone=None):
        """
        Returns a 2D data (horizontal dimensions) reshaped from 1D,
        according to geometry.

        :param data: the 1D data (or 3D with a T and Z dimensions,
          or 2D with either a T/Z dimension, to be specified),
          of dimension concording with geometry. In case data is 3D, T must be
          first dimension and Z the second.
        :param first_dimension: in case data is 2D, specify what is the first
          dimension of data among ('T', 'Z')
        :param subzone: optional, among ('C', 'CI'), for LAM grids only, informes that
          data is resp. on the C or C+I zone off the C+I(+E) zone.
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
               "vertical dimension of data must be 1 or self.vcoordinate.levels=" \
               + str(self.vcoordinate.levels)

        if d4 or 1 not in (nb_validities, nb_levels):  # data as 4D or truly 4D
            shp = self.get_datashape(dimT=nb_validities, force_dimZ=nb_levels,
                                     d4=True, subzone=subzone)
        else:
            shp = self.get_datashape(dimT=nb_validities, force_dimZ=nb_levels,
                                     subzone=subzone)
        return data.reshape(shp)

    def horizontally_flattened(self, data):
        """
        Returns a copy of *data* with horizontal dimensions flattened.
        *data* must be 4D for simplicity reasons.
        """
        assert len(data.shape) == 4
        data3D = numpy.empty(tuple(list(data.shape[:2]) + [self.gridpoints_number]))
        for t in range(data.shape[0]):
            for k in range(data.shape[1]):
                data3D[t, k, :] = data[t, k, :, :].flatten()
        return data3D

    def fill_maskedvalues(self, data, fill_value=None):
        """
        Returns a copy of *data* with masked values filled with *fill_value*.
        """
        assert isinstance(data, numpy.ma.masked_array)
        return data.filled(fill_value)

    def gimme_corners_ij(self, subzone=None):
        """
        Returns the indices (i, j) of the four corners of a rectangular grid,
        as a dict(corner=(i, j)) with corner in: \n
        ll = lower-left / lr = lower-right / ur = upper-right / ul = upper-left.

        (0, 0) is always the lower-left corner of the grid.

        :param subzone: for LAM fields, returns the corners of the subzone.
        """
        if 'LAMzone' not in self.grid or self.grid['LAMzone'] is None:
            ll = (0, 0)
            lr = (self.dimensions['X'] - 1, 0)
            ul = (0, self.dimensions['Y'] - 1)
            ur = (self.dimensions['X'] - 1, self.dimensions['Y'] - 1)
        elif self.grid['LAMzone'] == 'CIE':
            if subzone in (None, 'CIE'):
                ll = (0, 0)
                lr = (self.dimensions['X'] - 1, 0)
                ul = (0, self.dimensions['Y'] - 1)
                ur = (self.dimensions['X'] - 1, self.dimensions['Y'] - 1)
            elif subzone == 'CI':
                ll = (self.dimensions['X_CIoffset'],
                      self.dimensions['Y_CIoffset'])
                lr = (self.dimensions['X_CIoffset'] + self.dimensions['X_CIzone'] - 1,
                      self.dimensions['Y_CIoffset'])
                ul = (self.dimensions['X_CIoffset'],
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_CIzone'] - 1)
                ur = (self.dimensions['X_CIoffset'] + self.dimensions['X_CIzone'] - 1,
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_CIzone'] - 1)
            elif subzone == 'C':
                ll = (self.dimensions['X_CIoffset'] + self.dimensions['X_Iwidth'],
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_Iwidth'])
                lr = (self.dimensions['X_CIoffset'] + self.dimensions['X_Iwidth'] + self.dimensions['X_Czone'] - 1,
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_Iwidth'])
                ul = (self.dimensions['X_CIoffset'] + self.dimensions['X_Iwidth'],
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_Iwidth'] + self.dimensions['Y_Czone'] - 1)
                ur = (self.dimensions['X_CIoffset'] + self.dimensions['X_Iwidth'] + self.dimensions['X_Czone'] - 1,
                      self.dimensions['Y_CIoffset'] + self.dimensions['Y_Iwidth'] + self.dimensions['Y_Czone'] - 1)
        elif self.grid['LAMzone'] == 'CI':
            if subzone in (None, 'CI'):
                ll = (0, 0)
                lr = (self.dimensions['X_CIzone'] - 1, 0)
                ul = (0, self.dimensions['Y_CIzone'] - 1)
                ur = (self.dimensions['X_CIzone'] - 1, self.dimensions['Y_CIzone'] - 1)
            elif subzone == 'C':
                ll = (self.dimensions['X_Iwidth'],
                      self.dimensions['Y_Iwidth'])
                lr = (self.dimensions['X_Iwidth'] + self.dimensions['X_Czone'] - 1,
                      self.dimensions['Y_Iwidth'])
                ul = (self.dimensions['X_Iwidth'],
                      self.dimensions['Y_Iwidth'] + self.dimensions['Y_Czone'] - 1)
                ur = (self.dimensions['X_Iwidth'] + self.dimensions['X_Czone'] - 1,
                      self.dimensions['Y_Iwidth'] + self.dimensions['Y_Czone'] - 1)
        return {'ll':ll, 'lr':lr, 'ul':ul, 'ur':ur}

    def gimme_corners_ll(self, subzone=None, position=None):
        """
        Returns the lon/lat of the four corners of a rectangular grid,
        as a dict(corner=(lon, lat)) with corner in: \n
        ll = lower-left / lr = lower-right / ur = upper-right / ul = upper-left.

        :param subzone: for LAM grids, returns the corners of the subzone.
        :param position: position of corners with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        corners = self.gimme_corners_ij(subzone=subzone)
        for c in corners.keys():
            i = corners[c][0]
            j = corners[c][1]
            corners[c] = self.ij2ll(i, j, position)
        return corners

    def minmax_ll(self, subzone=None):
        """Return min/max of lon/lat."""
        (imin, jmin) = self.gimme_corners_ij(subzone)['ll']
        (imax, jmax) = self.gimme_corners_ij(subzone)['ur']
        border = [(imin, j) for j in range(jmin, jmax + 1)] + \
                 [(imax, j) for j in range(jmin, jmax + 1)] + \
                 [(i, jmin) for i in range(imin, imax + 1)] + \
                 [(i, jmax) for i in range(imin, imax + 1)]
        ilist, jlist = list(zip(*border))
        (lons, lats) = self.ij2ll(numpy.array(ilist),
                                  numpy.array(jlist))
        lonmin, lonmax = lons.min(), lons.max()
        latmin, latmax = lats.min(), lats.max()
        if numpy.ma.masked in (lonmin, lonmax, latmin, latmax):
            #space-view geometry for example
            lons, lats = self.get_lonlat_grid(subzone=subzone)
            lonmin, lonmax = lons.min(), lons.max()
            latmin, latmax = lats.min(), lats.max()
        return {'lonmin':lonmin, 'lonmax':lonmax,
                'latmin':latmin, 'latmax':latmax}

    def point_is_inside_domain_ll(self, lon, lat,
                                  margin=-0.1,
                                  subzone=None,
                                  position=None):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.

        :param lon: longitude of point(s) in degrees.
        :param lat: latitude of point(s) in degrees.
        :param margin: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        :param subzone: considers only a subzone among ('C', 'CI') of the domain.
        :param position: position of the grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """
        (Xmin, Ymin) = self.gimme_corners_ij(subzone)['ll']
        (Xmax, Ymax) = self.gimme_corners_ij(subzone)['ur']
        lon, lat = as_numpy_array(lon), as_numpy_array(lat)
        p = self.ll2ij(lon, lat, position)
        inside = numpy.logical_and(numpy.logical_and(Xmin + margin <= p[0],
                                                     p[0] <= Xmax - margin),
                                   numpy.logical_and(Ymin + margin <= p[1],
                                                     p[1] <= Ymax - margin))
        return inside.squeeze()

    def point_is_inside_domain_ij(self,
                                  i=None,
                                  j=None,
                                  margin=-0.1,
                                  subzone=None):
        """
        Returns True if the point(s) of i/j coordinates is(are) inside the
        field.

        :param i: X index of point
        :param j: Y index of point.
        :param margin: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        :param subzone: considers only a subzone among ('C', 'CI') of the domain.
        """
        if self.datashape['j'] and j is None:
            raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
        if self.datashape['i'] and i is None:
            raise epygramError("*i* is mandatory when field has one horizontal dimension")

        (Xmin, Ymin) = self.gimme_corners_ij(subzone)['ll']
        (Xmax, Ymax) = self.gimme_corners_ij(subzone)['ur']

        i = None if i is None else as_numpy_array(i)
        j = None if j is None else as_numpy_array(j)

        if i is not None:
            inside_i = numpy.logical_and(Xmin + margin <= i,
                                         i <= Xmax - margin)
        if j is not None:
            inside_j = numpy.logical_and(Ymin + margin <= j,
                                         j <= Ymax - margin)
        if i is j is None:
            inside = True
        elif i is None:
            inside = inside_j
        elif j is None:
            inside = inside_i
        else:
            inside = numpy.logical_and(inside_i, inside_j)
        return inside

    def nearest_points(self, lon, lat, request,
                       position=None,
                       external_distance=None,
                       squeeze=True):
        """
        Returns the (i, j) positions of the nearest points.

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
                Only valid with request={'n':'1'}
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

        In case of a {'n':'2*2'} request, order of points obtained on a rectangular grid is
        (and must be) such that first and second point share the i position (hence
        third and forth point also share the i position) and first and third share the j
        position (hence the second and forth also share the j position).
        """
        if not numpy.all(self.point_is_inside_domain_ll(lon, lat, position=position)):
            raise ValueError("point (" + str(lon) + ", " + str(lat) +
                             ") is out of field domain.")

        lon, lat = as_numpy_array(lon).flatten(), as_numpy_array(lat).flatten()

        (i0, j0) = self.ll2ij(lon, lat, position)
        i_int = numpy.floor(i0).astype('int')
        j_int = numpy.floor(j0).astype('int')
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

        # compute points position
        if request == {'n':'1'} and not external_distance:
            result = [(numpy.rint(i0).astype('int'),
                       numpy.rint(j0).astype('int'))]
            result = moveaxis(numpy.array(result), -1, 0)
        else:
            # square: size
            if external_distance:
                n = 1
            elif request.get('radius'):
                resolution = self.resolution_ll(lon, lat)
                n = max(numpy.around(float(request['radius']) / float(resolution)).astype('int') * 2,
                        1)
            elif nsquare_match:
                n = int(nsquare_match.group('n'))
            else:
                raise epygramError("unrecognized **request**: " + str(request))
            # square: indexes
            (ii, jj) = _increments(n)
            ii = [i_int + di for di in ii]
            jj = [j_int + dj for dj in jj]
            result = [(i, j) for i in ii for j in jj]
            result = moveaxis(numpy.array(result), -1, 0)

            # filter: if external distance
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

        # filter: if radius
        if request.get('radius'):
            result = list(result)  # We transform into list to be able to modify the length
            for ipt, points in enumerate(result):
                if request.get('shape', 'circle') == 'circle':
                    result[ipt] = [(i, j)
                                   for (i, j) in points
                                   if self.distance((lon[ipt], lat[ipt]), self.ij2ll(i, j)) <= request['radius']]
                elif request.get('shape') == 'square':
                    result[ipt] = [(i, j)
                                   for (i, j) in points
                                   if all(abs(numpy.array(self.ll2xy(lon[ipt], lat[ipt])) - numpy.array(self.ij2xy(i, j))) <= request['radius'])]
                assert len(result[ipt]) > 0, "no points found: radius may be too small."
                if not numpy.all(self.point_is_inside_domain_ij(*zip(*result[ipt]))):
                    raise epygramError("one point (" + str(lon) + ", " + str(lat) +
                                       ") too close to field domain borders.")
            if squeeze and len(result) == 1:
                result = result[0]
        else:
            # check all points in domain
            for point in moveaxis(numpy.array(result), 0, -1):
                if not numpy.all(self.point_is_inside_domain_ij(*point)):
                    raise epygramError("one point (" + str(lon) + ", " + str(lat) +
                                       ") too close to field domain borders.")
            if squeeze:
                result = result.squeeze()
        return result

    def _what_grid_dimensions(self, out=sys.stdout,
                              arpifs_var_names=False,
                              spectral_geometry=None):
        """
        Writes in file a summary of the grid & dimensions of the field.

        :param out: the output open file-like object
        :param arpifs_var_names: if True, prints the equivalent 'arpifs' variable
          names.
        :param spectral_geometry: an optional dict containing the spectral
          truncatures {'in_X':, 'in_Y':}.
        """
        varname = ''
        grid = self.grid
        dimensions = self.dimensions
        if 'LAMzone' in grid:
            write_formatted(out, "Zone", grid['LAMzone'])
            if arpifs_var_names:
                varname = ' (NDLON)'
            write_formatted(out, "Total points in X" + varname,
                            dimensions['X'])
            if arpifs_var_names:
                varname = ' (NDGL)'
            write_formatted(out, "Total points in Y" + varname,
                            dimensions['Y'])
            if grid['LAMzone'] == 'CIE':
                write_formatted(out, "Points of C+I in X",
                                dimensions['X_CIzone'])
                write_formatted(out, "Points of C+I in Y",
                                dimensions['Y_CIzone'])
                if arpifs_var_names:
                    varname = ' (NDLUN-1)'
                write_formatted(out,
                                "Low-left X offset for I zone" + varname,
                                dimensions['X_CIoffset'])
                if arpifs_var_names:
                    varname = ' (NDGUN-1)'
                write_formatted(out,
                                "Low-left Y offset for I zone" + varname,
                                dimensions['Y_CIoffset'])
                write_formatted(out, "Width of I strip in X",
                                dimensions['X_Iwidth'])
                write_formatted(out, "Width of I strip in Y",
                                dimensions['Y_Iwidth'])
            elif grid['LAMzone'] == 'CI':
                write_formatted(out, "Width of I strip in X",
                                dimensions['X_Iwidth'])
                write_formatted(out, "Width of I strip in Y",
                                dimensions['Y_Iwidth'])
            if spectral_geometry is not None:
                if arpifs_var_names:
                    varname = ' (NMSMAX)'
                write_formatted(out, "Truncation in X" + varname,
                                spectral_geometry['in_X'])
                if arpifs_var_names:
                    varname = ' (NSMAX)'
                write_formatted(out, "Truncation in Y" + varname,
                                spectral_geometry['in_Y'])
        else:
            if arpifs_var_names:
                varname = ' (NDLON)'
            write_formatted(out, "Total points in X" + varname,
                            dimensions['X'])
            write_formatted(out, "Total points in Y" + varname,
                            dimensions['Y'])
        out.write(separation_line)

    #FIXME: cleanme def __eq__(self, other):
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
        # known issue __eq__/must be defined both or none, else inheritance is broken
        return super(D3RectangularGridGeometry, self).__hash__()"""


class D3UnstructuredGeometry(D3RectangularGridGeometry):
    """Handles the geometry for an unstructured 3-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['unstructured'])),
            position_on_horizontal_grid=dict(
                default='center',
                values=set(['center'])),
        )
    )

    def __init__(self, *args, **kwargs):
        super(D3UnstructuredGeometry, self).__init__(*args, **kwargs)

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


class D3AcademicGeometry(D3RectangularGridGeometry):
    """Handles the geometry for an academic 3-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['academic'])),
            projection=dict(
                type=FPDict,
                info="Handles projection information."),
            geoid=dict(
                type=FPDict,
                optional=True,
                default=FPDict({}),
                info="Of no meaning in this geometry.")
        )
    )

    def __init__(self, *args, **kwargs):
        super(D3RectangularGridGeometry, self).__init__(*args, **kwargs)
        if self.grid['input_position'] != (0, 0):
            raise NotImplementedError("For now, only input_position = (0, 0) is allowed for academic geometries.")
        self._center_lon = (self.dimensions['X'] - 1) / 2.
        self._center_lat = (self.dimensions['Y'] - 1) / 2.

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
        Replaces the method defined in D3RectangularGridGeometry to deal with
        1D or 2D simulations.
        """
        offset = super(D3AcademicGeometry, self)._getoffset(position)
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
        Returns a academic V2DGeometry.

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
        vcoordinate = fpx.geometry(structure='V',
                                   typeoffirstfixedsurface=255,
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

        kwargs_geom = dict(structure='V2D',
                           name=self.name,
                           grid=FPDict(grid),
                           projection=projection,
                           dimensions=FPDict(dimensions),
                           position_on_horizontal_grid='center' if position is None else position,
                           vcoordinate=vcoordinate
                           )
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return fpx.geometry(**kwargs_geom)


class D3RegLLGeometry(D3RectangularGridGeometry):
    """
    Handles the geometry for a Regular Lon/Lat 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['regular_lonlat'])),
        )
    )

    @property
    def isglobal(self):
        """
        :return: True if geometry is global
        """
        return self.dimensions['X'] * self.grid['X_resolution'].get('degrees') >= 360. and \
               self.dimensions['Y'] * self.grid['Y_resolution'].get('degrees') >= 180.

    def __init__(self, *args, **kwargs):
        super(D3RegLLGeometry, self).__init__(*args, **kwargs)

        if self.grid['input_position'] == ((float(self.dimensions['X']) - 1) / 2.,
                                           (float(self.dimensions['Y']) - 1) / 2.):
            self._center_lon = self.grid['input_lon']
            self._center_lat = self.grid['input_lat']
        elif self.grid['input_position'] == (0, 0):
            self._center_lon = Angle(round(self.grid['input_lon'].get('degrees') +
                                           self.grid['X_resolution'].get('degrees') *
                                           (self.dimensions['X'] - 1) / 2, _rd),
                                     'degrees')
            self._center_lat = Angle(round(self.grid['input_lat'].get('degrees') +
                                           self.grid['Y_resolution'].get('degrees') *
                                           (self.dimensions['Y'] - 1) / 2, _rd),
                                     'degrees')
        elif self.grid['input_position'] == (0, self.dimensions['Y'] - 1):
            self._center_lon = Angle(round(self.grid['input_lon'].get('degrees') +
                                           self.grid['X_resolution'].get('degrees') *
                                           (self.dimensions['X'] - 1) / 2, _rd),
                                     'degrees')
            self._center_lat = Angle(round(self.grid['input_lat'].get('degrees') -
                                           self.grid['Y_resolution'].get('degrees') *
                                           (self.dimensions['Y'] - 1) / 2, _rd),
                                     'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1, 0):
            self._center_lon = Angle(round(self.grid['input_lon'].get('degrees') -
                                           self.grid['X_resolution'].get('degrees') *
                                           (self.dimensions['X'] - 1) / 2, _rd),
                                     'degrees')
            self._center_lat = Angle(round(self.grid['input_lat'].get('degrees') +
                                           self.grid['Y_resolution'].get('degrees') *
                                           (self.dimensions['Y'] - 1) / 2, _rd),
                                     'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1,
                                             self.dimensions['Y'] - 1):
            self._center_lon = Angle(round(self.grid['input_lon'].get('degrees') -
                                           self.grid['X_resolution'].get('degrees') *
                                           (self.dimensions['X'] - 1) / 2, _rd),
                                     'degrees')
            self._center_lat = Angle(round(self.grid['input_lat'].get('degrees') -
                                           self.grid['Y_resolution'].get('degrees') *
                                           (self.dimensions['Y'] - 1) / 2, _rd),
                                     'degrees')
        else:
            raise NotImplementedError("this 'input_position': " +
                                      str(self.grid['input_position']))

        # earth-round grids: wrap TODO:
        corners = self.gimme_corners_ll()
        if abs(abs(degrees_nearest_mod(corners['ul'][0], corners['ur'][0]) -
                   corners['ur'][0]) -
               self.grid['X_resolution'].get('degrees')) <= config.epsilon:
            self._earthround = True
        else:
            self._earthround = False

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


    def global_shift_center(self, longitude_shift):
        """
        Shifts the center of the geometry by *longitude_shift* (in degrees).
        *longitude_shift* has to be a multiple of the grid's resolution in
        longitude.
        """
        corners = self.gimme_corners_ll()
        zip_width = abs(degrees_nearest_mod(corners['ur'][0] - corners['ul'][0], 0.))
        zip_minus_resolution = round(zip_width - self.grid['X_resolution'].get('degrees'),
                                     _rd)
        if abs(zip_minus_resolution) < config.epsilon:
            as_int = 1e6  # decimal error
            if abs((longitude_shift * as_int) %
                   (self.grid['X_resolution'].get('degrees') * as_int)) > config.epsilon:
                raise epygramError(("*longitude_shift* ({}) has to be a multiple" +
                                    " of the grid's resolution in longitude ({}).").
                                   format(longitude_shift, self.grid['X_resolution'].get('degrees')))
            self._center_lon = Angle(self._center_lon.get('degrees') + longitude_shift,
                                     'degrees')
            self.grid['input_lon'] = Angle(self.grid['input_lon'].get('degrees') + longitude_shift,
                                           'degrees')
        else:
            raise epygramError("unable to shift center if " +
                               "lon_max - lon_min != X_resolution")

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
        return [self.xy2ll(*xy) for xy in xy_linspace]

    @_need_pyproj_geod
    def distance(self, end1, end2):
        """
        Computes the distance between two points along a straight line in the
        geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.

        Warning: requires the :mod:`pyproj` module.
        """
        plast = end1
        distance = 0
        for p in self.linspace(end1, end2, 1000)[1:]:
            distance += self._pyproj_geod.inv(plast[0], plast[1], *p)[2]
            plast = p
        return distance

    def make_zoom_geometry(self, zoom, extra_10th=False):
        """
        Returns a new geometry with the points contained in *zoom*.

        :param zoom: a dict(lonmin=, lonmax=, latmin=, latmax=).
        :param extra_10th: if True, add 1/10th of the X/Y extension of the zoom
                           (only usefull for the regular_lonlat grid implementation).
        """
        kwargs_zoomgeom = {'structure':self.structure,
                           'vcoordinate':self.vcoordinate.deepcopy(),
                           'position_on_horizontal_grid':self.position_on_horizontal_grid,
                           'geoid':self.geoid}
        if extra_10th:
            dx = (degrees_nearest_mod(zoom['lonmax'], 0.) -
                  degrees_nearest_mod(zoom['lonmin'], 0.)) / 10.
            dy = (zoom['latmax'] - zoom['latmin']) / 10.
            zoom = {'lonmin':zoom['lonmin'] - dx,
                    'lonmax':zoom['lonmax'] + dx,
                    'latmin':zoom['latmin'] - dy,
                    'latmax':zoom['latmax'] + dy}
        imin, jmin = self.ll2ij(zoom['lonmin'], zoom['latmin'])
        imax, jmax = self.ll2ij(zoom['lonmax'], zoom['latmax'])
        if imin > imax:
            gridmin = self.gimme_corners_ll()['ll'][0]
            diff_lonmin = (gridmin - degrees_nearest_mod(zoom['lonmin'],
                                                         gridmin))
            Xres = self.grid['X_resolution'].get('degrees')
            shift = (diff_lonmin // Xres + 1) * Xres
            shifted_self = self.deepcopy()
            shifted_self.global_shift_center(-shift)
            return shifted_self.make_zoom_geometry(zoom, extra_10th=False)  # zoom already includes the extra part
        elif imin == imax:  # means 360deg wide
            imin = 0
            imax = self.dimensions['X'] - 1
        imin = max(int(numpy.ceil(imin)),
                   0)
        imax = min(int(numpy.floor(imax)),
                   self.dimensions['X'] - 1)
        jmin = max(int(numpy.ceil(jmin)),
                   0)
        jmax = min(int(numpy.floor(jmax)),
                   self.dimensions['Y'] - 1)
        kwargs_zoomgeom['dimensions'] = {'X':imax - imin + 1,
                                         'Y':jmax - jmin + 1}
        lonmin, latmin = self.ij2ll(imin, jmin)
        kwargs_zoomgeom['name'] = self.name
        kwargs_zoomgeom['grid'] = {'input_position':(0, 0),
                                   'input_lon':Angle(lonmin, 'degrees'),
                                   'input_lat':Angle(latmin, 'degrees'),
                                   'X_resolution':self.grid['X_resolution'],
                                   'Y_resolution':self.grid['Y_resolution']}
        return fpx.geometry(**kwargs_zoomgeom)

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
        return (numpy.degrees(numpy.arctan2(x2 - x1, y2 - y1)) + 180.) % 360. - 180.

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


class D3RotLLGeometry(D3RegLLGeometry):
    """
    Handles the geometry for a Rotated Lon/Lat 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['rotated_lonlat']))
        )
    )

    def __init__(self, *args, **kwargs):
        super(D3RegLLGeometry, self).__init__(*args, **kwargs)
        if self.grid['input_position'] == ((float(self.dimensions['X']) - 1) / 2.,
                                           (float(self.dimensions['Y']) - 1) / 2.):
            self._center_rlon = self.grid['input_lon']
            self._center_rlat = self.grid['input_lat']
        elif self.grid['input_position'] == (0, 0):
            self._center_rlon = Angle(round(self.grid['input_lon'].get('degrees') +
                                            self.grid['X_resolution'].get('degrees') *
                                            (self.dimensions['X'] - 1) / 2, _rd),
                                      'degrees')
            self._center_rlat = Angle(round(self.grid['input_lat'].get('degrees') +
                                            self.grid['Y_resolution'].get('degrees') *
                                            (self.dimensions['Y'] - 1) / 2, _rd),
                                      'degrees')
        elif self.grid['input_position'] == (0, self.dimensions['Y'] - 1):
            self._center_rlon = Angle(round(self.grid['input_lon'].get('degrees') +
                                            self.grid['X_resolution'].get('degrees') *
                                            (self.dimensions['X'] - 1) / 2, _rd),
                                      'degrees')
            self._center_rlat = Angle(round(self.grid['input_lat'].get('degrees') -
                                            self.grid['Y_resolution'].get('degrees') *
                                            (self.dimensions['Y'] - 1) / 2, _rd),
                                      'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1, 0):
            self._center_rlon = Angle(round(self.grid['input_lon'].get('degrees') -
                                            self.grid['X_resolution'].get('degrees') *
                                            (self.dimensions['X'] - 1) / 2, _rd),
                                      'degrees')
            self._center_rlat = Angle(round(self.grid['input_lat'].get('degrees') +
                                            self.grid['Y_resolution'].get('degrees') *
                                            (self.dimensions['Y'] - 1) / 2, _rd),
                                      'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1,
                                             self.dimensions['Y'] - 1):
            self._center_rlon = Angle(round(self.grid['input_lon'].get('degrees') -
                                            self.grid['X_resolution'].get('degrees') *
                                            (self.dimensions['X'] - 1) / 2, _rd),
                                      'degrees')
            self._center_rlat = Angle(round(self.grid['input_lat'].get('degrees') -
                                            self.grid['Y_resolution'].get('degrees') *
                                            (self.dimensions['Y'] - 1) / 2, _rd),
                                      'degrees')
        else:
            raise NotImplementedError("this 'input_position': " +
                                      str(self.grid['input_position']))
        if self.grid.get('rotation', Angle(0., 'degrees')).get('degrees') != 0.:
            raise NotImplementedError("rotation != Angle(0.)")

        lon, lat = self.xy2ll(self._center_rlon.get('degrees'),
                              self._center_rlat.get('degrees'))
        self._center_lon = Angle(lon, 'degrees')
        self._center_lat = Angle(lat, 'degrees')

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


class D3ProjectedGeometry(D3RectangularGridGeometry):
    """
    Handles the geometry for a Projected 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic',
                            'space_view'])),
            projection=dict(
                type=FPDict,
                info="Handles projection information."),
        )
    )

    _ghost_attributes = D3RectangularGridGeometry._ghost_attributes + ['_proj']

    @property
    def secant_projection(self):
        """ Is the projection secant to the sphere ? (or tangent)"""
        return ('secant_lat' in self.projection or
                'secant_lat1' in self.projection)

    def __init__(self, *args, **kwargs):
        super(D3ProjectedGeometry, self).__init__(*args, **kwargs)
        import pyproj

        def compute_center_proj(p, center):
            if center == self.grid['input_position']:
                self._center_lon = self.grid['input_lon']
                self._center_lat = self.grid['input_lat']
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
                self._center_lon = Angle(center_lon, 'degrees')
                self._center_lat = Angle(center_lat, 'degrees')

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
                self._K = (math.log(m1) - math.log(m2)) / \
                          (math.log(t1) - math.log(t2))
            else:
                lat_1 = self.projection['reference_lat'].get('degrees')
                lat_2 = self.projection['reference_lat'].get('degrees')
                self._K = abs(self.projection['reference_lat'].get('cos_sin')[1])
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_1=lat_1, lat_2=lat_2,
                            **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = pyproj.Proj(proj=proj,
                                     lon_0=self.projection['reference_lon'].get('degrees'),
                                     lat_1=lat_1, lat_2=lat_2,
                                     x_0=-x0, y_0=-y0,
                                     **self.geoid)
        elif self.name == 'mercator':
            if self.secant_projection:
                lat_ts = self.projection['secant_lat'].get('degrees')
            else:
                lat_ts = 0.
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_ts=lat_ts,
                            **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = pyproj.Proj(proj=proj,
                                     lon_0=self.projection['reference_lon'].get('degrees'),
                                     lat_ts=lat_ts,
                                     x_0=-x0, y_0=-y0,
                                     **self.geoid)
        elif self.name == 'polar_stereographic':
            lat_0 = self.projection['reference_lat'].get('degrees')
            if self.secant_projection:
                lat_ts = self.projection['secant_lat'].get('degrees')
            else:
                lat_ts = self.projection['reference_lat'].get('degrees')
            p = pyproj.Proj(proj=proj,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            lat_0=lat_0, lat_ts=lat_ts,
                            **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = pyproj.Proj(proj=proj,
                                     lon_0=self.projection['reference_lon'].get('degrees'),
                                     lat_0=lat_0, lat_ts=lat_ts,
                                     x_0=-x0, y_0=-y0,
                                     **self.geoid)
        elif self.name == 'space_view':
            latSat = self.projection['satellite_lat'].get('degrees')
            lonSat = self.projection['satellite_lon'].get('degrees')
            height = self.projection['satellite_height']  # Height above ellipsoid
            if latSat != 0:
                raise epygramError("Only space views with satellite_lat=0 are allowed")
            p = pyproj.Proj(proj=proj,
                            h=height,
                            lon_0=lonSat,
                            **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = pyproj.Proj(proj=proj,
                                     h=height,
                                     lon_0=lonSat,
                                     x_0=-x0, y_0=-y0,
                                     **self.geoid)
        else:
            raise NotImplementedError("projection: " + self.name)

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
            geom_kwargs = copy.copy(self._attributes)
            geom_kwargs['grid'] = copy.copy(self.grid)
            geom_kwargs['dimensions'] = copy.copy(self.dimensions)
            if subzone == 'CI':
                io = self.dimensions.get('X_CIoffset', 0)
                jo = self.dimensions.get('Y_CIoffset', 0)
                centerPoint = (io + (float(self.dimensions['X_CIzone']) - 1) / 2.,
                               jo + (float(self.dimensions['Y_CIzone']) - 1) / 2.)  # Coordinates of center point
                geom_kwargs['grid']['LAMzone'] = subzone
                geom_kwargs['dimensions'] = {'X':self.dimensions['X_CIzone'],
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
                geom_kwargs['grid']['LAMzone'] = None
                geom_kwargs['dimensions'] = {'X':self.dimensions['X_Czone'],
                                             'Y':self.dimensions['Y_Czone']}
            geom_kwargs['grid']['input_lon'] = self._center_lon
            geom_kwargs['grid']['input_lat'] = self._center_lat
            geom_kwargs['grid']['input_position'] = centerPoint
            new_geom = fpx.geometry(**geom_kwargs)
        else:
            new_geom = self

        return new_geom

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
        kwargs_vcoord = {'structure': 'V',
                         'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [0]}
        vcoordinate = fpx.geometry(**kwargs_vcoord)
        geometry = fpx.geometrys.almost_clone(self,
                                              structure='H2D',
                                              vcoordinate=vcoordinate)
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
        Returns a projected V2DGeometry.

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
        vcoordinate = fpx.geometry(structure='V',
                                   typeoffirstfixedsurface=255,
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
        kwargs_geom = dict(structure='V2D',
                           name=self.name,
                           grid=FPDict(grid),
                           dimensions=FPDict(dimensions),
                           projection=FPDict(projection),
                           geoid=self.geoid,
                           position_on_horizontal_grid='center' if position is None else position,
                           vcoordinate=vcoordinate)
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return fpx.geometry(**kwargs_geom)

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
        theta = k * (lons - self.projection['reference_lon'].get('degrees')) - \
                self.projection.get('rotation', 0.).get('degrees')
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
        return super(D3ProjectedGeometry, self).__hash__()"""


class D3GaussGeometry(D3Geometry):
    """
    Handles the geometry for a Global Gauss grid 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['rotated_reduced_gauss', 'reduced_gauss', 'regular_gauss'])),
        )
    )

    _ghost_attributes = D3Geometry._ghost_attributes + ['_buffered_gauss_grid']

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
        kwargs_vcoord = {'structure': 'V',
                         'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [0]}
        vcoordinate = fpx.geometry(**kwargs_vcoord)
        geometry = fpx.geometrys.almost_clone(self,
                                              structure='H2D',
                                              vcoordinate=vcoordinate)
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

footprints.collectors.get(tag='geometrys').fasttrack = ('structure', 'name')
