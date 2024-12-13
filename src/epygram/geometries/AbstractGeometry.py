#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for 3D geometries of fields.
"""

import numpy
import math
import copy
import sys
import re

import footprints
from footprints import FPDict, proxy as fpx

from epygram import epygramError, config
from epygram.config import rounding_decimal as _rd
from epygram.util import (RecursiveObject, degrees_nearest_mod, Angle,
                          positive_longitudes, longitudes_between_minus180_180,
                          separation_line, write_formatted,
                          as_numpy_array, moveaxis,
                          is_scalar, Comparator)

from .VGeometry import VGeometry
from . import _need_pyproj_geod

epylog = footprints.loggers.getLogger(__name__)
_re_nearest_sq = re.compile(r'(?P<n>\d+)\*(?P<m>\d+)')


class Geometry(RecursiveObject):
    """
    Handles the geometry for a 3-Dimensions Field.
    Abstract mother class.
    """

    # ghost attributes are ignored when comparing 2 objects between them
    _ghost_attributes = RecursiveObject._ghost_attributes + ['_puredict', '_observer']  # footprints special attributes

    def __init__(self, grid, dimensions, vcoordinate,
                 position_on_horizontal_grid='__unknown__', geoid=None):
        """
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
        self.add_attr_dict('grid')
        self.add_attr_dict('dimensions')
        self.add_attr_class('vcoordinate', VGeometry)
        self.add_attr_inlist('position_on_horizontal_grid', ['upper-right', 'upper-left',
                                                             'lower-left', 'lower-right',
                                                             'center-left', 'center-right',
                                                             'lower-center', 'upper-center',
                                                             'center', '__unknown__'])
        self.add_attr_dict('geoid')

        self.grid = grid
        self.dimensions = dimensions
        self.vcoordinate = vcoordinate
        self.position_on_horizontal_grid = position_on_horizontal_grid
        self.geoid = config.default_geoid if geoid is None else geoid

        # Checks !
        self._consistency_check()

    def _consistency_check(self):
        # implemented in child classes
        pass

    def _GRIB2_sample(self, prefix):
        """Build GRIB2 sample name."""
        from epygram.extra.griberies.tables import typeoffixedsurface2sample as levels
        if self.structure == 'H2D':
            return '_'.join([prefix,
                             levels.get(self.vcoordinate.typeoffirstfixedsurface, 'sfc'),
                             'grib2'])
        else:
            raise NotImplementedError()

    @property
    def structure(self):
        """
        Returns the structure of the grid which depends on the prsent dimensions
        """
        has_k = len(self.vcoordinate.levels) > 1
        has_j = self.dimensions.get('Y', 1) != 1 or self.dimensions.get('lat_number', 1) != 1
        has_i = self.dimensions.get('X', 1) != 1 or self.dimensions.get('max_lon_number', 1) != 1

        return {(True, True, True): '3D',
                (True, True, False): 'H2D',
                (True, False, False): 'H1D',
                (False, False, False): 'Point',
                (True, False, True): 'V2D',
                (False, False, True): 'V1D',
               }[(has_i, has_j, has_k)]
        
    @property
    def rectangular_grid(self):
        """ Is the grid rectangular ? """
        return isinstance(self, RectangularGridGeometry)

    @property
    def projected_geometry(self):
        """ Is the geometry a projection ? """
        from .ProjectedGeometry import ProjectedGeometry #Imported here to prevent circular import
        return isinstance(self, ProjectedGeometry)

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
        geometry.vcoordinate.typeOfFirstFixedSurface = surface_type
        geometry.vcoordinate.levels = list(range(len(self.vcoordinate.levels))) if levels is None else levels
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
        """Returns a Geometry at coordinates *(lon,lat)* in degrees."""
        from .UnstructuredGeometry import UnstructuredGeometry #Imported here to prevent circular import
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[0])
        return UnstructuredGeometry(name='unstructured',
                                    vcoordinate=vcoordinate,
                                    dimensions={'X':1, 'Y':1},
                                    grid={'longitudes':[lon],
                                          'latitudes':[lat]},
                                    position_on_horizontal_grid='center'
                                    )

    def make_profile_geometry(self, lon, lat):
        """Returns a V1DGeometry at coordinates *(lon,lat)* in degrees."""
        from .UnstructuredGeometry import UnstructuredGeometry #Imported here to prevent circular import
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[])
        return UnstructuredGeometry(name='unstructured',
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
        Returns a Geometry.

        :param end1: must be a tuple (lon, lat) in degrees.
        :param end2: must be a tuple (lon, lat) in degrees.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param position: defines the position of data in the grid (for projected grids only)
        """
        from .UnstructuredGeometry import UnstructuredGeometry #Imported here to prevent circular import
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
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[])
        kwargs_geom = {'name':'unstructured',
                       'vcoordinate':vcoordinate,
                       'dimensions':{'X':len(transect), 'Y':1},
                       'grid':{'longitudes':[p[0] for p in transect],
                               'latitudes':[p[1] for p in transect]},
                       'position_on_horizontal_grid':'center' if position is None else position}
        if self.geoid:
            kwargs_geom['geoid'] = self.geoid
        return UnstructuredGeometry(**kwargs_geom)

    def make_physicallevels_geometry(self):
        """
        Returns a new geometry excluding levels with no physical meaning.

        :param getdata: if False returns a field without data
        """
        geometry = self.deepcopy()
        if self.vcoordinate.typeoffirstfixedsurface == 118:
            # build vertical geometry
            kwargs_vcoord = {'typeoffirstfixedsurface': self.vcoordinate.typeoffirstfixedsurface,
                             'position_on_grid': self.vcoordinate.position_on_grid,
                             'grid': copy.copy(self.vcoordinate.grid),
                             'levels': copy.copy(self.vcoordinate.levels)}
            # Suppression of levels above or under physical domain
            for level in list(kwargs_vcoord['levels']):
                if level < 1 or level > len(self.vcoordinate.grid['gridlevels']) - 1:
                    kwargs_vcoord['levels'].remove(level)
            # build geometry
            geometry.vcoordinate = VGeometry(**kwargs_vcoord)
        return geometry

    def make_zoom_geometry(self, zoom, extra_10th=False):
        """
        Returns an unstructured geometry with the points contained in *zoom*.

        :param zoom: a dict(lonmin=, lonmax=, latmin=, latmax=).
        :param extra_10th: if True, add 1/10th of the X/Y extension of the zoom
                           (only usefull for the regular_lonlat grid implementation).
        """
        from .UnstructuredGeometry import UnstructuredGeometry #Imported here to prevent circular import
        (lons, lats) = self.get_lonlat_grid()
        kwargs_zoomgeom = {'vcoordinate':self.vcoordinate.deepcopy(),
                           'position_on_horizontal_grid':self.position_on_horizontal_grid,
                           'geoid':self.geoid}
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
        kwargs_zoomgeom['name'] = 'unstructured'
        kwargs_zoomgeom['grid'] = {'longitudes':zoomlons,
                                   'latitudes':zoomlats}
        return UnstructuredGeometry(**kwargs_zoomgeom)

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
        assert isinstance(other, Geometry), "Other must be a geometry object"
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


class RectangularGridGeometry(Geometry):
    """
    Handles the geometry for a rectangular 3-Dimensions Field.
    Abstract.
    """

    @property
    def isglobal(self):
        """
        :return: True if geometry is global
        """
        return False  # Apart for global lon-lat

    def suggested_GRIB2_sample(self, spectral=False):
        if self.structure == 'H2D':
            if not spectral:
                # return self._GRIB2_sample(regular_ll')
                return 'GRIB2'
            else:
                return self._GRIB2_sample('sh')
        else:
            raise NotImplementedError()

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

        newgeom = self.deepcopy()
        if 'LAMzone' in newgeom.grid:
            newgeom.grid['LAMzone'] = None
        if 'input_position' in newgeom.grid:
            coords_00 = self.ij2ll(first_i, first_j, position='center')
            newgeom.grid['input_position'] = (0, 0)
            input_lon = Angle(coords_00[0], 'degrees') if isinstance(self.grid['input_lon'], Angle) else coords_00[0]
            input_lat = Angle(coords_00[1], 'degrees') if isinstance(self.grid['input_lat'], Angle) else coords_00[1]
            newgeom.grid['input_lon'] = input_lon
            newgeom.grid['input_lat'] = input_lat
        newgeom.dimensions = {'X':last_i - first_i + 1, 'Y':last_j - first_j + 1}
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
        newgeom = self.deepcopy()
        if 'LAMzone' in self.grid:
            newgeom.grid['LAMzone'] = None
        if 'input_position' in newgeom.grid:
            coords_00 = self.ij2ll(0, 0, position='center')
            newgeom.grid['input_position'] = (0, 0)
            input_lon = Angle(coords_00[0], 'degrees') if isinstance(self.grid['input_lon'], Angle) else coords_00[0]
            input_lat = Angle(coords_00[1], 'degrees') if isinstance(self.grid['input_lat'], Angle) else coords_00[1]
            newgeom.grid['input_lon'] = input_lon
            newgeom.grid['input_lat'] = input_lat
        newgeom.dimensions = {'X':int(math.ceil(newgeom.dimensions['X'] / float(sample_x))),
                              'Y':int(math.ceil(newgeom.dimensions['Y'] / float(sample_y)))}
        newgeom.grid['X_resolution'] = newgeom.grid['X_resolution'] * sample_x
        newgeom.grid['Y_resolution'] = newgeom.grid['Y_resolution'] * sample_y
        newgeom.vcoordinate.levels = newgeom.vcoordinate.levels[::sample_z]
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
        return super(RectangularGridGeometry, self).__hash__()"""


class LLGeometry(RectangularGridGeometry):
    """
    Handles the geometry for a Regular Lon/Lat 3-Dimensions Field.
    Abstract.

    TODO: Is this class really necessary. Maybe we could make only one class
          with LLGeometry, RegLLGeometry and RotLLGEometry?
          Is a RegLLGeometry equivalent to a RotLLGeometry with a null rotation angle?
    """

    @property
    def isglobal(self):
        """
        :return: True if geometry is global
        """
        return self.dimensions['X'] * self.grid['X_resolution'].get('degrees') >= 360. and \
               self.dimensions['Y'] * self.grid['Y_resolution'].get('degrees') >= 180.

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
        kwargs_zoomgeom = {'vcoordinate':self.vcoordinate.deepcopy(),
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
        return RegLLGeometry(**kwargs_zoomgeom) #FIXME why necessarily a regular and not a rotated?

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
