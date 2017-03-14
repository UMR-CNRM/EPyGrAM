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
from footprints import FootprintBase, FPDict, proxy as fpx

from epygram import epygramError, config
from epygram.util import RecursiveObject, degrees_nearest_mod, Angle, \
                         separation_line, write_formatted, stretch_array, \
                         nearlyEqual, set_figax, set_map_up
from .VGeometry import VGeometry

epylog = footprints.loggers.getLogger(__name__)
_re_nearest_sq = re.compile('(?P<n>\d+)\*(?P<m>\d+)')


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
                            'regular_lonlat',
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
                type=str,
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
                default=FPDict({}),
                info="To specify geoid shape; actually used in projected" +
                     " geometries only.")
        )
    )

    @property
    def rectangular_grid(self):
        """ Is the grid rectangular ? """
        return isinstance(self, D3RectangularGridGeometry)

    @property
    def projected_geometry(self):
        """ Is the geometry a projection ? """
        return 'projection' in self._attributes

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(D3Geometry, self).__init__(*args, **kwargs)

        # Checks !
        self._consistency_check()

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
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
                d4=True requires nb_validities > 0
        - *nb_validities* is the number of validities represented in data values
        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, extracts
          the levels resp. from the C or C+I zone off the C+I(+E) zone.

        levels are internally stored with the vertical dimension first whereas this
        method puts the time in first dimension.
        """

        if d4:
            if nb_validities < 1:
                raise ValueError("nb_validities must be >=1 when d4==True")

        levels = numpy.array(self.vcoordinate.levels)

        # We add the horizontal axis
        h_shape2D = self.get_datashape(force_dimZ=1, d4=True, subzone=subzone)[-2:]
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
            h_shape = self.get_datashape(force_dimZ=1)
            if len(levels.shape) == 1 + len(h_shape):
                original_has_time = False
            elif len(levels.shape) == 2 + len(h_shape):
                original_has_time = True
            else:
                raise epygramError("Wrong number of dimensions")
            if levels.shape[-len(h_shape):] != h_shape:
                raise epygramError("Shape of self.vcoordinate.levels does not agree with horizontal dimensions")
            if subzone is not None:
                levels = self.extract_subzone(levels, subzone)
            if d4 and ((not self.datashape['i']) or (not self.datashape['j'])):
                shape = levels.shape[:-len(h_shape)]  # shape without the horizontal dimensions
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

    def _getoffset(self, position=None):
        """Returns the offset to use for this position."""

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

    def distance(self, end1, end2):
        """
        Computes the distance between two points along a Great Circle.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).

        Warning: requires the :mod:`pyproj` module.
        """
        import pyproj

        g = pyproj.Geod(ellps='sphere')
        distance = g.inv(end1[0], end1[1], end2[0], end2[1])[2]

        return distance

    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up along a Great Circle.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        *num* is the number of points, including point1 and point2.

        Warning: requires the :mod:`pyproj` module.
        """
        import pyproj

        if num < 2:
            raise epygramError("'num' must be at least 2.")
        g = pyproj.Geod(ellps='sphere')
        transect = g.npts(end1[0], end1[1], end2[0], end2[1], num - 2)
        transect.insert(0, end1)
        transect.append(end2)

        return transect

    def azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points following a Great Circle.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).

        Warning: requires the :mod:`pyproj` module.
        """
        import pyproj

        g = pyproj.Geod(ellps='sphere')
        return g.inv(end1[0], end1[1], end2[0], end2[1])[0]

    def make_point_geometry(self, lon, lat):
        """
        Returns a PointGeometry.
        """
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
        """
        Returns a V1DGeometry.
        """
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

        Args: \n
        - *end1* must be a tuple (lon, lat).
        - *end2* must be a tuple (lon, lat).
        - *points_number* defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        - *resolution* defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        - *position* defines the position of data in the grid (for projected grids only)
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
        return fpx.geometry(structure='V2D',
                            name='unstructured',
                            vcoordinate=vcoordinate,
                            dimensions={'X':len(transect), 'Y':1},
                            grid={'longitudes':[p[0] for p in transect],
                                  'latitudes':[p[1] for p in transect]},
                            position_on_horizontal_grid='center' if position is None else position)

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

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def plotgeometry(self,
                     color='blue',
                     borderonly=True,
                     **kwargs):
        """
        Makes a simple plot of the geometry, with a number of options.

        Requires :mod:`matplotlib`

        Options: \n
        - *color*: color of the plotting. \n
        - *borderonly*: if True, only plot the border of the grid, else the
          whole grid. Ignored for global geometries.

        For other options, cf. :class:`epygram.fields.H2DField`.

        This method uses (hence requires) 'matplotlib' and 'basemap' libraries.
        """
        import matplotlib.pyplot as plt
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)
        fig, ax = set_figax(*kwargs.get('over', (None, None)))
        if self.name == 'academic':
            raise epygramError("We cannot plot lon/lat of an academic grid.")
        if kwargs.get('use_basemap') is None:
            bm_args = {k:kwargs[k]
                       for k in ('gisquality', 'subzone', 'specificproj', 'zoom')
                       if k in kwargs}
            bm = self.make_basemap(**bm_args)
        else:
            bm = kwargs.get('use_basemap')
        map_args = {k:kwargs[k]
                    for k in ('drawrivers', 'drawcoastlines', 'drawcountries',
                              'meridians', 'parallels',
                              'departments', 'boundariescolor',
                              'bluemarble', 'background')
                    if k in kwargs}
        set_map_up(bm, ax, **map_args)
        (lons, lats) = self.get_lonlat_grid(subzone=kwargs.get('subzone'))
        if borderonly and 'gauss' not in self.name:
            lons = numpy.array(list(lons[0, :]) + list(lons[-1, :]) +
                               list(lons[1:-1, 0]) + list(lons[1:-1, -1]))
            lats = numpy.array(list(lats[0, :]) + list(lats[-1, :]) +
                               list(lats[1:-1, 0]) + list(lats[1:-1, -1]))
        x, y = bm(lons, lats)
        xf = x.flatten()
        yf = y.flatten()
        bm.scatter(xf, yf,
                   s=kwargs.get('pointsize', 20),
                   marker=',',
                   color=color,
                   linewidths=0,
                   ax=ax)
        if kwargs.get('title') is None:
            ax.set_title(str(self.name))
        else:
            ax.set_title(kwargs.get('title'))

        return fig, ax

    def what(self, out=sys.stdout,
             vertical_geometry=True,
             arpifs_var_names=False,
             spectral_geometry=None):
        """
        Writes in file a summary of the geometry.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *vertical_geometry*: if True, writes the vertical geometry of the
          field.
        - *arpifs_var_names*: if True, prints the equivalent 'arpifs' variable
          names.
        - *spectral_geometry*: an optional dict containing the spectral
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
                            'regular_lonlat', 'academic', 'unstructured']))
        )
    )

    def _get_grid(self, indextype, subzone=None, position=None):
        """
        Returns a tuple of two tables containing the two indexes of each
        point, with 2D shape.

        - *indextype*: either 'ij', 'xy' or 'll' to get
          i,j indexes, x,y coordinates or lon,lat coordinates
        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        - *position*: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        Jmax = self.dimensions['Y']
        Imax = self.dimensions['X']
        igrid = numpy.zeros(Jmax * Imax)
        jgrid = numpy.zeros(Jmax * Imax)
        for j in range(Jmax):
            for i in range(Imax):
                igrid[j * Imax + i] = i
                jgrid[j * Imax + i] = j
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
        """Returns the number of gridpoints of the grid."""
        shp = self.get_datashape(dimT=1, force_dimZ=1, subzone=subzone)
        return shp[0] * shp[1]

    def get_lonlat_grid(self, subzone=None, position=None, d4=False, nb_validities=0):
        """
        Returns a tuple of two tables containing one the longitude of each
        point, the other the latitude, with 2D shape.

        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        - *position*: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
                d4=True requires nb_validities > 0
        - *nb_validities* is the number of validities represented in data values

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

        return (lons, lats)

    def extract_subzone(self, data, subzone):
        """
        Extracts the subzone C or CI from a LAM field.

        Args: \n
        - *data*: the data values with shape concording with geometry.
        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, extracts
          the data resp. from the C or C+I zone off the C+I(+E) zone.
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

        geom_kwargs = copy.deepcopy(self._attributes)
        geom_kwargs.pop('dimensions')
        if 'LAMzone' in geom_kwargs['grid']:
            geom_kwargs['grid']['LAMzone'] = None
        if 'input_position' in geom_kwargs['grid']:
            coords_00 = self.ij2ll(first_i, first_j)
            geom_kwargs['grid']['input_position'] = (0, 0)
            geom_kwargs['grid']['input_lon'] = Angle(coords_00[0], 'degrees')
            geom_kwargs['grid']['input_lat'] = Angle(coords_00[1], 'degrees')
        geom_kwargs['dimensions'] = {'X':last_i - first_i, 'Y':last_j - first_j}
        newgeom = fpx.geometry(**geom_kwargs)  # create new geometry object

        return newgeom

    def get_datashape(self, dimT=1, force_dimZ=None,
                      d4=False,
                      subzone=None):
        """
        Returns the data shape according to the geometry.
        - *force_dimZ*: if supplied, force the Z dimension instead of that
          of the vertical geometry
        - *dimT* if supplied, is the time dimension to be added to the
          data shape
        - *d4*: if True,  shape is 4D
                if False, shape has only those > 1
        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, informes that
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

    def reshape_data(self, data, first_dimension=None, d4=False,
                     subzone=None):
        """
        Returns a 2D data (horizontal dimensions) reshaped from 1D,
        according to geometry.

        - *data*: the 1D data (or 3D with a T and Z dimensions,
          or 2D with either a T/Z dimension, to be specified),
          of dimension concording with geometry. In case data is 3D, T must be
          first dimension and Z the second.
        - *first_dimension*: in case data is 2D, specify what is the first
          dimension of data among ('T', 'Z')
        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, informes that
          data is resp. on the C or C+I zone off the C+I(+E) zone.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
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

    def gimme_corners_ij(self, subzone=None):
        """
        Returns the indices (i, j) of the four corners of a rectangular grid,
        as a dict(corner=(i, j)) with corner in: \n
        ll = lower-left / lr = lower-right / ur = upper-right / ul = upper-left.

        (0, 0) is the lower-left corner of the grid.

        - *subzone*: for LAM fields, returns the corners of the subzone.
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

        - *subzone*: for LAM grids, returns the corners of the subzone.
        - *position*: position of corners with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        corners = self.gimme_corners_ij(subzone=subzone)
        for c in corners.keys():
            i = corners[c][0]
            j = corners[c][1]
            corners[c] = self.ij2ll(i, j, position)

        return corners

    def point_is_inside_domain_ll(self, lon, lat,
                                  margin=-0.1,
                                  subzone=None,
                                  position=None):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.

        Args: \n
        - *lon*: longitude of point(s) in degrees.
        - *lat*: latitude of point(s) in degrees.
        - *margin*: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        - *subzone*: considers only a subzone among ('C', 'CI') of the domain.
        - *position*: position of the grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        (Xmin, Ymin) = self.gimme_corners_ij(subzone)['ll']
        (Xmax, Ymax) = self.gimme_corners_ij(subzone)['ur']
        try:
            N = len(lon)
        except Exception:
            N = 1
        if N == 1:
            inside = Xmin + margin <= self.ll2ij(lon, lat, position)[0] <= Xmax - margin and \
                     Ymin + margin <= self.ll2ij(lon, lat, position)[1] <= Ymax - margin
        else:
            inside = []
            for i in range(N):
                p = self.ll2ij(lon[i], lat[i], position)
                inside.append(Xmin + margin <= p[0] <= Xmax - margin and
                              Ymin + margin <= p[1] <= Ymax - margin)

        return inside

    def point_is_inside_domain_ij(self,
                                  i=None,
                                  j=None,
                                  margin=-0.1,
                                  subzone=None):
        """
        Returns True if the point(s) of i/j coordinates is(are) inside the
        field.

        Args: \n
        - *i* and *j*: indexes of point.
        - *margin*: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        - *subzone*: considers only a subzone among ('C', 'CI') of the domain.
        """

        if self.datashape['j'] and j is None:
            raise epygramError("*j* is mandatory when field has a two horizontal dimensions")
        if self.datashape['i'] and j is None:
            raise epygramError("*i* is mandatory when field has one horizontal dimension")

        (Xmin, Ymin) = self.gimme_corners_ij(subzone)['ll']
        (Xmax, Ymax) = self.gimme_corners_ij(subzone)['ur']
        try:
            N = len(i)
        except Exception:
            N = 1
        if N == 1:
            inside_i = (i is None) or (Xmin + margin <= i <= Xmax - margin)
            inside_j = (j is None) or (Ymin + margin <= j <= Ymax - margin)
            inside = inside_i and inside_j
        else:
            inside = []
            for n in range(N):
                inside_i = (i is None) or (Xmin + margin <= i[n] <= Xmax - margin)
                inside_j = (j is None) or (Ymin + margin <= j[n] <= Ymax - margin)
                inside.append(inside_i and inside_j)

        return inside

    def nearest_points(self, lon, lat, request,
                       position=None,
                       external_distance=None):
        """
        Returns the (i, j) positions of the nearest points.

        :rtype list of (lon, lat) tuples, or (lon, lat) if one only point.

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
        """

        if not self.point_is_inside_domain_ll(lon, lat, position=position):
            raise ValueError("point (" + str(lon) + ", " + str(lat) +
                             ") is out of field domain.")

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
            points = [(numpy.rint(i0).astype('int'),
                       numpy.rint(j0).astype('int'))]
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
            points = [(i, j) for i in ii for j in jj]
            # filter: if external distance
            if external_distance:
                mindistance = None
                for p in points:
                    dist = abs(external_distance['external_field'].getvalue_ij(*p, one=True) - external_distance['target_value'])
                    if mindistance is None or dist < mindistance:
                        points = [p]
                        mindistance = dist
            # filter: if radius
            if request.get('radius'):
                if request.get('shape', 'circle') == 'circle':
                    points = [(i, j)
                              for (i, j) in points
                              if self.distance((lon, lat), self.ij2ll(i, j)) <= request['radius']]
                elif request.get('shape') == 'square':
                    points = [(i, j)
                              for (i, j) in points
                              if all(abs(numpy.array(self.ll2xy(lon, lat)) - numpy.array(self.ij2xy(i, j))) <= request['radius'])]
                assert len(points) > 0, "no points found: radius may be too small."

        # check all points in domain
        for point in points:
            if not self.point_is_inside_domain_ij(*point):
                raise epygramError("point (" + str(lon) + ", " + str(lat) +
                                   ") too close to field domain borders.")

        return points[0] if request.get('n') == '1' else points

    def _what_grid_dimensions(self, out,
                              arpifs_var_names=False,
                              spectral_geometry=None):
        """
        Writes in file a summary of the grid & dimensions of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *arpifs_var_names*: if True, prints the equivalent 'arpifs' variable
          names.
        - *spectral_geometry*: an optional dict containing the spectral
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


class D3UnstructuredGeometry(D3RectangularGridGeometry):
    """Handles the geometry for an unstructured 3-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['unstructured']))
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(D3UnstructuredGeometry, self).__init__(*args, **kwargs)

    def _consistency_check(self):
        """Check that the geometry is consistent."""
        grid_keys = ['latitudes', 'longitudes']
        if set(self.grid.keys()) != set(grid_keys) and \
           set(self.grid.keys()) != set(['DDH_domain']):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys) + " or " +
                               str(['DDH_domain']))

    def get_lonlat_grid(self, subzone=None, position=None, d4=False, nb_validities=0):
        """
        Returns a tuple of two tables containing one the longitude of each
        point, the other the latitude, with 2D shape.

        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        - *position*: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
                d4=True requires nb_validities > 0
        - *nb_validities* is the number of validities represented in data values


        Shape of 2D data on Rectangular grids: \n
          - grid[0,0] is SW, grid[-1,-1] is NE \n
          - grid[0,-1] is SE, grid[-1,0] is NW
        """

        if self._getoffset(position) != (0., 0.):
            raise epygramError('We can only retrieve latitude and longitude of mass point on an unstructured grid')
        lons = numpy.array(self.grid['longitudes'])
        lats = numpy.array(self.grid['latitudes'])
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

        return (lons, lats)

    def ij2ll(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.
        - *position*: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        (lons, lats) = self.get_lonlat_grid(position=position)
        return (lons[j, i], lats[j, i])

    def ll2ij(self, lon, lat, position=None):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints.
        - *position*: lat lon position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float.
        """

        (lons, lats) = self.get_lonlat_grid(position=position)
        result = None
        for j in range(lons.shape[0]):
            for i in range(lons.shape[1]):
                if (lons[j, i], lats[j, i]) == (lon, lat):
                    if result is None:
                        result = (i, j)
                    else:
                        raise epygramError("Several points have the same coordinates.")
        if result is None:
            raise epygramError("No point found with these coordinates.")
        return result

    def nearest_points(self, lon, lat, interpolation,
                       position=None,
                       external_distance=None):
        """
        Returns the (i, j) position of the points needed to perform an
        interpolation. This is a list of (lon, lat) tuples.

        Args: \n
        - *lon*: longitude of point in degrees.
        - *lat*: latitude of point in degrees.
        - *interpolation* can be:
          - 'nearest' to get the nearest point only
          - 'linear' to get the 2*2 points bordering the (*lon*, *lat*) position
          - 'cubic' to get the 4*4 points bordering the (*lon*, *lat*) position
        - *position*: position in the model cell of the lat lon position.
          Defaults to self.position_on_horizontal_grid.
        - *external_distance* can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':a_3DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = |target_value - external_field.data|
        """

        if interpolation != 'nearest':
            raise NotImplementedError("*interpolation* != 'nearest' for UnstructuredGeometry.")
        if external_distance is not None:
            raise NotImplementedError("*external_distance* is not None for UnstructuredGeometry.")

        (lons, lats) = self.get_lonlat_grid(position=position)

        if lons.ndim == 2:
            dist0 = (lon - lons[0, 0]) ** 2 + (lat - lats[0, 0]) ** 2
        else:
            dist0 = (lon - lons[0]) ** 2 + (lat - lats[0]) ** 2
        i0 = 0
        j0 = 0
        if lons.ndim == 2:
            for j in range(0, lons.shape[0]):
                for i in range(0, lons.shape[1]):
                    dist = (lon - lons[j, i]) ** 2 + (lat - lats[j, i]) ** 2
                    if dist < dist0:
                        dist0 = dist
                        i0 = i
                        j0 = j
        else:
            for k in range(lons.shape[0]):
                dist = (lon - lons[k]) ** 2 + (lat - lats[k]) ** 2
                if dist < dist0:
                    dist0 = dist
                    i0 = k
                    j0 = 0

        return (i0, j0)

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.
        *point* must be a tuple (lon, lat).
        """
        return self.resolution_ij(*self.nearest_points(lon, lat, 'nearest'))

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        """

        (lons, lats) = self.get_lonlat_grid()
        result = None
        for jp in range(lons.shape[0]):
            for ip in range(lons.shape[1]):
                if (ip, jp) != (i, j):
                    dist = self.distance((lons[j, i], lats[j, i]), (lons[jp, ip], lats[jp, ip]))
                    if result is None or dist < result:
                        result = dist
        if result is None:
            raise epygramError("Resolution cannot be computed on one point grid.")
        return result

    def _what_grid(self, out):
        """
        Writes in file a summary of the grid of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).

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

    def make_basemap(self, gisquality='i', specificproj=None,
                     zoom=None, ax=None, **kwargs):
        """
        Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
        projection (if available). This is designed to avoid explicit handling
        of deep horizontal geometry attributes.

        Args: \n
        - *gisquality*: defines the quality of GIS contours, cf. Basemap doc. \n
          Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
        - *specificproj*: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        - *zoom*: specifies the lon/lat borders of the map, implying hereby a
          'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=). \n
          Overwrites *specificproj*.
        - *ax*: a matplotlib ax on which to plot; if None, plots will be done
          on matplotlib.pyplot.gca()
        """
        from mpl_toolkits.basemap import Basemap

        # corners
        if self.dimensions['Y'] == 1:
            (lon, lat) = self.get_lonlat_grid()
            llcrnrlon = numpy.amin(lon)
            urcrnrlon = numpy.amax(lon)
            llcrnrlat = numpy.amin(lat)
            urcrnrlat = numpy.amax(lat)
        else:
            (llcrnrlon, llcrnrlat) = self.ij2ll(*self.gimme_corners_ij()['ll'])
            (urcrnrlon, urcrnrlat) = self.ij2ll(*self.gimme_corners_ij()['ur'])

        # make basemap
        if zoom is not None:
            # zoom case
            specificproj = 'cyl'
            llcrnrlon = zoom['lonmin']
            llcrnrlat = zoom['latmin']
            urcrnrlon = zoom['lonmax']
            urcrnrlat = zoom['latmax']
        if specificproj is None:
            # defaults
            if llcrnrlat <= -89.0 or \
               urcrnrlat >= 89.0:
                proj = 'cyl'
            else:
                proj = 'merc'
            (lons, lats) = self.get_lonlat_grid()
            if lons.ndim == 1:
                lonmax = lons[:].max()
                lonmin = lons[:].min()
            else:
                lonmax = lons[:, -1].max()
                lonmin = lons[:, 0].min()
            if lats.ndim == 1:
                latmax = lats[:].max()
                latmin = lats[:].min()
            else:
                latmax = lats[-1, :].max()
                latmin = lats[0, :].min()
            b = Basemap(resolution=gisquality, projection=proj,
                        llcrnrlon=lonmin,
                        llcrnrlat=latmin,
                        urcrnrlon=lonmax,
                        urcrnrlat=latmax,
                        ax=ax)
        else:
            # specificproj
            if hasattr(self, '_center_lon') and hasattr(self, '_center_lat'):
                lon0 = self._center_lon.get('degrees')
                lat0 = self._center_lat.get('degrees')
            else:
                lon0 = lat0 = None
            if specificproj == 'kav7':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'ortho':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            lat_0=lat0,
                            ax=ax)
            elif specificproj == 'cyl':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'moll':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif isinstance(specificproj, tuple) and \
                 specificproj[0] == 'nsper' and \
                 isinstance(specificproj[1], dict):
                sat_height = specificproj[1].get('sat_height', 3000) * 1000.
                b = Basemap(resolution=gisquality,
                            projection=specificproj[0],
                            lon_0=specificproj[1].get('lon', lon0),
                            lat_0=specificproj[1].get('lat', lat0),
                            satellite_height=sat_height,
                            ax=ax)

        return b


class D3AcademicGeometry(D3RectangularGridGeometry):
    """Handles the geometry for an academic 3-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['academic']))
        )
    )

    def _consistency_check(self):
        """Check that the geometry is consistent."""

        grid_keys = ['LAMzone', 'X_resolution', 'Y_resolution']
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
            offset = (offset[1], 0)

        return offset

    def ij2xy(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of gridpoints, \n
        (*x, y*) being the coordinates in the projection.
        - *position*: position to return with respect to model cell. Defaults
          to self.position_on_horizontal_grid.

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
        (*x, y*) being the coordinates in the projection, \n
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        - *position*: position represented by (x,y) with respect to model cell.
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
        (*i, j*) being the coordinates in the 2D matrix of gridpoints, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.
        - *position*: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        !!! This routine has no sense for this geometry, it is identity.
        """

        if isinstance(i, list) or isinstance(i, tuple):
            i = numpy.array(i)
        if isinstance(j, list) or isinstance(j, tuple):
            j = numpy.array(j)
        (oi, oj) = self._getoffset(position)

        return (i + 1 + oi, j + 1 + oj)

    def ll2ij(self, lon, lat, position=None):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        - *position*: lat lon position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        !!! This routine has no sense for this geometry, it is identity.
        """

        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)
        (oi, oj) = self._getoffset(position)

        return (lon - 1 - oi, lat - 1 - oj)

    def ll2xy(self, lon, lat):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*x, y*) being the coordinates in the projection.

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
        (*x, y*) being the coordinates in the projection, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """

        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)

        return self.ij2ll(*self.xy2ij(x, y))

    def distance(self, end1, end2):
        """
        Computes the distance between two points along a straight line in the
        geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        """

        return numpy.sqrt(((end1[0] - end2[0]) * self.grid['X_resolution']) ** 2 +
                          ((end1[1] - end2[1]) * self.grid['Y_resolution']) ** 2)

    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up in the geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        *num* is the number of points, including point1 and point2.
        """

        if num < 2:
            raise epygramError("'num' must be at least 2.")
        return list(zip(numpy.linspace(end1[0], end2[0], num=num),
                        numpy.linspace(end1[1], end2[1], num=num)))

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.
        *point* must be a tuple (lon, lat).
        """
        return min(self.grid['X_resolution'], self.grid['Y_resolution'])

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        """
        return min(self.grid['X_resolution'], self.grid['Y_resolution'])

    def azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points in geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        """

        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)

        return (numpy.degrees(numpy.arctan2(x2 - x1, y2 - y1)) + 180.) % 360. - 180.

    def _what_position(self, out):
        """
        Writes in file a summary of the position of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        """
        write_formatted(out, "Kind of Geometry", 'Academic')
        if 'latitude' in self.grid:
            write_formatted(out, "Latitude", self.grid['latitude'])
        if 'longitude' in self.grid:
            write_formatted(out, "Longitude", self.grid['longitude'])


class D3RegLLGeometry(D3RectangularGridGeometry):
    """
    Handles the geometry for a Regular Lon/Lat 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['regular_lonlat']))
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(D3RegLLGeometry, self).__init__(*args, **kwargs)

        if self.grid['input_position'] == ((float(self.dimensions['X']) - 1) / 2.,
                                           (float(self.dimensions['Y']) - 1) / 2.):
            self._center_lon = self.grid['input_lon']
            self._center_lat = self.grid['input_lat']
        elif self.grid['input_position'] == (0, 0):
            self._center_lon = Angle(self.grid['input_lon'].get('degrees') +
                                     self.grid['X_resolution'].get('degrees') *
                                     (self.dimensions['X'] - 1) / 2,
                                     'degrees')
            self._center_lat = Angle(self.grid['input_lat'].get('degrees') +
                                     self.grid['Y_resolution'].get('degrees') *
                                     (self.dimensions['Y'] - 1) / 2,
                                     'degrees')
        elif self.grid['input_position'] == (0, self.dimensions['Y'] - 1):
            self._center_lon = Angle(self.grid['input_lon'].get('degrees') +
                                     self.grid['X_resolution'].get('degrees') *
                                     (self.dimensions['X'] - 1) / 2,
                                     'degrees')
            self._center_lat = Angle(self.grid['input_lat'].get('degrees') -
                                     self.grid['Y_resolution'].get('degrees') *
                                     (self.dimensions['Y'] - 1) / 2,
                                     'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1, 0):
            self._center_lon = Angle(self.grid['input_lon'].get('degrees') -
                                     self.grid['X_resolution'].get('degrees') *
                                     (self.dimensions['X'] - 1) / 2,
                                     'degrees')
            self._center_lat = Angle(self.grid['input_lat'].get('degrees') +
                                     self.grid['Y_resolution'].get('degrees') *
                                     (self.dimensions['Y'] - 1) / 2,
                                     'degrees')
        elif self.grid['input_position'] == (self.dimensions['X'] - 1,
                                             self.dimensions['Y'] - 1):
            self._center_lon = Angle(self.grid['input_lon'].get('degrees') -
                                     self.grid['X_resolution'].get('degrees') *
                                     (self.dimensions['X'] - 1) / 2,
                                     'degrees')
            self._center_lat = Angle(self.grid['input_lat'].get('degrees') -
                                     self.grid['Y_resolution'].get('degrees') *
                                     (self.dimensions['Y'] - 1) / 2,
                                     'degrees')
        else:
            raise NotImplementedError("this 'input_position': " +
                                      str(self.grid['input_position']))

    def getcenter(self):
        """
        Returns the coordinate of the grid center as a tuple (center_lon, center_lat).
        """

        return (self._center_lon, self._center_lat)

    def _consistency_check(self):
        """
        Check that the geometry is consistent.
        """

        if self.name == 'regular_lonlat':
            grid_keys = ['input_lon', 'input_lat', 'input_position',
                         'X_resolution', 'Y_resolution']
            if set(self.grid.keys()) != set(grid_keys):
                raise epygramError("grid attribute must consist in keys: " +
                                   str(grid_keys))
            dimensions_keys = ['X', 'Y']
            if set(self.dimensions.keys()) != set(dimensions_keys):
                raise epygramError("dimensions attribute must consist in keys: " + str(dimensions_keys))

    def ij2xy(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of gridpoints, \n
        (*x, y*) being the coordinates in the projection.
        - *position*: position to return with respect to model cell.
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
        x = Xorigin + (i - i0 + oi) * Xresolution
        y = Yorigin + (j - j0 + oj) * Yresolution

        return (x, y)

    def xy2ij(self, x, y, position=None):
        """
        (*x, y*) being the coordinates in the projection, \n
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        - *position*: position represented by (x,y) with respect to model cell.
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

        return (i, j)

    def ij2ll(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of gridpoints, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.
        - *position*: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        return self.ij2xy(i, j, position)

    def ll2ij(self, lon, lat, position=None):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        - *position*: lat lon position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.

        Caution: (*i,j*) are float.
        """

        lon = degrees_nearest_mod(lon, self._center_lon.get('degrees'))
        return self.xy2ij(lon, lat, position)

    def ll2xy(self, lon, lat):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*x, y*) being the coordinates in the projection.

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
        (*x, y*) being the coordinates in the projection, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.

        Note that origin of coordinates in projection is the center of the C+I
        domain.
        """

        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)

        return self.ij2ll(*self.xy2ij(x, y))

    def make_basemap(self, gisquality='i', subzone=None, specificproj=None,
                     zoom=None, ax=None):
        """
        Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
        projection (if available). This is designed to avoid explicit handling
        of deep horizontal geometry attributes.

        Args: \n
        - *gisquality*: defines the quality of GIS contours, cf. Basemap doc. \n
          Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
          among: 'C', 'CI'.
        - *specificproj*: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        - *zoom*: specifies the lon/lat borders of the map, implying hereby
          a 'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=).\n
          Overwrites *specificproj*.
        - *ax*: a matplotlib ax on which to plot; if None, plots will be done
          on matplotlib.pyplot.gca()
        """
        from mpl_toolkits.basemap import Basemap

        # corners
        (llcrnrlon, llcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ll'])
        (urcrnrlon, urcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ur'])

        # make basemap
        if zoom is not None:
            # zoom case
            specificproj = 'cyl'
            llcrnrlon = zoom['lonmin']
            llcrnrlat = zoom['latmin']
            urcrnrlon = zoom['lonmax']
            urcrnrlat = zoom['latmax']
        if specificproj is None:
            # defaults
            if llcrnrlat <= -89.0 or \
               urcrnrlat >= 89.0:
                proj = 'cyl'
            else:
                proj = 'merc'
            b = Basemap(resolution=gisquality, projection=proj,
                        llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                        urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                        ax=ax)
        else:
            # specificproj
            lon0 = self._center_lon.get('degrees')
            lat0 = self._center_lat.get('degrees')
            if specificproj == 'kav7':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'ortho':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            lat_0=lat0,
                            ax=ax)
            elif specificproj == 'cyl':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'moll':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif isinstance(specificproj, tuple) and \
                 specificproj[0] == 'nsper' and \
                 isinstance(specificproj[1], dict):
                sat_height = specificproj[1].get('sat_height', 3000) * 1000.
                b = Basemap(resolution=gisquality,
                            projection=specificproj[0],
                            lon_0=specificproj[1].get('lon', lon0),
                            lat_0=specificproj[1].get('lat', lat0),
                            satellite_height=sat_height,
                            ax=ax)

        return b

    def global_shift_center(self, longitude_shift):
        """
        Shifts the center of the geometry by *longitude_shift* (in degrees).
        *longitude_shift* has to be a multiple of the grid's resolution in
        longitude.
        """

        corners = self.gimme_corners_ll()
        if abs(abs(degrees_nearest_mod(corners['ur'][0] - corners['ul'][0], 0.)) - self.grid['X_resolution'].get('degrees')) \
           < config.epsilon:
            if abs(longitude_shift % self.grid['X_resolution'].get('degrees')) > config.epsilon:
                raise epygramError("*longitude_shift* has to be a multiple" +
                                   " of the grid's resolution in longitude.")
            self._center_lon = Angle(self._center_lon.get('degrees') + longitude_shift,
                                     'degrees')
            self.grid['input_lon'] = Angle(self.grid['input_lon'].get('degrees') + longitude_shift,
                                           'degrees')
        else:
            raise epygramError("unable to shift center if " +
                               "lon_max - lon_min != X_resolution.")

    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up in the geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        *num* is the number of points, including point1 and point2.
        """

        if num < 2:
            raise epygramError("'num' must be at least 2.")
        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        xy_linspace = list(zip(numpy.linspace(x1, x2, num=num),
                               numpy.linspace(y1, y2, num=num)))

        return [self.xy2ll(*xy) for xy in xy_linspace]

    def distance(self, end1, end2):
        """
        Computes the distance between two points along a straight line in the
        geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).

        Warning: requires the :mod:`pyproj` module.
        """
        import pyproj

        plast = end1
        distance = 0
        g = pyproj.Geod(ellps='sphere')
        for p in self.linspace(end1, end2, 1000)[1:]:
            distance += g.inv(plast[0], plast[1], *p)[2]
            plast = p

        return distance

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.
        *point* must be a tuple (lon, lat).
        """
        return self.resolution_ij(*self.ll2ij(lon, lat))

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
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

    def azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points in geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        """

        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)

        return (numpy.degrees(numpy.arctan2(x2 - x1, y2 - y1)) + 180.) % 360. - 180.

    def _what_grid(self, out, arpifs_var_names=False):
        """
        Writes in file a summary of the grid of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *arpifs_var_names*: if True, prints the equivalent 'arpifs' variable
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
                        grid['X_resolution'].get('degrees') * dimensions['X'])
        if arpifs_var_names:
            varname = ' (ELY)'
        write_formatted(out, "Domain width in Y, in deg" + varname,
                        grid['Y_resolution'].get('degrees') * dimensions['Y'])
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

    def __eq__(self, other):
        """Test of equality by recursion on the object's attributes."""
        if self.__class__ == other.__class__ and \
           set(self.__dict__.keys()) == set(other.__dict__.keys()):
            selfcp = self.copy()
            othercp = other.copy()
            for obj in [selfcp, othercp]:
                obj.grid.pop('input_lon', None)
                obj.grid.pop('input_lat', None)
                obj.grid.pop('input_position', None)
                obj._proj = None
            ok = super(D3RegLLGeometry, selfcp).__eq__(othercp)
        else:
            ok = False
        return ok


class D3ProjectedGeometry(D3RectangularGridGeometry):
    """
    Handles the geometry for a Projected 3-Dimensions Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic', 'space_view'])),
            projection=dict(
                type=FPDict,
                info="Handles projection information."),
            projtool=dict(
                optional=True,
                default=config.default_projtool,
                info="To use pyproj or epygram.myproj."),
            geoid=dict(
                type=FPDict,
                optional=True,
                default=FPDict(config.default_geoid),
                info="Geoid definition in projections.")
        )
    )

    @property
    def secant_projection(self):
        """ Is the projection secant to the sphere ? (or tangent)"""
        return 'secant_lat' in self.projection \
            or 'secant_lat1' in self.projection

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(D3ProjectedGeometry, self).__init__(*args, **kwargs)

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
                    (center[0] - self.grid['input_position'][0]) * self.grid['X_resolution'],
                    (center[1] - self.grid['input_position'][1]) * self.grid['Y_resolution'],
                    'xy2ll')
                # xc, yc: coordinates of center point in non rotated proj
                xc = x1 + dx
                yc = y1 + dy
                center_lon, center_lat = p(xc, yc, inverse=True)
                self._center_lon = Angle(center_lon, 'degrees')
                self._center_lat = Angle(center_lat, 'degrees')

        if self.projtool == 'pyproj':
            import pyproj
            projtool = pyproj
            projdict = {'lambert':'lcc',
                        'mercator':'merc',
                        'polar_stereographic':'stere',
                        'space_view':'geos'}
        elif self.projtool == 'myproj':
            raise Warning("use of 'myproj' projtool is DEPRECATED ! Should rather use 'pyproj' instead, check config.")
            from epygram import myproj
            projtool = myproj
            projdict = {'lambert':'lambert',
                        'mercator':'mercator',
                        'polar_stereographic':'polar_stereographic',
                        'space_view':'space_view'}
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
            p = projtool.Proj(proj=proj,
                              lon_0=self.projection['reference_lon'].get('degrees'),
                              lat_1=lat_1, lat_2=lat_2,
                              **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = projtool.Proj(proj=proj,
                                       lon_0=self.projection['reference_lon'].get('degrees'),
                                       lat_1=lat_1, lat_2=lat_2,
                                       x_0=-x0, y_0=-y0,
                                       **self.geoid)
        elif self.name == 'mercator':
            if self.secant_projection:
                lat_ts = self.projection['secant_lat'].get('degrees')
            else:
                lat_ts = 0.
            p = projtool.Proj(proj=proj,
                              lon_0=self.projection['reference_lon'].get('degrees'),
                              lat_ts=lat_ts,
                              **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = projtool.Proj(proj=proj,
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
            p = projtool.Proj(proj=proj,
                              lon_0=self.projection['reference_lon'].get('degrees'),
                              lat_0=lat_0, lat_ts=lat_ts,
                              **self.geoid)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = projtool.Proj(proj=proj,
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
            p = projtool.Proj(proj=proj,
                              h=height,
                              lon_0=lonSat)
            compute_center_proj(p, centerPoint)
            x0, y0 = p(self._center_lon.get('degrees'),
                       self._center_lat.get('degrees'))
            self._proj = projtool.Proj(proj=proj,
                                       h=height,
                                       lon_0=lonSat,
                                       x_0=-x0, y_0=-y0)
        else:
            raise NotImplementedError("projection: " + self.name)

    def _rotate_axis(self, x, y, direction):
        """
        Internal method used to rotate x/y coordinates to handle rotated geometries.
        *dircetion*, if evals to 'xy2ll', direction used when converting (x, y) to (lat, lon).
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
        """
        Check that the geometry is consistent.
        """
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
        grid_keys = ['LAMzone', 'X_resolution', 'Y_resolution',
                     'input_lon', 'input_lat', 'input_position']
        if set(self.grid.keys()) != set(grid_keys):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys))
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

    def select_subzone(self, subzone):
        """
        If a LAMzone defines the geometry, select only the *subzone* from it
        and return a new geometry object.
        *subzone* among ('C', 'CI').
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
                centerPoint = ((float(self.dimensions['X']) - 1) / 2.,
                               (float(self.dimensions['Y']) - 1) / 2.)  # Coordinates of center point
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
        Returns the coordinate of the grid center as a tuple (center_lon, center_lat).
        """

        return (self._center_lon, self._center_lat)

    def ij2xy(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints, \n
        (*x, y*) being the coordinates in the projection.
        - *position*: position to return with respect to model cell.
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
        (*x, y*) being the coordinates in the projection, \n
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints.
        - *position*: position represented by (x,y) with respect to model cell.
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
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*x, y*) being the coordinates in the projection.

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
        (*x, y*) being the coordinates in the projection, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.

        Note that origin of coordinates in projection is the center of the
        C+I domain.
        """

        if isinstance(x, list) or isinstance(x, tuple):
            x = numpy.array(x)
        if isinstance(y, list) or isinstance(y, tuple):
            y = numpy.array(y)

        (a, b) = self._rotate_axis(x, y, direction='xy2ll')
        ll = self._proj(a, b, inverse=True)

        return ll

    def ij2ll(self, i, j, position=None):
        """
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.
        - *position*: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        return self.xy2ll(*self.ij2xy(i, j, position))

    def ll2ij(self, lon, lat, position=None):
        """
        (*lon, lat*) being the lon/lat coordinates in degrees, \n
        (*i, j*) being the coordinates in the 2D matrix of CIE gridpoints.
        - *position*: lat lon position with respect to the model cell.
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
        - *position*: grid position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        f = fpx.field(structure='H2D', geometry=self, fid={'geometry':'Map Factor'})
        lats = self.get_lonlat_grid(position=position)[1]
        data = self.map_factor(lats)
        f.setdata(data)

        return f

    def make_basemap(self, gisquality='i', subzone=None, specificproj=None,
                     zoom=None):
        """
        Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
        projection (if available). This is designed to avoid explicit handling
        of deep horizontal geometry attributes.

        Args: \n
        - *gisquality*: defines the quality of GIS contours, cf. Basemap doc. \n
          Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
        - *subzone*: defines the LAM subzone to be included, in LAM case,
          among: 'C', 'CI'.
        - *specificproj*: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        - *zoom*: specifies the lon/lat borders of the map, implying hereby a
          'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=). \n
          Overwrites *specificproj*.
        """
        from mpl_toolkits.basemap import Basemap

        if zoom is not None:
            if specificproj not in [None, 'cyl']:
                raise epygramError("projection can only be cyl in zoom mode.")
            specificproj = 'cyl'

        if specificproj is None:
            # corners
            if self.projection['rotation'].get('degrees') == 0.:
                (llcrnrlon, llcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ll'])
                (urcrnrlon, urcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ur'])
            else:
                (imin, jmin) = self.gimme_corners_ij(subzone)['ll']
                (imax, jmax) = self.gimme_corners_ij(subzone)['ur']
                border = [(imin, j) for j in range(jmin, jmax + 1)] + \
                         [(imax, j) for j in range(jmin, jmax + 1)] + \
                         [(i, jmin) for i in range(imin, imax + 1)] + \
                         [(i, jmax) for i in range(imin, imax + 1)]
                ilist, jlist = list(zip(*border))
                (x, y) = self.ij2xy(numpy.array(ilist), numpy.array(jlist))  # in model coordinates
                (x, y) = self._rotate_axis(x, y, direction='xy2ll')  # non-rotated coordinates
                (llcrnrlon, llcrnrlat) = self.xy2ll(*self._rotate_axis(x.min(), y.min(), direction='ll2xy'))
                (urcrnrlon, urcrnrlat) = self.xy2ll(*self._rotate_axis(x.max(), y.max(), direction='ll2xy'))
            # defaults
            if self.name == 'lambert':
                if self.secant_projection:
                    lat_0 = (self.projection['secant_lat1'].get('degrees') +
                             self.projection['secant_lat2'].get('degrees')) / 2.
                    b = Basemap(resolution=gisquality, projection='lcc',
                                lat_1=self.projection['secant_lat1'].get('degrees'),
                                lat_2=self.projection['secant_lat2'].get('degrees'),
                                lat_0=lat_0,
                                lon_0=self.projection['reference_lon'].get('degrees'),
                                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
                else:
                    b = Basemap(resolution=gisquality, projection='lcc',
                                lat_0=self.projection['reference_lat'].get('degrees'),
                                lon_0=self.projection['reference_lon'].get('degrees'),
                                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif self.name == 'mercator':
                if self.secant_projection:
                    lat = 'secant_lat'
                else:
                    lat = 'reference_lat'
                b = Basemap(resolution=gisquality, projection='merc',
                            lat_ts=self.projection[lat].get('degrees'),
                            lon_0=self._center_lon.get('degrees'),
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif self.name == 'polar_stereographic':
                if self.secant_projection:
                    lat = 'secant_lat'
                else:
                    lat = 'reference_lat'
                b = Basemap(resolution=gisquality, projection='stere',
                            lat_ts=self.projection[lat].get('degrees'),
                            lat_0=numpy.copysign(90., self.projection[lat].get('degrees')),
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif self.name == 'space_view':
                b = Basemap(resolution=gisquality, projection='geos',
                            lon_0=self.projection['satellite_lon'].get('degrees'))
            else:
                raise epygramError("Projection name unknown.")
        else:
            # corners
            if zoom:
                llcrnrlon = zoom['lonmin']
                llcrnrlat = zoom['latmin']
                urcrnrlon = zoom['lonmax']
                urcrnrlat = zoom['latmax']
            else:
                (imin, jmin) = self.gimme_corners_ij(subzone)['ll']
                (imax, jmax) = self.gimme_corners_ij(subzone)['ur']
                border = [(imin, j) for j in range(jmin, jmax + 1)] + \
                         [(imax, j) for j in range(jmin, jmax + 1)] + \
                         [(i, jmin) for i in range(imin, imax + 1)] + \
                         [(i, jmax) for i in range(imin, imax + 1)]
                ilist, jlist = list(zip(*border))
                (lons, lats) = self.ij2ll(numpy.array(ilist), numpy.array(jlist))
                llcrnrlon, urcrnrlon = lons.min(), lons.max()
                llcrnrlat, urcrnrlat = lats.min(), lats.max()
            # specificproj
            lon0 = self._center_lon.get('degrees')
            lat0 = self._center_lat.get('degrees')
            if specificproj == 'kav7':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif specificproj == 'ortho':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            lat_0=lat0)
            elif specificproj == 'cyl':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif specificproj == 'moll':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            elif isinstance(specificproj, tuple) and \
                 specificproj[0] == 'nsper' and \
                 isinstance(specificproj[1], dict):
                sat_height = specificproj[1].get('sat_height', 3000) * 1000.
                b = Basemap(resolution=gisquality,
                            projection=specificproj[0],
                            lon_0=specificproj[1].get('lon', lon0),
                            lat_0=specificproj[1].get('lat', lat0),
                            satellite_height=sat_height)
            else:
                raise epygramError('unknown **specificproj**: ' + str(specificproj))

        return b

    def linspace(self, end1, end2, num):
        """
        Returns evenly spaced points over the specified interval.
        Points are lined up in the geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        *num* is the number of points, including point1 and point2.
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
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
        """

        (x1, y1) = self.ll2xy(*end1)
        (x2, y2) = self.ll2xy(*end2)
        # distance on the projected surface
        distance = numpy.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        # map factor computed as the mean over 100 points along the line
        # between point1 and point2
        mean_map_factor = numpy.array([self.map_factor(lat) for (_, lat) in self.linspace(end1, end2, 100)]).mean()

        return distance / mean_map_factor

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.
        *point* must be a tuple (lon, lat).
        """
        return self.resolution_ij(*self.ll2ij(lon, lat))

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
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

    def azimuth(self, end1, end2):
        """
        Initial bearing from *end1* to *end2* points in geometry.
        *end1* must be a tuple (lon, lat).
        *end2* must be a tuple (lon, lat).
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

        Args: \n
        - *end1* must be a tuple (lon, lat).
        - *end2* must be a tuple (lon, lat).
        - *points_number* defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        - *resolution* defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        - *position* defines the position of data in the grid (defaults to 'center')
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

        return fpx.geometry(**kwargs_geom)

    def compass_grid(self, subzone=None, position=None):
        """
        Get the compass grid, i.e. the angle between Y-axis and North for each
        gridpoint.

        - *subzone*: optional, among ('C', 'CI'), for LAM grids only, returns
          the grid resp. for the C or C+I zone off the C+I+E zone. \n
          Default is no subzone, i.e. the whole field.
        - *position*: position of lonlat grid with respect to the model cell.
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
                self.projection.get('rotation', 0.).get('degrees')  # TOBECHECKED: rotation

        return theta

    def reproject_wind_on_lonlat(self, u, v, lon, lat,
                                 map_factor_correction=True,
                                 reverse=False):
        """
        Reprojects a wind vector (u, v) on the grid onto real
        axes, i.e. with components on true zonal/meridian axes.
        lon/lat are the point(s) coordinates.

        If *map_factor_correction*, applies a
        correction of magnitude due to map factor.

        If *reverse*, apply the reverse reprojection.
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

    def _what_projection(self, out, arpifs_var_names=False):
        """
        Writes in file a summary of the projection of the field's grid.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *arpifs_var_names*: if True, prints the equivalent 'arpifs' variable
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
        if self.grid.get('LAMzone', False):
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
            varname = ' (ELAT0)'
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
                        grid['X_resolution'] * dimX)
        if arpifs_var_names:
            varname = ' (ELY)'
        write_formatted(out,
                        "Domain width (of C+I) in Y, in metres" + varname,
                        grid['Y_resolution'] * dimY)
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

    def __eq__(self, other):
        """Test of equality by recursion on the object's attributes."""
        if self.__class__ == other.__class__ and \
           set(self.__dict__.keys()) == set(other.__dict__.keys()):
            selfcp = self.deepcopy()
            othercp = other.deepcopy()
            for obj in [selfcp, othercp]:
                obj.grid.pop('input_lon', None)
                obj.grid.pop('input_lat', None)
                obj.grid.pop('input_position', None)
                obj._proj = None
            ok = super(D3ProjectedGeometry, selfcp).__eq__(othercp)
            ok = ok and (self.getcenter() == other.getcenter())
        else:
            ok = False
        return ok


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

    def _consistency_check(self):
        """Check that the geometry is consistent."""

        grid_keys = ['dilatation_coef', 'latitudes']
        if self.name == 'rotated_reduced_gauss':
            grid_keys.extend(['pole_lon', 'pole_lat'])
        if set(self.grid.keys()) != set(grid_keys):
            raise epygramError("grid attribute must consist in keys: " +
                               str(grid_keys))
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
        (*i, j*) being the coordinates in the 2D matrix of gridpoints, \n
        (*lon, lat*) being the lon/lat coordinates in degrees.
        - *position*: lat lon position to return with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        if self._getoffset(position) != (0., 0.):
            raise NotImplementedError("horizontal staggered grids for reduced\
                                       gauss grid are not implemented.")

        if isinstance(i, list) or isinstance(i, tuple) or\
           isinstance(i, numpy.ndarray):
            lat = [self.grid['latitudes'][n].get('degrees') for n in j]
            lon = [(numpy.pi * 2 * i[n]) / self.dimensions['lon_number_by_lat'][j[n]] for n in range(len(j))]
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

        (*lon, lat*) being the lon/lat coordinates in degrees, \n

        - *position*: position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        if isinstance(lon, list) or isinstance(lon, tuple):
            lon = numpy.array(lon)
        if isinstance(lat, list) or isinstance(lat, tuple):
            lat = numpy.array(lat)

        return self.nearest_points(lon, lat, {'n': '1'}, position)

    def _allocate_colocation_grid(self, compressed=False, as_float=False):
        """
        Creates the array for lonlat grid.
        Just a trick to avoid recomputing the array for several fields that
        share their geometry.
        If *compressed*, return 1D arrays, else 2D masked arrays.
        If *as_float*, return arrays with dtype float64, else int64.
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
    def gridpoints_number(self, **useless):
        """Returns the number of gridpoints of the grid."""
        return sum(self.dimensions['lon_number_by_lat'])

    def get_lonlat_grid(self, position=None, d4=False, nb_validities=0, **useless):
        """
        Returns a tuple of two tables containing one the longitude of each
        point, the other the latitude, with 2D shape.

        Shape of 2D data in Gauss grids: \n
          - grid[0, 0:Nj] is first (Northern) band of latitude, masked after
            Nj = number of longitudes for latitude j \n
          - grid[-1, 0:Nj] is last (Southern) band of latitude (idem).

        Args:\n
        - *position*: position of lonlat grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
                d4=True requires nb_validities > 0
        - *nb_validities* is the number of validities represented in data values
        """
        # !!! **useless enables the method to receive arguments specific to
        #     other geometries but useless here ! Do not remove.

        if hasattr(self, '_buffered_gauss_grid') and self._buffered_gauss_grid.get('filled'):
            lons = self._buffered_gauss_grid['lons']
            lats = self._buffered_gauss_grid['lats']
        else:
            (igrid, jgrid) = self._allocate_colocation_grid(compressed=True, as_float=False)
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
                    self._buffered_gauss_grid['lons'][...] = lons[...]
                    self._buffered_gauss_grid['lats'][...] = lats[...]
                self._buffered_gauss_grid['filled'] = True

        if d4:
            lons, lats = self._reshape_lonlat_4d(lons, lats, nb_validities)
        elif not d4 and nb_validities != 0:
            raise ValueError("*nb_validities* must be 0 when d4==False")

        return (lons, lats)

    def get_datashape(self, force_dimZ=None, dimT=None, d4=False, **useless):
        """
        Returns the data shape according to the geometry.

        - *force_dimZ*: if supplied, force the Z dimension instead of that
          of the vertical geometry
        - *dimT* if supplied, is the time dimension to be added to the
          data shape
        - *d4*: if True,  shape is 4D (need to specify *dimT*)
                if False, shape is 3D if dimZ > 1 else 2D
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

        - *data*: the 1D data (or 3D with a T and Z dimensions,
          or 2D with either a T/Z dimension, to be specified),
          of dimension concording with geometry. In case data is 3D, T must be
          first dimension and Z the second.
        - *first_dimension*: in case data is 2D, specify what is the first
          dimension of data among ('T', 'Z')
        - *d4*: if True,  returned values are shaped in a 4 dimensions array
                if False, shape of returned values is determined with respect to geometry
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

        if fill_value is None:
            fill_value = data.fill_value
        data_copy = data.copy()
        for t in range(data_copy.shape[0]):
            for k in range(data_copy.shape[1]):
                for j in range(self.dimensions['lat_number']):
                    for i in range(self.dimensions['lon_number_by_lat'][j]):
                        if data_copy.mask[t, k, j, i]:
                            data_copy[t, k, j, i] = fill_value

        return data_copy

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

    def make_basemap(self,
                     gisquality='i',
                     specificproj=None,
                     zoom=None,
                     ax=None, **kwargs):
        """
        Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
        projection (if available). This is designed to avoid explicit handling
        of deep horizontal geometry attributes.

        Args: \n
        - *gisquality*: defines the quality of GIS contours, cf. Basemap doc. \n
          Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
        - *specificproj*: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        - *zoom*: specifies the lon/lat borders of the map, implying hereby a
          'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=). \n
          Overwrites *specificproj*.
        - *ax*: a matplotlib ax on which to plot; if None, plots will be done
          on matplotlib.pyplot.gca()
        """
        # !!! **kwargs enables the method to receive arguments specific to
        #     other geometries, useless here ! Do not remove.
        from mpl_toolkits.basemap import Basemap

        gisquality = 'l'  # forced for Gauss, for time-consumption reasons...

        # corners
        llcrnrlon = -180
        llcrnrlat = -90
        urcrnrlon = 180
        urcrnrlat = 90

        # make basemap
        if zoom is not None:
            # zoom case
            specificproj = 'cyl'
            llcrnrlon = zoom['lonmin']
            llcrnrlat = zoom['latmin']
            urcrnrlon = zoom['lonmax']
            urcrnrlat = zoom['latmax']
        if specificproj is None:
            # defaults
            if 'rotated' in self.name:
                lon0 = self.grid['pole_lon'].get('degrees')
            else:
                lon0 = 0.0
            b = Basemap(resolution=gisquality, projection='moll',
                        lon_0=lon0,
                        ax=ax)
        else:
            # specificproj
            if 'rotated' in self.name:
                lon0 = self.grid['pole_lon'].get('degrees')
                lat0 = self.grid['pole_lat'].get('degrees')
            else:
                lon0 = 0.0
                lat0 = 0.0
            if specificproj == 'kav7':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'ortho':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            lat_0=lat0,
                            ax=ax)
            elif specificproj == 'cyl':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif specificproj == 'moll':
                b = Basemap(resolution=gisquality, projection=specificproj,
                            lon_0=lon0,
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            elif isinstance(specificproj, tuple) and \
                 specificproj[0] == 'nsper' and \
                 isinstance(specificproj[1], dict):
                sat_height = specificproj[1].get('sat_height', 3000) * 1000.
                b = Basemap(resolution=gisquality,
                            projection=specificproj[0],
                            lon_0=specificproj[1].get('lon', lon0),
                            lat_0=specificproj[1].get('lat', lat0),
                            satellite_height=sat_height,
                            ax=ax)

        return b

    def resolution_ll(self, lon, lat):
        """
        Returns the local resolution at the nearest point of lon/lat.
        It's the distance between this point and its closest neighbour.
        *point* must be a tuple (lon, lat).
        """
        return self.resolution_ij(*self.ll2ij(lon, lat))

    def resolution_ij(self, i, j):
        """
        Returns the distance to the nearest point of (i,j) point.
        (*i, j*) being the coordinates in the 2D matrix of gridpoints.
        """

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

    def point_is_inside_domain_ll(self, lon, lat,
                                  margin=-0.1,
                                  position=None):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.

        Args: \n
        - *lon*: longitude of point(s) in degrees.
        - *lat*: latitude of point(s) in degrees.
        - *margin*: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        - *position*: position of the grid with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        try:
            N = len(lon)
        except Exception:
            N = 1
        if N == 1:
            inside = True
        else:
            inside = [True for _ in range(N)]

        return inside

    def point_is_inside_domain_ij(self, i, j, margin=-0.1):
        """
        Returns True if the point(s) of lon/lat coordinates is(are) inside the
        field.

        Args: \n
        - *lon*: longitude of point(s) in degrees.
        - *lat*: latitude of point(s) in degrees.
        - *margin*: considers the point inside if at least 'margin' points far
          from the border. The -0.1 default is a safety for precision errors.
        """

        try:
            N = len(i)
        except Exception:
            N = 1
        if N == 1:
            inside = True
            if j >= self.dimensions['lat_number'] or j < 0:
                inside = False
            if i >= self.dimensions['lon_number_by_lat'][j] or i < 0:
                inside = False
        else:
            inside = [self.point_is_inside_domain_ij(i[n], j[n], margin=margin)
                      for n in range(N)]

        return inside

    def _rotate_stretch(self, lon, lat, reverse=False):
        """
        Internal method used to transform true lon/lat into rotated and stretched lon/lat
        (*lon, lat*) being the lon/lat coordinates in degrees.

        If *reverse*, do the reverse transform.

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
            lon = degrees_nearest_mod(lon, self.grid.get('pole_lon', Angle(90., 'degrees')).get('degrees'))
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
        (*lon*, *lat*) point in degrees.
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
        - *position*: grid position with respect to the model cell.
          Defaults to self.position_on_horizontal_grid.
        """

        f = fpx.field(structure='H2D', geometry=self, fid={'geometry':'Map Factor'})
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
        lon/lat are the point(s) coordinates on real sphere.

        If *map_factor_correction*, applies a
        correction of magnitude due to map factor.

        If *reverse*, apply the reverse reprojection.
        """

        pc = self.grid['dilatation_coef']
        plac = self.grid['pole_lat'].get('degrees')
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
                epylog.warning('*map_factor_correction* should be **True** !')
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
                       external_distance=None):
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
        """

        if self._getoffset(position) != (0., 0.):
            raise NotImplementedError("horizontal staggered grids for " +
                                      "reduced gauss grid are not implemented.")

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
            """

            if not numpy.all(numpy.array(latitudes[1:]) <= numpy.array(latitudes[:-1])):
                raise ValueError('latitudes must be in descending order')
            distmin = None
            nearest = None
            for n in range(0, len(latitudes)):
                dist = latrs - latitudes[n]
                if distmin is None or abs(dist) < abs(distmin):
                    distmin = dist
                    nearest = n
            if num == 2:
                lats = [nearest,
                        (nearest - numpy.copysign(1, distmin)).astype('int')]
            elif num == 3:
                lats = [nearest - 1, nearest, nearest + 1]
            elif num % 2:  # odd
                lats = [nearest - k for k in reversed(range(1, num // 2 + 1))] + \
                       [nearest, ] + \
                       [nearest + k for k in range(1, num // 2 + 1)]
            elif not num % 2:  # even
                lats = [nearest - k for k in reversed(range(1, num // 2))] + \
                       [nearest, ] + \
                       [nearest + k for k in range(1, num // 2 + 1)]
            return lats

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
            """

            lonrs = numpy.radians(lonrs)
            if latnum < 0:
                # We are near the pole, have to look-up on the first latitudes
                # circles, symmetrically with regards to the pole
                j = abs(latnum) - 1
                lonnummax = self.dimensions['lon_number_by_lat'][j]
                i = ((lonrs - numpy.pi) % (numpy.pi * 2)) * \
                    lonnummax / (numpy.pi * 2)
            elif latnum >= len(latitudes):
                # j = len(latitudes) - 1
                j = len(latitudes) - (latnum - len(latitudes)) - 1  # TOBECHECKED: next circle past the pole
                lonnummax = self.dimensions['lon_number_by_lat'][j]
                i = ((lonrs - numpy.pi) % (numpy.pi * 2)) * \
                    lonnummax / (numpy.pi * 2)
            else:
                j = latnum
                lonnummax = self.dimensions['lon_number_by_lat'][latnum]
                i = lonrs * lonnummax / (numpy.pi * 2)
            if num == 1:
                lons = (numpy.rint(i).astype('int') % lonnummax, j)
            elif num == 2:
                lons = [(numpy.floor(i).astype('int') % lonnummax, j),
                        ((numpy.floor(i).astype('int') + 1) % lonnummax, j)]
            elif num % 2:  # odd
                raise NotImplementedError("but is it necessary ?")
            elif not num % 2:  # even
                ii = numpy.floor(i).astype('int')
                lons = [((ii - k) % lonnummax, j) for k in reversed(range(1, num // 2))] + \
                       [(ii % lonnummax, j)] + \
                       [((ii + k) % lonnummax, j) for k in reversed(range(1, num // 2 + 1))]

            return lons

        def nearest(lon, lat, lonrs, latrs):
            """
            Internal method used to find the nearest point.
            lon/lat are the true coordinate, lonrs/latrs are rotated and
            streched coordinates.
            """

            distance = None
            nearpoint = None
            for latnum in nearest_lats(latrs, 3):
                point = nearest_lons(lonrs, latnum, 1)
                dist = self.distance((lon, lat), self.ij2ll(*point))
                if distance is None or dist < distance:
                    nearpoint = point
                    distance = dist
            return nearpoint

        def nearests(lonrs, latrs, n):
            """
            Internal methods used to find the n*n points surrunding the
            point lonrs/latrs, lonrs/latrs are rotated and stretched
            coordinates.
            """

            points = []
            for j in nearest_lats(latrs, n):
                points.extend(nearest_lons(lonrs, j, n))
            return points

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
            if isinstance(lat, numpy.ndarray) and lat.shape != ():
                j = numpy.zeros(len(lat), dtype=numpy.int)
                i = numpy.zeros(len(lat), dtype=numpy.int)
                for k in range(0, len(lat)):
                    (i[k], j[k]) = nearest(lon[k], lat[k], lonrs[k], latrs[k])
            else:
                (i, j) = nearest(lon, lat, lonrs, latrs)
            points = [(i, j)]
        # 2.: several points are needed
        else:
            # 2.1: how many ?
            if external_distance:
                n = 1
            elif request.get('radius'):
                if isinstance(lat, numpy.ndarray):
                    raise NotImplementedError("request:{'radius':..., 'shape':'radius'} with several points.")
                resolution = self.resolution_ll(lon, lat)
                n = max(numpy.around(float(request['radius']) / float(resolution)).astype('int') * 2,
                        1)
            elif nsquare_match:
                n = int(nsquare_match.group('n'))
            else:
                raise epygramError("unrecognized **request**: " + str(request))
            # 2.2: get points
            if isinstance(lat, numpy.ndarray):
                i = numpy.zeros((len(lat), n * n), dtype=numpy.int)
                j = numpy.zeros((len(lat), n * n), dtype=numpy.int)
                for k in range(len(lat)):
                    points = nearests(lonrs[k], latrs[k], n)
                    for p in range(n * n):
                        i[k, p], j[k, p] = points[p]
                points = [(i[:, p], j[:, p]) for p in range(n * n)]
            else:
                points = nearests(lonrs, latrs, n)
            # 2.3: only select the nearest with regards to external_distance
            if external_distance:
                distance = None
                for p in points:
                    dist = abs(external_distance['external_field'].getvalue_ij(*p, one=True) - external_distance['target_value'])
                    if distance is None or dist < distance:
                        nearpoint = p
                        distance = dist
                points = [nearpoint]
            elif request.get('shape') == 'circle':
                points = [(i, j)
                          for (i, j) in points
                          if self.distance((lon, lat), self.ij2ll(i, j)) <= request['radius']]

        return points[0] if request == {'n':'1'} else points

    def _what_grid_dimensions(self, out, spectral_geometry=None):
        """
        Writes in file a summary of the grid & dimensions of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *spectral_geometry*: an optional dict containing the spectral
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

    def __eq__(self, other):
        """Test of equality by recursion on the object's attributes."""
        if self.__class__ == other.__class__ and \
           set(self.__dict__.keys()).discard('_buffered_gauss_grid') == \
           set(other.__dict__.keys()).discard('_buffered_gauss_grid'):
            selfcp = self.deepcopy()
            selfcp._clear_buffered_gauss_grid()
            othercp = other.deepcopy()
            othercp._clear_buffered_gauss_grid()
            ok = super(D3GaussGeometry, selfcp).__eq__(othercp)
        else:
            ok = False
        return ok

footprints.collectors.get(tag='geometrys').fasttrack = ('structure', 'name')
