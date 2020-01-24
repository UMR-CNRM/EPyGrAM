#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the Basemap-interfaced facilities to plot field.

.. deprecated:: 1.3.11
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import numpy
from mpl_toolkits.basemap import Basemap

import footprints

from epygram import epygramError

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([make_basemap,],
                               plugin_name)
    from epygram.geometries.D3Geometry import D3ProjectedGeometry
    D3ProjectedGeometry.make_basemap = make_basemap


def make_basemap(self,
                 gisquality='i',
                 subzone=None,
                 specificproj=None,
                 zoom=None,
                 ax=None):
    """
    Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
    projection (if available). This is designed to avoid explicit handling
    of deep horizontal geometry attributes.

    .. deprecated:: 1.3.9

    :param gisquality: defines the quality of GIS contours, cf. Basemap doc.
      Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
    :param subzone: defines the LAM subzone to be included, in LAM case,
      among: 'C', 'CI'.
    :param specificproj: enables to make basemap on the specified projection,
      among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
      In 'nsper' case, the {} may contain:\n
      - 'sat_height' = satellite height in km;
      - 'lon' = longitude of nadir in degrees;
      - 'lat' = latitude of nadir in degrees. \n
      Overwritten by *zoom*.
    :param zoom: specifies the lon/lat borders of the map, implying hereby a
      'cyl' projection.
      Must be a dict(lonmin=, lonmax=, latmin=, latmax=).
      Overwrites *specificproj*.
    :param ax: a matplotlib ax on which to plot; if None, plots will be done
      on matplotlib.pyplot.gca()
    """
    epylog.info("The 'make_basemap' method is deprecated, please use 'default_cartopy_CRS'")
    if zoom is not None:
        if specificproj not in [None, 'cyl', 'merc']:
            raise epygramError("projection can only be cyl/merc in zoom mode.")
        specificproj = True

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
                b = Basemap(resolution=gisquality,
                            projection='lcc',
                            lat_1=self.projection['secant_lat1'].get('degrees'),
                            lat_2=self.projection['secant_lat2'].get('degrees'),
                            lat_0=lat_0,
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
            else:
                b = Basemap(resolution=gisquality,
                            projection='lcc',
                            lat_0=self.projection['reference_lat'].get('degrees'),
                            lon_0=self.projection['reference_lon'].get('degrees'),
                            llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                            ax=ax)
        elif self.name == 'mercator':
            if self.secant_projection:
                lat = 'secant_lat'
            else:
                lat = 'reference_lat'
            b = Basemap(resolution=gisquality,
                        projection='merc',
                        lat_ts=self.projection[lat].get('degrees'),
                        lon_0=self._center_lon.get('degrees'),
                        llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                        urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                        ax=ax)
        elif self.name == 'polar_stereographic':
            if self.secant_projection:
                lat = 'secant_lat'
            else:
                lat = 'reference_lat'
            b = Basemap(resolution=gisquality,
                        projection='stere',
                        lat_ts=self.projection[lat].get('degrees'),
                        lat_0=numpy.copysign(90., self.projection[lat].get('degrees')),
                        lon_0=self.projection['reference_lon'].get('degrees'),
                        llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                        urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                        ax=ax)
        elif self.name == 'space_view':
            b = Basemap(resolution=gisquality,
                        projection='geos',
                        lon_0=self.projection['satellite_lon'].get('degrees'),
                        ax=ax)
        else:
            raise epygramError("Projection name unknown.")
    else:
        # corners
        if zoom:
            llcrnrlon = zoom['lonmin']
            llcrnrlat = zoom['latmin']
            urcrnrlon = zoom['lonmax']
            urcrnrlat = zoom['latmax']
            if llcrnrlat <= -89.0 or \
               urcrnrlat >= 89.0:
                specificproj = 'cyl'
            else:
                specificproj = 'merc'
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
            b = Basemap(resolution=gisquality,
                        projection=specificproj,
                        lon_0=lon0,
                        llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                        urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                        ax=ax)
        elif specificproj == 'ortho':
            b = Basemap(resolution=gisquality,
                        projection=specificproj,
                        lon_0=lon0,
                        lat_0=lat0,
                        ax=ax)
        elif specificproj in ('cyl', 'merc'):
            b = Basemap(resolution=gisquality,
                        projection=specificproj,
                        llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                        urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                        ax=ax)
        elif specificproj == 'moll':
            b = Basemap(resolution=gisquality,
                        projection=specificproj,
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
        else:
            raise epygramError('unknown **specificproj**: ' + str(specificproj))

    return b
