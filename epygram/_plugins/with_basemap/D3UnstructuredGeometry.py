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


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([make_basemap,],
                               plugin_name)
    from epygram.geometries.D3Geometry import D3UnstructuredGeometry
    D3UnstructuredGeometry.make_basemap = make_basemap


def make_basemap(self,
                 gisquality='i',
                 specificproj=None,
                 zoom=None,
                 ax=None,
                 **_):
    """
    Returns a :class:`matplotlib.basemap.Basemap` object of the 'ad hoc'
    projection (if available). This is designed to avoid explicit handling
    of deep horizontal geometry attributes.

    :param gisquality: defines the quality of GIS contours, cf. Basemap doc. \n
      Possible values (by increasing quality): 'c', 'l', 'i', 'h', 'f'.
    :param specificproj: enables to make basemap on the specified projection,
      among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc).
      In 'nsper' case, the {} may contain:\n
      - 'sat_height' = satellite height in km;
      - 'lon' = longitude of nadir in degrees;
      - 'lat' = latitude of nadir in degrees. \n
      Overwritten by *zoom*.
    :param zoom: specifies the lon/lat borders of the map, implying hereby a
      'cyl' projection.
      Must be a dict(lonmin=, lonmax=, latmin=, latmax=). \n
      Overwrites *specificproj*.
    :param ax: a matplotlib ax on which to plot; if None, plots will be done
      on matplotlib.pyplot.gca()
    """
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
