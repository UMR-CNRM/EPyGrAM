#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend Geometry with plotting methods using cartopy.
"""

import cartopy.crs as ccrs

import footprints

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([default_cartopy_CRS],
                               plugin_name)
    from epygram.geometries import RectangularGridGeometry
    # defaults arguments for cartopy plots
    RectangularGridGeometry.default_cartopy_CRS = default_cartopy_CRS


def default_cartopy_CRS(self):
    """
    Create a cartopy.crs appropriate to the Geometry.

    By default, a PlateCarree (if the domain gets close to a pole)
    or a Miller projection is returned.
    """
    (lons, lats) = self.get_lonlat_grid()
    if lons.ndim == 1:
        lonmax = lons[:].max()
        lonmin = lons[:].min()
        latmax = lats[:].max()
        latmin = lats[:].min()
    else:
        lonmax = lons[:, -1].max()
        lonmin = lons[:, 0].min()
        latmax = lats[-1, :].max()
        latmin = lats[0, :].min()
    center_lon = (lonmax + lonmin) / 2.
    if latmin <= -80.0 or latmax >= 84.0:
        crs = ccrs.PlateCarree(center_lon)
    else:
        if 'a' in self.geoid and 'b' in self.geoid:
            globe = ccrs.Globe(semimajor_axis=self.geoid['a'], semiminor_axis=self.geoid['b'])
        elif 'ellps' in self.geoid:
            globe = ccrs.Globe(ellipse=self.geoid['ellps'])
        else:
            raise NotImplementedError("self.geoid={}".format(str(self.geoid)))
        crs = ccrs.Mercator(center_lon, globe=globe)
    return crs
