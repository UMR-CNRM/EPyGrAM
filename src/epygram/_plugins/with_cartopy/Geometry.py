#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend Geometry with plotting methods using cartopy.
"""

import numpy

import cartopy.crs as ccrs

import footprints

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([cartopy_CRS_reproject],
                               plugin_name)
    from epygram.geometries import Geometry
    # defaults arguments for cartopy plots
    Geometry.cartopy_CRS_reproject = cartopy_CRS_reproject


def cartopy_CRS_reproject(self, lons, lats, projection=None):
    """Reproject lons/lats onto a cartopy CRS projection coordinates."""
    if projection is None:
        projection = self.default_cartopy_CRS()
    if isinstance(lons, (float, int)):
        lons = [lons]
        lats = [lats]
    tmp = projection.transform_points(ccrs.PlateCarree(),
                                      numpy.array(lons),
                                      numpy.array(lats))
    return tmp[..., 0], tmp[..., 1]
