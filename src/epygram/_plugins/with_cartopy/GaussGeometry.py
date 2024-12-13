#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend Geometry with plotting methods using cartopy.
"""

from cartopy import crs as ccrs

import footprints


epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([default_cartopy_CRS],
                               plugin_name)
    from epygram.geometries import GaussGeometry
    # defaults arguments for cartopy plots
    GaussGeometry.default_cartopy_CRS = default_cartopy_CRS


def default_cartopy_CRS(self):
        """
        Create a cartopy.crs appropriate to the Geometry.
        """
        if 'rotated' in self.name:
            lon_0 = self.grid['pole_lon'].get('degrees')
        else:
            lon_0 = 0.
        #Mollweide projection does not handle elliptical globes
        #Thus, we do not set the globe attribute from the geometry geoid 
        crs = ccrs.Mollweide(central_longitude=lon_0)
        return crs
