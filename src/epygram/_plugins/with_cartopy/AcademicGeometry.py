#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend Geometry with plotting methods using cartopy.
"""

import footprints

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([default_cartopy_CRS],
                               plugin_name)
    from epygram.geometries import AcademicGeometry
    # defaults arguments for cartopy plots
    AcademicGeometry.default_cartopy_CRS = default_cartopy_CRS


def default_cartopy_CRS(self):  # TODO: externalize to cartopy plugin
        """
        Create a cartopy.crs appropriate to the Geometry.

        In this case, cartopy is not used but raw matplotlib,
        so returned CRS is None.
        """
        return None
