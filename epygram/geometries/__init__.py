#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all geometries classes.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from .Geometry import (Geometry, RectangularGridGeometry, ProjectedGeometry,
                       UnstructuredGeometry, AcademicGeometry, GaussGeometry,
                       gauss_latitudes, RegLLGeometry, RotLLGeometry,
                       _need_pyproj_geod)
from .VGeometry import VGeometry
from .SpectralGeometry import SpectralGeometry, truncation_from_gridpoint_dims
from . import domain_making

Pressure = 100
Altitude = 102
Height = 103
PV = 109
HybridH = 118
HybridP = 119

#: Mapping of vertical coordinates
vertical_coordinates = {'Pressure':Pressure,
                        'Altitude':Altitude,
                        'Height':Height,
                        'PV':PV,
                        'HybridH':HybridH,
                        'HybridP':HybridP}


def build_surf_VGeometry():
    """Build a surface vertical geometry."""
    return VGeometry(levels=[0], typeoffirstfixedsurface=1)
