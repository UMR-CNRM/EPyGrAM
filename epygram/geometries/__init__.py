#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all geometries classes.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from footprints import proxy as fpx

from .H2DGeometry import H2DGeometry
from .D3Geometry import D3Geometry, D3AcademicGeometry, _need_pyproj_geod
from .VGeometry import VGeometry
from .V1DGeometry import V1DGeometry
from .V2DGeometry import V2DGeometry
from .H1DGeometry import H1DGeometry
from .PointGeometry import PointGeometry
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
    return fpx.geometry(levels=[0], structure='V', typeoffirstfixedsurface=1)


def build_geometry(help_on=None, **kwargs):
    """
    Proxy to build geometry from scratch.
    All arguments are to be given to build geometry.

    :param help_on: name of the attribute to have some help on.
    """
    g = fpx.geometry(**kwargs)
    if g is None and help_on is not None:
        print('H' * 80)
        print('Help_on:', help_on)
        attr_map = fpx.geometrys.build_attrmap()[help_on]
        attr_map = {am['name']:am['values'] for am in attr_map}
        for c in sorted(attr_map.keys()):
            print('* Class', c, ':', attr_map[c])
    return g
