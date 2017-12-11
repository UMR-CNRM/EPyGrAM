#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities for building a LAM domain.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from footprints import proxy as fpx

# parameters
#: minimum width of Extension zone
Ezone_minimum_width = 11
#: threshold in degrees towards Equator for the domain center, to choose lambert/mercator
threshold_mercator_lambert = 1.
#: threshold in degrees towards Pole for the min/max latitude, to choose lambert/polar_stereographic
threshold_pole_distance_lambert = 1.
maxdims_security_barrier = 10000
vkw = {'structure': 'V',
       'typeoffirstfixedsurface': 1,
       'levels': [1]}
vgeom = fpx.geometry(**vkw)
projections_s2g = {'L':'lambert', 'M':'mercator', 'PS':'polar_stereographic'}
projections_g2p = {'lambert':'Lambert (conformal conic)', 'mercator':'Mercator', 'polar_stereographic':'Polar Stereographic'}
projections_s2p = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}
projections_g2s = {v:k for k, v in projections_s2g.items()}


def default_Iwidth(resolution):
    """
    Return default Iwidth depending on the resolution.
    """
    return 16 if resolution < 2000. else 8  # FIXME: go beyond for smaller resolutions ?
