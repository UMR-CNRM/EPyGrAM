#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities for building a LAM domain.
"""

import math

from epygram.geometries import VGeometry

# parameters
#: minimum width of Extension zone
Ezone_minimum_width = 11
#: threshold in degrees towards Equator for the domain center, to choose lambert/mercator
threshold_mercator_lambert = 5.
#: threshold in degrees towards Pole for the min/max latitude, to choose lambert/polar_stereographic
threshold_pole_distance_lambert = 1.
#: default width of Izone, to compute width in number of points according to resolution
default_Izone_width_in_m = 20000.

maxdims_security_barrier = 10000
vgeom = VGeometry(typeoffirstfixedsurface= 1, levels=[1])
projections_s2g = {'L':'lambert', 'M':'mercator', 'PS':'polar_stereographic'}
projections_g2p = {'lambert':'Lambert (conformal conic)', 'mercator':'Mercator', 'polar_stereographic':'Polar Stereographic'}
projections_s2p = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}
projections_g2s = {v:k for k, v in projections_s2g.items()}


def default_Iwidth(resolution, Izone_width_in_m=default_Izone_width_in_m):
    """
    Return default Iwidth depending on the resolution.

    Algo:
    make it at least **Izone_width_in_m** wide, and not less than 8 points.
    """
    n = int(math.ceil(float(Izone_width_in_m) / resolution))
    if n % 2 != 0:  # make it even
        n += 1
    n = max(n, 8)  # not less than 8 points
    return n
