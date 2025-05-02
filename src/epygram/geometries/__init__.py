#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all geometries classes.
"""

import numpy

from bronx.syntax.decorators import nicedeco

@nicedeco
def _need_pyproj_geod(mtd):
    """
    Decorator for Geometry object: if the method needs a pyproj.Geod
    object to be set.
    """
    def with_geod(self, *args, **kwargs):
        if not hasattr(self, '_pyproj_geod'):
            self._set_geoid()
        return mtd(self, *args, **kwargs)
    return with_geod

from .AbstractGeometry import Geometry, RectangularGridGeometry

from .ProjectedGeometry import ProjectedGeometry
from .UnstructuredGeometry import UnstructuredGeometry
from .AcademicGeometry import AcademicGeometry
from .GaussGeometry import GaussGeometry
from .RegLLGeometry import RegLLGeometry
from .RotLLGeometry import RotLLGeometry

from .VGeometry import VGeometry
from .SpectralGeometry import SpectralGeometry, truncation_from_gridpoint_dims
from . import domain_making

Pressure = 100
Altitude = 102
Height = 103
PV = 109
HybridH = 118
HybridP = 119
GeneralVerticalHeightCoordinate = 150

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

def gauss_latitudes(nlat):
    """Compute the Gauss latitudes for a **nlat** points grid."""
    x, _ = numpy.polynomial.legendre.leggauss(nlat)
    return numpy.degrees(numpy.arcsin(x[::-1]))
