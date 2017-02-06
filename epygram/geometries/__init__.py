#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all geometries classes.
"""

from .H2DGeometry import H2DGeometry, H2DUnstructuredGeometry
from .D3Geometry import D3Geometry
from .VGeometry import VGeometry
from .V1DGeometry import V1DGeometry
from .V2DGeometry import V2DGeometry
from .PointGeometry import PointGeometry
from .SpectralGeometry import SpectralGeometry

Pressure = 100
Altitude = 102
Height = 103
PV = 109
HybridH = 118
HybridP = 119
