#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all classes needed to plot fields using cartopy.
"""
from __future__ import print_function, absolute_import, unicode_literals, division

from . import (H2DField, H2DVectorField, D3AcademicGeometry, D3GaussGeometry,
               D3Geometry, D3ProjectedGeometry, D3RectangularGeometry)


def activate():
    """Activate plugin."""
    H2DField.activate()
    H2DVectorField.activate()
    D3AcademicGeometry.activate()
    D3GaussGeometry.activate()
    D3Geometry.activate()
    D3ProjectedGeometry.activate()
    D3RectangularGeometry.activate()
