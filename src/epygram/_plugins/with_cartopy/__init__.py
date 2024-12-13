#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all classes needed to plot fields using cartopy.
"""

from . import (H2DField, H2DVectorField, AcademicGeometry, GaussGeometry,
               Geometry, ProjectedGeometry, RectangularGeometry)


def activate():
    """Activate plugin."""
    H2DField.activate()
    H2DVectorField.activate()
    AcademicGeometry.activate()
    GaussGeometry.activate()
    Geometry.activate()
    ProjectedGeometry.activate()
    RectangularGeometry.activate()
