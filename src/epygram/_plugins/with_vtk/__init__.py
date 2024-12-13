#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all classes needed to plot 3D fields through vtk.
"""

from . import _D3CommonField, D3VectorField, Geometry


def activate():
    """Activate plugin."""
    _D3CommonField.activate()
    D3VectorField.activate()
    Geometry.activate()
