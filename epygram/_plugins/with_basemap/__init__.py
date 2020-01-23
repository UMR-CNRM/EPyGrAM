#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all classes needed to plot fields using Basemap.
"""
from __future__ import print_function, absolute_import, unicode_literals, division

from . import (H2DField, H2DVectorField,
               D3UnstructuredGeometry, D3GaussGeometry, D3ProjectedGeometry,
               D3RegLLGeometry, D3RotLLGeometry)

def activate():
    """Activate plugin."""
    H2DField.activate()
    H2DVectorField.activate()
    D3UnstructuredGeometry.activate()
    D3GaussGeometry.activate()
    D3ProjectedGeometry.activate()
    D3RegLLGeometry.activate()
    D3RotLLGeometry.activate()
