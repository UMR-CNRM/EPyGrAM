#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for Point geometry of fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from .V1DGeometry import V1DGeometry
from .H2DGeometry import H2DUnstructuredGeometry


class PointGeometry(V1DGeometry, H2DUnstructuredGeometry):
    """
    Handles the geometry for a Geographical Point Field.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['Point'])),
            name=dict(
                values=set(['DDH:point', 'DDH:ij_point', 'DDH:quadrilateral',
                            'DDH:rectangle', 'DDH:globe', 'DDH:zonal_bands',
                            'unstructured']),
                optional=True),
            position_on_horizontal_grid=dict(
                values=['center'],
                default='center')
        )
    )
