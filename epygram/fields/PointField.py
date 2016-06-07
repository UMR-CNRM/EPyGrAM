#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to for Point (0D == 1 value) fields.
"""


from .D3Field import D3Field
from epygram.geometries import PointGeometry



class PointField(D3Field):
    """
    0-Dimension (point) field class.
    A field is defined by its identifier 'fid',
    its data, its geometry, and its validity.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['Point'])),
            geometry=dict(
                type=PointGeometry),
        )
    )


