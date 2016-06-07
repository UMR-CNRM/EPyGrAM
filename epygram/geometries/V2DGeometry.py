#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for Horizontal 2D geometries of fields.
"""

from .D3Geometry import D3Geometry, D3RectangularGridGeometry, \
                        D3UnstructuredGeometry, D3ProjectedGeometry
 
from epygram import epygramError



class V2DGeometry(D3Geometry):
    """
    Handles the geometry for a Vertical 2-Dimensions Field.
    Abstract mother class.
    """

    _abstract = True
    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),
        )
    )
    
    
    def _consistency_check(self):
        """Check that the geometry is consistent."""
        if self.dimensions['Y'] != 1:
            raise epygramError("V2DGeometry must have only one point in y-direction.")
        super(V2DGeometry, self)._consistency_check()



class V2DRectangularGridGeometry(V2DGeometry, D3RectangularGridGeometry):
    """
    Handles the geometry for a Verical 2-Dimensions Field for which the surface points
    come from a rectangular grid.
    Abstract.
    """

    _abstract = True
    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  #inheritance priority problem
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic', 'space_view', 'unstructured']))
        )
    )

class V2DUnstructuredGeometry(V2DRectangularGridGeometry, D3UnstructuredGeometry):
    """Handles the geometry for an unstructured Horizontal 2-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  #inheritance priority problem
            name=dict(
                values=set(['unstructured']))
        )
    )

class V2DProjectedGeometry(V2DRectangularGridGeometry, D3ProjectedGeometry):
    """Handles the geometry for a Projected Horizontal 2-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  #inheritance priority problem
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic', 'space_view']))
        )
    )