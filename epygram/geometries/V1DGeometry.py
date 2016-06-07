#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for Horizontal 2D geometries of fields.
"""

from .V2DGeometry import V2DUnstructuredGeometry

from epygram import epygramError



class V1DGeometry(V2DUnstructuredGeometry):
    """
    Handles the geometry for a Vertical 1-Dimension Field.
    Abstract mother class.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V1D'])),
            name=dict(
                values=set(['unstructured',
                            'DDH:point', 'DDH:ij_point', 'DDH:quadrilateral',
                            'DDH:rectangle', 'DDH:globe', 'DDH:zonal_bands']),
                optional=True)
        )
    )


    def _consistency_check(self):
        """Check that the geometry is consistent."""
        if self.dimensions['X'] != 1:
            raise epygramError("V1DGeometry must have only one point in x-direction.")
        super(V1DGeometry, self)._consistency_check()

    def distance(self, end1, end2):
        """No meaning with a one-point geometry."""
        raise epygramError("this method has no meaning with a one-point geometry.")

    def azimuth(self, end1, end2):
        """No meaning with a one-point geometry."""
        raise epygramError("this method has no meaning with a one-point geometry.")
