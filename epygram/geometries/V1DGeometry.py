#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes for Horizontal 2D geometries of fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import sys

from epygram import epygramError
from epygram.util import write_formatted
from .V2DGeometry import V2DUnstructuredGeometry


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

    def _what_grid(self, out=sys.stdout):
        """
        Writes in file a summary of the grid of the field.

        :param out: the output open file-like object
        """
        (lons, lats) = (self.grid['longitudes'], self.grid['latitudes'])
        write_formatted(out, "Kind of Geometry", 'Unstructured')
        write_formatted(out, "Longitude in deg", lons[0])
        write_formatted(out, "Latitude in deg", lats[0])

    def __eq__(self, other):
        """Test of equality by recursion on the object's attributes."""
        if self.__class__ == other.__class__ and \
           set(self._attributes.keys()) == set(other._attributes.keys()):
            for attr in self._attributes.keys():
                ok = self._attributes[attr] == other._attributes[attr]
                if not ok:
                    break
        else:
            ok = False
        return ok

    def __hash__(self):
        # known issue __eq__/must be defined both or none, else inheritance is broken
        return super(V1DGeometry, self).__hash__()
