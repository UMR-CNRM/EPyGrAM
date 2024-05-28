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
from epygram.util import separation_line, write_formatted
from .D3Geometry import (D3Geometry, D3RectangularGridGeometry,
                         D3UnstructuredGeometry, D3ProjectedGeometry,
                         D3AcademicGeometry)


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

    def what(self, out=sys.stdout,
             vertical_geometry=True,
             **_):
        """
        Writes in file a summary of the geometry.

        :param out: the output open file-like object
        :param vertical_geometry: if True, writes the vertical geometry of the
          field.
        """
        out.write("###########################\n")
        out.write("### HORIZONTAL GEOMETRY ###\n")
        out.write("###########################\n")
        self._what_grid_dimensions(out)
        self._what_grid(out)
        out.write(separation_line)
        out.write("\n")
        if vertical_geometry:
            self.vcoordinate.what(out)

    def _what_grid(self, out):
        """
        Writes in file a summary of the grid of the field.

        :param out: the output open file-like object
        """
        (lons, lats) = self.get_lonlat_grid()
        write_formatted(out, "Kind of Geometry", 'Unstructured')
        write_formatted(out, "Max Longitude in deg", lons.max())
        write_formatted(out, "Min Longitude in deg", lons.min())
        write_formatted(out, "Max Latitude in deg", lats.max())
        write_formatted(out, "Min Latitude in deg", lats.min())


class V2DRectangularGridGeometry(V2DGeometry, D3RectangularGridGeometry):
    """
    Handles the geometry for a Vertical 2-Dimensions Field for which the surface points
    come from a rectangular grid.
    Abstract.
    """

    _abstract = True
    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  # inheritance priority problem
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
                values=set(['V2D'])),  # inheritance priority problem
            name=dict(
                values=set(['unstructured'])),
            position_on_horizontal_grid=dict(
                default='center',
                values=set(['center'])),
        )
    )


class V2DAcademicGeometry(D3AcademicGeometry, V2DRectangularGridGeometry):
    """Handles the geometry for an academic Horizontal 2-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  # inheritance priority problem
            name=dict(
                values=set(['academic']))
        )
    )


class V2DProjectedGeometry(V2DRectangularGridGeometry, D3ProjectedGeometry):
    """Handles the geometry for a Projected Horizontal 2-Dimensions Field."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),  # inheritance priority problem
            name=dict(
                values=set(['lambert', 'mercator', 'polar_stereographic', 'space_view']))
        )
    )
