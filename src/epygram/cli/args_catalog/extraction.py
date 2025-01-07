#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with extraction stuff
"""

from . import _defaults

d = {
    'point_coordinates':[
        '-c', '--coordinates',
        dict(help="lon/lat coordinates of the point.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'section_starting_point':[
        '-s', '--starting_point',
        dict(help="lon/lat coordinate of starting point of the section.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'section_ending_point':[
        '-e', '--ending_point',
        dict(help="lon/lat coordinate of ending point of the section.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'horizontal_interpolation':[
        '-i', '--interpolation',
        dict(help="interpolation mode from field grid to point/section\
                   coordinates. Among ('nearest', 'linear', 'cubic',\
                   'bilinear'). Defaults to 'nearest'.",
             choices=['nearest', 'linear', 'cubic', 'bilinear'],
             default='nearest')],
    'external_distance':[
        '-z', '--external_distance',
        dict(help="for 'nearest' interpolation mode, the nearest is chosen\
                   among the 4 nearest points regarding the point whose value\
                   of field *EXT* is the closest to *VALUE*, with syntax:\
                   *external_distance* = 'VALUE; EXT'. EXT being a fid.",
             default=None)],
    'section_transect_points_number':[
        '-p', '--points_number',
        dict(help="number of points from starting point to ending point\
                   (included). Defaults to the number of points computed from\
                   the fields resolution, or from the resolution given via\
                   option -r.",
             type=int,
             default=None)],
    'section_transect_resolution':[
        '-r', '--resolution',
        dict(help="resolution of the section. Defaults to the fields\
                   resolution, or computed from the number of points given via\
                   option -p.",
             type=float,
             default=None)],
    'verticalcoord2pressure':[
        '-P', '--verticalcoord2pressure',
        dict(action='store_const',
             dest='Yconvert',
             const='pressure',
             help="compute pressure as vertical coordinate, if possible.",
             default=None)],
    'verticalcoord2height':[
        '-H', '--verticalcoord2height',
        dict(action='store_const',
             dest='Yconvert',
             const='height',
             help="compute height as vertical coordinate, if possible.",
             default=None)],
    'verticalcoord2altitude':[
        '-A', '--verticalcoord2altitude',
        dict(action='store_const',
             dest='Yconvert',
             const='altitude',
             help="compute altitude as vertical coordinate, if possible.",
             default=None)],
    'no_cheap_height_conversion':[
        '--no_cheap_height',
        dict(action='store_false',
             dest='cheap_height',
             help="for the computation of heights (-A/-H) to be done \
                   taking hydrometeors into account (in R \
                   computation) and NH Pressure departure \
                   (Non-Hydrostatic data). Slower but more accurate.",
             default=True)]
                      }

