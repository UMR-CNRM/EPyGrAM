#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Miscellaneous arguments
"""

from . import _defaults

d = {
    'LAMzone':[
        '-z', '--zone',
        dict(help="zone of the domain, in LAM case:\
                   'C' for zone C, 'CI' for zone C+I, 'CIE' for zone C+I+E.\
                   Default is 'CI'.",
             choices=['C', 'CI', 'CIE'],
             default='CI')],
    'array_flattening_order':[
        '--order',
        dict(help="for LAM arrays, whether to flatten in C (row-major) or\
                   Fortran (column-major) order 2D arrays. Default = 'C'.",
             default='C')],
    'flatten_horizontal_grids':[
        '--flatten',
        dict(help="for netCDF, flatten 2D horizontal grids to 1D.",
             action='store_true',
             default=False)],
    'operation_on_field':[
        '-x', '--operation',
        dict(help="do the requested operation on field right after reading it. \
                   Syntax: '-,273.15' (e.g. for K => C) or 'exp' to take the \
                   exponential of the field. \
                   The operand must be among (+,-,*,/) or 'normalize' or any\
                   numpy function.\
                   For '-' operand, use short-name option -x without spacetab\
                   between option and argument.",
             default=None)],
    'diffoperation_on_field':[
        '-X', '--diffoperation',
        dict(help="do the requested operation on difference field right after \
                   computing it. \
                   Syntax: '-,273.15' (e.g. for K => C) or 'exp' to take the \
                   exponential of the field. \
                   The operand must be among (+,-,*,/) or 'normalize' or any\
                   numpy function.\
                   For '-' operand, use short-name option -X without spacetab\
                   between option and argument.",
             default=None)],
    'pressure_unit_hpa':[
        '-u', '--pressure_unit_hpa',
        dict(action='store_true',
             help="converts pressure in Pa to hPa.",
             default=False)],
    'composition_with_field':[
        '--compose_with',
        dict(help="compose a transformed field with another field. \
                   Syntax: 'otherfield, composition, [file=otherfile], [preset='norm'|'ceil']'. \
                   The *composition* must be among (+,-,*,/).\
                   If *file* is filled (optional), the *otherfield* is read in it.\
                   If *preset* is filled (optional), one or several operation can be done \
                   on the *otherfield* before operation.\
                   Ex: 'SFX.FRAC_WATER, *, file=mypgd.fa, preset=ceil' \
                   will set non-water points of the plotted field (e.g. SFX.TS_WATER) to 0.",
             default=None)],
    'mask_threshold':[
        '--mt', '--mask_threshold',
        dict(help="set a threshold to mask values. E.g. 'min=0.0' will\
                   mask negative values. 'min=0.0,max=1e8' will mask values\
                   outside these boundaries.",
             default=None,
             dest='mask_threshold')],
    'femars.diff_to_avg':[
        '--avg', '--diff_to_avg',
        dict(help="instead of computing member to member differences, \
                   compute member-to-average-of-members differences.",
             action='store_true',
             default=False,
             dest='diff_to_avg')],
    'map_factor_correction':[
        '--mpf', '--map_factor_correction',
        dict(help="apply a map factor correction to wind.",
             action='store_true',
             default=False,
             dest='map_factor_correction')],
    'wind_components_are_projected_on':[
        '--wpo', '--wind_projected_on',
        dict(help="specify on which coordinates the wind components are projected.",
             default=None,
             choices=['lonlat', 'grid'],
             dest='wind_components_are_projected_on')],
                }

