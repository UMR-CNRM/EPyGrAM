#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments for domain_maker
"""

from . import _defaults

d = {
    'mode':[
        '-l',
        dict(dest='mode',
             action='store_const',
             help="mode: define domain by specifying a lon/lat domain that must\
                   be included inside.",
             const='lonlat_included',
             default='center_dims')],
    'no_display':[
        '-n', '--no_display',
        dict(action='store_true',
             help="run without displaying the domain.",
             default=False)],
    'maximize_CI_in_E':[
        '-m', '--maximize_CI_in_E',
        dict(action='store_true',
             help="forces the C+I zone to be the greatest possible inside a\
                   given (discrete) C+I+E size. \
                   In other words, with a given (discrete) C+I+E size, forces\
                   the E-zone to be the smallest possible.",
             default=False)],
    'truncation':[
        '-t', '--truncation',
        dict(help="(for namelists) the kind of truncation of spectral geometry\
                   to generate, among ('linear', 'quadratic', 'cubic').",
             default='linear')],
    'orography_subtruncation':[
        '-o', '--orography_subtruncation',
        dict(help="(for namelists) additional subtruncation for orography to be generated.",
             default='quadratic')],
                        }

