#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with output
"""

from . import _defaults

d = {
    'output':[
        '-o', '--output',
        dict(help="store graphical output in file in with specified format,\
                   among ('png', pdf). Pdf is kind of disadvised, for it is\
                   very slow and produces much too big files. 'X' to open GUI.",
             choices=['png', 'pdf', 'X'],
             default=_defaults.get('default_graphical_output'))],
    'outputfmt':[
        '-o', '--outputfmt',
        dict(help="specify format for output file,\
                   among ('png', pdf). Pdf is kind of disadvised, for it is\
                   very slow and produces much too big files. 'X' to open GUI.",
             choices=['png', 'pdf', 'X'],
             default=_defaults.get('default_graphical_output'))],
    'savefig':[
        '--sf', '--savefig',
        dict(help="save figure in file.",
             action='store_true',
             dest='savefig',
             default=False)],
    'outputfilename':[
        '-O', '--outputfilename',
        dict(help="store output in the specified filename (without format for\
                   graphical output, to be completed by -o/--output).",
             default=False)],
    'one_pdf':[
        '-p', '--pdf',
        dict(action='store_true',
             help="store output (eventually multiple plots) in one .pdf file.",
             default=False)],
    'noplot':[
        '-n', '--noplot',
        dict(action='store_true',
             help="disable plot. Profile/spectrum will be computed\
                   and written as text output but not plotted.",
             default=False)],
    'GeoPoints_llv':[
        '--llv',
        dict(action='store_true',
             help="simplify GeoPoints format to LAT, LON, VALUE only.",
             default=False)],
    'precision':[
        '-e', '--precision',
        dict(help="precision on values (number of decimals in scientific\
                   format). Default = 4.",
             default=4)],
    'lonlat_precision':[
        '-E', '--lonlat_precision',
        dict(help="precision on longitudes/latitudes (number of decimals). \
                   Default = 4.",
             default=4)],
    'GeoPoints_precision':[
        '-e', '--precision',
        dict(help="precision on GeoPoints' values (number of decimals in \
                   scientific format). Default = " +
                   str(_defaults.get('GeoPoints_precision', 4)) + '.',
             default=_defaults.get('GeoPoints_precision', 4))],
    'GeoPoints_lonlat_precision':[
        '-E', '--lonlat_precision',
        dict(help="precision on GeoPoints' longitudes/latitudes (number of \
                   decimals). Default = " +
                   str(_defaults.get('GeoPoints_lonlat_precision', 4)) + '.',
             default=_defaults.get('GeoPoints_precision', 4))],
    'get_field_compression':[
        '-c', '--compression',
        dict(action='store_true',
             help="get compression options of each field.",
             default=False)],
    'get_field_details':[
        '-d', '--details',
        dict(help="get some details about each field.\
                   E.g. 'spectral' to inquire the spectralness of fields, or\
                   'compression' to inquire compression parameters, 'grid' or \
                   'comment' for LFI.",
             default=None)],
    'GRIB_short_fid':[
        '-s', '--grib_short_fid',
        dict(help="condense GRIB fid.",
             action='store_true',
             default=False)],
    'GRIB_other_options':[
        '--GRIB_other_options',
        dict(help="a set of key:value pairs, separated by commas, to be given to\
                   the GRIB message when writing.\
                   Ex: 'typeOfGeneratingProcess:12,productionStatusOfProcessedData:2'.",
             default=None)],
    'stdout':[
        '-o', '--stdout',
        dict(action='store_true',
             help="redirects output to standard output (rather than file).",
             default=False)],
    'output_format':[
        '-o', '--output_format',
        dict(help="format of conversion output, among:\
                   'grb' (GRIB2), 'geo' (GeoPoints), 'nc' (netCDF4).",
             choices=['grb', 'geo', 'nc'],
             required=True)],
    'reproject_wind':[
        '--rw', '--reproject_wind',
        dict(action='store_true',
             help="reprojects a wind vector (u, v) onto true\
                   zonal/meridian axes (assuming it is projected onto grid axes\
                   in fields).",
             dest='reproject_wind',
             default=False)],
    'only_maxdiff':[
        '--oxd', '--only_maxdiff',
        dict(action='store_true',
             help='if diffonly and activated, only print the maximum absolute\
                   difference (and min/max values of fields for magnitude',
             dest='only_maxdiff',
             default=False)],
    'GeoPoints_columns':[
        '--geopoints_cols',
        dict(dest='geopoints_cols',
             help="set manually the columns of GeoPoints, for instance 'LAT,LON,VALUE'.",
             default=None)],
    'GeoPoints_noheader':[
        '--geopoints_noheader',
        dict(action='store_true',
             dest='geopoints_noheader',
             help="discard the header of GeoPoints.",
             default=False)],
                  }

