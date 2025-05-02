#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with fields
"""

from epygram.extra import griberies
from . import _defaults

d = {
    'field':[
        '-f', '-F', '--field',
        dict(help="*fid* = *field identifier* of of the field(s) to be\
                   processed. Syntax depends on format:\
                   GRIB: handgrip, e.g. 'shortName:t,level:850'.\
                   FA: name, e.g. 'S050TEMPERATURE'; regular expressions may\
                   be used, such as 'S00[2-6]WIND.[U-V].PHYS',\
                   'SURFALBEDO*' or 'SURF?.OF.OZONE'.\
                   To obtain the list of fields in file, use the 'epy_what'\
                   tool.",
             default=None)],
    'windfield':[
        '-w', '--computewind',
        dict(help="to process wind as a vector.\
                   Syntax depends on format:\
                   FA: 'S*WIND' or 'P*VENT', 'H*VENT', 'V*VENT', 'CLSVENT.*',\
                   or 'CLS*.RAF'; * designing the level: combination of\
                   (*, ?, digits).\
                   GRIB: handgrip, e.g. 'shortName:u+v,level:850' or\
                   'indicatorOfParameter:33+34,level:850' or \
                   'parameterCategory:2,parameterNumber:2+3,level:850'.\
                   Not implemented in difference mode.",
             default=None)],
    'windfieldU':[
        '--wU', '--Ucomponentofwind',
        dict(help="to process wind as a vector. U component of wind.\
                   (same syntax as -f).\
                   Not implemented in difference mode.",
             dest="Ucomponentofwind",
             default=None)],
    'windfieldV':[
        '--wV', '--Vcomponentofwind',
        dict(help="to process wind as a vector. V component of wind.\
                   (same syntax as -f).\
                   Not implemented in difference mode.",
             dest="Vcomponentofwind",
             default=None)],
    'FA_field':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   To obtain the list of fields in file, use the 'fa_what'\
                   tool.",
             default=None)],
    'FA_windfield':[
        '-w', '--computewind',
        dict(help="to process wind module, using the following syntax:\
                   'S*WIND' or 'P*VENT', 'H*VENT', 'V*VENT', 'CLSVENT.*', 'CLS*.RAF',\
                   * designing the level: combination of (*, ?, digits).\
                   Not implemented in difference mode.",
             default=None)],
    'FA_multiple_fields':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   Regular expressions can be used, such as\
                   'S00[2-6]WIND.[U-V].PHYS' to process meridian and zonal winds of levels 2 to 6,\
                   'SURFALBEDO*' to process all kinds of albedo,\
                   'SURF?.OF.OZONE' to process all kinds of ozone (A,B,C).\
                   In case this option is not used, no more than '-l' option,\
                   all fields present in file are processed.\
                   Fields not found in file are ignored.\
                   To obtain the list of fields in file, use the 'fa_what'\
                   tool.",
             default=None)],
    'vertical_field':[
        '-F', '-f', '--field',
        dict(help="designation of the fields from which to extract profile.\
                   Syntax depends on format:\
                   GRIB: handgrip designing the parameter and the type of\
                   levels.\
                   FA: name of the upper-air field to be processed,\
                   with a * for level.\
                   To obtain the list of fields in file, use the 'epy_what'\
                   tool.",
             required=True,
             default=None)],
    'FA_vertical_field':[
        '-F', '-f', '--field',
        dict(help="name of the upper-air field to be processed, with a * for\
                   level. To obtain the list of fields in file, use the\
                   'fa_what' tool.",
             required=True,
             default=None)],
    'external_vertical_coord':[
        '-Z', '--externaleVCoord',
        dict(help="fid kind: fid is the field identifier to be used\
                   as vertical coordinate and kind is the vertical coordinate\
                   code (eg. 100 for pressure)",
             nargs=2,
             dest='Yconvert',
             default=None)],
    'DDHLFA_multiple_fields':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   Regular expressions can be used, such as\
                   'F*FLUVERTDYN' or 'VKK[0-1]'.\
                   In case this option is not used, no more than '-l' option,\
                   all fields present in file are processed.\
                   Fields not found in file are ignored.\
                   To obtain the list of fields in file, use the 'ddhlfa_what'\
                   tool.",
             default=None)],
    'DDHLFA_domain':[
        '-n', '--domain_number',
        dict(help="number of the domain to plot. \
                   Multiple domains must be given within quotes, separated by comma. \
                   To obtain the list of domains in file, \
                   use the 'ddhlfa_what' tool.",
             required=True)],
    'GRIB_field':[
        '-F', '-f', '--field',
        dict(help="(possibly partial) handgrip of the field to be processed.\
                   To obtain the list of fields in file (and their handgrip), \
                   use EPyGrAM's 'grib_what' tool or GRIB_API's 'grib_dump'.\
                   Ex: -f 'shortName:t,level:850'.",
             default=None)],
    'GRIB_secondfield':[
        '-2', '--secondfield',
        dict(help="(possibly partial) handgrip of the second field to be\
                   processed. To obtain the list of fields in file (and their\
                   handgrip), use EPyGrAM's 'grib_what' tool or GRIB_API's\
                   'grib_dump'. Ex: -f 'shortName:t,level:850'.",
             default=None)],
    '2fields_mode':[
        '-y', '--2fields_mode',
        dict(help="mode for 2-fields joining:\
                   - 'vectorize': sets a vector with each field as a component\
                   - 'superpose': superpose second field over the first with\
                     contourlines\
                   - '-': computes and plots the difference field\
                   - '+': computes and plots the sum field\
                   - '*': computes and plots the product field\
                   - '/': computes and plots the division field",
             default=None)],
    'list_of_fields':[
        '-l', '--listoffields',
        dict(help="name of an external file containing the list of fields\
                   to be processed, with the same specifications as for -f\
                   argument.",
             default=None)],
    'reverse_fields_selection':[
        '-r', '--reverse',
        dict(action='store_true',
             help="reverse fields selection: all but those requested.",
             default=False)],
    'sort_fields':[
        '-s', '--sortfields',
        dict(action='store_true',
             help="not for GRIB: sort fields with regards to\
                   their name and type.",
             default=False)],
    'GRIB_sort':[
        '-k', '--sortbykey',
        dict(help="GRIB only: sort fields with regards to the given fid key.\
                   Only for *mode* == 'one+list' or 'fid_list'.",
             dest='sortfields',
             default=False)],
    'GRIB_what_mode':[
        '-m', '--mode',
        dict(help="information display mode (GRIB only): 'one+list' gives the \
                   validity/geometry of the first field in GRIB, plus \
                   the list of fid; 'fid_list' gives only the fid of \
                   each field in GRIB; 'what' gives the values of the \
                   keys from each GRIB message that are used to generate \
                   an **epygram** field from the message (slower); \
                   'ls' gives the values of the 'ls' keys from each GRIB \
                   message; 'mars' gives the values of the 'mars' keys \
                   from each GRIB message.",
             choices=('one+list', 'fid_list', 'what', 'ls', 'mars'),
             default='one+list')],
    'FA_set_compression':[
        '-k', '--kompression',
        dict(help="bits number for FA gridpoint compression : KNBITS/KNBPDG.\
                   Defaults to " +
                   str(_defaults.get('FA_default_compression',
                                     {'KNBPDG':('undefined')}).get('KNBPDG', 'undefined')) + ".",
             type=int,
             default=None)],
    'GRIB2_packing':[
        '--GRIB2_packing',
        dict(help="packing for writing GRIB2s.\
                   Defaults to " +
                   str(griberies.defaults.GRIB2_keyvalue[5]) + ".",
             default=None)],
    'netCDF_compression':[
        '--nc_comp',
        dict(type=int,
             help="compression level to compress fields in netCDF. \
                   Ranges from 0 to 9. 0 is no compression, 1 is low-but-fast \
                   compression, 9 is high-but-slow compression. \
                   Default to " + str(_defaults.get('netCDF_default_compression',
                                                    "4")),
             choices=list(range(0, 10)),
             default=str(_defaults.get('netCDF_default_compression', 4)))],
                    }

