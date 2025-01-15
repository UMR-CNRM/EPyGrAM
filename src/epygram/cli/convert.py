#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for converting file formats. Spectral fields are converted into gridpoints."""

import argparse

import footprints
from bronx.syntax.parsing import str2dict

import epygram
from epygram import epylog as logger
from epygram.extra import griberies
from epygram.formats.conversion.functions import batch_convert
from epygram.formats.conversion import converters
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           runtime_args,
                           operational_args,
                           output_args)
from . import epilog

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    batch_convert(args.filenames,
                  args.output_format,
                  # technical
                  threads_number=args.threads_number,
                  progressmode=args.progressmode,
                  # instructions common to formats
                  fieldseed=args.fieldseed,
                  subzone=args.subzone,
                  grib_short_fid=args.grib_short_fid,
                  # GRIB specifics
                  suite=args.suite,
                  typeofgeneratingprocess=args.typeOfGeneratingProcess,
                  default_packing=args.default_packing,
                  other_grib_options=args.other_GRIB_options,
                  # netCDF specifics
                  compression=args.nc_comp,
                  flatten_horizontal_grids=args.flatten,
                  # GeoPoints specifics
                  order=args.order,
                  lonlat_precision=args.lonlat_precision,
                  precision=args.precision,
                  llv=args.llv,
                  columns=args.geopoints_cols,
                  no_header=args.geopoints_noheader
                  )


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)

    add_arg_to_parser(parser, files_args['several_files'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_args['field'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    add_arg_to_parser(parser, misc_args['LAMzone'])
    add_arg_to_parser(parser, output_args['output_format'])
    # compression/precision options
    add_arg_to_parser(parser, fields_args['netCDF_compression'])
    add_arg_to_parser(parser, fields_args['GRIB2_packing'])
    add_arg_to_parser(parser, output_args['GeoPoints_lonlat_precision'])
    add_arg_to_parser(parser, output_args['GeoPoints_precision'])
    # GRIB specifics
    add_arg_to_parser(parser, operational_args['suite'])
    add_arg_to_parser(parser, operational_args['typeOfGeneratingProcess'])
    add_arg_to_parser(parser, operational_args['numod'])
    add_arg_to_parser(parser, output_args['GRIB_other_options'])
    # GeoPoints specific
    add_arg_to_parser(parser, misc_args['array_flattening_order'])
    cols = parser.add_mutually_exclusive_group()
    add_arg_to_parser(cols, output_args['GeoPoints_llv'])
    add_arg_to_parser(cols, output_args['GeoPoints_columns'])
    add_arg_to_parser(parser, output_args['GeoPoints_noheader'])
    # netCDF
    add_arg_to_parser(parser, misc_args['flatten_horizontal_grids'])
    # others
    add_arg_to_parser(parser, output_args['GRIB_short_fid'])
    add_arg_to_parser(parser, runtime_args['threads_number'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_args['verbose'])
    add_arg_to_parser(status, runtime_args['percentage'])

    args = parser.parse_args()

    # 2. Initializations
    ####################

    # 2.1 options
    if args.zone in ('C', 'CI'):
        args.subzone = args.zone
    else:
        args.subzone = None
    if args.GRIB2_packing is not None:
        args.default_packing = str2dict(args.GRIB2_packing, try_convert=int)
    else:
        args.default_packing = griberies.defaults.GRIB2_keyvalue[5]
    if args.GRIB_other_options is not None:
        args.other_GRIB_options = str2dict(args.GRIB_other_options, try_convert=int)
    else:
        args.other_GRIB_options = {}
    if args.numod is not None:
        args.other_GRIB_options['generatingProcessIdentifier'] = args.numod
    if args.geopoints_cols is not None:
        args.args.geopoints_cols = [c.strip() for c in args.geopoints_cols.split(',')]
    assert args.filenames != [], \
           "must supply one or several filenames."
    args.threads_number = min(args.threads_number, len(args.filenames))
    args.progressmode = None
    if args.verbose:
        args.progressmode = 'verbose'
    elif args.percentage:
        if threads_number > 1:
            args.progressmode = None
        else:
            args.progressmode = 'percentage'

    # 2.2 list of fields to be processed
    if args.field is not None:
        args.fieldseed = [args.field]
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as listfile:
            args.fieldseed = [line.replace('\n', '').strip() for line in listfile.readlines()]
    else:
        args.fieldseed = None
    return args

