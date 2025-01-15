#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for extracting from a resource the value of meteorological fields on a geographical point."""

import argparse
import sys
import os

from footprints import FPDict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog as logger
from epygram.util import write_formatted_table
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           runtime_args,
                           extraction_args,
                           output_args)

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    point(
         args.filename,
         args.fieldseed,
         args.coordinates,
         refname=args.refname,
         diffonly=args.diffonly,
         pressure_unit_hpa=args.pressure_unit_hpa,
         operation=args.operation,
         interpolation=args.interpolation,
         precision=args.precision,
         lonlat_precision=args.lonlat_precision,
         stdoutput=args.stdout,
         grib_short_fid=args.grib_short_fid,
         external_distance=args.external_distance,
         outputfilename=args.outputfilename)


def point(
         filename,
         fieldseed,
         coordinates,
         refname=None,
         diffonly=False,
         pressure_unit_hpa=False,
         operation=None,
         interpolation='nearest',
         precision=epygram.config.GeoPoints_precision,
         lonlat_precision=epygram.config.GeoPoints_lonlat_precision,
         stdoutput=False,
         grib_short_fid=False,
         external_distance=None,
         outputfilename=None):
    """
    Extract point.

    :param filename: name of the file to be processed.
    :param fieldseed: either a string or a list of strings, used as a seed for
                   generating the list of fields to be processed.
    :param coordinates: (lon, lat) coordinates (in °) of the point to be extracted.
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only prints the difference.
    :param pressure_unit_hpa: converts FA log(pressure) fields to hPa.
    :param operation: makes the requested operation
                   (e.g. {'operation':'-','scalar':273.15} or
                   {'operation':'exp'}) on the field.
    :param interpolation: kind of interpolation from grid to coordinates, among
                       ('nearest', 'linear', 'cubic').
    :param precision: number of digits for output fields values.
    :param lonlat_precision: number of digits for output lon/lat values.
    :param stdoutput: if True, output is redirected to stdout.
    :param grib_short_fid: if True, condense GRIB fid.
    :param external_distance: can be a dict containing the target point value
                           and an additional field fid, to which field the distance
                           for computing nearest point
                           is computed within the 4 horizontally nearest points; e.g.
                           {'target_value':4810, 'external_field':an_additional_fid}.
                           If so, the nearest point is selected with
                           distance = |target_value - external_field.data|
    :param outputfilename: specify an output filename.
    """
    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if external_distance:
        external_distance['external_field'] = resource.readfield(external_distance['external_field'])
        if external_distance['external_field'].spectral:
            external_distance['external_field'].sp2gp()

    if not diffmode:
        values = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(d) for d in fidlist]
        firstfield = True
        for f in fidlist:
            logger.info(f)
            field = resource.readfield(f)
            if field.spectral:
                field.sp2gp()
            if operation is not None:
                field.operation(**operation)
            if (pressure_unit_hpa and
                (field.fid['generic'].get('discipline') == 0 and
                 field.fid['generic'].get('parameterCategory') == 3 and
                 field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                field.scalar_operation('/', 100.)
            if firstfield and interpolation == 'nearest':
                (value, neighborinfo) = field.getvalue_ll(*coordinates,
                                                          interpolation=interpolation,
                                                          neighborinfo=True,
                                                          one=True,
                                                          external_distance=external_distance)
                firstfield = False
            else:
                value = field.getvalue_ll(*coordinates,
                                          interpolation=interpolation,
                                          one=True,
                                          external_distance=external_distance)
            values[f] = value
    else:
        values = {}
        refvalues = {}
        diffvalues = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        reffidlist = reference.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(d) for d in fidlist]
            reffidlist = [FPDict(d) for d in reffidlist]
        unionfidlist = list(set(fidlist).union(set(reffidlist)))
        intersectionfidlist = list(set(fidlist).intersection(set(reffidlist)))
        unionfidlist.sort()
        firstfield = True
        for f in unionfidlist:
            logger.info(f)
            if f in fidlist:
                field = resource.readfield(f)
                if field.spectral:
                    field.sp2gp()
                if operation is not None:
                    field.operation(**operation)
                if (pressure_unit_hpa and
                    (field.fid['generic'].get('discipline') == 0 and
                     field.fid['generic'].get('parameterCategory') == 3 and
                     field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                    field.scalar_operation('/', 100.)
                if firstfield and interpolation == 'nearest':
                    (value, neighborinfo) = field.getvalue_ll(*coordinates,
                                                              interpolation=interpolation,
                                                              neighborinfo=True,
                                                              one=True,
                                                              external_distance=external_distance)
                    firstfield = False
                else:
                    value = field.getvalue_ll(*coordinates,
                                              interpolation=interpolation,
                                              one=True,
                                              external_distance=external_distance)
                values[f] = value
            if f in reffidlist:
                field = reference.readfield(f)
                if field.spectral:
                    field.sp2gp()
                if operation is not None:
                    field.operation(**operation)
                if (pressure_unit_hpa and
                    (field.fid['generic'].get('discipline') == 0 and
                     field.fid['generic'].get('parameterCategory') == 3 and
                     field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                    field.scalar_operation('/', 100.)
                if firstfield and interpolation == 'nearest':
                    (value, neighborinfo) = field.getvalue_ll(*coordinates,
                                                              interpolation=interpolation,
                                                              neighborinfo=True,
                                                              one=True,
                                                              external_distance=external_distance)
                    firstfield = False
                else:
                    value = field.getvalue_ll(*coordinates,
                                              interpolation=interpolation,
                                              one=True,
                                              external_distance=external_distance)
                refvalues[f] = value
            if f in intersectionfidlist:
                diffvalues[f] = values[f] - refvalues[f]

    # Output
    coordstr = "(" + '{:.{precision}{type}}'.format(coordinates[0], type='F', precision=lonlat_precision) + ", " + \
                     '{:.{precision}{type}}'.format(coordinates[1], type='F', precision=lonlat_precision) + ")"
    if interpolation == 'nearest':
        gridpointstr = "(" + '{:.{precision}{type}}'.format(neighborinfo[0], type='F', precision=lonlat_precision) + ", " + \
                             '{:.{precision}{type}}'.format(neighborinfo[1], type='F', precision=lonlat_precision) + ")"
    suffix = "point.out"
    if not diffmode:
        fidlist.sort()
        noutfields = len(fidlist)
        if outputfilename:
            filename = outputfilename
        else:
            if noutfields == 1:
                parameter = smooth_string(fidlist[0])
                filename = '.'.join([resource.container.abspath,
                                     parameter,
                                     str(coordinates[0]) + "E" + str(coordinates[1]) + "N",
                                     suffix])
            else:
                filename = '.'.join([resource.container.abspath,
                                     str(noutfields) + "fields",
                                     str(coordinates[0]) + "E" + str(coordinates[1]) + "N",
                                     suffix])
        output = [['# field', 'value']]
        for f in fidlist:
            if resource.format == 'GRIB' and grib_short_fid:
                fid = f.get('GRIB1', f)
                fid = '/'.join([str(fid[k]) for k in sorted(fid.keys())])
            else:
                fid = str(f).strip()
            output.append([fid, values[f]])
    else:
        noutfields = len(unionfidlist)
        if outputfilename:
            filename = outputfilename
        else:
            if noutfields == 1:
                parameter = smooth_string(unionfidlist[0])
                filename = resource.container.absdir + os.path.sep + \
                    '.'.join(['diff',
                              resource.container.basename + '-' + reference.container.basename,
                              parameter,
                              str(coordinates[0]) + "E" + str(coordinates[1]) + "N",
                              suffix])
            else:
                filename = resource.container.absdir + os.path.sep + \
                    '.'.join(['diff',
                              resource.container.basename + '-' + reference.container.basename,
                              str(noutfields) + "fields",
                              str(coordinates[0]) + "E" + str(coordinates[1]) + "N",
                              suffix])
        if diffonly:
            output = [['# field', resource.container.basename + '-' + reference.container.basename]]
        else:
            output = [['# field', resource.container.basename, reference.container.basename, '-diff-']]
        for f in unionfidlist:
            if resource.format == 'GRIB' and grib_short_fid:
                fid = f.get('GRIB1', f)
                fid = '/'.join([str(fid[k]) for k in sorted(fid.keys())])
            else:
                fid = str(f).strip()
            if diffonly:
                output.append([fid, diffvalues[f]])
            else:
                output.append([fid, values[f], refvalues[f], diffvalues[f]])

    if stdoutput:
        out = sys.stdout
    else:
        out = open(filename, 'w')
    if interpolation == 'nearest':
        head = "# asked coord: " + coordstr + " / gridpoint used: " + gridpointstr
    else:
        head = "# asked coord: " + coordstr + " / " + interpolation + " interpolation"
    out.write(head + '\n')
    out.write('-' * len(head) + '\n')
    write_formatted_table(out, output, precision=precision)


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_args['field'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    add_arg_to_parser(parser, extraction_args['point_coordinates'])
    add_arg_to_parser(parser, extraction_args['horizontal_interpolation'])
    add_arg_to_parser(parser, extraction_args['external_distance'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, misc_args['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_args['operation_on_field'])
    add_arg_to_parser(parser, output_args['precision'])
    add_arg_to_parser(parser, output_args['lonlat_precision'])
    add_arg_to_parser(parser, output_args['GRIB_short_fid'])
    add_arg_to_parser(parser, output_args['stdout'])
    add_arg_to_parser(parser, output_args['outputfilename'])
    add_arg_to_parser(parser, runtime_args['verbose'])
    args = parser.parse_args()

    # 2. Initializations
    ####################
    # 2.1 options
    if args.Drefname is not None:
        args.refname = args.Drefname
        args.diffonly = True
    else:
        args.refname = args.refname
        args.diffonly = False
    args.coordinates = args.coordinates.split(',')
    args.coordinates = tuple([float(i) for i in args.coordinates])
    if args.external_distance is not None:
        external_distance = args.external_distance.split(';')
        if ':' in external_distance[1] or '=' in external_distance[1]:  # GRIB
            external_distance[1] = {i.replace('=', ':').split(':')[0].strip():i.replace('=', ':').split(':')[1] for i in args.field.split(',')}
            for k, v in external_distance[1].items():
                try:
                    external_distance[1][k] = int(v)
                except ValueError:
                    pass
        else:
            external_distance[1] = external_distance[1]
        external_distance = {'target_value':float(external_distance[0]),
                             'external_field':external_distance[1]}
        args.external_distance = external_distance
    if args.operation is not None:
        _operation = args.operation.split(',')
        args.operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            args.operation['operand'] = float(_operation.pop(0).strip())
    # 2.2 list of fields to be processed
    if args.field is not None:
        args.fieldseed = args.field
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as lf:
            args.fieldseed = [line.replace('\n', '').strip() for line in lf.readlines()]
    else:
        args.fieldseed = None

    return args

