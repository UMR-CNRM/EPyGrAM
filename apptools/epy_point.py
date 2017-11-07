#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import argparse
import sys

from footprints import FPDict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog
from epygram.util import write_formatted_table
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  misc_options, runtime_options,
                                  extraction_options, output_options)


def main(filename,
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
        epylog.warning(" ".join(["tool NOT TESTED with format",
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
            epylog.info(f)
            field = resource.readfield(f)
            if field.spectral:
                field.sp2gp()
            if operation is not None:
                field.operation(**operation)
            if pressure_unit_hpa and \
               (field.fid['generic'].get('discipline') == 0 and
                field.fid['generic'].get('parameterCategory') == 3 and
                field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
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
            epylog.info(f)
            if f in fidlist:
                field = resource.readfield(f)
                if field.spectral:
                    field.sp2gp()
                if operation is not None:
                    field.operation(**operation)
                if pressure_unit_hpa and \
                   (field.fid['generic'].get('discipline') == 0 and
                    field.fid['generic'].get('parameterCategory') == 3 and
                    field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
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
                if pressure_unit_hpa and \
                   (field.fid['generic'].get('discipline') == 0 and
                    field.fid['generic'].get('parameterCategory') == 3 and
                    field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
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
                filename = resource.container.absdir + \
                           '.'.join(['diff',
                                     resource.container.basename + '-' + reference.container.basename,
                                     parameter,
                                     str(coordinates[0]) + "E" + str(coordinates[1]) + "N",
                                     suffix])
            else:
                filename = resource.container.absdir + \
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
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description="An EPyGrAM tool for extracting from a resource the value of meteorological fields\
                                                  on a geographical point.",
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['field'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, extraction_options['point_coordinates'])
    add_arg_to_parser(parser, extraction_options['horizontal_interpolation'])
    add_arg_to_parser(parser, extraction_options['external_distance'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, misc_options['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_options['operation_on_field'])
    add_arg_to_parser(parser, output_options['precision'])
    add_arg_to_parser(parser, output_options['lonlat_precision'])
    add_arg_to_parser(parser, output_options['GRIB_short_fid'])
    add_arg_to_parser(parser, output_options['stdout'])
    add_arg_to_parser(parser, output_options['outputfilename'])
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    # 2. Initializations
    ####################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.Drefname is not None:
        refname = args.Drefname
        diffonly = True
    else:
        refname = args.refname
        diffonly = False
    coordinates = args.coordinates.split(',')
    coordinates = tuple([float(i) for i in coordinates])
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
    else:
        external_distance = False
    if args.operation is not None:
        _operation = args.operation.split(',')
        operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            operation['operand'] = float(_operation.pop(0).strip())
    else:
        operation = None

    # 2.2 list of fields to be processed
    if args.field is not None:
        fieldseed = args.field
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            fieldseed = l.readlines()
        for n in range(len(fieldseed)):
            fieldseed[n] = fieldseed[n].replace('\n', '').strip()
    else:
        fieldseed = None

    # 3. Main
    #########
    main(six.u(args.filename),
         fieldseed,
         coordinates,
         refname=refname,
         diffonly=diffonly,
         pressure_unit_hpa=args.pressure_unit_hpa,
         operation=operation,
         interpolation=args.interpolation,
         precision=args.precision,
         lonlat_precision=args.lonlat_precision,
         stdoutput=args.stdout,
         grib_short_fid=args.grib_short_fid,
         external_distance=external_distance,
         outputfilename=args.outputfilename)
