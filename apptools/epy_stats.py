#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import argparse
import sys
import os

from footprints import FPDict

import epygram
from epygram import epylog
from epygram.util import printstatus, write_formatted_table
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  misc_options, runtime_options,
                                  output_options)


def main(filename,
         fieldseed,
         refname=None,
         diffonly=False,
         only_maxdiff=False,
         subzone=None,
         pressure_unit_hpa=False,
         operation=None,
         diffoperation=None,
         precision=4,
         progressmode=None,
         stdoutput=False,
         outputfilename=None,
         grib_short_fid=False):
    """
    Computes stats over fields.

    :param filename: name of the file to be processed.
    :param fieldseed: either a fid or a list of fid, used as a seed for
                   generating the list of fields to be processed.
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only prints the difference statistics.
    :param only_maxdiff: if diffonly and True, only print the maximum absolute
                      difference (and min/max values of fields for magnitude
    :param subzone: LAM zone among ('C', 'CI', 'CIE').
    :param pressure_unit_hpa: for FA, converts log(pressure) fields to hPa.
    :param operation: makes the requested operation
                   (e.g. {'operation':'-','scalar':273.15} or
                   {'operation':'exp'}) on the field.
    :param diffoperation: makes the requested operation
                       (e.g. {'operation':'-','scalar':273.15} or
                       {'operation':'exp'}) on the difference field.
    :param precision: number of decimals for writing floats.
    :param progressmode: among ('verbose', 'percentage', None).
    :param stdoutput: if True, output is redirected to stdout.
    :param outputfilename: specify an output filename.
    :param grib_short_fid: if True, condense GRIB fid.
    """
    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI', 'TIFFMF'):
        epylog.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if not diffmode:
        stats = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(d) for d in fidlist]
        numfields = len(fidlist)
        n = 1
        for f in fidlist:
            if progressmode == 'verbose':
                epylog.info(str(f))
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            field = resource.readfield(f)
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            if field.spectral:
                field.sp2gp()
            if operation is not None:
                field.operation(**operation)
            if (pressure_unit_hpa and
                (field.fid['generic'].get('discipline') == 0 and
                 field.fid['generic'].get('parameterCategory') == 3 and
                 field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                field.scalar_operation('/', 100.)
            stats[f] = field.stats(subzone=subzone)
    else:
        stats = {}
        refstats = {}
        diffstats = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        reffidlist = reference.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(d) for d in fidlist]
            reffidlist = [FPDict(d) for d in reffidlist]
        unionfidlist = list(set(fidlist).union(set(reffidlist)))
        intersectionfidlist = list(set(fidlist).intersection(set(reffidlist)))
        unionfidlist.sort()
        numfields = len(unionfidlist)
        n = 1
        for f in unionfidlist:
            if progressmode == 'verbose':
                epylog.info(f)
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if f in fidlist:
                field = resource.readfield(f)
                if not field.geometry.grid.get('LAMzone', False):
                    subzone = None
                if field.spectral:
                    field.sp2gp()
                if operation is not None:
                    field.operation(**operation)
                if (pressure_unit_hpa and
                    (field.fid['generic'].get('discipline') == 0 and
                     field.fid['generic'].get('parameterCategory') == 3 and
                     field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                    field.scalar_operation('/', 100.)
                stats[f] = field.stats(subzone=subzone)
            if f in reffidlist:
                reffield = reference.readfield(f)
                if not reffield.geometry.grid.get('LAMzone', False):
                    subzone = None
                if reffield.spectral:
                    reffield.sp2gp()
                if operation is not None:
                    reffield.operation(**operation)
                if (pressure_unit_hpa and
                    (reffield.fid['generic'].get('discipline') == 0 and
                     reffield.fid['generic'].get('parameterCategory') == 3 and
                     reffield.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25))):
                    reffield.scalar_operation('/', 100.)
                refstats[f] = reffield.stats(subzone=subzone)
            if f in intersectionfidlist:
                diff = field - reffield
                if diffoperation is not None:
                    field.operation(**diffoperation)
                diffstats[f] = diff.stats(subzone=subzone)

    # Output
    single_params = ['min', 'max', 'mean', 'std', 'quadmean', 'nonzero']
    diffonly_params = ['min', 'max', 'bias', 'rmsd', 'nonzero']
    diff_params = ['min', 'max', 'mean', 'std', 'bias', 'rmsd', 'nonzero']
    singlemap = {'min':'min', 'max':'max', 'mean':'mean', 'std':'std',
                 'bias':None, 'rmsd':None,
                 'quadmean':'quadmean', 'nonzero':'nonzero'}
    diffmap = {'min':'min', 'max':'max', 'mean':None, 'std':None,
               'bias':'mean', 'rmsd':'quadmean', 'nonzero':'nonzero'}
    suffix = "stats.out"
    if not diffmode:
        printlist = fidlist
        printlist.sort()
        noutfields = len(printlist)
        parameter = epygram.util.soft_string(printlist[0])
        if outputfilename:
            filename = outputfilename
        else:
            if noutfields == 1:
                filename = '.'.join([resource.container.abspath,
                                     parameter,
                                     suffix])
            else:
                filename = '.'.join([resource.container.abspath,
                                     str(noutfields) + "fields",
                                     suffix])
        head = resource.container.basename
    else:
        if diffonly:
            printlist = intersectionfidlist
        else:
            printlist = unionfidlist
        printlist.sort()
        noutfields = len(printlist)
        parameter = epygram.util.soft_string(printlist[0])
        if outputfilename:
            filename = outputfilename
        else:
            if noutfields == 1:
                filename = os.path.join(resource.container.absdir,
                                        '.'.join(['diff',
                                                  resource.container.basename + '-' + reference.container.basename,
                                                  parameter,
                                                  suffix]))
            else:
                filename = os.path.join(resource.container.absdir,
                                        '.'.join(['diff',
                                                  resource.container.basename + '-' + reference.container.basename,
                                                  str(noutfields) + "fields",
                                                  suffix]))
        head = "-diff-"
    if stdoutput:
        out = sys.stdout
    else:
        out = open(filename, 'w')
    output = ["# " + head]
    if not diffmode:
        output.extend(single_params)
        output = [output]
        for f in printlist:
            if resource.format == 'GRIB' and grib_short_fid:
                fid = f.get('GRIB1', f)
                fid = '/'.join([str(fid[k]) for k in sorted(fid.keys())])
            else:
                fid = str(f).strip()
            output.append([fid] + [stats[f].get(singlemap[k], None) for k in single_params])
        write_formatted_table(out, output, precision=precision)
    elif diffonly:
        if only_maxdiff:
            output.extend(['min_flds', 'max_flds', 'abs_max_diff', 'magnitude'])
        else:
            output.extend(diffonly_params)
        output = [output]
        for f in printlist:
            if resource.format == 'GRIB' and grib_short_fid:
                fid = f.get('GRIB1', f)
                fid = '/'.join([str(fid[k]) for k in sorted(fid.keys())])
            else:
                fid = str(f).strip()
            if only_maxdiff:
                mini = min(stats[f].get('min'),
                           refstats[f].get('min'))
                maxi = max(stats[f].get('max'),
                           refstats[f].get('max'))
                diffmax = max(abs(diffstats[f].get('max')),
                              abs(diffstats[f].get('min')))
                magnitude = (diffmax /  # diffmax / max value
                             max([abs(mini), abs(maxi), epygram.config.epsilon]))
                magnitude = '{:.1f}%'.format(magnitude * 100.)
                output.append([fid] + [mini, maxi, diffmax, magnitude])
            else:
                output.append([fid] + [diffstats[f].get(diffmap[k], None)
                                       for k in diffonly_params])
        write_formatted_table(out, output, precision=precision,
                              alignments=['<', '>'] if only_maxdiff else ['<', '^'])
    else:
        for f in printlist:
            if resource.format == 'GRIB' and grib_short_fid:
                fid = f.get('GRIB1', f)
                fid = '/'.join([str(fid[k]) for k in sorted(fid.keys())])
            else:
                fid = str(f).strip()
            output = ['# ' + fid]
            output.extend(diff_params)
            output = [output]
            if f in stats:
                output.append([resource.container.basename] + [stats[f].get(singlemap[k], None) for k in diff_params])
            else:
                output.append([resource.container.basename] + ['-' for k in diff_params])
            if f in refstats:
                output.append([reference.container.basename] + [refstats[f].get(singlemap[k], None) for k in diff_params])
            else:
                output.append([reference.container.basename] + ['-' for k in diff_params])
            if f in diffstats:
                output.append([head] + [diffstats[f].get(diffmap[k], None) for k in diff_params])
            else:
                output.append([head] + ['-' for k in diff_params])
            write_formatted_table(out, output, precision=precision)
    out.close()
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description='An EPyGrAM tool for computing basic statistics on \
                                                  meteorological fields.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['FA_multiple_fields'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_options['only_maxdiff'])
    add_arg_to_parser(parser, misc_options['LAMzone'])
    add_arg_to_parser(parser, misc_options['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_options['operation_on_field'])
    add_arg_to_parser(parser, misc_options['diffoperation_on_field'])
    add_arg_to_parser(parser, output_options['precision'])
    add_arg_to_parser(parser, output_options['GRIB_short_fid'])
    add_arg_to_parser(parser, output_options['stdout'])
    add_arg_to_parser(parser, output_options['outputfilename'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_options['verbose'])
    add_arg_to_parser(status, runtime_options['percentage'])

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
    if args.zone in ('C', 'CI'):
        subzone = args.zone
    else:
        subzone = None
    if args.verbose:
        progressmode = 'verbose'
    elif args.percentage:
        progressmode = 'percentage'
    else:
        progressmode = None
    if args.operation is not None:
        _operation = args.operation.split(',')
        operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            operation['operand'] = float(_operation.pop(0).strip())
    else:
        operation = None
    if args.diffoperation is not None:
        _diffoperation = args.diffoperation.split(',')
        diffoperation = {'operation':_diffoperation.pop(0).strip()}
        if len(_diffoperation) > 0:
            diffoperation['operand'] = float(_diffoperation.pop(0).strip())
    else:
        diffoperation = None

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
         refname=refname,
         diffonly=diffonly,
         only_maxdiff=args.only_maxdiff,
         subzone=subzone,
         pressure_unit_hpa=args.pressure_unit_hpa,
         operation=operation,
         diffoperation=diffoperation,
         precision=args.precision,
         progressmode=progressmode,
         stdoutput=args.stdout,
         outputfilename=args.outputfilename,
         grib_short_fid=args.grib_short_fid)
