#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division

import six
import argparse

import epygram
from epygram import epylog, epygramError
from epygram.util import printstatus
from epygram.args_catalog import add_arg_to_parser, \
                                 files_management, fields_management, \
                                 runtime_options


def main(filename, refname, fieldseed,
         progressmode=None,
         operation=False,
         in_place=False):
    """
    Args:
        filename: name of the file to be processed.
        refname: name of the source file from which to add fields.
        fieldseed: either a fid or a list of fid, used as a seed for
                   generating the list of fields to be processed.
        progressmode: among ('verbose', 'percentage', None).
        operation: among ('diff', 'reversediff', 'add', 'multiply') for
                   operation to be done between fields of file and source.
        in_place: if True, the field(s) is(are) replaced "in place",
                  not in a new file.
    """

    target = epygram.formats.resource(filename, openmode='a')
    if target.format not in ('GRIB', 'FA', 'LFI'):
        epylog.warning(" ".join(["tool NOT TESTED with format",
                                 target.format, "!"]))
    source = epygram.formats.resource(refname, openmode='r')
    fidlist = source.find_fields_in_resource(seed=fieldseed)

    def treat_one_field(f, s, nt, operation=False):
        """
        Unitary task to do on one field.
        f = fid
        s = source,
        nt = newtarget
        """
        # set reading options
        if s.format == 'GRIB':
            options = [k[0] for k in epygram.config.GRIB_default_packing[1]] + \
                      [k[0] for k in epygram.config.GRIB_default_packing[2]]
            options.extend(list(epygram.config.GRIB_default_ordering.keys()))
            options.append('centre')
            options = {'get_info_as_json':options}
        else:
            options = {}
        # read source field
        field = s.readfield(f, **options)
        # read target field and perform operation
        if operation and f in target.listfields():
            originfield = target.readfield(f)
            if operation == 'diff':
                field.operation('*', -1.)
                field.operation('+', originfield)
            elif operation == 'reversediff':
                field.operation('-', originfield)  # filename - source
            elif operation == 'add':
                field.operation('+', originfield)
            elif operation == 'multiply':
                field.operation('*', originfield)
        # set writing options
        if s.format == 'FA':
            options = {'compression':source.fieldscompression.get(f, None)}
        elif s.format == 'GRIB':
            grib_edition = field.fid.get('GRIB1', field.fid)['editionNumber']
            import json
            options = json.loads(field.comment)
            for k, v in options.items():
                if isinstance(k, unicode):
                    del options[k]
                    k = str(k)
                if isinstance(v, unicode):
                    options[k] = str(v)
                else:
                    options[k] = v
            options = {'sample':'file:' + s.container.abspath,
                       'packing':{k:options.pop(k) for k in
                                  [k for k in epygram.config.GRIB_default_packing[grib_edition]]
                                  if k in options},
                       'ordering':{k:options.pop(k) for k in
                                   epygram.config.GRIB_default_ordering.keys()},
                       'other_GRIB_options':options
                       }
        # write field
        nt.writefield(field, **options)

    # process fields
    if in_place:
        numfields = len(fidlist)
        n = 1
        for f in fidlist:
            if progressmode == 'verbose':
                epylog.info(f)
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if source.format == 'GRIB' and f in target.listfields():
                raise epygramError('field ' + str(f) + ' is already in' +
                                   ' target resource: unable to modify GRIB' +
                                   ' messages "in place".')
            treat_one_field(f, source, target, operation=operation)
    else:
        # open new target
        tobecopied = {'fmt':target.format}
        if source.format == 'FA':
            tobecopied.update({'headername':source.headername,
                               'validity':source.validity,
                               'cdiden':source.cdiden,
                               'default_compression':source.default_compression,
                               'processtype':source.processtype})
        elif source.format == 'LFI':
            tobecopied.update({'compressed':source.compressed})
        newtarget = epygram.formats.resource(filename + '.movefield.out',
                                             openmode='w',
                                             **tobecopied)
        targetfields = target.listfields()
        to_add_fields = [f for f in fidlist if f not in targetfields]
        numfields = len(targetfields) + len(to_add_fields)
        n = 1
        for f in targetfields:
            if progressmode == 'verbose':
                epylog.info(f)
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if f in fidlist:
                treat_one_field(f, source, newtarget, operation=operation)  # replace
            else:
                treat_one_field(f, target, newtarget)  # copy (untouched)
        for f in to_add_fields:
            treat_one_field(f, source, newtarget)  # add
# end of main() ###############################################################


if __name__ == '__main__':

    # ## 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description="An EPyGrAM tool for moving field(s) from a resource to another\
                                                  (fields are not removed from the source file).",
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    add_arg_to_parser(parser, files_management['source_file'])
    add_arg_to_parser(parser, files_management['in_place'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['field'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    diff = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diff, files_management['replace_by_diff'])
    add_arg_to_parser(diff, files_management['replace_by_reversediff'])
    add_arg_to_parser(diff, files_management['replace_by_addition'])
    add_arg_to_parser(diff, files_management['replace_by_product'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_options['verbose'])
    add_arg_to_parser(status, runtime_options['percentage'])

    args = parser.parse_args()

    # ## 2. Initializations
    ######################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.verbose:
        progressmode = 'verbose'
    elif args.percentage:
        progressmode = 'percentage'
    else:
        progressmode = None

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

    # ## 3. Main
    ###########
    main(six.u(args.filename),
         six.u(args.refname),
         fieldseed,
         progressmode=progressmode,
         operation=args.replace_op,
         in_place=args.in_place)

###########
### END ###
###########
