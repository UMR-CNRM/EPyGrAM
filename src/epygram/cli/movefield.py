#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for moving field(s) from a resource to another (fields are not removed from the source file)."""

import argparse

from bronx.fancies.display import printstatus

import epygram
from epygram import epylog as logger
from epygram import epygramError
from epygram.extra import griberies
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           runtime_args)

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    movefield(args.filename,
              args.refname,
              args.fieldseed,
              progressmode=args.progressmode,
              operation=args.replace_op,
              in_place=args.in_place)


def movefield(filename,
              refname,
              fieldseed,
              progressmode=None,
              operation=False,
              in_place=False):
    """
    Move fields.

    :param filename: name of the file to be processed.
    :param refname: name of the source file from which to add fields.
    :param fieldseed: either a fid or a list of fid, used as a seed for
                      generating the list of fields to be processed.
    :param progressmode: among ('verbose', 'percentage', None).
    :param operation: among ('diff', 'reversediff', 'add', 'multiply') for
                      operation to be done between fields of file and source.
    :param in_place: if True, the field(s) is(are) replaced "in place",
                     not in a new file.
    """
    target = epygram.formats.resource(filename, openmode='a')
    if target.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
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
            read_misc_metadata = set(griberies.defaults.GRIB1_packing.keys())
            read_misc_metadata.update(set(griberies.defaults.GRIB2_keyvalue[5].keys()))
            read_misc_metadata.update(set(griberies.defaults.GRIB1_ordering.keys()))
            read_misc_metadata.add('centre')
            options = {'read_misc_metadata':list(read_misc_metadata)}
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
            options = {'sample':'file:' + s.container.abspath}
        # write field
        nt.writefield(field, **options)

    # process fields
    if in_place:
        numfields = len(fidlist)
        n = 1
        for f in fidlist:
            if progressmode == 'verbose':
                logger.info(f)
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
            tobecopied.update({'headername':target.headername,
                               'validity':target.validity,
                               'cdiden':target.cdiden,
                               'default_compression':target.default_compression,
                               'processtype':target.processtype})
        elif source.format == 'LFI':
            tobecopied.update({'compressed':target.compressed})
        newtarget = epygram.formats.resource(filename + '.movefield.out',
                                             openmode='w',
                                             **tobecopied)
        targetfields = target.listfields()
        to_add_fields = [f for f in fidlist if f not in targetfields]
        numfields = len(targetfields) + len(to_add_fields)
        n = 1
        for f in targetfields:
            if progressmode == 'verbose':
                logger.info(f)
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if f in fidlist:
                treat_one_field(f, source, newtarget, operation=operation)  # replace
            else:
                treat_one_field(f, target, newtarget)  # copy (untouched)
        for f in to_add_fields:
            treat_one_field(f, source, newtarget)  # add


def get_args():

    # ## 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, files_args['source_file'])
    add_arg_to_parser(parser, files_args['in_place'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_args['field'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    diff = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diff, files_args['replace_by_diff'])
    add_arg_to_parser(diff, files_args['replace_by_reversediff'])
    add_arg_to_parser(diff, files_args['replace_by_addition'])
    add_arg_to_parser(diff, files_args['replace_by_product'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_args['verbose'])
    add_arg_to_parser(status, runtime_args['percentage'])
    args = parser.parse_args()

    # 2. Initializations
    ####################
    # 2.1 options
    if args.verbose:
        args.progressmode = 'verbose'
    elif args.percentage:
        args.progressmode = 'percentage'
    else:
        args.progressmode = None
    # 2.2 list of fields to be processed
    if args.field is not None:
        args.fieldseed = args.field
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            args.fieldseed = [line.replace('\n', '').strip() for line in l.readlines()]
    else:
        args.fieldseed = None

    return args

