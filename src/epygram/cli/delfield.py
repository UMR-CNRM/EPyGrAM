#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for removing field(s) from a resource."""

import argparse

from bronx.fancies.display import printstatus

import epygram
from epygram import epygramError
from epygram import epylog as logger
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           runtime_args)
from epygram.extra import griberies

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    delfield(args.filename,
             args.fieldseed,
             reverse=args.reverse,
             progressmode=args.progressmode,
             in_place=args.in_place)


def delfield(filename,
             fieldseed,
             reverse=False,
             progressmode=None,
             in_place=False):
    """
    Delete fields.

    :param filename: name of the file to be processed.
    :param fieldseed: either a fid or a list of fid, used as a seed for
                      generating the list of fields to be processed.
    :param reverse: if True, reverse the field selection: deletes all but
                    selected fields.
    :param progressmode: among ('verbose', 'percentage', None)
    :param in_place: if True, the field(s) is(are) deleted "in place",
                     not in a new file.
    """
    source = epygram.formats.resource(filename, openmode='a')
    fidlist = source.find_fields_in_resource(seed=fieldseed)
    if source.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
                                 source.format, "!"]))

    if in_place:
        if not (hasattr(source, "delfield") and callable(getattr(source, "delfield"))):
            raise epygramError('unable to remove messages from this format' +
                               ' resource "in place".')
        if reverse:
            for f in source.listfields():
                if f not in fidlist:
                    source.delfield(f)
        else:
            for f in fidlist:
                source.delfield(f)
    else:
        tobecopied = {'fmt':source.format}
        if source.format == 'FA':
            tobecopied.update({'headername':source.headername,
                               'validity':source.validity,
                               'cdiden':source.cdiden,
                               'default_compression':source.default_compression,
                               'processtype':source.processtype})
        elif source.format == 'LFI':
            tobecopied.update({'compressed':source.compressed})
        output = epygram.formats.resource(filename + '.delfield.out',
                                          openmode='w',
                                          **tobecopied)
        outputfidlist = source.listfields()
        if reverse:
            outputfidlist = fidlist
        else:
            for f in fidlist:
                if f in outputfidlist:
                    outputfidlist.remove(f)
        numfields = len(outputfidlist)
        n = 1
        for f in outputfidlist:
            if progressmode == 'verbose':
                logger.info(str(f))
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if source.format == 'GRIB':
                read_misc_metadata = set(griberies.defaults.GRIB1_packing.keys())
                read_misc_metadata.update(set(griberies.defaults.GRIB2_keyvalue[5].keys()))
                read_misc_metadata.update(set(griberies.defaults.GRIB1_ordering.keys()))
                read_misc_metadata.add('centre')
                options = {'read_misc_metadata':list(read_misc_metadata)}
            else:
                options = {}
            field = source.readfield(f, **options)
            if source.format == 'FA':
                options = {'compression':source.fieldscompression.get(f, None)}
            elif source.format == 'GRIB':
                options = {'sample':'file:' + source.container.abspath}
            output.writefield(field, **options)


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, files_args['in_place'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_args['field'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    add_arg_to_parser(parser, fields_args['reverse_fields_selection'])
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
        with open(listfile.abspath, 'r') as lf:
            args.fieldseed = [line.replace('\n', '').strip() for line in lf.readlines()]
    else:
        args.fieldseed = None
    
    return args

