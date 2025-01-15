#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for transforming spectral fields from a FA to gridpoint."""

import argparse

from bronx.fancies.display import printstatus

import epygram
from epygram import epylog as logger
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
    fa_sp2gp(args.filename,
             progressmode=args.progressmode,
             compression=args.kompression,
             in_place=args.in_place)


def fa_sp2gp(filename,
             progressmode=None,
             compression=None,
             in_place=False):
    """
    Converts spectral fields of a FA to gridpoint.

    :param filename: name of the file to be processed.
    :param progressmode: among ('verbose', 'percentage', None).
    :param compression: number of bits for gridpoint-converted fields to be
                     encoded on.
    """
    source = epygram.formats.resource(filename, openmode='a', fmt='FA')
    if in_place:
        output = source
    else:
        output = epygram.formats.resource(filename + '.allgp', openmode='w', fmt='FA',
                                          headername=source.headername,
                                          validity=source.validity,
                                          cdiden=source.cdiden,
                                          default_compression=source.default_compression,
                                          processtype=source.processtype)

    fidlist = source.listfields()
    numfields = len(fidlist)
    n = 1
    for f in fidlist:
        if progressmode == 'verbose':
            logger.info(f)
        elif progressmode == 'percentage':
            printstatus(n, numfields)
            n += 1
        field = source.readfield(f)
        fieldcompression = None
        if isinstance(field, epygram.fields.H2DField):
            if field.spectral:
                field.sp2gp()
                if compression is not None:
                    fieldcompression = {'KNBPDG':compression}
                else:
                    fieldcompression = {'KNBPDG':source.fieldscompression[f]['KNBCSP']}
            else:
                if in_place:
                    continue
                fieldcompression = source.fieldscompression[f]
            if fieldcompression.get('KNBPDG') == 0:
                fieldcompression['KNGRIB'] = 0
        output.writefield(field, fieldcompression)


def get_args():

    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, files_args['in_place'])
    add_arg_to_parser(parser, fields_args['FA_set_compression'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_args['verbose'])
    add_arg_to_parser(status, runtime_args['percentage'])
    args = parser.parse_args()

    if args.verbose:
        args.progressmode = 'verbose'
    elif args.percentage:
        args.progressmode = 'percentage'
    else:
        args.progressmode = None

    return args
