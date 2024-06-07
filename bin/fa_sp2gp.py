#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import os
import sys
import argparse

from bronx.fancies.display import printstatus

# Automatically set the python path
package_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_path, 'src'))
import epygram
from epygram import epylog
from epygram.args_catalog import (add_arg_to_parser, files_management,
                                  fields_management, runtime_options)


def main(filename,
         progressmode=None,
         compression=None,
         in_place=False):
    """
    Converts spectral fields to gridpoint.

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
            epylog.info(f)
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
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description="An EPyGrAM tool for transforming spectral fields from a FA to gridpoint.",
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    add_arg_to_parser(parser, files_management['in_place'])
    add_arg_to_parser(parser, fields_management['FA_set_compression'])
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
    if args.verbose:
        progressmode = 'verbose'
    elif args.percentage:
        progressmode = 'percentage'
    else:
        progressmode = None

    # 3. Main
    #########
    main(six.u(args.filename),
         progressmode=progressmode,
         compression=args.kompression,
         in_place=args.in_place)
