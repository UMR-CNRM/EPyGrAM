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
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  runtime_options)


def main(filename, fieldseed,
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
        epylog.warning(" ".join(["tool NOT TESTED with format",
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
                epylog.info(str(f))
            elif progressmode == 'percentage':
                printstatus(n, numfields)
                n += 1
            if source.format == 'GRIB':
                options = [k[0] for k in epygram.config.GRIB_default_packing[1]] + \
                          [k[0] for k in epygram.config.GRIB_default_packing[2]]
                options.extend(list(epygram.config.GRIB_default_ordering.keys()))
                options.append('centre')
                options = {'get_info_as_json':options}
            else:
                options = {}
            field = source.readfield(f, **options)
            if source.format == 'FA':
                options = {'compression':source.fieldscompression.get(f, None)}
            elif source.format == 'GRIB':
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
                options = {'sample':'file:' + source.container.abspath,
                           'packing':{k:options.pop(k, v) for (k, v) in
                                      epygram.config.GRIB_default_packing[grib_edition].items()},
                           'ordering':{k:options.pop(k, v) for (k, v) in
                                       epygram.config.GRIB_default_ordering.items()},
                           'other_GRIB_options':options
                           }
            output.writefield(field, **options)
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description='An EPyGrAM tool for removing field(s) from a resource.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    add_arg_to_parser(parser, files_management['in_place'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['field'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, fields_management['reverse_fields_selection'])
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
    main(six.u(args.filename), fieldseed,
         reverse=args.reverse,
         progressmode=progressmode,
         in_place=args.in_place)
