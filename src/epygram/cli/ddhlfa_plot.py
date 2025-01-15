#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for simple plots of meteorological fields from a DDHLFA."""

import argparse

import epygram
from epygram import epylog as logger
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           graphical_args,
                           runtime_args)
from epygram.fields.V1DField import plotprofiles

import matplotlib.pyplot as plt

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    ddhlfa_plot(args.filename,
                args.fieldseed,
                args.domains,
                refname=args.refname,
                legend=args.legend,
                unit=args.unit,
                logscale=args.logscale)

def ddhlfa_plot(filename,
                fieldseed,
                domains,
                refname=None,
                legend=None,
                unit='SI',
                logscale=False):
    """
    Plot DDHLFA.

    :param filename: name of the file to be processed.
    :param fieldseed: either a string or a list of strings, used as a seed for
                   generating the list of fields to be processed.
    :param domains: list of domains number within the DDHLFA.
    :param refname: name of the reference file to be compared to.
    :param legend: legend of plot.
    :param unit: scientifical unit of plot.
    :param logscale: if True, sets vertical logscale.
    """
    diffmode = refname is not None
    resource = epygram.formats.resource(filename, openmode='r', fmt='DDHLFA')
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r', fmt='DDHLFA')

    if not diffmode:
        plots = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='V1D')
        for f in fidlist:
            logger.info(f)
            fields = resource.readfield(f)
            fields = epygram.base.FieldSet([fields[i - 1] for i in domains])
            if legend is not None:
                title = legend
            else:
                title = f + "\n" + str(fields[0].validity.get())
            plots[f] = plotprofiles(fields,
                                    labels=['Domain ' + str(d) for d in domains],
                                    unit=unit,
                                    title=title,
                                    logscale=logscale)
            plt.show()
    else:
        plots = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='V1D')
        reffidlist = reference.find_fields_in_resource(seed=fieldseed, fieldtype='V1D')
        intersectionfidlist = list(set(fidlist).intersection(set(reffidlist)))
        intersectionfidlist.sort()
        for f in intersectionfidlist:
            logger.info(f)
            fields = resource.readfield(f)
            fields = epygram.base.FieldSet([fields[i - 1] for i in domains])
            reffields = reference.readfield(f)
            reffields = epygram.base.FieldSet([reffields[i - 1] for i in domains])
            plots[f] = {}
            for d in range(len(domains)):
                if legend is not None:
                    title = legend
                else:
                    title = f + ' on domain ' + str(domains[d]) + "\n" + str(fields[d].validity.get())
                profs = epygram.base.FieldSet((fields[d], reffields[d]))
                labels = [resource.container.basename,
                          reference.container.basename]
                plots[f][str(domains[d])] = plotprofiles(profs,
                                                         labels=labels,
                                                         unit=unit,
                                                         title=title,
                                                         logscale=logscale)
                plt.show()
# end of main() ###############################################################


def get_args():

    # ## 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)

    add_arg_to_parser(parser, files_args['principal_file'])
    flds = parser.add_mutually_exclusive_group(required=True)
    add_arg_to_parser(flds, fields_args['DDHLFA_multiple_fields'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    add_arg_to_parser(parser, fields_args['DDHLFA_domain'])
    add_arg_to_parser(parser, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(parser, graphical_args['legend'])
    add_arg_to_parser(parser, graphical_args['scientifical_unit'])
    add_arg_to_parser(parser, graphical_args['vertical_logscale'])
    add_arg_to_parser(parser, runtime_args['verbose'])

    args = parser.parse_args()

    # ## 2. Initializations
    ######################

    # 2.1 arguments
    if isinstance(args.domain_number, int):
        args.domains = [args.domain_number]
    elif isinstance(args.domain_number, str):
        args.domains = [int(i) for i in args.domain_number.split(',')]

    # 2.2 list of fields to be processed
    if args.field is not None:
        fieldseed = args.field
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            args.fieldseed = [line.replace('\n', '').strip() for line in l.readlines()]
    else:
        args.fieldseed = None
    return args

