#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division

import argparse

import epygram
from epygram import epylog
from epygram.args_catalog import add_arg_to_parser, files_management, \
                                 fields_management, graphical_options, \
                                 runtime_options
from epygram.fields.V1DField import plotprofiles

import matplotlib.pyplot as plt


def main(filename,
         fieldseed,
         domains,
         refname=None,
         legend=None,
         unit='SI',
         logscale=False):
    """
    Args:
        filename: name of the file to be processed.
        fieldseed: either a string or a list of strings, used as a seed for
                   generating the list of fields to be processed.
        domains: list of domains number within the DDHLFA.
        refname: name of the reference file to be compared to.
        legend: legend of plot.
        unit: scientifical unit of plot.
        logscale: if True, sets vertical logscale.
    """

    diffmode = refname is not None
    resource = epygram.formats.resource(filename, openmode='r', fmt='DDHLFA')
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r', fmt='DDHLFA')

    if not diffmode:
        plots = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='V1D')
        for f in fidlist:
            epylog.info(f)
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
            epylog.info(f)
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


if __name__ == '__main__':
    ### 1. Parse arguments
    ######################
    epygramstr = 'EPyGrAM'
    parser = argparse.ArgumentParser(description='An ' + epygramstr + ' tool for simple\
                                                  plots of meteorological fields from a DDHLFA.',
                                     epilog='End of help for: %(prog)s (' + epygramstr + ' v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    flds = parser.add_mutually_exclusive_group(required=True)
    add_arg_to_parser(flds, fields_management['DDHLFA_multiple_fields'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, fields_management['DDHLFA_domain'])
    add_arg_to_parser(parser, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(parser, graphical_options['legend'])
    add_arg_to_parser(parser, graphical_options['scientifical_unit'])
    add_arg_to_parser(parser, graphical_options['vertical_logscale'])
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    ### 2. Initializations
    ######################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 arguments
    if isinstance(args.domain_number, int):
        domains = [args.domain_number]
    elif isinstance(args.domain_number, str):
        domains = args.domain_number.split(',')
        domains = [int(i) for i in domains]

    # 2.2 list of fields to be processed
    if args.field is not None:
        fieldseed = args.field
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            fieldseed = l.readlines()
        for i in range(0, len(fieldseed)):
            fieldseed[i] = fieldseed[i].replace('\n', '').strip()
    else:
        fieldseed = None

    ### 3. Main
    ###########
    main(args.filename,
         fieldseed,
         domains,
         refname=args.refname,
         legend=args.legend,
         unit=args.unit,
         logscale=args.logscale)



###########
### END ###
###########
