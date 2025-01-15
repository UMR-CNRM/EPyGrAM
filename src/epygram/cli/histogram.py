#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for making histograms of meteorological fields from a resource."""

import argparse
import numpy

from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog as logger
from epygram import epygramError
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           output_args,
                           runtime_args,
                           graphical_args)

import matplotlib.pyplot as plt

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    histogram(args.filename,
              args.fieldseed,
              refname=args.refname,
              diffonly=args.diffonly,
              computewind=args.computewind,
              subzone=args.subzone,
              operation=args.operation,
              diffoperation=args.diffoperation,
              legend=args.legend,
              output=args.output if args.output != 'X' else None,
              outputfilename=args.outputfilename,
              minmax=args.minmax,
              diffminmax=args.diffminmax,
              bins=args.bins,
              diffbins=args.diffbins,
              center_hist_on_0=args.center_hist_on_0,
              diff_center_hist_on_0=args.diffcenter_hist_on_0,
              zoom=args.zoom,
              figures_dpi=args.figures_dpi,
              mask_threshold=args.mask_threshold)


def histogram(filename,
              fieldseed,
              refname=None,
              diffonly=False,
              computewind=False,
              subzone=None,
              operation=None,
              diffoperation=None,
              legend=None,
              output=False,
              outputfilename=None,
              minmax=None,
              diffminmax=None,
              bins=50,
              diffbins=50,
              center_hist_on_0=False,
              diff_center_hist_on_0=True,
              zoom=None,
              figures_dpi=epygram.config.default_figures_dpi,
              mask_threshold=None
              ):
    """
    Build histograms.

    :param filename: name of the file to be processed.
    :param fieldseed: field identifier.
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only plots the difference histogram.
    :param computewind: from fieldseed, gets U and V components of wind, and
                     computes the module; plots barbs and module together.
    :param subzone: LAM zone among ('C', 'CI', None).
    :param operation: makes the requested operation
                   (e.g. {'operation':'-','operand':273.15} or
                   {'operation':'exp'}) on the field before hist.
    :param diffoperation: makes the requested operation
                       (e.g. {'operation':'-','operand':273.15} or
                       {'operation':'exp'}) on the difference field before hist.
    :param legend: legend to be written over plot.
    :param output: output format, among ('png', 'pdf', False).
    :param outputfilename: specify an output filename for the plot
                        (completed by output format).
    :param minmax: tuple giving (or not) min and max fields values to be selected.
    :param diffminmax: idem for difference fields.
    :param bins: number of bins or bins edges.
    :param diffbins: idem for difference fields.
    :param zoom: a dict(lonmin, lonmax, latmin, latmax) on which to build the hist.
    :param figures_dpi: quality of saved figures.
    :param mask_threshold: dict with min and/or max value(s) to mask outside.
    """
    if outputfilename and not output:
        raise epygramError('*output* format must be defined if outputfilename is supplied.')

    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if not diffmode:
        if not computewind:
            field = resource.readfield(fieldseed)
            assert isinstance(field, epygram.fields.H2DField), \
                   ' '.join(['Oops ! Looks like',
                             str(fieldseed),
                             'is not known as a horizontal 2D Field by epygram.',
                             'Add it to ~/.epygram/user_Field_Dict_FA.csv ?'])
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            if field.spectral:
                field.sp2gp()
            if operation is not None:
                field.operation(**operation)
            if legend is not None:
                title = legend
            else:
                title = str(fieldseed) + "\n" + str(field.validity.get())
            if zoom is not None:
                field = field.extract_zoom(zoom)
            plot, _ = field.histogram(subzone=subzone,
                                      title=title,
                                      mask_threshold=mask_threshold,
                                      range=minmax,
                                      bins=bins,
                                      center_hist_on_0=center_hist_on_0)
            if not output:
                plt.show()
        else:
            assert len(fieldseed) == 2
            (Ufid, Vfid) = fieldseed
            U = resource.readfield(Ufid)
            field = U
            assert isinstance(field, epygram.fields.H2DField), \
                   ' '.join(['Oops ! Looks like',
                             str(fieldseed),
                             'is not known as a horizontal 2D Field by epygram.',
                             'Add it to ~/.epygram/user_Field_Dict_FA.csv ?'])
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            if U.spectral:
                U.sp2gp()
            if operation is not None:
                U.operation(**operation)
            V = resource.readfield(Vfid)
            if V.spectral:
                V.sp2gp()
            if operation is not None:
                V.operation(**operation)
            vectwind = epygram.fields.make_vector_field(U, V)
            field = vectwind.to_module()
            if legend is not None:
                title = legend
            else:
                title = 'module ' + str(fieldseed) + "\n" + str(U.validity.get())
            if zoom is not None:
                field = field.extract_zoom(zoom)
            plot, _ = field.histogram(subzone=subzone,
                                      title=title,
                                      mask_threshold=mask_threshold,
                                      range=minmax,
                                      bins=bins,
                                      center_hist_on_0=center_hist_on_0)
            if not output:
                plt.show()
    else:
        field = resource.readfield(fieldseed)
        assert isinstance(field, epygram.fields.H2DField), \
               ' '.join(['Oops ! Looks like',
                         str(fieldseed),
                         'is not known as a horizontal 2D Field by epygram.',
                         'Add it to ~/.epygram/user_Field_Dict_FA.csv ?'])
        if not field.geometry.grid.get('LAMzone', False):
            subzone = None
        if field.spectral:
            field.sp2gp()
        if operation is not None:
            field.operation(**operation)
        if not diffonly:
            if legend is not None:
                title = legend
            else:
                title = resource.container.basename + " : " + str(fieldseed) + "\n" + str(field.validity.get())
            if zoom is not None:
                field = field.extract_zoom(zoom)
        reffield = reference.readfield(fieldseed)
        if reffield.spectral:
            reffield.sp2gp()
        if operation is not None:
            reffield.operation(**operation)
        if not diffonly:
            if legend is not None:
                title = legend
            else:
                title = str(fieldseed)
            if zoom is not None:
                reffield = reffield.extract_zoom(zoom)
            label = [resource.container.basename, reference.container.basename]
            plot, _ = field.histogram(subzone=subzone,
                                      title=title,
                                      mask_threshold=mask_threshold,
                                      together_with=reffield,
                                      center_hist_on_0=center_hist_on_0,
                                      range=minmax,
                                      bins=bins,
                                      label=label)
        if legend is not None:
            title = legend
        else:
            title = resource.container.basename + " - " + reference.container.basename + " : " + str(fieldseed)
        diff = field - reffield
        if diffoperation is not None:
            diff.operation(**diffoperation)
        diffplot, _ = diff.histogram(subzone=subzone,
                                     title=title,
                                     mask_threshold=mask_threshold,
                                     range=diffminmax,
                                     bins=diffbins,
                                     center_hist_on_0=diff_center_hist_on_0)
        if not output:
            plt.show()

    # Output
    if output:
        logger.info("save plots...")
        suffix = '.'.join(['hist', output])
        parameter = smooth_string(fieldseed)
        # main resource
        if not diffonly:
            if not diffmode and outputfilename:
                outputfile = '.'.join([outputfilename, output])
            else:
                outputfile = '.'.join([resource.container.abspath,
                                       parameter,
                                       suffix])
            plot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)
        # diff
        if diffmode:
            if not outputfilename:
                outputfile = resource.container.absdir + \
                             '.'.join(['diff',
                                       resource.container.basename + '-' + reference.container.basename,
                                       parameter,
                                       suffix])
            else:
                outputfile = '.'.join([outputfilename, output])
            diffplot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, fields_args['field'])
    add_arg_to_parser(parser, fields_args['windfieldU'])
    add_arg_to_parser(parser, fields_args['windfieldV'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_args['output'])
    add_arg_to_parser(parser, output_args['outputfilename'])
    add_arg_to_parser(parser, misc_args['LAMzone'])
    add_arg_to_parser(parser, graphical_args['bins'])
    add_arg_to_parser(parser, graphical_args['diffbins'])
    add_arg_to_parser(parser, graphical_args['minmax'])
    add_arg_to_parser(parser, graphical_args['diffminmax'])
    add_arg_to_parser(parser, graphical_args['center_hist_on_0'])
    add_arg_to_parser(parser, graphical_args['diff_center_hist_on_0'])
    add_arg_to_parser(parser, graphical_args['legend'])
    add_arg_to_parser(parser, graphical_args['lonlat_zoom'])
    add_arg_to_parser(parser, graphical_args['figures_dpi'])
    add_arg_to_parser(parser, misc_args['operation_on_field'])
    add_arg_to_parser(parser, misc_args['diffoperation_on_field'])
    add_arg_to_parser(parser, misc_args['mask_threshold'])
    add_arg_to_parser(parser, runtime_args['verbose'])
    args = parser.parse_args()

    # 2. Initializations
    ####################
    # 2.1 options
    if args.Drefname is not None:
        args.refname = args.Drefname
        args.diffonly = True
    else:
        args.refname = args.refname
        args.diffonly = False
    args.diffmode = refname is not None
    if args.zone in ('C', 'CI'):
        args.subzone = args.zone
    elif args.zone == 'CIE':
        args.subzone = None
    if args.bins is not None:
        if 'range' in args.bins:
            args.bins = list(numpy.arange(
                *[float(b) for b in args.bins.strip('range').strip('(').strip(')').split(',')]))
        else:
            args.bins = [float(b) for b in args.bins.split(',')]
        if len(args.bins) == 1:
            args.bins = int(args.bins[0])
    else:
        args.bins = None
    if args.diffbins is not None:
        if 'range' in args.diffbins:
            args.diffbins = list(numpy.arange(
                *[float(b) for b in args.diffbins.strip('range').strip('(').strip(')').split(',')]))
        else:
            args.diffbins = [float(b) for b in args.diffbins.split(',')]
        if len(args.diffbins) == 1:
            args.diffbins = int(args.diffbins[0])
    else:
        args.diffbins = None
    if args.minmax is not None:
        args.minmax = args.minmax.split(',')
    else:
        args.minmax = None
    if args.diffminmax is not None:
        diffminmax = args.diffminmax.split(',')
    else:
        args.diffminmax = None
    if args.zoom is not None:
        args.zoom = str2dict(args.zoom, float)
    else:
        args.zoom = None
    if args.operation is not None:
        _operation = args.operation.split(',')
        args.operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            args.operation['operand'] = float(_operation.pop(0).strip())
    else:
        args.operation = None
    if args.diffoperation is not None:
        _diffoperation = args.diffoperation.split(',')
        args.diffoperation = {'operation':_diffoperation.pop(0).strip()}
        if len(_diffoperation) > 0:
            args.diffoperation['operand'] = float(_diffoperation.pop(0).strip())
    else:
        args.diffoperation = None
    if args.mask_threshold is not None:
        args.mask_threshold = str2dict(args.mask_threshold, float)
    else:
        args.mask_threshold = None
    # 2.2 field to be processed
    args.computewind = False
    if args.field is not None:
        args.fieldseed = args.field
    elif args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        args.fieldseed = (args.Ucomponentofwind, args.Vcomponentofwind)
        if None in args.fieldseed:
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        args.computewind = True
        if args.diffmode:
            raise NotImplementedError("diffmode (-d/D) AND wind mode (--wU/wV) options together.")
    else:
        raise epygramError("Need to specify a field (-f) or two wind fields (--wU/--wV).")

    return args
