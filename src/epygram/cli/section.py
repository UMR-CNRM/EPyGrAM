#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for extracting (and plotting) vertical sections of meteorological fields from a resource."""

import argparse

from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog as logger
from epygram import epygramError

import matplotlib.pyplot as plt
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           output_args,
                           runtime_args,
                           graphical_args,
                           extraction_args)

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    section(
        args.filename,
        args.fieldseed,
        args.starting_point,
        args.ending_point,
        refname=args.refname,
        diffonly=args.diffonly,
        points_number=args.points_number,
        resolution=args.resolution,
        Yconvert=args.Yconvert,
        interpolation=args.interpolation,
        cheap_height=args.cheap_height,
        operation=args.operation,
        diffoperation=args.diffoperation,
        graphicmode=args.graphicmode,
        minmax=args.minmax,
        diffminmax=args.diffminmax,
        levelsnumber=args.levelsnumber,
        difflevelsnumber=args.difflevelsnumber,
        colormap=args.colormap,
        diffcolormap=args.diffcolormap,
        center_cmap_on_0=args.center_cmap_on_0,
        diffcenter_cmap_on_0=args.diffcenter_cmap_on_0,
        global_shift_center=args.global_shift_center,
        zoom=args.zoom,
        legend=args.legend,
        output=args.output if args.output != 'X' else None,
        outputfilename=args.outputfilename,
        figures_dpi=args.figures_dpi,
        section_abscissa=args.section_abscissa,
        mask_threshold=args.mask_threshold)


def section(
    filename,
    fieldseed,
    starting_point,
    ending_point,
    refname=None,
    diffonly=False,
    points_number=None,
    resolution=None,
    Yconvert=None,
    interpolation='nearest',
    cheap_height=True,
    operation=None,
    diffoperation=None,
    graphicmode='colorshades',
    minmax=None,
    diffminmax=None,
    levelsnumber=50,
    difflevelsnumber=50,
    colormap='jet',
    diffcolormap='RdBu_r',
    center_cmap_on_0=False,
    diffcenter_cmap_on_0=True,
    global_shift_center=None,
    zoom=None,
    legend=None,
    output=False,
    outputfilename=None,
    figures_dpi=epygram.config.default_figures_dpi,
    section_abscissa='distance',
    mask_threshold=None):
    """
    Extract and plot section.

    :param filename: name of the file to be processed.
    :param fieldseed: field identifier.
    :param starting_point: (lon, lat) coordinates of the starting point of section.
    :param ending_point: (lon, lat) coordinates of the ending point of section.
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only plots the difference field.
    :param points_number: specific number of points on the transect.
    :param resolution: resolution of the transect.
    :param Yconvert: among ('pressure', 'height', 'altitude'), or [fid, kind],
                  to convert the vertical coordinate.
                  For height/altitude, implies the read of T and q
                  profiles and optionally pressure departure, hydrometeors.
                  If [fid, kind] is profided, the vertical coordinate is read using the fid
                  and is considered to be of the kind specified (eg. 100 for pressure).
    :param interpolation: kind of interpolation from grid to coordinates, among
                       ('nearest', 'linear', 'cubic').
    :param cheap_height: if True, do not take hydrometeors nor pressure departure
                      into account for computing height from pressure.
    :param operation: makes the requested operation
                   (e.g. {'operation':'-','scalar':273.15} or
                   {'operation':'exp'}) on the field before plot.
    :param diffoperation: makes the requested operation
                       (e.g. {'operation':'-','scalar':273.15} or
                       {'operation':'exp'}) on the difference field before plot.
    :param graphicmode: graphical display, among ('colorshades', 'contourlines',
                     'points').
    :param minmax: tuple giving (or not) min and max fields values to be plotted.
    :param diffminmax: idem for difference fields.
    :param levelsnumber: number of color discretization/isolines for fields plots.
    :param difflevelsnumber: idem for difference fields.
    :param colormap: name of the colormap for fields plots.
    :param diffcolormap: idem for difference fields.
    :param center_cmap_on_0: to center the colormap on 0.
    :param diffcenter_cmap_on_0: NOT to center the diffcolormap on 0.
    :param global_shift_center: for global lon/lat grids, shift the center by the
        requested angle (in degrees). Enables a [0,360] grid
        to be shifted to a [-180,180] grid, for instance (with -180 argument).
    :param zoom: a dict(ymin, ymax) to restrain the plot.
    :param legend: legend to be written over plot.
    :param output: output format, among ('png', 'pdf', False).
    :param outputfilename: specify an output filename for the section (completed
                        by the requested output format).
    :param figures_dpi: quality of saved figures.
    :param section_abscissa: abscissa of section, among ('distance', 'lon', 'lat').
    :param mask_threshold: dict with min and/or max value(s) to mask outside.
    """
    if outputfilename and not output:
        raise epygramError('*output* format must be defined if outputfilename is supplied.')

    map_vcoord = {'pressure': 100,
                  'altitude': 102,
                  'height': 103}
    if not isinstance(Yconvert, list):
        Yconvert = map_vcoord.get(Yconvert, Yconvert)

    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')
    section = resource.extractsection(fieldseed,
                                      starting_point, ending_point,
                                      points_number=points_number,
                                      resolution=resolution,
                                      vertical_coordinate=Yconvert if isinstance(Yconvert, int) else None,
                                      interpolation=interpolation,
                                      cheap_height=cheap_height,
                                      global_shift_center=global_shift_center)
    if Yconvert is not None and not isinstance(Yconvert, int):
        # args.Yconvert is a list containing the fid of the field to use as a vertical coordinate
        # and the GRIB code corresponding to the kind of this vertical level coordinate
        sectionVCoord = resource.extractsection(Yconvert[0],
                                                starting_point, ending_point,
                                                points_number=points_number,
                                                resolution=resolution,
                                                vertical_coordinate=None,
                                                interpolation=interpolation,
                                                global_shift_center=global_shift_center)
        section.use_field_as_vcoord(sectionVCoord, force_kind=int(Yconvert[1]))
    if operation is not None:
        section.operation(**operation)
    if diffmode:
        refsection = reference.extractsection(fieldseed,
                                              starting_point, ending_point,
                                              points_number=points_number,
                                              resolution=resolution,
                                              vertical_coordinate=Yconvert if isinstance(Yconvert, int) else None,
                                              interpolation=interpolation,
                                              cheap_height=cheap_height,
                                              global_shift_center=global_shift_center)
        if Yconvert is not None and not isinstance(Yconvert, int):
            # args.Yconvert is a list containing the fid of the field to use as a vertical coordinate
            # and the GRIB code corresponding to the kind of this vertical level coordinate
            sectionVCoord = reference.extractsection(Yconvert[0],
                                                     starting_point, ending_point,
                                                     points_number=points_number,
                                                     resolution=resolution,
                                                     vertical_coordinate=None,
                                                     interpolation=interpolation,
                                                     global_shift_center=global_shift_center)
            refsection.use_field_as_vcoord(sectionVCoord, force_kind=int(Yconvert[1]))

        if operation is not None:
            refsection.operation(**operation)
        diff = section - refsection
        if diffoperation is not None:
            diff.operation(**diffoperation)
    # make plots
    if not diffmode:
        sectionplot, _ = section.plotfield(graphicmode=graphicmode,
                                           minmax=minmax,
                                           levelsnumber=levelsnumber,
                                           center_cmap_on_0=center_cmap_on_0,
                                           colormap=colormap,
                                           zoom=zoom,
                                           title=legend,
                                           x_is=section_abscissa,
                                           mask_threshold=mask_threshold)
    else:
        if not diffonly:
            legend = resource.container.basename + " : " + \
                     str(section.fid.get(resource.format, section.fid)) + "\n" + \
                     str(section.validity.get())
            sectionplot, _ = section.plotfield(graphicmode=graphicmode,
                                               minmax=minmax,
                                               levelsnumber=levelsnumber,
                                               center_cmap_on_0=center_cmap_on_0,
                                               colormap=colormap,
                                               zoom=zoom,
                                               title=legend,
                                               x_is=section_abscissa,
                                               mask_threshold=mask_threshold)
            legend = reference.container.basename + " : " + \
                     str(refsection.fid.get(reference.format, refsection.fid)) + "\n" + \
                     str(refsection.validity.get())
            refplot, _ = refsection.plotfield(graphicmode=graphicmode,
                                              minmax=minmax,
                                              levelsnumber=levelsnumber,
                                              center_cmap_on_0=center_cmap_on_0,
                                              colormap=colormap,
                                              zoom=zoom,
                                              title=legend,
                                              x_is=section_abscissa,
                                              mask_threshold=mask_threshold)
        if legend is None:
            legend = resource.container.basename + " - " + \
                     reference.container.basename + " : " + \
                     str(section.fid.get(resource.format, section.fid))
        diffplot, _ = diff.plotfield(graphicmode=graphicmode,
                                     minmax=diffminmax,
                                     levelsnumber=difflevelsnumber,
                                     center_cmap_on_0=diffcenter_cmap_on_0,
                                     colormap=diffcolormap,
                                     zoom=zoom,
                                     title=legend,
                                     x_is=section_abscissa,
                                     mask_threshold=mask_threshold)

    if not output:
        plt.show()

    # Output
    if output:
        logger.info("save plots...")
        parameter = smooth_string(section.fid.get(resource.format, section.fid))
        suffix = '.'.join(['_'.join([parameter,
                                     str(starting_point[0]) + "E" + str(starting_point[1]) + "N",
                                     str(ending_point[0]) + "E" + str(ending_point[1]) + "N"]),
                           "section",
                           output])
        # main resource
        if not diffonly:
            if not diffmode and outputfilename:
                outputfile = '.'.join([outputfilename, output])
            else:
                outputfile = '.'.join([resource.container.abspath, suffix])
            sectionplot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)
        # reference
        if diffmode and not diffonly:
            outputfile = '.'.join([reference.container.abspath, suffix])
            refplot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)
        # diff
        if diffmode:
            if not outputfilename:
                outputfile = resource.container.absdir + '.'.join(['diff',
                                                                   '-'.join([resource.container.basename,
                                                                             reference.container.basename]),
                                                                   suffix])
            else:
                outputfile = '.'.join([outputfilename, output])
            diffplot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, fields_args['vertical_field'])
    add_arg_to_parser(parser, extraction_args['section_starting_point'])
    add_arg_to_parser(parser, extraction_args['section_ending_point'])
    add_arg_to_parser(parser, extraction_args['horizontal_interpolation'],
                      default='linear',
                      help="interpolation mode from field grid to point/section\
                            coordinates. Among ('nearest', 'linear', 'cubic').\
                            Defaults to 'linear'.")
    t = parser.add_mutually_exclusive_group()
    add_arg_to_parser(t, extraction_args['section_transect_points_number'])
    add_arg_to_parser(t, extraction_args['section_transect_resolution'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_args['output'])
    add_arg_to_parser(parser, output_args['outputfilename'])
    Y = parser.add_mutually_exclusive_group()
    add_arg_to_parser(Y, extraction_args['verticalcoord2pressure'])
    add_arg_to_parser(Y, extraction_args['verticalcoord2height'])
    add_arg_to_parser(Y, extraction_args['verticalcoord2altitude'])
    add_arg_to_parser(Y, fields_args['external_vertical_coord'])
    add_arg_to_parser(parser, extraction_args['no_cheap_height_conversion'])
    add_arg_to_parser(parser, graphical_args['legend'])
    add_arg_to_parser(parser, graphical_args['graphicmode'])
    add_arg_to_parser(parser, graphical_args['minmax'])
    add_arg_to_parser(parser, graphical_args['diffminmax'])
    add_arg_to_parser(parser, graphical_args['levels_number'])
    add_arg_to_parser(parser, graphical_args['diff_levels_number'])
    add_arg_to_parser(parser, graphical_args['colormap'])
    add_arg_to_parser(parser, graphical_args['diffcolormap'])
    add_arg_to_parser(parser, graphical_args['center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_args['diff_center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_args['vertical_zoom'])
    add_arg_to_parser(parser, graphical_args['figures_dpi'])
    add_arg_to_parser(parser, graphical_args['section_abscissa'])
    add_arg_to_parser(parser, graphical_args['global_shift_center'])
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
    starting_point = args.starting_point.split(',')
    args.starting_point = tuple([float(i) for i in starting_point])
    ending_point = args.ending_point.split(',')
    args.ending_point = tuple([float(i) for i in ending_point])
    if args.zoom is not None:
        args.zoom = str2dict(args.zoom, float)
    if args.minmax is not None:
        args.minmax = args.minmax.split(',')
    if args.diffminmax is not None:
        args.diffminmax = args.diffminmax.split(',')
    if args.operation is not None:
        _operation = args.operation.split(',')
        args.operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            args.operation['operand'] = float(_operation.pop(0).strip())
    if args.diffoperation is not None:
        _diffoperation = args.diffoperation.split(',')
        args.diffoperation = {'operation':_diffoperation.pop(0).strip()}
        if len(_diffoperation) > 0:
            args.diffoperation['operand'] = float(_diffoperation.pop(0).strip())
    # FIXME: problem in plot with mask_threshold
    if args.mask_threshold is not None:
        args.mask_threshold = str2dict(args.mask_threshold, float)
    # 2.2 field to be processed
    args.fieldseed = args.field

    return args

