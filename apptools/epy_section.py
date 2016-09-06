#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import argparse

import epygram
from epygram import epylog, epygramError

import matplotlib.pyplot as plt
from epygram.args_catalog import add_arg_to_parser, \
                                 files_management, fields_management, \
                                 misc_options, output_options, \
                                 runtime_options, graphical_options, \
                                 extraction_options



def main(filename,
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
         diffcolormap='seismic',
         center_cmap_on_0=False,
         diffcenter_cmap_on_0=True,
         zoom=None,
         legend=None,
         output=False,
         outputfilename=None,
         figures_dpi=epygram.config.default_figures_dpi,
         section_abscissa='distance',
         mask_threshold=None):
    """
    Args:
        filename: name of the file to be processed.
        fieldseed: field identifier.
        starting_point: (lon, lat) coordinates of the starting point of section.
        ending_point: (lon, lat) coordinates of the ending point of section.
        refname: name of the reference file to be compared to.
        diffonly: if True, only plots the difference field.
        points_number: specific number of points on the transect.
        resolution: resolution of the transect.
        Yconvert: among ('pressure', 'height', 'altitude'),
                  to convert the vertical coordinate.
                  For height/altitude, implies the read of T and q 
                  profiles and optionally pressure departure, hydrometeors.
        interpolation: kind of interpolation from grid to coordinates, among
                       ('nearest', 'linear', 'cubic').
        cheap_height: if True, do not take hydrometeors nor pressure departure
                      into account for computing height from pressure.
        operation: makes the requested operation
                   (e.g. {'operation':'-','scalar':273.15} or
                   {'operation':'exp'}) on the field before plot.
        diffoperation: makes the requested operation
                       (e.g. {'operation':'-','scalar':273.15} or
                       {'operation':'exp'}) on the difference field before plot.
        graphicmode: graphical display, among ('colorshades', 'contourlines',
                     'points').
        minmax: tuple giving (or not) min and max fields values to be plotted.
        diffminmax: idem for difference fields.
        levelsnumber: number of color discretization/isolines for fields plots.
        difflevelsnumber: idem for difference fields.
        colormap: name of the colormap for fields plots.
        diffcolormap: idem for difference fields.
        center_cmap_on_0: to center the colormap on 0.
        diffcenter_cmap_on_0: NOT to center the diffcolormap on 0.
        zoom: a dict(ymin, ymax) to restrain the plot.
        legend: legend to be written over plot.
        output: output format, among ('png', 'pdf', False).
        outputfilename: specify an output filename for the section (completed
                        by the requested output format).
        figures_dpi: quality of saved figures.
        section_abscissa: abscissa of section, among ('distance', 'lon', 'lat').
        mask_threshold: dict with min and/or max value(s) to mask outside.
    """

    if outputfilename and not output:
        raise epygramError('*output* format must be defined if outputfilename is supplied.')

    map_vcoord = {'pressure':100,
                  'altitude':102,
                  'height':103}
    Yconvert = map_vcoord.get(Yconvert, Yconvert)

    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI'):
        epylog.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname != None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')
    section = resource.extractsection(fieldseed,
                                      starting_point, ending_point,
                                      points_number=points_number,
                                      resolution=resolution,
                                      vertical_coordinate=Yconvert,
                                      interpolation=interpolation,
                                      cheap_height=cheap_height)
    if operation is not None:
        section.operation(**operation)
    if diffmode:
        refsection = reference.extractsection(fieldseed,
                                              starting_point, ending_point,
                                              points_number=points_number,
                                              resolution=resolution,
                                              vertical_coordinate=Yconvert,
                                              interpolation=interpolation,
                                              cheap_height=cheap_height)
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
        if legend == None:
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
        epylog.info("save plots...")
        parameter = epygram.util.linearize2str(section.fid.get(resource.format, section.fid))
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

# end of main() ###############################################################



if __name__ == '__main__':

    ### 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description='An EPyGrAM tool for extracting (and plotting) vertical sections \
                                                  of meteorological fields from a resource.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    add_arg_to_parser(parser, fields_management['vertical_field'])
    add_arg_to_parser(parser, extraction_options['section_starting_point'])
    add_arg_to_parser(parser, extraction_options['section_ending_point'])
    add_arg_to_parser(parser, extraction_options['horizontal_interpolation'],
                      default='linear',
                      help="interpolation mode from field grid to point/section\
                            coordinates. Among ('nearest', 'linear', 'cubic').\
                            Defaults to 'linear'.")
    t = parser.add_mutually_exclusive_group()
    add_arg_to_parser(t, extraction_options['section_transect_points_number'])
    add_arg_to_parser(t, extraction_options['section_transect_resolution'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_options['output'])
    add_arg_to_parser(parser, output_options['outputfilename'])
    Y = parser.add_mutually_exclusive_group()
    add_arg_to_parser(Y, extraction_options['verticalcoord2pressure'])
    add_arg_to_parser(Y, extraction_options['verticalcoord2height'])
    add_arg_to_parser(Y, extraction_options['verticalcoord2altitude'])
    add_arg_to_parser(parser, extraction_options['no_cheap_height_conversion'])
    add_arg_to_parser(parser, graphical_options['legend'])
    add_arg_to_parser(parser, graphical_options['graphicmode'])
    add_arg_to_parser(parser, graphical_options['minmax'])
    add_arg_to_parser(parser, graphical_options['diffminmax'])
    add_arg_to_parser(parser, graphical_options['levels_number'])
    add_arg_to_parser(parser, graphical_options['diff_levels_number'])
    add_arg_to_parser(parser, graphical_options['colormap'])
    add_arg_to_parser(parser, graphical_options['diffcolormap'])
    add_arg_to_parser(parser, graphical_options['center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['diff_center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['vertical_zoom'])
    add_arg_to_parser(parser, graphical_options['figures_dpi'])
    add_arg_to_parser(parser, graphical_options['section_abscissa'])
    add_arg_to_parser(parser, misc_options['operation_on_field'])
    add_arg_to_parser(parser, misc_options['diffoperation_on_field'])
    add_arg_to_parser(parser, misc_options['mask_threshold'])
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    ### 2. Initializations
    ######################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.Drefname != None:
        refname = args.Drefname
        diffonly = True
    else:
        refname = args.refname
        diffonly = False
    starting_point = args.starting_point.split(',')
    starting_point = tuple([float(i) for i in starting_point])
    ending_point = args.ending_point.split(',')
    ending_point = tuple([float(i) for i in ending_point])
    if args.zoom != None:
        zoom = epygram.util.parse_str2dict(args.zoom, float)
    else:
        zoom = None
    if args.minmax != None:
        minmax = args.minmax.split(',')
    else:
        minmax = args.minmax
    if args.diffminmax != None:
        diffminmax = args.diffminmax.split(',')
    else:
        diffminmax = args.diffminmax
    if args.operation != None:
        _operation = args.operation.split(',')
        operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            operation['operand'] = float(_operation.pop(0).strip())
    else:
        operation = None
    if args.diffoperation != None:
        _diffoperation = args.diffoperation.split(',')
        diffoperation = {'operation':_diffoperation.pop(0).strip()}
        if len(_diffoperation) > 0:
            diffoperation['operand'] = float(_diffoperation.pop(0).strip())
    else:
        diffoperation = None
    #FIXME: problem in plot with mask_threshold
    if args.mask_threshold is not None:
        mask_threshold = epygram.util.parse_str2dict(args.mask_threshold, float)
    else:
        mask_threshold = None

    # 2.2 field to be processed
    fieldseed = args.field

    ### 3. Main
    ###########
    main(args.filename,
         fieldseed,
         starting_point,
         ending_point,
         refname=refname,
         diffonly=diffonly,
         points_number=args.points_number,
         resolution=args.resolution,
         Yconvert=args.Yconvert,
         interpolation=args.interpolation,
         cheap_height=args.cheap_height,
         operation=operation,
         diffoperation=diffoperation,
         graphicmode=args.graphicmode,
         minmax=minmax,
         diffminmax=diffminmax,
         levelsnumber=args.levelsnumber,
         difflevelsnumber=args.difflevelsnumber,
         colormap=args.colormap,
         diffcolormap=args.diffcolormap,
         center_cmap_on_0=args.center_cmap_on_0,
         diffcenter_cmap_on_0=args.diffcenter_cmap_on_0,
         zoom=zoom,
         legend=args.legend,
         output=args.output,
         outputfilename=args.outputfilename,
         figures_dpi=args.figures_dpi,
         section_abscissa=args.section_abscissa,
         mask_threshold=mask_threshold)

###########
### END ###
###########
