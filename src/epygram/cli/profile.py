#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for extracting (and plotting) vertical profiles of meteorological fields from a resource."""

import argparse
import numpy

from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog as logger
from epygram.fields.V1DField import plotprofiles
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           output_args,
                           runtime_args,
                           graphical_args,
                           extraction_args)

import matplotlib.pyplot as plt

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    profile(
         args.filename,
         args.fieldseed,
         args.coordinates,
         interpolation=args.interpolation,
         refname=args.refname,
         diffonly=args.diffonly,
         operation=args.operation,
         diffoperation=args.diffoperation,
         legend=args.legend,
         output=args.output if args.output != 'X' else None,
         outputfilename=args.outputfilename,
         zoom=args.zoom,
         Yconvert=args.Yconvert,
         cheap_height=args.cheap_height,
         noplot=args.noplot,
         unit=args.unit,
         logscale=args.logscale,
         emagramlike=args.emagramlike,
         external_distance=args.external_distance,
         figures_dpi=args.figures_dpi)


def profile(
         filename,
         fieldseed,
         coordinates,
         interpolation='nearest',
         refname=None,
         diffonly=False,
         operation=None,
         diffoperation=None,
         legend=None,
         output=False,
         outputfilename=None,
         zoom=None,
         Yconvert=None,
         cheap_height=True,
         noplot=False,
         unit='SI',
         logscale=False,
         emagramlike=False,
         external_distance=None,
         figures_dpi=epygram.config.default_figures_dpi):
    """
    Extract and plot profiles.

    :param filename: name of the file to be processed.
    :param fieldseed: either a string or a list of strings, used as a seed for
                   generating the list of fields to be processed.
    :param coordinates: (lon, lat) coordinates (in °) of the point to be extracted.
                     If None, make horizontally-averaged profile.
    :param interpolation: kind of interpolation from grid to coordinates, among
                       ('nearest', 'linear', 'cubic').
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only plots the difference field.
    :param operation: makes the requested operation
                   (e.g. {'operation':'-','operand':273.15} or
                   {'operation':'exp'}) on the field before plot.
    :param diffoperation: makes the requested operation
                       (e.g. {'operation':'-','operand':273.15} or
                       {'operation':'exp'}) on the difference field before plot.
    :param legend: legend of plot.
    :param output: output format for graphics, among ('png', 'pdf', False).
    :param outputfilename: specify an output filename for the profile.
                        The graphical output, if requested, will be completed
                        by the requested output format.
    :param zoom: a dict(ymin, ymax) to restrain the plot.
    :param Yconvert: among ('pressure', 'height', 'altitude'), or [fid, kind],
                  to convert the vertical coordinate.
                  For height/altitude, implies the read of T and q
                  profiles and optionally pressure departure, hydrometeors.
                  If [fid, kind] is profided, the vertical coordinate is read using the fid
                  and is considered to be of the kind specified (eg. 100 for pressure).
    :param cheap_height: if True, do not take hydrometeors nor pressure departure
                      into account for computing height from pressure.
    :param noplot: if True, disable the plot of profiles.
    :param unit: scientifical unit of plot.
    :param logscale: if True, sets vertical logscale.
    :param emagramlike: if True, tilts the Temperature axis "as an emagram".
    :param external_distance: can be a dict containing the target point value
                           and an additional field fid, to which field the distance
                           for computing nearest point
                           is computed within the 4 horizontally nearest points; e.g.
                           {'target_value':4810, 'external_field':an_additional_fid}.
                           If so, the nearest point is selected with
                           distance = |target_value - external_field.data|
    :param figures_dpi: quality of saved figures.
    """
    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI'):
        logger.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    nosave = False
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if external_distance:
        external_distance['external_field'] = resource.readfield(external_distance['external_field'])
        if external_distance['external_field'].spectral:
            external_distance['external_field'].sp2gp()

    def make_profile_from(r, fieldseed, coordinates,
                          vertical_coordinate,
                          interpolation,
                          cheap_height,
                          external_distance):
        """Encapsulation of profile extracting from any resource 'r'."""

        map_vcoord = {'pressure':100,
                      'altitude':102,
                      'height':103}
        if not isinstance(vertical_coordinate, list):
            vertical_coordinate = map_vcoord.get(vertical_coordinate, vertical_coordinate)
        profile = r.extractprofile(fieldseed,
                                   *coordinates,
                                   vertical_coordinate=vertical_coordinate if isinstance(vertical_coordinate, int) else None,
                                   interpolation=interpolation,
                                   cheap_height=cheap_height,
                                   external_distance=external_distance)
        if vertical_coordinate is not None and not isinstance(vertical_coordinate, int):
            # vertical_coordinate is a list containing the fid of the field to use as a vertical coordinate
            # and the GRIB code corresponding to the kind of this vertical level coordinate
            profileVCoord = r.extractprofile(vertical_coordinate[0],
                                             *coordinates,
                                             vertical_coordinate=None,
                                             interpolation=interpolation)
            profile.use_field_as_vcoord(profileVCoord, force_kind=int(vertical_coordinate[1]))

        if operation is not None:
            profile.operation(**operation)

        return profile
    # end of make_profile_from() ##############################################

    profile = make_profile_from(resource, fieldseed, coordinates,
                                Yconvert,
                                interpolation,
                                cheap_height,
                                external_distance)
    if not diffmode:
        if not noplot:
            if emagramlike and profile.geometry.vcoordinate.typeoffirstfixedsurface != 100:
                emagramlike = False
            if legend is not None:
                title = legend
            else:
                title = profile.comment + '\n' + str(profile.validity.get())
            plot, _ = plotprofiles(profile,
                                   fidkey=resource.format,
                                   unit=unit,
                                   logscale=logscale,
                                   title=title,
                                   ema=emagramlike,
                                   zoom=zoom)
            if not output:
                plt.show()
    else:
        refprofile = make_profile_from(reference, fieldseed, coordinates,
                                       Yconvert,
                                       interpolation,
                                       cheap_height,
                                       external_distance)
        if not diffonly:
            labels = [resource.container.basename,
                      reference.container.basename,
                      '- diff -']
        else:
            labels = [resource.container.basename + ' - ' + reference.container.basename]
        toplot = epygram.base.FieldSet()
        if not diffonly:
            toplot.append(profile)
            toplot.append(refprofile)
        if profile.geometry == refprofile.geometry:
            diffprofile = profile - refprofile
            if diffoperation is not None:
                diffprofile.operation(**diffoperation)
            toplot.append(diffprofile)
            same_Z = True
        else:
            epygram.logger.warning('profiles vertical grids differ (surface pressure ?): cannot compute difference.')
            same_Z = False
            nosave = True
            # raise epygram.epygramError("unable to compute profiles difference because of vertical grids differ: \
            #                             surely, because surface pressure differ.")
        if not noplot:
            if emagramlike and profile.geometry.vcoordinate.typeoffirstfixedsurface != 100:
                emagramlike = False
            if legend is not None:
                title = legend
            else:
                title = str(profile.fid.get(resource.format, profile.fid)) + '\n' + profile.comment
            if same_Z:
                plot, _ = plotprofiles(toplot,
                                       labels=labels,
                                       fidkey=resource.format,
                                       unit=unit,
                                       logscale=logscale,
                                       title=title,
                                       ema=emagramlike,
                                       zoom=zoom)
            else:
                plot, ax = plt.subplots()
                for i in range(len(toplot)):
                    plot, _ = plotprofiles(toplot[i],
                                           over=(plot, ax),
                                           labels=[labels[i]],
                                           fidkey=resource.format,
                                           unit=unit,
                                           logscale=logscale,
                                           title=title,
                                           ema=emagramlike,
                                           zoom=zoom)

            if not output:
                plt.show()

    # Text Output
    parameter = smooth_string(profile.fid.get(resource.format, profile.fid))
    if None in coordinates:
        position = 'mean'
    else:
        position = str(coordinates[0]) + "E" + str(coordinates[1]) + "N"
    suffix = "profile.out"
    if not diffmode:
        filename = '.'.join([resource.container.abspath,
                             parameter, position, suffix])
        profiletosave = profile
    elif not nosave:
        filename = '.'.join([resource.container.absdir + 'diff',
                             resource.container.basename + '-' + reference.container.basename,
                             parameter, position, suffix])
        profiletosave = diffprofile
    if outputfilename:
        filename = outputfilename
    if not nosave:
        L = len(profile.geometry.vcoordinate.levels)
        precision = 4
        fieldnamelength = 16
        length_Z = len(str(profile.geometry.vcoordinate.typeoffirstfixedsurface))
        Z = numpy.array(profile.geometry.vcoordinate.levels).flatten()

        flds = '{:<{width}}'.format("Z: " + str(profile.geometry.vcoordinate.typeoffirstfixedsurface),
                                    width=length_Z) \
               + '{:^{width}}'.format(parameter, width=fieldnamelength + 2)
        with open(filename, 'w') as o:
            o.write(flds + "\n")
            for k in range(0, L):
                line = '{:^{width}}'.format(str(Z[k]), width=length_Z)
                valstr = '{:.{precision}{type}}'.format(profiletosave.getdata()[k],
                                                        type='E',
                                                        precision=precision)
                line += '{:^{width}}'.format(valstr, width=fieldnamelength + 2)
                o.write(line + "\n")

    # Graphical Output
    if output:
        logger.info("save plot...")
        suffix = 'profile.' + output
        if outputfilename:
            filename = '.'.join([outputfilename, output])
        else:
            if not diffmode:
                filename = '.'.join([resource.container.abspath,
                                     parameter, position, suffix])
            else:
                filename = '.'.join([resource.container.absdir + 'diff',
                                     resource.container.basename + '-' + reference.container.basename,
                                     parameter, position, suffix])

        plot.savefig(filename, bbox_inches='tight', dpi=figures_dpi)


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    add_arg_to_parser(parser, fields_args['vertical_field'])
    add_arg_to_parser(parser, extraction_args['point_coordinates'],
                      required=False, default=None,
                      help=(extraction_args['point_coordinates'][-1]['help'] +
                            'If not given, make horizontally-averaged profile.'))
    add_arg_to_parser(parser, extraction_args['horizontal_interpolation'])
    add_arg_to_parser(parser, extraction_args['external_distance'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_args['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_args['output'])
    add_arg_to_parser(parser, output_args['outputfilename'])
    add_arg_to_parser(parser, output_args['noplot'])
    Y = parser.add_mutually_exclusive_group()
    add_arg_to_parser(Y, extraction_args['verticalcoord2pressure'])
    add_arg_to_parser(Y, extraction_args['verticalcoord2height'])
    add_arg_to_parser(Y, extraction_args['verticalcoord2altitude'])
    add_arg_to_parser(Y, fields_args['external_vertical_coord'])
    add_arg_to_parser(parser, extraction_args['no_cheap_height_conversion'])
    add_arg_to_parser(parser, graphical_args['legend'])
    add_arg_to_parser(parser, graphical_args['emagram_like_profiles'])
    add_arg_to_parser(parser, graphical_args['vertical_logscale'])
    add_arg_to_parser(parser, graphical_args['scientifical_unit'])
    add_arg_to_parser(parser, graphical_args['vertical_zoom'])
    add_arg_to_parser(parser, graphical_args['figures_dpi'])
    add_arg_to_parser(parser, misc_args['operation_on_field'])
    add_arg_to_parser(parser, misc_args['diffoperation_on_field'])
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
    if args.coordinates is not None:
        args.coordinates = args.coordinates.split(',')
        args.coordinates = tuple([float(i) for i in args.coordinates])
    else:
        args.coordinates = (None, None)
    if args.zoom is not None:
        zoom = dict()
        for limit in args.zoom.split(','):
            l, v = limit.split('=')
            zoom[l.strip()] = float(v)
        args.zoom = zoom
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
    if args.external_distance is not None:
        external_distance = args.external_distance.split(';')
        assert len(external_distance) == 2, "syntax must be: 'VALUE; FIELD_ID'"
        external_distance = {'target_value':float(external_distance[0]),
                             'external_field':external_distance[1]}
        args.external_distance = external_distance
    # 2.2 field to be processed
    args.fieldseed = args.field

    return args

