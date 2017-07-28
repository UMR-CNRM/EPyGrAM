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
from epygram.args_catalog import add_arg_to_parser, \
                                 files_management, fields_management, \
                                 misc_options, output_options, \
                                 runtime_options, graphical_options

import matplotlib.pyplot as plt


def main(filename,
         fieldseed,
         refname=None,
         diffonly=False,
         computewind=False,
         subzone=None,
         operation=None,
         diffoperation=None,
         pressure_unit_hpa=False,
         legend=None,
         output=False,
         outputfilename=None,
         specificproj=None,
         gisquality='i',
         graphicmode='colorshades',
         minmax=None,
         diffminmax=None,
         levelsnumber=50,
         difflevelsnumber=50,
         colormap='jet',
         diffcolormap='seismic',
         center_cmap_on_0=False,
         diffcenter_cmap_on_0=True,
         parallels='auto',
         meridians='auto',
         french_depts=False,
         drawrivers=False,
         vectors_subsampling=20,
         zoom=None,
         pointsize=20,
         symbol='barbs',
         figures_dpi=epygram.config.default_figures_dpi,
         wind_components_are_projected_on='grid',
         bluemarble=0.,
         background=False,
         composition=None,
         mask_threshold=None,
         quiverkey=None,
         global_shift_center=None,
         ):
    """
    Args:
        filename: name of the file to be processed.
        fieldseed: field identifier.
        refname: name of the reference file to be compared to.
        diffonly: if True, only plots the difference field.
        computewind: from fieldseed, gets U and V components of wind, and
                     computes the module; plots barbs and module together.
        subzone: LAM zone among ('C', 'CI', None).
        operation: makes the requested operation
                   (e.g. {'operation':'-','operand':273.15} or
                   {'operation':'exp'}) on the field before plot.
        diffoperation: makes the requested operation
                       (e.g. {'operation':'-','operand':273.15} or
                       {'operation':'exp'}) on the difference field before plot.
        pressure_unit_hpa: converts FA log(pressure) fields to hPa.
        legend: legend to be written over plot.
        output: output format, among ('png', 'pdf', False). Overwritten by
                outputfilename.
        outputfilename: specify an output filename for the plot
                        (completed by output format).
        specificproj: specific projection name (cf.
                      :mod:`epygram.fields.H2DField`.plotfield() doc).
        gisquality: quality of coastlines and countries boundaries.
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
        meridians and parallels: enable to fine-tune the choice of lines to
                                 plot, with either:
                                 - 'auto': automatic scaling to the basemap extents
                                 - 'default': range(0,360,10) and range(-90,90,10)
                                 - a list of values
                                 - a grid step, e.g. 5 to plot each 5 degree.
                                 - None: no one is plot
                                 - *meridian* == 'greenwich' // 'datechange' // 'greenwich+datechange'
                                   *parallel* == 'equator' // 'polarcircles' // 'tropics' or any
                                   combination (+) will plot only these.
        french_depts: draws french departments instead of countries boundaries.
        drawrivers: draw rivers.
        vectors_subsampling: subsampling ratio of vectors plots.
        zoom: a dict(lonmin, lonmax, latmin, latmax) on which to build the plot.
        pointsize: size of the point, case graphicmode='points'.
        symbol: among ('barbs', 'arrows', 'stream') for vector plots.
        figures_dpi: quality of saved figures.
        wind_components_are_projected_on: inform the plot on which axes the
                                          vector components are projected on
                                          ('grid' or 'lonlat').
        bluemarble: if >0., displays NASA's "blue marble" as background with
        given transparency.
        background: if True, set a background color to continents and oceans.
        composition: dict containing info for a composition of fields.(cf. code for doc)
        mask_threshold: dict with min and/or max value(s) to mask outside.
        quiverkey: options to be passed to plotfield to activate a quiver key (cf. pyplot.quiverkey).
        global_shift_center: for global lon/lat grids, shift the center by the
                             requested angle (in degrees). Enables a [0,360] grid
                             to be shifted to a [-180,180] grid, for instance (with -180 argument).
    """

    quiverkey = epygram.util.ifNone_emptydict(quiverkey)
    if outputfilename and not output:
        raise epygramError('*output* format must be defined if outputfilename is supplied.')

    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format not in ('GRIB', 'FA', 'LFI', 'TIFFMF'):
        if resource.format == 'DDHLFA':
            raise epygramError('use ddhlfa_plot.py tool for DDHLFA files.')
        epylog.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    bm = None
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
            if global_shift_center is not None:
                field.global_shift_center(global_shift_center)
            if composition is not None:
                if 'file' in composition:
                    toreadin = epygram.formats.resource(composition['file'], 'r')
                else:
                    toreadin = resource
                composefield = toreadin.readfield(composition['fid'])
                assert isinstance(composefield, epygram.fields.H2DField), \
                       ' '.join(['Oops ! Looks like',
                                 str(fieldseed),
                                 'is not known as a horizontal 2D Field by epygram.',
                                 'Add it to ~/.epygram/user_Field_Dict_FA.csv ?'])
                if 'file' in composition:
                    del toreadin
                if 'preset' in composition:
                    composefield.operation(composition['preset'])
                field.operation(composition['operation'], operand=composefield)
            if pressure_unit_hpa and \
               (field.fid['generic'].get('discipline') == 0 and
                field.fid['generic'].get('parameterCategory') == 3 and
                field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
                field.operation('/', 100.)
            if legend is not None:
                title = legend
            else:
                title = str(fieldseed) + "\n" + str(field.validity.get())
            if field.geometry.name != 'academic':
                bm = field.geometry.make_basemap(subzone=subzone, specificproj=specificproj,
                                                 zoom=zoom, gisquality=gisquality)
            plot, _ = field.plotfield(subzone=subzone,
                                      specificproj=specificproj,
                                      title=title,
                                      use_basemap=bm,
                                      minmax=minmax,
                                      graphicmode=graphicmode,
                                      levelsnumber=levelsnumber,
                                      drawrivers=drawrivers,
                                      parallels=parallels,
                                      meridians=meridians,
                                      colormap=colormap,
                                      center_cmap_on_0=center_cmap_on_0,
                                      departments=french_depts,
                                      pointsize=pointsize,
                                      bluemarble=bluemarble,
                                      background=background,
                                      mask_threshold=mask_threshold,
                                      zoom=zoom)
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
            if global_shift_center is not None:
                U.global_shift_center(global_shift_center)
                V.global_shift_center(global_shift_center)
            vectwind = epygram.fields.make_vector_field(U, V)
            if resource.format == 'FA' and 'WIND.U.PHYS' in U.fid['FA']:
                map_factor_correction = False
            else:
                map_factor_correction = True
            if legend is not None:
                title = legend
            else:
                title = str(fieldseed) + "\n" + str(U.validity.get())
            if field.geometry.name != 'academic':
                bm = U.geometry.make_basemap(subzone=subzone,
                                             specificproj=specificproj,
                                             zoom=zoom,
                                             gisquality=gisquality)
            fig, _ = vectwind.plotfield(subzone=subzone,
                                        use_basemap=bm,
                                        title=title,
                                        parallels=parallels,
                                        meridians=meridians,
                                        subsampling=vectors_subsampling,
                                        symbol=symbol,
                                        departments=french_depts,
                                        bluemarble=bluemarble,
                                        background=background,
                                        plot_module=not bluemarble,
                                        quiverkey=quiverkey,
                                        plot_module_options=dict(minmax=minmax,
                                                                 graphicmode=graphicmode,
                                                                 levelsnumber=levelsnumber,
                                                                 colormap=colormap,
                                                                 # departments=french_depts,
                                                                 pointsize=pointsize,
                                                                 # bluemarble=bluemarble,
                                                                 # background=background,
                                                                 mask_threshold=mask_threshold),
                                        components_are_projected_on=wind_components_are_projected_on,
                                        map_factor_correction=map_factor_correction)
            plot = fig
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
        if global_shift_center is not None:
            field.global_shift_center(global_shift_center)
        if pressure_unit_hpa and \
           (field.fid['generic'].get('discipline') == 0 and
            field.fid['generic'].get('parameterCategory') == 3 and
            field.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
            field.operation('/', 100.)
        if not diffonly:
            if legend is not None:
                title = legend
            else:
                title = resource.container.basename + " : " + str(fieldseed) + "\n" + str(field.validity.get())
            if field.geometry.name != 'academic':
                bm = field.geometry.make_basemap(subzone=subzone, specificproj=specificproj,
                                                 gisquality=gisquality, zoom=zoom)
            plot, _ = field.plotfield(subzone=subzone,
                                      specificproj=specificproj,
                                      title=title,
                                      use_basemap=bm,
                                      minmax=minmax,
                                      graphicmode=graphicmode,
                                      levelsnumber=levelsnumber,
                                      drawrivers=drawrivers,
                                      parallels=parallels,
                                      meridians=meridians,
                                      colormap=colormap,
                                      center_cmap_on_0=center_cmap_on_0,
                                      departments=french_depts,
                                      pointsize=pointsize,
                                      bluemarble=bluemarble,
                                      background=background,
                                      mask_threshold=mask_threshold)
        else:
            bm = None
        reffield = reference.readfield(fieldseed)
        if reffield.spectral:
            reffield.sp2gp()
        if operation is not None:
            reffield.operation(**operation)
        if global_shift_center is not None:
            field.global_shift_center(global_shift_center)
        if pressure_unit_hpa and \
           (reffield.fid['generic'].get('discipline') == 0 and
            reffield.fid['generic'].get('parameterCategory') == 3 and
            reffield.fid['generic'].get('parameterNumber') in (0, 1, 2, 8, 11, 25)):
            reffield.operation('/', 100.)
        if not diffonly:
            if legend is not None:
                title = legend
            else:
                title = reference.container.basename + " : " + str(fieldseed) + "\n" + str(reffield.validity.get())
            refplot, _ = reffield.plotfield(subzone=subzone,
                                            specificproj=specificproj,
                                            title=title,
                                            use_basemap=bm,
                                            minmax=minmax,
                                            graphicmode=graphicmode,
                                            levelsnumber=levelsnumber,
                                            drawrivers=drawrivers,
                                            parallels=parallels,
                                            meridians=meridians,
                                            colormap=colormap,
                                            center_cmap_on_0=center_cmap_on_0,
                                            departments=french_depts,
                                            pointsize=pointsize,
                                            bluemarble=bluemarble,
                                            background=background,
                                            mask_threshold=mask_threshold)
        if legend is not None:
            title = legend
        else:
            title = resource.container.basename + " - " + reference.container.basename + " : " + str(fieldseed)
        diff = field - reffield
        if bm is None and field.geometry.name != 'academic':
            bm = diff.geometry.make_basemap(subzone=subzone,
                                            specificproj=specificproj,
                                            zoom=zoom)
        if diffoperation is not None:
            diff.operation(**diffoperation)
        diffplot, _ = diff.plotfield(subzone=subzone,
                                     specificproj=specificproj,
                                     title=title,
                                     use_basemap=bm,
                                     minmax=diffminmax,
                                     graphicmode=graphicmode,
                                     levelsnumber=difflevelsnumber,
                                     drawrivers=drawrivers,
                                     parallels=parallels,
                                     meridians=meridians,
                                     colormap=diffcolormap,
                                     center_cmap_on_0=diffcenter_cmap_on_0,
                                     departments=french_depts,
                                     pointsize=pointsize,
                                     bluemarble=bluemarble,
                                     background=background,
                                     mask_threshold=mask_threshold)
        if not output:
            plt.show()

    # Output
    if output:
        epylog.info("save plots...")
        suffix = '.'.join(['plot', output])
        parameter = epygram.util.linearize2str(fieldseed)
        # main resource
        if not diffonly:
            if not diffmode and outputfilename:
                outputfile = '.'.join([outputfilename, output])
            else:
                outputfile = '.'.join([resource.container.abspath,
                                       parameter,
                                       suffix])
            plot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)
        # reference
        if diffmode and not diffonly:
            outputfile = '.'.join([reference.container.abspath,
                                   parameter,
                                   suffix])
            refplot.savefig(outputfile, bbox_inches='tight', dpi=figures_dpi)
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
# end of main() ###############################################################


if __name__ == '__main__':

    # ## 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description='An EPyGrAM tool for simple plots\
                                                  of meteorological fields from a resource.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    add_arg_to_parser(parser, fields_management['field'])
    add_arg_to_parser(parser, fields_management['windfieldU'])
    add_arg_to_parser(parser, fields_management['windfieldV'])
    diffmodes = parser.add_mutually_exclusive_group()
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(diffmodes, files_management['file_to_refer_in_diffonly'])
    add_arg_to_parser(parser, output_options['output'])
    add_arg_to_parser(parser, output_options['outputfilename'])
    add_arg_to_parser(parser, misc_options['LAMzone'])
    add_arg_to_parser(parser, graphical_options['specific_map_projection'])
    add_arg_to_parser(parser, graphical_options['graphicmode'])
    add_arg_to_parser(parser, graphical_options['minmax'])
    add_arg_to_parser(parser, graphical_options['diffminmax'])
    add_arg_to_parser(parser, graphical_options['levels_number'])
    add_arg_to_parser(parser, graphical_options['diff_levels_number'])
    add_arg_to_parser(parser, graphical_options['colormap'])
    add_arg_to_parser(parser, graphical_options['diffcolormap'])
    add_arg_to_parser(parser, graphical_options['center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['diff_center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['legend'])
    add_arg_to_parser(parser, graphical_options['gis_quality'])
    add_arg_to_parser(parser, graphical_options['draw_rivers'])
    add_arg_to_parser(parser, graphical_options['french_departments'])
    add_arg_to_parser(parser, graphical_options['parallels'])
    add_arg_to_parser(parser, graphical_options['meridians'])
    add_arg_to_parser(parser, graphical_options['vectors_subsampling'])
    add_arg_to_parser(parser, graphical_options['points_size'])
    add_arg_to_parser(parser, graphical_options['lonlat_zoom'])
    add_arg_to_parser(parser, graphical_options['vector_symbol'])
    add_arg_to_parser(parser, graphical_options['quiverkey'])
    add_arg_to_parser(parser, graphical_options['figures_dpi'])
    add_arg_to_parser(parser, graphical_options['bluemarble'])
    add_arg_to_parser(parser, graphical_options['background'])
    add_arg_to_parser(parser, graphical_options['global_shift_center'])
    add_arg_to_parser(parser, misc_options['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_options['operation_on_field'])
    add_arg_to_parser(parser, misc_options['diffoperation_on_field'])
    add_arg_to_parser(parser, misc_options['composition_with_field'])
    add_arg_to_parser(parser, misc_options['mask_threshold'])
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    # ## 2. Initializations
    ######################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.Drefname is not None:
        refname = args.Drefname
        diffonly = True
    else:
        refname = args.refname
        diffonly = False
    diffmode = refname is not None
    if args.zone in ('C', 'CI'):
        subzone = args.zone
    elif args.zone == 'CIE':
        subzone = None
    if args.minmax is not None:
        minmax = args.minmax.split(',')
    else:
        minmax = None
    if args.diffminmax is not None:
        diffminmax = args.diffminmax.split(',')
    else:
        diffminmax = None
    if args.zoom is not None:
        zoom = epygram.util.parse_str2dict(args.zoom, float)
    else:
        zoom = None
    if args.operation is not None:
        _operation = args.operation.split(',')
        operation = {'operation':_operation.pop(0).strip()}
        if len(_operation) > 0:
            operation['operand'] = float(_operation.pop(0).strip())
    else:
        operation = None
    if args.diffoperation is not None:
        _diffoperation = args.diffoperation.split(',')
        diffoperation = {'operation':_diffoperation.pop(0).strip()}
        if len(_diffoperation) > 0:
            diffoperation['operand'] = float(_diffoperation.pop(0).strip())
    else:
        diffoperation = None
    if args.compose_with is not None:
        _composition = args.compose_with.split(',')
        composition = {'fid':_composition.pop(0).strip(),
                       'operation':_composition.pop(0).strip()}
        while len(_composition) > 0:
            composition.update(epygram.util.parse_str2dict(_composition.pop(0).strip()))
    else:
        composition = None
    if args.projection is not None and 'nsper' in args.projection:
        specificproj = ('nsper', {})
        for item in args.projection.split(',')[1:]:
            k, v = item.replace('=', ':').split(':')
            specificproj[1][k.strip()] = float(v)
    else:
        specificproj = args.projection
    if args.parallels == 'None':
        parallels = None
    else:
        if ',' in args.parallels:
            parallels = [p.strip() for p in args.parallels.split(',')]
        else:
            parallels = args.parallels
    if args.meridians == 'None':
        meridians = None
    else:
        if ',' in args.meridians:
            meridians = [m.strip() for m in args.meridians.split(',')]
        else:
            meridians = args.meridians
    if args.mask_threshold is not None:
        mask_threshold = epygram.util.parse_str2dict(args.mask_threshold, float)
    else:
        mask_threshold = None
    if args.quiverkey is not None:
        quiverkey = epygram.util.parse_str2dict(args.quiverkey, float)
    else:
        quiverkey = {}

    # 2.2 field to be processed
    computewind = False
    if args.field is not None:
        fieldseed = args.field
    elif args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        fieldseed = (args.Ucomponentofwind, args.Vcomponentofwind)
        if None in fieldseed:
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        computewind = True
        if diffmode:
            raise NotImplementedError("diffmode (-d/D) AND wind mode (--wU/wV) options together.")
    else:
        raise epygramError("Need to specify a field (-f) or two wind fields (--wU/--wV).")

    # ## 3. Main
    ###########
    main(six.u(args.filename),
         fieldseed,
         refname=refname,
         diffonly=diffonly,
         computewind=computewind,
         subzone=subzone,
         operation=operation,
         diffoperation=diffoperation,
         pressure_unit_hpa=args.pressure_unit_hpa,
         legend=args.legend,
         output=args.output,
         outputfilename=args.outputfilename,
         specificproj=specificproj,
         gisquality=args.gisquality,
         graphicmode=args.graphicmode,
         minmax=minmax,
         diffminmax=diffminmax,
         levelsnumber=args.levelsnumber,
         difflevelsnumber=args.difflevelsnumber,
         colormap=args.colormap,
         diffcolormap=args.diffcolormap,
         center_cmap_on_0=args.center_cmap_on_0,
         diffcenter_cmap_on_0=args.diffcenter_cmap_on_0,
         parallels=parallels,
         meridians=meridians,
         french_depts=args.depts,
         drawrivers=args.drawrivers,
         vectors_subsampling=args.vectors_subsampling,
         zoom=zoom,
         pointsize=args.pointsize,
         symbol=args.vector_symbol,
         figures_dpi=args.figures_dpi,
         bluemarble=args.bluemarble,
         background=args.background,
         composition=composition,
         mask_threshold=mask_threshold,
         quiverkey=quiverkey,
         global_shift_center=args.global_shift_center)

###########
### END ###
###########
