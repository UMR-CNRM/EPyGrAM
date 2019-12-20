#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import argparse

from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog, epygramError
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  misc_options, output_options,
                                  runtime_options, graphical_options)

import matplotlib.pyplot as plt
import cartopy.feature as cf

CFEATURES = [f for f in dir(cf) if all([c.isupper() for c in f])]


def read_and_preprocess(resource,
                        fid,
                        operation,
                        global_shift_center,
                        pressure_pa2hpa,
                        zoom):
    """Read field in resource, and preprocess if requested."""
    field = resource.readfield(fid)
    assert isinstance(field, epygram.fields.H2DField), \
        ' '.join(['Oops ! Looks like {} is not known as a horizontal 2D Field by epygram.',
                  'Add it to ~/.epygram/user_Field_Dict_FA.csv ?']).format(fid)
    if field.spectral:
        field.sp2gp()
    if zoom is not None:
        field = field.extract_zoom(zoom)
    if operation is not None:
        field.operation(**operation)
    if global_shift_center is not None:
        field.global_shift_center(global_shift_center)
    if pressure_pa2hpa:
        field.operation('/', 100.)
    return field


def main(filename,
         fid=None,
         Ufid=None,
         Vfid=None,
         refname=None,
         diffonly=False,
         # pre-processing
         operation=None,
         diffoperation=None,
         pressure_pa2hpa=False,
         global_shift_center=None,
         zoom=None,
         # figure
         title=None,
         difftitle=None,
         # geometry
         subzone=None,
         # graphical settings
         plot_method='pcolormesh',
         minmax=None,
         diffminmax=None,
         mask_threshold=None,
         colorsnumber=50,
         diffcolorsnumber=50,
         colormap='plasma',
         diffcolormap='RdBu_r',
         center_cmap_on_0=False,
         diffcenter_cmap_on_0=True,
         scatter_kw=None,
         # cartography
         parallels='auto',
         meridians='auto',
         french_depts=False,
         cartopy_features=[],
         # wind/vectors
         vectors_subsampling=20,
         vector_plot_method='quiver',
         wind_components_are_projected_on='grid',
         quiverkey=None,
         map_factor_correction=False,
         # output
         savefig=False,
         outputfilename=None,
         figures_dpi=epygram.config.default_figures_dpi,
         outputfmt='png'):
    """
    Plot fields.

    :param filename: name of the file to be processed.
    :param fid: field identifier.
    :param Ufid: U-component of wind field identifier.
    :param Vfid: V-component of wind field identifier.
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only plots the difference field.

    Pre-processing:

    :param operation: makes the requested operation
        (e.g. {'operation':'-','operand':273.15} or
        {'operation':'exp'}) on the field before plot.
    :param diffoperation: makes the requested operation
        (e.g. {'operation':'-','operand':273.15} or
        {'operation':'exp'}) on the difference field before plot.
    :param pressure_pa2hpa: converts pressure fields to hPa.
    :param global_shift_center: for global lon/lat grids, shift the center by the
        requested angle (in degrees). Enables a [0,360] grid
        to be shifted to a [-180,180] grid, for instance (with -180 argument).
    :param zoom: a dict(lonmin, lonmax, latmin, latmax) on which to build the plot.

    Figure:

    :param title: title to be written over plot.
    :param difftitle: title to be written over diff plot.

    Geometry:

    :param subzone: LAM zone among ('C', 'CI', None).

    Graphical settings:

    :param plot_method: matplotlib plotting method to be used, among
        ('pcolormesh', 'contourf', 'contour', 'scatter').
    :param minmax: tuple giving (or not) min and max fields values to be plotted.
    :param diffminmax: idem for difference fields.
    :param colorsnumber: number of color discretization/isolines for fields plots.
    :param diffcolorsnumber: idem for difference fields.
    :param colormap: name of the colormap for fields plots.
    :param diffcolormap: idem for difference fields.
    :param center_cmap_on_0: to center the colormap on 0.
    :param diffcenter_cmap_on_0: NOT to center the diffcolormap on 0.
    :param mask_threshold: dict with min and/or max value(s) to mask outside.
    :param scatter_kw: kwargs to be passed to matplotlib's ax.scatter().
        Only for plot_method = 'scatter'.

    Cartography:
    :param meridians and parallels: enable to fine-tune the choice of lines to
        plot, with either:
        - 'auto': automatic scaling to the map extents
        - 'default': range(0,360,10) and range(-90,90,10)
        - a list of values
        - a grid step, e.g. 5 to plot each 5 degree.
        - None: no one is plot
    :param french_depts: draws french departments instead of countries boundaries.
    :param cartopy_features: list of cartopy.feature.??? features.

    Vector plots:

    :param vectors_subsampling: subsampling ratio of vectors plots.
    :param vector_plot_method: among ('quiver', 'barbs', 'streamplot') for vector plots.
    :param wind_components_are_projected_on: inform the plot on which axes the
        vector components are projected on
        ('grid' or 'lonlat').
    :param quiverkey: options to be passed to plotfield to activate a quiver key
        (cf. pyplot.quiverkey).
    :param map_factor_correction: if True, applies a correction of magnitude
        to vector due to map factor.

    Output:

    :param savefig: save figures to file, instead of interactive plot
    :param outputfmt: output format, among ('png', 'pdf', ...).
        Overwritten by outputfilename.
    :param outputfilename: specify an output filename for the plot,
        including format as extension.
    :param figures_dpi: quality of saved figures.
    """
    # 0/ checks, determine mode, initializations
    # checks
    assert not all([f is None for f in (fid, Ufid, Vfid)]), "Mandatory arguments: *fid* OR *Ufid/Vfid*."
    if fid is not None:
        assert Ufid is Vfid is None, "Exclusive arguments: *fid* OR *Ufid/Vfid*."
    if Ufid is not None:
        assert Vfid is not None, "Arguments Ufid/Vfid got by pair."
    # mode
    windmode = fid is None
    diffmode = refname is not None
    assert not (diffmode and windmode), "Exclusive options: diff to reference file and wind plot."
    # plot options
    plot_kwargs = dict(
        # geometry
        subzone=subzone,
        # graphical settings
        plot_method=plot_method,
        minmax=minmax,
        mask_threshold=mask_threshold,
        scatter_kw=scatter_kw,
        # cartography
        parallels=parallels,
        meridians=meridians,
        epygram_departments=french_depts,
        cartopy_features=cartopy_features,
        # colormapping
        colormap=colormap,
        colorsnumber=colorsnumber,
        center_cmap_on_0=center_cmap_on_0)
    # pre-processing options
    preprocess_options = dict(
        operation=operation,
        global_shift_center=global_shift_center,
        pressure_pa2hpa=pressure_pa2hpa,
        zoom=zoom)

    # 1/ resource(s)
    resource = epygram.formats.resource(filename, openmode='r')
    if resource.format == 'DDHLFA':
        raise epygramError('use ddhlfa_plot.py tool for DDHLFA files.')
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if windmode:
        # 2.1/ wind
        u = read_and_preprocess(resource, Ufid,
                                **preprocess_options)
        v = read_and_preprocess(resource, Vfid,
                                **preprocess_options)
        field = epygram.fields.make_vector_field(u, v)
        if title is None:
            title = "\n".join([str(fid), str(field.validity.get())])
        fig, _ = field.cartoplot(map_factor_correction=map_factor_correction,
                                 subsampling=vectors_subsampling,
                                 components_are_projected_on=wind_components_are_projected_on,
                                 vector_plot_method=vector_plot_method,
                                 vector_plot_kwargs=None,
                                 quiverkey=quiverkey,
                                 # module_plot_kwargs
                                 title=title,
                                 **plot_kwargs)

    else:
        # 2.2/ scalar field
        if diffmode:
            # 2.2.1/ diff of scalar fields
            # read field and ref, compute diff
            field = read_and_preprocess(resource, fid,
                                        **preprocess_options)
            ref_field = read_and_preprocess(reference, fid,
                                            **preprocess_options)
            diff_field = field - ref_field
            if diffoperation is not None:
                diff_field.operation(**diffoperation)
            if title is None:
                title = "\n".join([str(fid), str(field.validity.get())])
            if not diffonly:
                # plot field and ref
                fig, _ = field.cartoplot(title="\n".join([resource.container.basename, title]),
                                         **plot_kwargs)
                ref_fig, _ = ref_field.cartoplot(title="\n".join([reference.container.basename, title]),
                                                 **plot_kwargs)
            # set diff specifics
            if difftitle is None:
                difftitle = "\n".join([resource.container.basename + ' - ' +
                                       reference.container.basename,
                                       title])
            plot_kwargs.update(
                minmax=diffminmax,
                colorsnumber=diffcolorsnumber,
                colormap=diffcolormap,
                center_cmap_on_0=diffcenter_cmap_on_0)
            # plot diff
            diff_fig, _ = diff_field.cartoplot(title=difftitle,
                                               **plot_kwargs)
        else:
            # 2.2.2/ plot single scalar fields
            field = read_and_preprocess(resource, fid,
                                        **preprocess_options)
            if title is None:
                title = "\n".join([str(fid), str(field.validity.get())])
            fig, _ = field.cartoplot(title=title, **plot_kwargs)
    # 3/ output
    if savefig:
        epylog.info("save plot(s)...")
        parameter = smooth_string(fid)
        save_kwargs = dict(bbox_inches='tight', dpi=figures_dpi)
        if diffmode:
            if not diffonly:
                fig.savefig('.'.join([resource.container.abspath,
                                      parameter,
                                      outputfmt]),
                            **save_kwargs)
                ref_fig.savefig('.'.join([reference.container.abspath,
                                          parameter,
                                          outputfmt]),
                                **save_kwargs)
            diff_fig.savefig(resource.container.absdir +
                             '.'.join(['diff',
                                       resource.container.basename + '-' + reference.container.basename,
                                       parameter,
                                       outputfmt]),
                             **save_kwargs)
        else:
            if outputfilename is None:
                if outputfmt is not None:
                    outputfilename = '.'.join([resource.container.abspath,
                                               parameter,
                                               outputfmt])
                else:
                    raise epygramError('*outputfmt* or *outputfilename* must be supplied if *savefig*.')
            fig.savefig(outputfilename, **save_kwargs)
    else:
        plt.show()
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
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

    add_arg_to_parser(parser, misc_options['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_options['operation_on_field'])
    add_arg_to_parser(parser, misc_options['mask_threshold'])
    add_arg_to_parser(parser, misc_options['wind_components_are_projected_on'])
    add_arg_to_parser(parser, misc_options['map_factor_correction'])
    add_arg_to_parser(parser, misc_options['LAMzone'])
    # graphics
    add_arg_to_parser(parser, graphical_options['plot_method'])
    add_arg_to_parser(parser, graphical_options['minmax'])
    add_arg_to_parser(parser, graphical_options['levels_number'])
    add_arg_to_parser(parser, graphical_options['colormap'], default='plasma')
    add_arg_to_parser(parser, graphical_options['center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['title'])
    add_arg_to_parser(parser, graphical_options['cartopy_features'],
                      help="cartopy features (cartopy.feature.*), separated by comma " +
                      str(CFEATURES))
    add_arg_to_parser(parser, graphical_options['french_departments'])
    add_arg_to_parser(parser, graphical_options['parallels'])
    add_arg_to_parser(parser, graphical_options['meridians'])
    add_arg_to_parser(parser, graphical_options['vectors_subsampling'])
    add_arg_to_parser(parser, graphical_options['scatter_kw'])
    add_arg_to_parser(parser, graphical_options['lonlat_zoom'])
    add_arg_to_parser(parser, graphical_options['vector_plot_method'])
    add_arg_to_parser(parser, graphical_options['quiverkey'])
    add_arg_to_parser(parser, graphical_options['figures_dpi'])
    add_arg_to_parser(parser, graphical_options['global_shift_center'])
    # diff
    add_arg_to_parser(parser, misc_options['diffoperation_on_field'])
    add_arg_to_parser(parser, graphical_options['diffminmax'])
    add_arg_to_parser(parser, graphical_options['diff_levels_number'])
    add_arg_to_parser(parser, graphical_options['diffcolormap'])
    add_arg_to_parser(parser, graphical_options['diff_center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_options['difftitle'])
    # output
    add_arg_to_parser(parser, output_options['outputfmt'])
    add_arg_to_parser(parser, output_options['outputfilename'], default=None)
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    # 2. Initializations
    ####################
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
        zoom = str2dict(args.zoom, float)
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
    if args.parallels == 'None':
        parallels = None
    elif ',' in args.parallels:
        parallels = [float(p.strip()) for p in args.parallels.split(',')]
    else:
        try:
            parallels = float(args.parallels)
        except ValueError:
            parallels = args.parallels
    if args.meridians == 'None':
        meridians = None
    elif ',' in args.meridians:
            meridians = [float(m.strip()) for m in args.meridians.split(',')]
    else:
        try:
            meridians = float(args.meridians)
        except ValueError:
            meridians = args.meridians
    if args.mask_threshold is not None:
        mask_threshold = str2dict(args.mask_threshold, float)
    else:
        mask_threshold = None
    if args.quiverkey is None or args.quiverkey == '':
        quiverkey = None
    else:
        quiverkey = str2dict(args.quiverkey, float)
    if args.scatter_kw is not None:
        scatter_kw = str2dict(args.scatter_kw, int)
    else:
        scatter_kw = None
    if args.cartopy_features is not None:
        cartopy_features = args.cartopy_features.split(',')
    else:
        cartopy_features = []
    if args.outputfilename or args.outputfmt:
        savefig = True
    else:
        savefig = False

    # 2.2 field to be processed
    if args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        if None in (args.Ucomponentofwind, args.Vcomponentofwind):
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        if diffmode:
            raise NotImplementedError("diffmode (-d/D) AND wind mode (--wU/wV) options together.")
    elif args.field is None:
        raise epygramError("Need to specify a field (-f) or two wind fields (--wU/--wV).")

    # 3. Main
    #########
    main(six.u(args.filename),
         fid=args.field,
         Ufid=args.Ucomponentofwind,
         Vfid=args.Vcomponentofwind,
         refname=refname,
         diffonly=diffonly,
         # pre-processing
         operation=operation,
         diffoperation=diffoperation,
         pressure_pa2hpa=args.pressure_unit_hpa,
         global_shift_center=args.global_shift_center,
         zoom=zoom,
         # figure
         title=args.title,
         difftitle=args.difftitle,
         # geometry
         subzone=subzone,
         # graphical settings
         plot_method=args.plot_method,
         minmax=minmax,
         diffminmax=diffminmax,
         mask_threshold=mask_threshold,
         colorsnumber=args.levelsnumber,
         diffcolorsnumber=args.difflevelsnumber,
         colormap=args.colormap,
         diffcolormap=args.diffcolormap,
         center_cmap_on_0=args.center_cmap_on_0,
         diffcenter_cmap_on_0=args.diffcenter_cmap_on_0,
         scatter_kw=scatter_kw,
         # cartography
         parallels=parallels,
         meridians=meridians,
         french_depts=args.depts,
         cartopy_features=cartopy_features,
         # wind/vectors
         vectors_subsampling=args.vectors_subsampling,
         vector_plot_method=args.vector_plot_method,
         wind_components_are_projected_on=args.wind_components_are_projected_on,
         quiverkey=quiverkey,
         map_factor_correction=args.map_factor_correction,
         # output
         savefig=savefig,
         outputfilename=args.outputfilename,
         figures_dpi=args.figures_dpi,
         outputfmt=args.outputfmt)
