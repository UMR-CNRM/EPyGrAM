#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for simple plots of meteorological fields from a resource."""

import argparse

from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epygramError
from epygram import epylog as logger
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           output_args,
                           runtime_args,
                           graphical_args)

import matplotlib.pyplot as plt
import cartopy.feature as cf

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')

    cartoplot(args.filename,
         fid=args.field,
         Ufid=args.Ucomponentofwind,
         Vfid=args.Vcomponentofwind,
         refname=args.refname,
         diffonly=args.diffonly,
         # pre-processing
         operation=args.operation,
         diffoperation=args.diffoperation,
         pressure_pa2hpa=args.pressure_unit_hpa,
         global_shift_center=args.global_shift_center,
         zoom=args.zoom,
         # figure
         title=args.title,
         difftitle=args.difftitle,
         # geometry
         subzone=args.subzone,
         # graphical settings
         plot_method=args.plot_method,
         minmax=args.minmax,
         diffminmax=args.diffminmax,
         mask_threshold=args.mask_threshold,
         colorsnumber=args.levelsnumber,
         diffcolorsnumber=args.difflevelsnumber,
         colormap=args.colormap,
         diffcolormap=args.diffcolormap,
         center_cmap_on_0=args.center_cmap_on_0,
         diffcenter_cmap_on_0=args.diffcenter_cmap_on_0,
         scatter_kw=args.scatter_kw,
         # cartography
         parallels=args.parallels,
         meridians=args.meridians,
         french_depts=args.depts,
         cartopy_features=args.cartopy_features,
         # wind/vectors
         vectors_subsampling=args.vectors_subsampling,
         vector_plot_method=args.vector_plot_method,
         wind_components_are_projected_on=args.wind_components_are_projected_on,
         quiverkey=args.quiverkey,
         map_factor_correction=args.map_factor_correction,
         # output
         savefig=args.savefig,
         outputfilename=args.outputfilename,
         figures_dpi=args.figures_dpi,
         outputfmt=args.outputfmt if args.outputfmt != 'X' else None)

def read_and_preprocess(resource,
                        fid,
                        operation,
                        global_shift_center,
                        pressure_pa2hpa,
                        zoom):
    """Read field in resource, and preprocess if requested."""
    field = resource.readfield(fid)
    assert isinstance(field, (epygram.fields.H2DField, epygram.fields.H2DVectorField)), \
        ' '.join(['Oops ! Looks like {} is not known as a horizontal 2D Field by epygram.',
                  'Add it to ~/.epygram/user_Field_Dict_{}.csv ?']).format(fid, resource.format)
    if field.spectral:
        field.sp2gp()
    if global_shift_center is not None:
        field.global_shift_center(global_shift_center)
    if zoom is not None:
        field = field.extract_zoom(zoom)
    if operation is not None:
        field.operation(**operation)
    if pressure_pa2hpa:
        field.operation('/', 100.)
    return field


def cartoplot(filename,
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
         wind_components_are_projected_on=None,
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
        vector components are projected on ('grid' or 'lonlat').
        If None (default), look for information in the field, or raise error.
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
            if isinstance(field, epygram.fields.H2DVectorField):
                # Scalar field in a vector field is a true color image
                fig, _ = field.cartoimage(title=title, **plot_kwargs)
            else:
                fig, _ = field.cartoplot(title=title, **plot_kwargs)
    # 3/ output
    if savefig:
        logger.info("save plot(s)...")
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
# end of cartoplot() ###############################################################


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

    add_arg_to_parser(parser, misc_args['pressure_unit_hpa'])
    add_arg_to_parser(parser, misc_args['operation_on_field'])
    add_arg_to_parser(parser, misc_args['mask_threshold'])
    add_arg_to_parser(parser, misc_args['wind_components_are_projected_on'])
    add_arg_to_parser(parser, misc_args['map_factor_correction'])
    add_arg_to_parser(parser, misc_args['LAMzone'])
    # graphics
    add_arg_to_parser(parser, graphical_args['plot_method'])
    add_arg_to_parser(parser, graphical_args['minmax'])
    add_arg_to_parser(parser, graphical_args['levels_number'])
    add_arg_to_parser(parser, graphical_args['colormap'], default='plasma')
    add_arg_to_parser(parser, graphical_args['center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_args['title'])
    add_arg_to_parser(parser, graphical_args['cartopy_features'],
                      help="cartopy features (cartopy.feature.*), separated by comma " +
                      str([f for f in dir(cf) if all([c.isupper() for c in f])]))
    add_arg_to_parser(parser, graphical_args['french_departments'])
    add_arg_to_parser(parser, graphical_args['parallels'])
    add_arg_to_parser(parser, graphical_args['meridians'])
    add_arg_to_parser(parser, graphical_args['vectors_subsampling'])
    add_arg_to_parser(parser, graphical_args['scatter_kw'])
    add_arg_to_parser(parser, graphical_args['lonlat_zoom'])
    add_arg_to_parser(parser, graphical_args['vector_plot_method'])
    add_arg_to_parser(parser, graphical_args['quiverkey'])
    add_arg_to_parser(parser, graphical_args['figures_dpi'])
    add_arg_to_parser(parser, graphical_args['global_shift_center'])
    # diff
    add_arg_to_parser(parser, misc_args['diffoperation_on_field'])
    add_arg_to_parser(parser, graphical_args['diffminmax'])
    add_arg_to_parser(parser, graphical_args['diff_levels_number'])
    add_arg_to_parser(parser, graphical_args['diffcolormap'])
    add_arg_to_parser(parser, graphical_args['diff_center_cmap_on_0'])
    add_arg_to_parser(parser, graphical_args['difftitle'])
    # output
    add_arg_to_parser(parser, output_args['outputfmt'])
    add_arg_to_parser(parser, output_args['outputfilename'], default=None)
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
    args.diffmode = args.refname is not None
    if args.zone in ('C', 'CI'):
        args.subzone = args.zone
    elif args.zone == 'CIE':
        args.subzone = None
    if args.minmax is not None:
        args.minmax = args.minmax.split(',')
    else:
        args.minmax = None
    if args.diffminmax is not None:
        args.diffminmax = args.diffminmax.split(',')
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
    if args.parallels == 'None':
        args.parallels = None
    elif ',' in args.parallels:
        args.parallels = [float(p.strip()) for p in args.parallels.split(',')]
    else:
        try:
            args.parallels = float(args.parallels)
        except ValueError:
            args.parallels = args.parallels
    if args.meridians == 'None':
        args.meridians = None
    elif ',' in args.meridians:
            args.meridians = [float(m.strip()) for m in args.meridians.split(',')]
    else:
        try:
            args.meridians = float(args.meridians)
        except ValueError:
            args.meridians = args.meridians
    if args.mask_threshold is not None:
        args.mask_threshold = str2dict(args.mask_threshold, float)
    else:
        args.mask_threshold = None
    if args.quiverkey is None or args.quiverkey == '':
        args.quiverkey = None
    else:
        args.quiverkey = str2dict(args.quiverkey, float)
    if args.scatter_kw is not None:
        args.scatter_kw = str2dict(args.scatter_kw, int)
    else:
        args.scatter_kw = None
    if args.cartopy_features is not None:
        args.cartopy_features = args.cartopy_features.split(',')
    else:
        args.cartopy_features = []
    if args.outputfilename or args.outputfmt != 'X':
        args.savefig = True
    else:
        args.savefig = False

    # 2.2 field to be processed
    if args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        if None in (args.Ucomponentofwind, args.Vcomponentofwind):
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        if args.diffmode:
            raise NotImplementedError("diffmode (-d/D) AND wind mode (--wU/wV) options together.")
    elif args.field is None:
        raise epygramError("Need to specify a field (-f) or two wind fields (--wU/--wV).")

    return args
