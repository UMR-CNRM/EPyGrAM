#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with graphical options
"""

from . import _defaults

d = {
    'legend':[
        '-L', '--legend',
        dict(help="legend to be written over field plot.",
             default=None)],
    'title':[
        '--title',
        dict(help="title to be written over field plot.",
             default=None)],
    'difftitle':[
        '--difftitle',
        dict(help="title to be written over diff field plot.",
             default=None)],
    'scientifical_unit':[
        '-u', '--unit',
        dict(help="optional unit for labeling plot axis. Defaults to 'SI'.",
             default='SI')],
    'vertical_logscale':[
        '-s', '--logscale',
        dict(action='store_true',
             help="plots with vertical logscale.",
             default=False)],
    'specific_map_projection':[
        '-j', '--projection',
        dict(help="specific graphical projection of plots, among\
                   ('ll', 'geoid').\
                   Default is the actual projection of fields,\
                   or Mollweide projection for Gauss geometry.\
                   If set to 'geoid' the field is plotted around\
                   the geoid.",
             default=None)],
    'graphicmode':[
        '-g', '--graphicmode',
        dict(help="graphical mode for plots, among\
                   ('colorshades', 'contourlines', 'points'). Default is\
                   'colorshades'. There is a known bug (yet unsolved) with\
                   Arpege & contourlines.",
             choices=['colorshades', 'contourlines', 'points'],
             default='colorshades')],
    'plot_method':[
        '--pm', '--plot_method',
        dict(help="plot method, among\
                   ('pcolormesh', 'contourf', 'contour', 'scatter'). Default is\
                   'pcolormesh' for rectangular grids, else 'contourf'.\
                   There is a known bug (yet unsolved) with\
                   Arpege & contourlines.",
             choices=['pcolormesh', 'contourf', 'contour', 'scatter'],
             dest='plot_method',
             default='__default__')],
    'plotmode':[
        'plotmode',
        dict(type=str,
             help="Kind of 3d plot, among ('contour', 'color', 'volume',\
                   'vector', 'stream')",
             choices=['contour', 'color', 'volume', 'vectors', 'streamlines', 'donothing'])],
    'minmax':[
        '-m', '--minmax',
        dict(help="min and max values for the plot.\
                   Syntax: 'min, max'. '0.0, max' also works.\
                   Default is the field min/max values.\
                   In diff mode, this is valuable for resource and reference\
                   only, (min, max) for difference plot should be defined\
                   with --diffminmax option. For negative values, use\
                   short-name option -m without spacetab between option and\
                   argument.",
             default=None)],
    'diffminmax':[
        '-M', '--diffminmax',
        dict(help="min and max values for the difference plot.\
                   Syntax: 'min, max'. '0.0, max' also works.\
                   Default is the difference field min/max values.\
                   For negative values, use short-name option -M without\
                   spacetab between option and argument.",
             default=None)],
    'levels_number':[
        '-n', '--levelsnumber',
        dict(help="number of levels for contours and shades.\
                   Default is 50.",
             type=int,
             default=50)],
    'diff_levels_number':[
        '-N', '--difflevelsnumber',
        dict(help="number of levels for difference contours and shades.\
                   Default is 50.",
             type=int,
             default=50)],
    'colorminmax':[
        '-c', '--colorminmax',
        dict(help="Colors associated to min and max values. Syntax: 'color1, color2'.",
             default=None,
             dest='colorminmax')],
    'diffcolorminmax':[
        '--diffcolorminax',
        dict(help="Colors associated to min and max values for diff plot. \
                   Syntax: 'color1, color2'.",
             default=None,
             dest='diffcolorminmax')],
    'alpha':[
        '--alpha',
        dict(help="Alpha value (for contour, color and vectors plots)",
             default=1.,
             type=float,
             dest='alpha')],
    'alphaminmax':[
        '--alphaminmax',
        dict(help="Alpha values associated to min and max values \
                   for volume plot. \
                   Syntax: 'alpha1, alpha2'.",
             default=None,
             dest='alphaminmax')],
    'diffalphaminmax':[
        '--diffalphaminax',
        dict(help="Alpha values associated to min and max values \
                   for diff volume plot. \
                   Syntax: 'alpha1, alpha2'.",
             default=None,
             dest='diffalphaminmax')],
    'colormap':[
        '-c', '--colormap',
        dict(help="name of the **matplotlib** colormap to use.\
                   Default is 'jet'\
                   (Cf. http://matplotlib.org/examples/color/colormaps_reference.html \
                   for standard matplotlib colormaps, or epygram.config.epygram_colormaps.keys() \
                   in a python interpreter for epygram's own colormaps).\
                   Custom colormaps can be defined (http://colormap.org or\
                   manually) and added in userconfig,\
                   in usercolormaps = {'my_cmap':'path_to_my_cmap'}.",
             default='jet')],
    'diffcolormap':[
        '-C', '--diffcolormap',
        dict(help="name of the **matplotlib** colormap to use for diff.\
                   Default is 'RdBu_r'\
                   (Cf. http://matplotlib.org/examples/color/colormaps_reference.html \
                   for standard matplotlib colormaps, or epygram.config.epygram_colormaps.keys() \
                   in a python interpreter for epygram's own colormaps).\
                   Custom colormaps can be defined (http://colormap.org or\
                   manually) and added in userconfig,\
                   in usercolormaps = {'my_cmap':'path_to_my_cmap'}.",
             default='RdBu_r')],
    'center_cmap_on_0':[
        '-t', '--center_cmap_on_0',
        dict(action='store_true',
             help="to center the colormap on the value 0. Can be useful\
                   for wind plots for instance.",
             default=False)],
    'diff_center_cmap_on_0':[
        '-T', '--diffcenter_cmap_on_0',
        dict(action='store_false',
             help="NOT to center the colormap of diff plots on the value 0.\
                   May be useful for fluxes decumulation.",
             default=True)],
    'gis_quality':[
        '-q', '--gisquality',
        dict(help="quality of the GIS used to draw coastlines, rivers and\
                   countries; among ('c', 'l', 'i', 'h', 'f'), by increasing\
                   quality.",
             choices=['c', 'l', 'i', 'h', 'f'],
             default='i')],
    'draw_rivers':[
        '-r', '--drawrivers',
        dict(action='store_true',
             help="draw rivers on map. (Much more slow)",
             default=False)],
    'french_departments':[
        '--depts',
        dict(action='store_true',
             help="draw french departments on map (instead of countries\
                   boundaries).",
             default=False)],
    'parallels':[
        '--parallels',
        dict(help="""tune the choice of lines to plot, among:
                     'auto': automatic scaling to the map extents (default) |
                     'default': range(0,360,10) and range(-90,90,10) |
                     a list of values |
                     a grid step, e.g. 5 to plot each 5 degree |
                     None: no one is plot""",
             type=str,
             default='auto')],
    'meridians':[
        '--meridians',
        dict(help="Same as parallels, cf. parallels doc.",
             type=str,
             default='auto')],
    'hide_axes':[
        '--hide_axes',
        dict(help="To hide axe arrows.",
             action='store_true',
             default=False)],
    'vectors_subsampling':[
        '-s', '--vectors_subsampling',
        dict(help="Subsampling factor for plotting vectors barbs\
                   (-w: computewind option). Defaults to 20.",
             type=int,
             default=20)],
    'vectors_verticalsubsampling':[
        '--vectors_verticalsubsampling',
        dict(help="Vertical subsampling factor for plotting vectors barbs.\
                   Defaults to 1.",
             type=int,
             default=1)],
    'vectors_scale_factor':[
        '-S', '--vectors_scale_factor',
        dict(help="Scale factor to apply on vectors",
             type=float,
             default=1.)],
    'diffvectors_scale_factor':[
        '--diffvectors_scale_factor',
        dict(help="Scale factor to apply on diff vectors",
             type=float,
             default=1.)],
    'streamlines_time':[
        '-t', '--stream_time',
        dict(help="Integration time for stream lines or tubes",
             default=1,
             type=float,
             dest='streamlines_time')],
    'diffstreamlines_time':[
        '-T', '--diffstream_time',
        dict(help="Integration time for diff stream lines or tubes",
             default=None,
             dest='diffstreamlines_time')],
    'points_size':[
        '-p', '--pointsize',
        dict(help="size of points for *graphicmode* == 'points'.\
                   Defaults to 20.",
             type=int,
             default=20)],
    'lonlat_zoom':[
        '--zoom',
        dict(help="optional zoom on the specified region of the plot.\
                   Forces to 'cyl' projection. Syntax: 'lonmin=-5, lonmax=1.2,\
                   latmin=40.8, latmax=51'. Overwrites 'projection' option.",
             default=None)],
    'vertical_zoom':[
        '--zoom',
        dict(help="optional zoom (vertical reduction) on the profile plot. \
                   Ex: 'ymax=150, ymin=850'. The unit must be that of the \
                   vertical coordinate requested (hPa, m, level number).",
             default=None)],
    'z_factor':[
        '--z_factor',
        dict(help="factor to apply on z values (to modify aspect ratio of the plot).",
             type=float,
             default=None)],
    'spectra_zoom':[
        '--zoom',
        dict(help="optional zoom on the spectra plot. Ex: 'xmax=10, ymin=1.0'.",
             default=None)],
    'emagram_like_profiles':[
        '-e', '--emagramlike',
        dict(action='store_true',
             help="plots profiles in a emagram-like style (only with\
                   -P/--hybridP2pressure). Should not be used for other than\
                   Temperature profiles.",
             default=False)],
    'superpose_spectra_plots':[
        '-s', '--superposeplots',
        dict(action='store_true',
             help='for superposing spectra of all requested fields in one\
                   plot.',
             default=False)],
    'spectra_slopes':[
        '-k', '--kindofslopes',
        dict(help="optional kind of slopes to be plotted, with syntax:\
                   sequence of (exp offset label) between simple quotes ('),\
                   with slope = offset * k**exp.\
                   offset is optional, with default = 1.0.\
                   Label is optional and denotes how do the exp appears in\
                   legend (by default, exp=0.5 will appear 1/2 => add a label\
                   0.5 for it to appear unchanged.). If label is provided,\
                   offset must be provided as well.\
                   Ex: 'exp1, exp2 offset2 label2, exp3 offset3'.\
                   or  '-3 1.0, -5./3. 1e-2 -5/3' (default).",
             default='-3 1.0 -3, -5./3. 1e-2 -5/3')],
    'vector_symbol':[
        '--vs', '--vector_symbol',
        dict(help="symbol to be used for vectors, among ('barbs', 'arrows', 'stream').",
             dest='vector_symbol',
             default=_defaults['vector_symbol'])],
    'vector_plot_method':[
        '--vpm', '--vector_plot_method',
        dict(help="symbol to be used for vectors, among ('quiver', 'barbs', 'streamplot').",
             dest='vector_plot_method',
             default='quiver')],
    'quiverkey':[
        '--qk', '--quiverkey',
        dict(help="arguments to be passed to pyplot.quiverkey(), in case\
                   *vector_symbol* == 'arrows'. E.g. X=1.05,Y=1.05,U=10.,label='10m/s'.",
             dest='quiverkey',
             default=None)],
    'figures_dpi':[
        '--fd', '--figures_dpi',
        dict(help="quality of saved figures.",
             dest='figures_dpi',
             type=int,
             default=_defaults['default_figures_dpi'])],
    'resolution_increase':[
        '--resolution_increase',
        dict(help="Resolution increase factor",
             type=int,
             default=1)],
    'window_size':[
        '--window_size',
        dict(help="window size. Syntax: 'width, height'",
             default=None)],
    'bluemarble':[
        '--bluemarble',
        dict(help="displays NASA's \"blue marble\" as background, with\
                   a transparency set to the given value [0.0, 1.0].",
             type=float,
             default=0.0)],
    'ground':[
        '--ground',
        dict(help="'bluemarble' to display NASA's \"blue marble\" on ground\
                   or a url (using '${z}', '${x}' and '${y}' as place holders)\
                   to plot maptiles (use single quote around the url in the\
                   shall to prevent '$' to be expanded)",
             type=str,
             default=None)],
    'background':[
        '--background',
        dict(help="sets a background color to continents and oceans.",
             action='store_true',
             default=False)],
    'background_color':[
        '--bg', '--background_color',
        dict(help="backgound color.",
             type=str,
             default='Black',
             dest='background_color')],
    'section_abscissa':[
        '--section_abscissa', '--sa',
        dict(help="abscissa of section, among ('distance', 'lon', 'lat').",
             default='distance')],
    'global_shift_center':[
        '--global_shift_center', '--gsc',
        dict(help="for global lon/lat grids, shift the center by the requested \
                   angle (in degrees). Enables a [0,360] grid to be shifted to \
                   a [-180,180] grid, for instance (with -180 argument).",
             type=float,
             default=None)],
    'bins':[
        '-b', '--bins',
        dict(help="number of bins or bins edges (separated by commas). \
                   'range(0,100,12.5)' also works. Default is 50.",
             type=str,
             default='50')],
    'diffbins':[
        '-B', '--diffbins',
        dict(help="number of bins or bins edges (separated by commas) \
                   for diff plot. \
                   'range(0,100,12.5)' also works. Default is 50.",
             type=str,
             default='50')],
    'center_hist_on_0':[
        '-t', '--center_hist_on_0',
        dict(action='store_true',
             help="to center the histogram on the value 0. Can be useful\
                   for wind plots for instance.",
             default=False)],
    'diff_center_hist_on_0':[
        '-T', '--diffcenter_hist_on_0',
        dict(action='store_false',
             help="NOT to center the histogram of diff on the value 0.\
                   May be useful for fluxes decumulation.",
             default=True)],
    'focal_point':[
        '--focal_point',
        dict(type=str,
             help="Focal point: 'x,y,z'.",
             default=None)],
    'camera':[
        '--camera',
        dict(type=str,
             help="Camera position: 'x, y, z'.",
             default=None)],
    'scatter_kw':[
        '--skw', '--scatter_kw',
        dict(help="arguments to be passed to pyplot.scatter(), in case\
                   *plot_method* == 'scatter'.",
             dest='scatter_kw',
             default=None)],
    'cartopy_features':[
        '--cpyf', '--cartopy_features',
        dict(help="cartopy features (cartopy.feature.*), separated by comma",
             dest='cartopy_features',
             default=None)],
    }

