#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains a catalog of reusable command-line arguments for argparse,
with a function *add_arg_to_parser* to import them to an *argparse* parser.

The arguments are classified by categories, each category being a dict.

Each argument is a list composed as follows:
['name', 'alternate_name1', 'alternate_name2', ..., dict(kw1, kw2, kw3, ...)]
where the keywords 'kwi' are argparse.ArgumentParser.add_argument() optional
arguments.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from epygram import config, util
from epygram.extra import griberies
import multiprocessing


def add_arg_to_parser(parser, arg, **flychanges):
    """
    Wrapper to add one item *arg* of the following dictionaries to a *parser*.

    *flychanges* enable to change argument options on the fly.
    """

    arg[-1].update(flychanges.items())
    parser.add_argument(*arg[:-1], **arg[-1])

_defaults = {}
_defaults.update(config.__dict__)

# : Arguments dealing with files
files_management = {
    'principal_file':[
        'filename',
        dict(type=str,
             help='name of the file to be processed.')],
    'several_files':[
        'filenames',
        dict(type=str,
             nargs='+',
             help='names of the files to be processed.')],
    'file_to_refer_in_diff':[
        '-d', '--diff',
        dict(type=str,
             dest='refname',
             help='name of the 2nd (reference) file to be processed, to which\
                   comparison is done.')],
    'file_to_refer_in_diffonly':[
        '-D', '--diffonly',
        dict(type=str,
             dest='Drefname',
             help='same as -d/--diff, but only fields difference is processed.')],
    'source_file':[
        '-s', '--source',
        dict(type=str,
             dest='refname',
             help='name of the 2nd file from which to extract the fields.',
             required=True)],
    'replace_by_diff':[
        '-d', '--diff',
        dict(action='store_const',
             const='diff',
             dest='replace_op',
             help='replaces the fields by the difference fields\
                   (filename - source).',
             default=False)],
    'replace_by_reversediff':[
        '-r', '--reversediff',
        dict(action='store_const',
             const='reversediff',
             dest='replace_op',
             help='replaces the fields by the reverse difference fields\
                   (source - filename).',
             default=False)],
    'replace_by_addition':[
        '-a', '--add',
        dict(action='store_const',
             const='add',
             dest='replace_op',
             help='replaces the fields by the addition fields\
                   (filename + source).',
             default=False)],
    'replace_by_product':[
        '-m', '--multiply',
        dict(action='store_const',
             const='multiply',
             dest='replace_op',
             help='replaces the fields by the product fields\
                   (filename * source).',
             default=False)],
    'in_place':[
        '-i', '--in_place',
        dict(action='store_true',
             help='the operation on fields is done "in place" on the file,\
                   not on a new file.',
             default=False)],
                    }

# : Arguments dealing with fields
fields_management = {
    'field':[
        '-f', '-F', '--field',
        dict(help="*fid* = *field identifier* of of the field(s) to be\
                   processed. Syntax depends on format:\
                   GRIB: handgrip, e.g. 'shortName:t,level:850'.\
                   FA: name, e.g. 'S050TEMPERATURE'; regular expressions may\
                   be used, such as 'S00[2-6]WIND.[U-V].PHYS',\
                   'SURFALBEDO*' or 'SURF?.OF.OZONE'.\
                   To obtain the list of fields in file, use the 'epy_what'\
                   tool.",
             default=None)],
    'windfield':[
        '-w', '--computewind',
        dict(help="to process wind as a vector.\
                   Syntax depends on format:\
                   FA: 'S*WIND' or 'P*VENT', 'H*VENT', 'V*VENT', 'CLSVENT.*',\
                   or 'CLS*.RAF'; * designing the level: combination of\
                   (*, ?, digits).\
                   GRIB: handgrip, e.g. 'shortName:u+v,level:850' or\
                   'indicatorOfParameter:33+34,level:850' or \
                   'parameterCategory:2,parameterNumber:2+3,level:850'.\
                   Not implemented in difference mode.",
             default=None)],
    'windfieldU':[
        '--wU', '--Ucomponentofwind',
        dict(help="to process wind as a vector. U component of wind.\
                   (same syntax as -f).\
                   Not implemented in difference mode.",
             dest="Ucomponentofwind",
             default=None)],
    'windfieldV':[
        '--wV', '--Vcomponentofwind',
        dict(help="to process wind as a vector. V component of wind.\
                   (same syntax as -f).\
                   Not implemented in difference mode.",
             dest="Vcomponentofwind",
             default=None)],
    'FA_field':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   To obtain the list of fields in file, use the 'fa_what'\
                   tool.",
             default=None)],
    'FA_windfield':[
        '-w', '--computewind',
        dict(help="to process wind module, using the following syntax:\
                   'S*WIND' or 'P*VENT', 'H*VENT', 'V*VENT', 'CLSVENT.*', 'CLS*.RAF',\
                   * designing the level: combination of (*, ?, digits).\
                   Not implemented in difference mode.",
             default=None)],
    'FA_multiple_fields':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   Regular expressions can be used, such as\
                   'S00[2-6]WIND.[U-V].PHYS' to process meridian and zonal winds of levels 2 to 6,\
                   'SURFALBEDO*' to process all kinds of albedo,\
                   'SURF?.OF.OZONE' to process all kinds of ozone (A,B,C).\
                   In case this option is not used, no more than '-l' option,\
                   all fields present in file are processed.\
                   Fields not found in file are ignored.\
                   To obtain the list of fields in file, use the 'fa_what'\
                   tool.",
             default=None)],
    'vertical_field':[
        '-F', '-f', '--field',
        dict(help="designation of the fields from which to extract profile.\
                   Syntax depends on format:\
                   GRIB: handgrip designing the parameter and the type of\
                   levels.\
                   FA: name of the upper-air field to be processed,\
                   with a * for level.\
                   To obtain the list of fields in file, use the 'epy_what'\
                   tool.",
             required=True,
             default=None)],
    'FA_vertical_field':[
        '-F', '-f', '--field',
        dict(help="name of the upper-air field to be processed, with a * for\
                   level. To obtain the list of fields in file, use the\
                   'fa_what' tool.",
             required=True,
             default=None)],
    'DDHLFA_multiple_fields':[
        '-f', '-F', '--field',
        dict(help="name of the field to be processed.\
                   Regular expressions can be used, such as\
                   'F*FLUVERTDYN' or 'VKK[0-1]'.\
                   In case this option is not used, no more than '-l' option,\
                   all fields present in file are processed.\
                   Fields not found in file are ignored.\
                   To obtain the list of fields in file, use the 'ddhlfa_what'\
                   tool.",
             default=None)],
    'DDHLFA_domain':[
        '-n', '--domain_number',
        dict(help="number of the domain to plot. \
                   Multiple domains must be given within quotes, separated by comma. \
                   To obtain the list of domains in file, \
                   use the 'ddhlfa_what' tool.",
             required=True)],
    'GRIB_field':[
        '-F', '-f', '--field',
        dict(help="(possibly partial) handgrip of the field to be processed.\
                   To obtain the list of fields in file (and their handgrip), \
                   use EPyGrAM's 'grib_what' tool or GRIB_API's 'grib_dump'.\
                   Ex: -f 'shortName:t,level:850'.",
             default=None)],
    'GRIB_secondfield':[
        '-2', '--secondfield',
        dict(help="(possibly partial) handgrip of the second field to be\
                   processed. To obtain the list of fields in file (and their\
                   handgrip), use EPyGrAM's 'grib_what' tool or GRIB_API's\
                   'grib_dump'. Ex: -f 'shortName:t,level:850'.",
             default=None)],
    '2fields_mode':[
        '-y', '--2fields_mode',
        dict(help="mode for 2-fields joining:\
                   - 'vectorize': sets a vector with each field as a component\
                   - 'superpose': superpose second field over the first with\
                     contourlines\
                   - '-': computes and plots the difference field\
                   - '+': computes and plots the sum field\
                   - '*': computes and plots the product field\
                   - '/': computes and plots the division field",
             default=None)],
    'list_of_fields':[
        '-l', '--listoffields',
        dict(help="name of an external file containing the list of fields\
                   to be processed, with the same specifications as for -f\
                   argument.",
             default=None)],
    'reverse_fields_selection':[
        '-r', '--reverse',
        dict(action='store_true',
             help="reverse fields selection: all but those requested.",
             default=False)],
    'sort_fields':[
        '-s', '--sortfields',
        dict(action='store_true',
             help="not for GRIB: sort fields with regards to\
                   their name and type.",
             default=False)],
    'GRIB_sort':[
        '-k', '--sortbykey',
        dict(help="GRIB only: sort fields with regards to the given fid key.\
                   Only for *mode* == 'one+list' or 'fid_list'.",
             dest='sortfields',
             default=False)],
    'GRIB_what_mode':[
        '-m', '--mode',
        dict(help="information display mode (GRIB only): 'one+list' gives the \
                   validity/geometry of the first field in GRIB, plus \
                   the list of fid; 'fid_list' gives only the fid of \
                   each field in GRIB; 'what' gives the values of the \
                   keys from each GRIB message that are used to generate \
                   an **epygram** field from the message (slower); \
                   'ls' gives the values of the 'ls' keys from each GRIB \
                   message; 'mars' gives the values of the 'mars' keys \
                   from each GRIB message.",
             choices=('one+list', 'fid_list', 'what', 'ls', 'mars'),
             default='one+list')],
    'FA_set_compression':[
        '-k', '--kompression',
        dict(help="bits number for FA gridpoint compression : KNBITS/KNBPDG.\
                   Defaults to " +
                   str(_defaults.get('FA_default_compression',
                                     {'KNBPDG':('undefined')}).get('KNBPDG', 'undefined')) + ".",
             type=int,
             default=None)],
    'GRIB2_packing':[
        '--GRIB2_packing',
        dict(help="packing for writing GRIB2s.\
                   Defaults to " +
                   str(griberies.defaults.GRIB2_keyvalue[5]) + ".",
             default=None)],
    'netCDF_compression':[
        '--nc_comp',
        dict(type=int,
             help="compression level to compress fields in netCDF. \
                   Ranges from 0 to 9. 0 is no compression, 1 is low-but-fast \
                   compression, 9 is high-but-slow compression. \
                   Default to " + str(_defaults.get('netCDF_default_compression',
                                                    "4")),
             choices=list(range(0, 10)),
             default=str(_defaults.get('netCDF_default_compression', 4)))],
                    }

# : Arguments dealing with output
output_options = {
    'output':[
        '-o', '--output',
        dict(help="store graphical output in file in with specified format,\
                   among ('png', pdf). Pdf is kind of disadvised, for it is\
                   very slow and produces much too big files. 'X' to open GUI.",
             choices=['png', 'pdf', 'X'],
             default=_defaults.get('default_graphical_output'))],
    'outputfmt':[
        '-o', '--outputfmt',
        dict(help="specify format for output file,\
                   among ('png', pdf). Pdf is kind of disadvised, for it is\
                   very slow and produces much too big files. 'X' to open GUI.",
             choices=['png', 'pdf', 'X'],
             default=_defaults.get('default_graphical_output'))],
    'savefig':[
        '--sf', '--savefig',
        dict(help="save figure in file.",
             action='store_true',
             dest='savefig',
             default=False)],
    'outputfilename':[
        '-O', '--outputfilename',
        dict(help="store output in the specified filename (without format for\
                   graphical output, to be completed by -o/--output).",
             default=False)],
    'one_pdf':[
        '-p', '--pdf',
        dict(action='store_true',
             help="store output (eventually multiple plots) in one .pdf file.",
             default=False)],
    'noplot':[
        '-n', '--noplot',
        dict(action='store_true',
             help="disable plot. Profile/spectrum will be computed\
                   and written as text output but not plotted.",
             default=False)],
    'GeoPoints_llv':[
        '--llv',
        dict(action='store_true',
             help="simplify GeoPoints format to LAT, LON, VALUE only.",
             default=False)],
    'precision':[
        '-e', '--precision',
        dict(help="precision on values (number of decimals in scientific\
                   format). Default = 4.",
             default=4)],
    'lonlat_precision':[
        '-E', '--lonlat_precision',
        dict(help="precision on longitudes/latitudes (number of decimals). \
                   Default = 4.",
             default=4)],
    'GeoPoints_precision':[
        '-e', '--precision',
        dict(help="precision on GeoPoints' values (number of decimals in \
                   scientific format). Default = " +
                   str(_defaults.get('GeoPoints_precision', 4)) + '.',
             default=_defaults.get('GeoPoints_precision', 4))],
    'GeoPoints_lonlat_precision':[
        '-E', '--lonlat_precision',
        dict(help="precision on GeoPoints' longitudes/latitudes (number of \
                   decimals). Default = " +
                   str(_defaults.get('GeoPoints_lonlat_precision', 4)) + '.',
             default=_defaults.get('GeoPoints_precision', 4))],
    'get_field_compression':[
        '-c', '--compression',
        dict(action='store_true',
             help="get compression options of each field.",
             default=False)],
    'get_field_details':[
        '-d', '--details',
        dict(help="get some details about each field.\
                   E.g. 'spectral' to inquire the spectralness of fields, or\
                   'compression' to inquire compression parameters, 'grid' or \
                   'comment' for LFI.",
             default=None)],
    'GRIB_short_fid':[
        '-s', '--grib_short_fid',
        dict(help="condense GRIB fid.",
             action='store_true',
             default=False)],
    'GRIB_other_options':[
        '--GRIB_other_options',
        dict(help="a set of key:value pairs, separated by commas, to be given to\
                   the GRIB message when writing.\
                   Ex: 'typeOfGeneratingProcess:12,productionStatusOfProcessedData:2'.",
             default=None)],
    'stdout':[
        '-o', '--stdout',
        dict(action='store_true',
             help="redirects output to standard output (rather than file).",
             default=False)],
    'output_format':[
        '-o', '--output_format',
        dict(help="format of conversion output, among:\
                   'grb' (GRIB2), 'geo' (GeoPoints), 'nc' (netCDF4).",
             choices=['grb', 'geo', 'nc'],
             required=True)],
    'reproject_wind':[
        '--rw', '--reproject_wind',
        dict(action='store_true',
             help="reprojects a wind vector (u, v) onto true\
                   zonal/meridian axes (assuming it is projected onto grid axes\
                   in fields).",
             dest='reproject_wind',
             default=False)],
    'only_maxdiff':[
        '--oxd', '--only_maxdiff',
        dict(action='store_true',
             help='if diffonly and activated, only print the maximum absolute\
                   difference (and min/max values of fields for magnitude',
             dest='only_maxdiff',
             default=False)],
    'GeoPoints_columns':[
        '--geopoints_cols',
        dict(dest='geopoints_cols',
             help="set manually the columns of GeoPoints, for instance 'LAT,LON,VALUE'.",
             default=None)],
    'GeoPoints_noheader':[
        '--geopoints_noheader',
        dict(action='store_true',
             dest='geopoints_noheader',
             help="discard the header of GeoPoints.",
             default=False)],
                  }

# : Miscellaneous arguments
misc_options = {
    'LAMzone':[
        '-z', '--zone',
        dict(help="zone of the domain, in LAM case:\
                   'C' for zone C, 'CI' for zone C+I, 'CIE' for zone C+I+E.\
                   Default is 'CI'.",
             choices=['C', 'CI', 'CIE'],
             default='CI')],
    'array_flattening_order':[
        '--order',
        dict(help="for LAM arrays, whether to flatten in C (row-major) or\
                   Fortran (column-major) order 2D arrays. Default = 'C'.",
             default='C')],
    'flatten_horizontal_grids':[
        '--flatten',
        dict(help="for netCDF, flatten 2D horizontal grids to 1D.",
             action='store_true',
             default=False)],
    'operation_on_field':[
        '-x', '--operation',
        dict(help="do the requested operation on field right after reading it. \
                   Syntax: '-,273.15' (e.g. for K => C) or 'exp' to take the \
                   exponential of the field. \
                   The operand must be among (+,-,*,/) or 'normalize' or any\
                   numpy function.\
                   For '-' operand, use short-name option -x without spacetab\
                   between option and argument.",
             default=None)],
    'diffoperation_on_field':[
        '-X', '--diffoperation',
        dict(help="do the requested operation on difference field right after \
                   computing it. \
                   Syntax: '-,273.15' (e.g. for K => C) or 'exp' to take the \
                   exponential of the field. \
                   The operand must be among (+,-,*,/) or 'normalize' or any\
                   numpy function.\
                   For '-' operand, use short-name option -X without spacetab\
                   between option and argument.",
             default=None)],
    'pressure_unit_hpa':[
        '-u', '--pressure_unit_hpa',
        dict(action='store_true',
             help="converts pressure in Pa to hPa.",
             default=False)],
    'composition_with_field':[
        '--compose_with',
        dict(help="compose a transformed field with another field. \
                   Syntax: 'otherfield, composition, [file=otherfile], [preset='norm'|'ceil']'. \
                   The *composition* must be among (+,-,*,/).\
                   If *file* is filled (optional), the *otherfield* is read in it.\
                   If *preset* is filled (optional), one or several operation can be done \
                   on the *otherfield* before operation.\
                   Ex: 'SFX.FRAC_WATER, *, file=mypgd.fa, preset=ceil' \
                   will set non-water points of the plotted field (e.g. SFX.TS_WATER) to 0.",
             default=None)],
    'mask_threshold':[
        '--mt', '--mask_threshold',
        dict(help="set a threshold to mask values. E.g. 'min=0.0' will\
                   mask negative values. 'min=0.0,max=1e8' will mask values\
                   outside these boundaries.",
             default=None,
             dest='mask_threshold')],
    'femars.diff_to_avg':[
        '--avg', '--diff_to_avg',
        dict(help="instead of computing member to member differences, \
                   compute member-to-average-of-members differences.",
             action='store_true',
             default=False,
             dest='diff_to_avg')],
    'map_factor_correction':[
        '--mpf', '--map_factor_correction',
        dict(help="apply a map factor correction to wind.",
             action='store_true',
             default=False,
             dest='map_factor_correction')],
    'wind_components_are_projected_on':[
        '--wpo', '--wind_projected_on',
        dict(help="specify on which coordinates the wind components are projected.",
             default=None,
             choices=['lonlat', 'grid'],
             dest='wind_components_are_projected_on')],
                }

# : Arguments dealing with graphical options
graphical_options = {
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

# : Arguments dealing with extraction stuff
extraction_options = {
    'point_coordinates':[
        '-c', '--coordinates',
        dict(help="lon/lat coordinates of the point.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'section_starting_point':[
        '-s', '--starting_point',
        dict(help="lon/lat coordinate of starting point of the section.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'section_ending_point':[
        '-e', '--ending_point',
        dict(help="lon/lat coordinate of ending point of the section.\
                   Syntax: 'lon, lat'.",
             required=True,
             default=None)],
    'horizontal_interpolation':[
        '-i', '--interpolation',
        dict(help="interpolation mode from field grid to point/section\
                   coordinates. Among ('nearest', 'linear', 'cubic',\
                   'bilinear'). Defaults to 'nearest'.",
             choices=['nearest', 'linear', 'cubic', 'bilinear'],
             default='nearest')],
    'external_distance':[
        '-z', '--external_distance',
        dict(help="for 'nearest' interpolation mode, the nearest is chosen\
                   among the 4 nearest points regarding the point whose value\
                   of field *EXT* is the closest to *VALUE*, with syntax:\
                   *external_distance* = 'VALUE; EXT'. EXT being a fid.",
             default=None)],
    'section_transect_points_number':[
        '-p', '--points_number',
        dict(help="number of points from starting point to ending point\
                   (included). Defaults to the number of points computed from\
                   the fields resolution, or from the resolution given via\
                   option -r.",
             type=int,
             default=None)],
    'section_transect_resolution':[
        '-r', '--resolution',
        dict(help="resolution of the section. Defaults to the fields\
                   resolution, or computed from the number of points given via\
                   option -p.",
             type=float,
             default=None)],
    'verticalcoord2pressure':[
        '-P', '--verticalcoord2pressure',
        dict(action='store_const',
             dest='Yconvert',
             const='pressure',
             help="compute pressure as vertical coordinate, if possible.",
             default=None)],
    'verticalcoord2height':[
        '-H', '--verticalcoord2height',
        dict(action='store_const',
             dest='Yconvert',
             const='height',
             help="compute height as vertical coordinate, if possible.",
             default=None)],
    'verticalcoord2altitude':[
        '-A', '--verticalcoord2altitude',
        dict(action='store_const',
             dest='Yconvert',
             const='altitude',
             help="compute altitude as vertical coordinate, if possible.",
             default=None)],
    'no_cheap_height_conversion':[
        '--no_cheap_height',
        dict(action='store_false',
             dest='cheap_height',
             help="for the computation of heights (-A/-H) to be done \
                   taking hydrometeors into account (in R \
                   computation) and NH Pressure departure \
                   (Non-Hydrostatic data). Slower but more accurate.",
             default=True)]
                      }

# : Arguments dealing with runtime options
runtime_options = {
    'verbose':[
        '-v', '--verbose',
        dict(action='store_true',
             help="run verbosely. Else, only messages of level Error will be\
                   displayed.",
             default=False)],
    'percentage':[
        '-p', '--percentage',
        dict(action='store_true',
             help="display the percentage done on the run.",
             default=False)],
    'threads_number':[
        '-t', '--threads_number',
        dict(help="number of threads to be run in parallel.",
             type=int,
             default=multiprocessing.cpu_count() // 2)],
                    }

# : Operational arguments
operational_options = {
    'suite':[
        '-S', '--suite',
        dict(help="name of the suite, among (oper, dble, research, test).\
                   Defaults to research. Used to build GRIB's \
                   *productionStatusOfProcessedData*.",
             default='research')],
    'typeOfGeneratingProcess':[
        '-g', '--typeOfGeneratingProcess',
        dict(help="GRIB's type of generating process.",
             choices=list(griberies.tables.typeOfGeneratingProcess_dict.keys()),
             default='Forecast')],
    'numod':[
        '-N', '--NUMOD',
        dict(help="model identifier (known as NUMOD at Meteo-France). \
                   A.k.a. 'generatingProcessIdentifier' in GRIB_API. \
                   Default is 255.",
             dest='numod',
             default=None,
             type=int)],
                       }

# : Arguments for domain_maker
domain_maker_options = {
    'mode':[
        '-l',
        dict(dest='mode',
             action='store_const',
             help="mode: define domain by specifying a lon/lat domain that must\
                   be included inside.",
             const='lonlat_included',
             default='center_dims')],
    'no_display':[
        '-n', '--no_display',
        dict(action='store_true',
             help="run without displaying the domain.",
             default=False)],
    'maximize_CI_in_E':[
        '-m', '--maximize_CI_in_E',
        dict(action='store_true',
             help="forces the C+I zone to be the greatest possible inside a\
                   given (discrete) C+I+E size. \
                   In other words, with a given (discrete) C+I+E size, forces\
                   the E-zone to be the smallest possible.",
             default=False)],
    'truncation':[
        '-t', '--truncation',
        dict(help="(for namelists) the kind of truncation of spectral geometry\
                   to generate, among ('linear', 'quadratic', 'cubic').",
             default='linear')],
    'orography_subtruncation':[
        '-o', '--orography_subtruncation',
        dict(help="(for namelists) additional subtruncation for orography to be generated.",
             default='quadratic')],
                        }
