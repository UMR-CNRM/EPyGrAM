#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import argparse
import time
import copy
from distutils.version import LooseVersion

import footprints
from footprints import FPDict, FPList
import taylorism
assert LooseVersion(taylorism.__version__) >= LooseVersion('1.0.7'), \
    'Version of package taylorism (here {}) need to be >= 1.0.7'.format(taylorism.__version__)
from bronx.fancies.display import printstatus
from bronx.syntax.parsing import str2dict

import epygram
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  misc_options, runtime_options,
                                  operational_options, output_options)
import griberies

epylog = footprints.loggers.getLogger(__name__)
fmt_dict = {'grb':'GRIB', 'nc':'netCDF', 'geo':'GeoPoints'}


class Converter(taylorism.Worker):

    _abstract = True
    _footprint = dict(
        info="Converts a file into another format.",
        attr=dict(
            output_format_suffix=dict(),
            filename=dict(
                info='name of the file to be processed',
                access='rwx'),
            fieldseed=dict(
                info="to select only a part of the fields in resource.",
                type=FPList,
                optional=True,
                default=None,
                access='rwx'),
            subzone=dict(
                info="to extract a subzone of LAM fields, in case, among ('C', 'CI').",
                optional=True,
                default=None,
                access='rwx'),
            progressmode=dict(
                info="to print either progress 'percentage' of job, \
                      or field fids ('verbose').",
                optional=True,
                default=None)
        )
    )


def convert(filename,
            output_format_suffix,
            get_write_kwargs=lambda *args:{},
            fieldseed=None,
            subzone=None,
            grib_short_fid=False,
            progressmode=None,
            **kwargs):
    """
    Conversion function.

    :param filename: name of the file to be processed.
    :param output_format_suffix: among 'grb' (GRIB2), 'nc' (netCDF4), 'geo' (GeoPoints)
    :param get_write_kwargs: function that gives the necessary write options
    :param fieldseed: either a fid or a list of fid, used as a seed for
                   generating the list of fields to be processed.
    :param subzone: LAM zone among ('C', 'CI', None).
    :param grib_short_fid: condense GRIB fid as string, in case converting a GRIB
                        file.
    :param progressmode: among ('verbose', 'percentage', None)

    Other kwargs are specific to output formats, and passed to
                     get_write_kwargs(), wherein they can be handled whatever for...
    """
    t0 = time.time()

    if output_format_suffix == 'geo':
        one_output_file = False
    else:
        one_output_file = True
    resource = epygram.formats.resource(filename, openmode='r')
    if one_output_file:
        output_resource = epygram.formats.resource('.'.join([filename, output_format_suffix]),
                                                   openmode='w',
                                                   fmt=fmt_dict[output_format_suffix])
    # warnings about formats
    if resource.format not in ('GRIB', 'FA'):
        epylog.warning(" ".join(["tool NOT TESTED with format",
                                 resource.format, "!"]))
    if output_format_suffix not in ('grb', 'nc', 'geo'):
        epylog.warning(" ".join(["tool NOT TESTED with output format",
                                 fmt_dict[output_format_suffix], "!"]))
    fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
    numfields = len(fidlist)
    if progressmode == 'percentage':
        n = 0
        printstatus(n, numfields)
        n += 1
    # loop over fields
    for f in fidlist:
        if progressmode == 'verbose':
            epylog.info(str(f))
        field = resource.readfield(f)
        if not field.geometry.grid.get('LAMzone', False):
            subzone = None
        # options for write
        write_kwargs = get_write_kwargs(resource, field, kwargs)
        if field.spectral:
            field.sp2gp()
        if subzone is not None:
            field.select_subzone(subzone)
        # output fid
        if resource.format == 'GRIB':
            grib_edition = field.fid[list(field.fid.keys())[0]]['editionNumber']
            fid = field.fid['GRIB' + str(grib_edition)]
        else:
            fid = field.fid[resource.format]
        if output_format_suffix == 'grb':
            if 'generic' not in field.fid:
                raise NotImplementedError("how to convert this fid to GRIB2 ?")
            field.fid['GRIB2'] = epygram.formats.fid_converter(field.fid['generic'], 'generic', 'GRIB2')
            assert all([k in field.fid['GRIB2']
                        for k in ('discipline', 'parameterCategory', 'parameterNumber')]), \
                   "missing key(s) among ('discipline', 'parameterCategory', 'parameterNumber')"
        else:
            try:
                field.fid[fmt_dict[output_format_suffix]] = epygram.formats.fid_converter(fid,
                                                                                          resource.format,
                                                                                          fmt_dict[output_format_suffix],
                                                                                          grib_short_fid=grib_short_fid)
            except NotImplementedError:
                field.fid[fmt_dict[output_format_suffix]] = str(field.fid[resource.format])
        # now write field !
        if not one_output_file:
            if output_format_suffix == 'geo':
                output_resource = epygram.formats.resource('.'.join([filename,
                                                                     field.fid[fmt_dict[output_format_suffix]],
                                                                     output_format_suffix]),
                                                           openmode='w',
                                                           fmt=fmt_dict[output_format_suffix],
                                                           other_attributes=({'FORMAT':'XYV'} if kwargs.get('llv') else None),
                                                           columns=kwargs.get('columns'),
                                                           no_header=kwargs.get('no_header', False)
                                                           )
            else:
                output_resource = epygram.formats.resource('.'.join([filename,
                                                                     field.fid[fmt_dict[output_format_suffix]],
                                                                     output_format_suffix]),
                                                           openmode='w',
                                                           fmt=fmt_dict[output_format_suffix])
        output_resource.writefield(field, **write_kwargs)
        # print progress if requested
        if progressmode == 'percentage':
            printstatus(n, numfields)
            n += 1

    t1 = time.time()
    return ' '.join([filename, ":",
                     str(numfields), "fields successfully converted from",
                     resource.format, "to", fmt_dict[output_format_suffix],
                     "in", str(t1 - t0), "s."])


# GRIB2 writer #
################
class Griber(Converter):

    _footprint = dict(
        info="Converts a file into GRIB2, with some MF operations rules to \
              define 'generatingProcessIdentifier'.",
        attr=dict(
            output_format_suffix=dict(
                values=['grb', ]),
            suite=dict(
                info="to specify grib_api's 'productionStatusOfProcessedData'. \
                      Among ('oper', 'dble', 'test', 'research', 'unknown').",
                optional=True,
                default='research',
                values=set(griberies.tables.productionStatusOfProcessedData_dict.keys()),
                access='rwx'),
            typeofgeneratingprocess=dict(
                info="cf. epygram.formats.griberies.tables.typeOfGeneratingProcess_dict.",
                optional=True,
                default='Forecast',
                values=set(griberies.tables.typeOfGeneratingProcess_dict.keys()),
                access='rwx'),
            default_packing=dict(
                info="to specify default packing of fields in GRIB target.",
                type=FPDict,
                optional=True,
                default=FPDict(griberies.defaults.GRIB2_keyvalue[5]),
                access='rwx'),
            specific_packing=dict(
                info="to specify packing of specific fields, with regards to \
                      their name.",
                type=FPDict,
                optional=True,
                default=FPDict({}),
                access='rwx'),
            other_grib_options=dict(
                info="to specify other GRIB options (key/value pairs).",
                type=FPDict,
                optional=True,
                default=FPDict({}),
                access='rwx'),
        )
    )

    def _task(self):
        def get_write_kwargs_for_GRIB(input_resource, field, kwargs):
            write_kwargs = {}
            packing = kwargs['specific_packing'].get(field.fid[input_resource.format],
                                                     copy.copy(kwargs['default_packing']))
            if input_resource.format == 'FA':
                if field.spectral:
                    packing['bitsPerValue'] = input_resource.fieldscompression[field.fid[input_resource.format]]['KNBCSP']
                else:
                    if input_resource.fieldscompression[field.fid[input_resource.format]]['KNGRIB'] != 0:
                        # compressed field
                        packing['bitsPerValue'] = input_resource.fieldscompression[field.fid[input_resource.format]]['KNBPDG']
                    else:
                        # not compressed field : no packing !
                        packing = {'packingType':'grid_ieee'}
                if packing.get('bitsPerValue') == 0:
                    packing.pop('bitsPerValue')
            write_kwargs['packing'] = packing
            write_kwargs['other_GRIB_options'] = kwargs.get('other_GRIB_options', {})
            write_kwargs['grib_edition'] = 2
            return write_kwargs

        # generate additional keys
        if not isinstance(self.typeofgeneratingprocess, int):
            typeOfGeneratingProcess = griberies.tables.typeOfGeneratingProcess_dict[self.typeofgeneratingprocess]
        other_GRIB_options = {'generatingProcessIdentifier':255,
                              'productionStatusOfProcessedData':griberies.tables.productionStatusOfProcessedData_dict[self.suite],
                              'typeOfGeneratingProcess':typeOfGeneratingProcess}
        other_GRIB_options.update(self.other_grib_options)
        # convert
        done = convert(self.filename,
                       self.output_format_suffix,
                       get_write_kwargs_for_GRIB,
                       fieldseed=self.fieldseed,
                       subzone=self.subzone,
                       progressmode=self.progressmode,
                       # format specific
                       default_packing=self.default_packing,
                       specific_packing=self.specific_packing,
                       other_GRIB_options=other_GRIB_options)

        return done


# netCDF writer #
#################
class NetCDFWriter(Converter):

    _footprint = dict(
        info="Converts a file into netCDF4.",
        attr=dict(
            output_format_suffix=dict(
                values=['nc', ]),
            grib_short_fid=dict(
                info="condense GRIB fid as string, in case converting a GRIB file.",
                optional=True,
                default=False,
                type=bool,
                access='rwx'),
            compression=dict(
                info="compression level to compress fields in netCDF. \
                      Ranges from 0 to 9. 0 is no compression, 1 is low-but-fast \
                      compression, 9 is high-but-slow compression. \
                      Default is 4.",
                type=int,
                optional=True,
                default=epygram.config.netCDF_default_compression,
                access='rwx'),
            flatten_horizontal_grids=dict(
                info="flatten 2D horizontal grids to 1D.",
                optional=True,
                type=bool,
                default=False),
        )
    )

    def _task(self):

        def get_write_kwargs_for_netCDF(input_resource, field, kwargs):
            write_kwargs = {}
            write_kwargs['compression'] = kwargs['compression']
            write_kwargs['adhoc_behaviour'] = {'flatten_horizontal_grids':kwargs['flatten_horizontal_grids']}
            return write_kwargs

        done = convert(self.filename,
                       self.output_format_suffix,
                       get_write_kwargs_for_netCDF,
                       fieldseed=self.fieldseed,
                       subzone=self.subzone,
                       progressmode=self.progressmode,
                       # format specific
                       grib_short_fid=self.grib_short_fid,
                       compression=self.compression,
                       flatten_horizontal_grids=self.flatten_horizontal_grids)

        return done


# GeoPoints writer #
####################
class GeoPointsWriter(Converter):

    _footprint = dict(
        info="Converts a file into GeoPoints.",
        attr=dict(
            output_format_suffix=dict(
                values=['geo', ]),
            grib_short_fid=dict(
                info="condense GRIB fid as string, in case converting a GRIB file.",
                optional=True,
                default=False,
                type=bool,
                access='rwx'),
            order=dict(
                info="for a rectangular H2DField, whether to flatten 2D arrays \
                      in 'C' (row-major) or 'F' (Fortran, column-major) order.",
                optional=True,
                default='C',
                access='rwx'),
            lonlat_precision=dict(
                info="number of digits for output lon/lat values.",
                optional=True,
                type=int,
                default=epygram.config.GeoPoints_lonlat_precision,
                access='rwx'),
            precision=dict(
                info="number of digits for output fields values.",
                optional=True,
                type=int,
                default=epygram.config.GeoPoints_precision,
                access='rwx'),
            llv=dict(
                info="lat/lon/value mode == simplified GeoPoints.",
                optional=True,
                type=bool,
                default=False,
                access='rwx'),
            columns=dict(
                info="simplified GeoPoints, specify columns to be written.",
                optional=True,
                type=FPList,
                default=None,
                access='rwx'),
            no_header=dict(
                info="If True, do not write header (openmode='w').",
                optional=True,
                type=bool,
                default=False)
        )
    )

    def _task(self):

        def get_write_kwargs_for_GeoPoints(input_resource, field, kwargs):
            write_kwargs = {}
            write_kwargs['order'] = kwargs.get('order', 'C')
            if 'lonlat_precision' in kwargs:
                write_kwargs['llprecision'] = kwargs['lonlat_precision']
            if 'precision' in kwargs:
                write_kwargs['precision'] = kwargs['precision']
            return write_kwargs

        done = convert(self.filename,
                       self.output_format_suffix,
                       get_write_kwargs_for_GeoPoints,
                       fieldseed=self.fieldseed,
                       subzone=self.subzone,
                       grib_short_fid=self.grib_short_fid,
                       progressmode=self.progressmode,
                       # format specific
                       order=self.order,
                       precision=self.precision,
                       lonlat_precision=self.lonlat_precision,
                       llv=self.llv,
                       columns=self.columns,
                       no_header=self.no_header)

        return done


################################################################################
# MAIN
def main(filenames,
         output_format_suffix,
         threads_number=1,
         progressmode=None,
         **kwargs):
    """
    Converts a series of files to *output_format_suffix*.

    Mandatory arguments:

    :param filenames: name(s) of the files to be processed
    :param output_format_suffix: among 'grb' (GRIB2), 'nc' (netCDF4), 'geo' (GeoPoints)

    Technical named (optional) arguments:

    :param threads_number: parallelisation of files processing
    :param progressmode: among ('verbose', 'percentage', None)

    Other named arguments depend on the output format, and are defined in the
    Workers footprints attributes !
    """
    # build a dummy Converter of the right type
    dummy_converter = footprints.proxy.worker(filename='', output_format_suffix=output_format_suffix)
    instructions_keys = dummy_converter.footprint_attributes
    del dummy_converter
    # build instructions
    common_instructions = {k:v for k, v in kwargs.items() if k in instructions_keys}
    # standard instructions for converters
    common_instructions['output_format_suffix'] = output_format_suffix
    if threads_number == 1:
        common_instructions['progressmode'] = progressmode
    if isinstance(filenames, six.string_types):
        filenames = [filenames]
    individual_instructions = {'filename':filenames}

    # run !
    taylorism.batch_main(common_instructions, individual_instructions,
                         scheduler=footprints.proxy.scheduler(limit='threads', max_threads=threads_number),
                         verbose=(progressmode == 'verbose'))
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description='An EPyGrAM tool for converting file formats. \
                                                  Spectral fields are converted into gridpoints.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['several_files'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['field'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, misc_options['LAMzone'])
    add_arg_to_parser(parser, output_options['output_format'])
    # compression/precision options
    add_arg_to_parser(parser, fields_management['netCDF_compression'])
    add_arg_to_parser(parser, fields_management['GRIB2_packing'])
    add_arg_to_parser(parser, output_options['GeoPoints_lonlat_precision'])
    add_arg_to_parser(parser, output_options['GeoPoints_precision'])
    # GRIB specifics
    add_arg_to_parser(parser, operational_options['suite'])
    add_arg_to_parser(parser, operational_options['typeOfGeneratingProcess'])
    add_arg_to_parser(parser, operational_options['numod'])
    add_arg_to_parser(parser, output_options['GRIB_other_options'])
    # GeoPoints specific
    add_arg_to_parser(parser, misc_options['array_flattening_order'])
    cols = parser.add_mutually_exclusive_group()
    add_arg_to_parser(cols, output_options['GeoPoints_llv'])
    add_arg_to_parser(cols, output_options['GeoPoints_columns'])
    add_arg_to_parser(parser, output_options['GeoPoints_noheader'])
    # netCDF
    add_arg_to_parser(parser, misc_options['flatten_horizontal_grids'])
    # others
    add_arg_to_parser(parser, output_options['GRIB_short_fid'])
    add_arg_to_parser(parser, runtime_options['threads_number'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_options['verbose'])
    add_arg_to_parser(status, runtime_options['percentage'])

    args = parser.parse_args()

    # 2. Initializations
    ####################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')
    epylog.info("Start.")

    # 2.1 options
    if args.zone in ('C', 'CI'):
        subzone = args.zone
    else:
        subzone = None
    if args.GRIB2_packing is not None:
        default_packing = str2dict(args.GRIB2_packing, try_convert=int)
    else:
        default_packing = griberies.defaults.GRIB2_keyvalue[5]
    if args.GRIB_other_options is not None:
        other_GRIB_options = str2dict(args.GRIB_other_options, try_convert=int)
    else:
        other_GRIB_options = {}
    if args.numod is not None:
        other_GRIB_options['generatingProcessIdentifier'] = args.numod
    if args.geopoints_cols is not None:
        args.geopoints_cols = [c.strip() for c in args.geopoints_cols.split(',')]
    assert args.filenames != [], \
           "must supply one or several filenames."
    threads_number = min(args.threads_number, len(args.filenames))
    progressmode = None
    if args.verbose:
        progressmode = 'verbose'
    elif args.percentage:
        if threads_number > 1:
            progressmode = None
        else:
            progressmode = 'percentage'

    # 2.2 list of fields to be processed
    if args.field is not None:
        fieldseed = [args.field]
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as listfile:
            fieldseed = listfile.readlines()
        for n in range(len(fieldseed)):
            fieldseed[n] = fieldseed[n].replace('\n', '').strip()
    else:
        fieldseed = None

    # 3. Main
    #########
    main([six.u(f) for f in args.filenames],
         args.output_format,
         # technical
         threads_number=threads_number,
         progressmode=progressmode,
         # instructions common to formats
         fieldseed=fieldseed,
         subzone=subzone,
         grib_short_fid=args.grib_short_fid,
         # GRIB specifics
         suite=args.suite,
         typeofgeneratingprocess=args.typeOfGeneratingProcess,
         default_packing=default_packing,
         other_grib_options=other_GRIB_options,
         # netCDF specifics
         compression=args.nc_comp,
         flatten_horizontal_grids=args.flatten,
         # GeoPoints specifics
         order=args.order,
         lonlat_precision=args.lonlat_precision,
         precision=args.precision,
         llv=args.llv,
         columns=args.geopoints_cols,
         no_header=args.geopoints_noheader
         )
