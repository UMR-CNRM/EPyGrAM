#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Converters: workers taylored to write in a dedicated format.
"""
import time
import copy

import footprints
from footprints import FPDict, FPList
import taylorism
from bronx.fancies.display import printstatus

from ... import config
from ...extra import griberies
from .functions import convert

__all__ = ['Griber', 'NetCDFWriter', 'GeoPointsWriter']


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
                info="cf. epygram.extra.griberies.tables.typeOfGeneratingProcess_dict.",
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
                default=config.netCDF_default_compression,
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
                       grib_short_fid=self.grib_short_fid,  # TODO: use griberies cfName ?
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
                default=config.GeoPoints_lonlat_precision,
                access='rwx'),
            precision=dict(
                info="number of digits for output fields values.",
                optional=True,
                type=int,
                default=config.GeoPoints_precision,
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

