#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Conversion functions (API).
"""
import time

import footprints
from bronx.fancies.display import printstatus

from .. import resource as eopen

epylog = footprints.loggers.getLogger(__name__)
fmt_dict = {'grb':'GRIB', 'nc':'netCDF', 'geo':'GeoPoints'}


def fid_converter(initial_fid,
                  initial_fmt,
                  target_fmt,
                  grib_short_fid=False):
    """
    Creates and returns the fid in format *target_fmt* from an *initial_fid* in
    *initial_fmt*.

    *grib_short_fid* condense GRIB fid in string.
    """
    if initial_fmt == 'generic' and target_fmt == 'GRIB2':
        target_fid = copy.copy(initial_fid)
    elif initial_fmt == 'GRIB' and target_fmt in ('netCDF', 'GeoPoints'):
        # TODO: ? this is very basic !
        if grib_short_fid:
            target_fid = '-'.join([str(initial_fid[k])
                                   for k in sorted(initial_fid.keys())])
        else:
            target_fid = str(initial_fid).replace(' ', '').replace("'", "").replace("{", "_")
            """
            FIXME: doesn't work
            try:
                from .GRIB import namesgribdef
                fid = copy.copy(initial_fid)
                fid.pop('name', None)
                fid.pop('shortName', None)
                fid.pop('editionNumber', None)
                fid.pop('tablesVersion', None)
                cfVarName = namesgribdef.cfVarName(fid,
                                                   'grib{}'.format(initial_fid['editionNumber']))
                if len(cfVarName) == 1:
                    target_fid = list(cfVarName.keys())[0]
            except Exception:
                pass"""
    else:
        raise NotImplementedError("this kind of conversion.")
    return target_fid


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
    resource = eopen(filename, openmode='r')
    if one_output_file:
        output_resource = eopen('.'.join([filename, output_format_suffix]),
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
            field.fid['GRIB2'] = fid_converter(field.fid['generic'], 'generic', 'GRIB2')
            assert all([k in field.fid['GRIB2']
                        for k in ('discipline', 'parameterCategory', 'parameterNumber')]), \
                   "missing key(s) among ('discipline', 'parameterCategory', 'parameterNumber')"
        else:
            try:
                field.fid[fmt_dict[output_format_suffix]] = fid_converter(fid,
                                                                          resource.format,
                                                                          fmt_dict[output_format_suffix],
                                                                          grib_short_fid=grib_short_fid)
            except NotImplementedError:
                field.fid[fmt_dict[output_format_suffix]] = str(field.fid[resource.format])
        # now write field !
        if not one_output_file:
            if output_format_suffix == 'geo':
                output_resource = eopen('.'.join([filename,
                                                  field.fid[fmt_dict[output_format_suffix]],
                                                  output_format_suffix]),
                                        openmode='w',
                                        fmt=fmt_dict[output_format_suffix],
                                        other_attributes=({'FORMAT':'XYV'} if kwargs.get('llv') else None),
                                        columns=kwargs.get('columns'),
                                        no_header=kwargs.get('no_header', False)
                                        )
            else:
                output_resource = eopen('.'.join([filename,
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


def batch_convert(filenames,
                  output_format_suffix,
                  threads_number=1,
                  progressmode=None,
                  **kwargs):
    """
    Converts a series of files to *output_format_suffix* in batch.

    Mandatory arguments:

    :param filenames: name(s) of the files to be processed
    :param output_format_suffix: among 'grb' (GRIB2), 'nc' (netCDF4), 'geo' (GeoPoints)

    Technical named (optional) arguments:

    :param threads_number: parallelisation of files processing
    :param progressmode: among ('verbose', 'percentage', None)

    Other named arguments depend on the output format, and are defined in the
    Workers footprints attributes !
    """
    import taylorism
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
    if isinstance(filenames, str):
        filenames = [filenames]
    individual_instructions = {'filename':filenames}

    # run !
    taylorism.batch_main(common_instructions, individual_instructions,
                         scheduler=footprints.proxy.scheduler(limit='threads', max_threads=threads_number),
                         verbose=(progressmode == 'verbose'))

