#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Module contains:

- all formats classes.
- utilities to play with resource formats:\n
  - guess the format of an existing resource in a given container;
  - open a Resource instance with a generic function (potentially guessing its format if in read mode);
  - FA fields utilities
  - conversion utilities
"""

import importlib
import copy
import os

import footprints
from footprints import proxy as fpx
from bronx.system.unistd import stderr_redirected

from .. import config, epygramError
from . import fafields

__all__ = []

epylog = footprints.loggers.getLogger(__name__)


def available_format(fmt):
    """Check availability of low-level libraries for activating format."""
    reason = None
    available = True
    try:
        if fmt in ('FA', 'LFI', 'DDHLFA', 'LFA'):
            import falfilfa4py
        elif fmt == 'GRIB':
            if config.GRIB_lowlevel_api.lower() in ('gribapi', 'grib_api'):
                import gribapi
            elif config.GRIB_lowlevel_api.lower() == 'eccodes':
                import eccodes
            else:
                available = False
                reason = 'Unknown config.GRIB_lowlevel_api={}'.format(config.GRIB_lowlevel_api)
        elif fmt in ('netCDF', 'netCDFMNH'):
            import netCDF4
        elif fmt == 'TIFFMF':
            from epygram.extra import pyexttiff
        elif fmt == 'HDF5SAF':
            import h5py
    except Exception as e:
        available = False
        reason = str(e)
    return available, reason


# Formats loading may have to follow an order,
# for common dynamic libraries of different versions. Order to be set in config.implemented_formats

#: list of formats actually available at runtime
runtime_available_formats = []
# import formats modules
for f in config.implemented_formats:
    _available, _reason = available_format(f)
    if _available:
        runtime_available_formats.append(f)
        importlib.import_module('.' + f, __name__)
    else:
        epylog.info(("Format: {} is deactivated at runtime (Error: {}). " +
                     "Please deactivate from config.implemented_formats or fix error.").format(f, _reason))


#############
# UTILITIES #
#############
def guess(filename):
    """
    Returns the name of the format of the resource located at a given
    **filename**, if succeeded.
    """

    formats_in_guess_order = copy.copy(runtime_available_formats)
    _guess_last_formats = ['DDHLFA', 'LFA', 'FA', 'LFI', ]  # because they may not be very clean at catching exceptions
    for glf in _guess_last_formats:
        if glf in formats_in_guess_order:
            formats_in_guess_order = [f for f in formats_in_guess_order
                                      if f != glf] + [glf]
    for f in formats_in_guess_order:
        try:
            if config.silent_guess_format:
                with stderr_redirected():
                    r = fpx.dataformat(filename=filename, openmode='r', format=f)
                    r.close()
            else:
                r = fpx.dataformat(filename=filename, openmode='r', format=f)
                r.close()
            fmt = f
            break
        except IOError:
            fmt = 'unknown'

    return fmt


def resource(filename, openmode, fmt=None, **kwargs):
    """
    Returns an instance of Resource of the requested **fmt** format,
    located at the given **filename**, open with the given **openmode**.

    :param filename: name (path) of the file to open
    :param openmode: opening mode ('r', 'a', 'w')
    :param fmt: format of the resource; with openmode 'r' or 'a', *fmt* is
                optional and can be guessed from the existing resource

    Other kwargs are passed to the resource constructor.
    """
    if openmode in ('r', 'a'):
        assert os.path.exists(filename), 'File not found: ' + filename
    if fmt is None and openmode in ('r', 'a'):
        fmt = guess(filename)
        if fmt == 'unknown':
            raise epygramError("unable to guess format of resource at: " +
                               filename)
    elif fmt is None and openmode == 'w':
        raise epygramError("must specify 'fmt' argument with\
                            'openmode' == 'w'.")

    return fpx.dataformat(filename=filename, openmode=openmode, format=fmt,
                          **kwargs)

from . import conversion

