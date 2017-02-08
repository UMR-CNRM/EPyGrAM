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
  - create a Resource instance with a generic function,
    eventually (if already existing) without knowing its format *a priori*;
  - some catalogs of GRIB numbers
  - ...
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import importlib
import copy
import os

from footprints import proxy as fpx
from epygram import config, epygramError, util

__all__ = []

_formats_in_loading_order = copy.copy(config.implemented_formats)
_loaded_first_formats = ['GRIB', 'FA', 'LFI', 'DDHLFA', 'LFA']  # because they might have been compiled with newer libraries versions
for lff in _loaded_first_formats[::-1]:
    if lff in _formats_in_loading_order:
        _formats_in_loading_order = [lff] + [f for f in _formats_in_loading_order if f != lff]

for f in _formats_in_loading_order:
    if f not in [m['name'] for m in config.usermodules]:
        importlib.import_module('.' + f, __name__)


#############
# UTILITIES #
#############
def guess(filename):
    """
    Returns the name of the format of the resource located at a given
    **filename**, if succeeded.
    """

    formats_in_guess_order = copy.copy(config.implemented_formats)
    _guess_last_formats = ['FA', 'LFI', ]  # because they might not be very clean at catching exceptions
    for glf in _guess_last_formats:
        if glf in formats_in_guess_order:
            formats_in_guess_order = [f for f in formats_in_guess_order
                                      if f != glf] + [glf]
    for f in formats_in_guess_order:
        try:
            r = fpx.dataformat(filename=filename, openmode='r', format=f)
            r.close()
            fmt = f
            break
        except IOError:
            fmt = 'unknown'

    return fmt


def resource(filename, openmode, fmt=None, **kwargs):
    """
    Returns an instance of Resource of the requested *fmt* format,
    located at the given *filename* (:class:`epygram.base.File`, ...),
    and with the given *openmode*.

    If *fmt* is not given, tries to guess it (only for *openmode* 'r' or 'a').

    Other *kwargs* are passed to the resource constructor.
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


def fid_converter(initial_fid, initial_fmt, target_fmt,
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
    else:
        raise NotImplementedError("this kind of conversion.")

    return target_fid
