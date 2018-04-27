#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
arpifs4py:

Contains the interface routines to arpifs code:
- LFI format
- FA format
- LFA format
- spectral transforms
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import os
import numpy as np

from . import ctypesForFortran

_libs4py = "libs4py.so"  # local name in the directory


# helpers
#########
def addReturnCode(func):
    def wrapper(*args, **kwargs):
        """
        This decorator adds an integer at the beginning of the "returned"
        signature of the Python function.
        """
        out = func(*args, **kwargs)
        out[1].insert(0, (np.int64, None, OUT))
        return out
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def treatReturnCode(func):
    def wrapper(*args):
        """
        This decorator raises a Python error if the integer returned by
        addReturnCode is different from 0.
        """
        result = func(*args)
        try:
            nout = len(result)
        except Exception:
            nout = 1
        if nout == 1:
            result = (result,)
        if result[0] != 0:
            raise RuntimeError("arpifs4py: Error code " + str(result[0]) + " was raised.")
        result = result[1:]
        if len(result) == 1:
            result = result[0]
        elif len(result) == 0:
            result = None
        return result
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def treatReturnCode_LFA(func):
    def wrapper(*args):
        """
        This decorator raises a Python error if the integer returned by
        addReturnCode is different from 0.
        The message raised is linked to the LFA library.
        """
        LFA_errors = {'-999':'Error while searching for an available logical\
                              unit.',
                      '-1':'Field not found in file.',
                      '1':'Field not found in file.',
                      '-6':'Article bigger than expected (argument kdimb).',
                      '-8':'Wrong data type (real, integer, char.).'}
        result = func(*args)
        try:
            nout = len(result)
        except Exception:
            nout = 1
        if nout == 1:
            result = (result,)
        if result[0] != 0:
            raise RuntimeError("Error code " + str(result[0]) +
                               " was raised : " +
                               LFA_errors[str(result[0])])
        result = result[1:]
        if len(result) == 1:
            result = result[0]
        elif len(result) == 0:
            result = None
        return result
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def complete_GRIB_samples_path_from_dynamic_gribapi(lib):
    import subprocess
    import re
    import grib_utilities
    ldd_out = subprocess.check_output(['ldd', lib])
    libs_grib_api = {}
    for apilib in ('libgrib_api', 'libeccodes'):
        for line in ldd_out.splitlines():
            _re = '.*\s*({}.*) => (.*)(/lib/{}.*\.so.*)\s\(0x.*'.format(apilib, apilib)
            match = re.match(_re, str(line))
            if match:
                libs_grib_api[match.group(1)] = match.group(2)
    for l in set(libs_grib_api.values()):
        grib_utilities.complete_grib_paths(l, 'grib_api', reset=False)


# common parameters and objects
###############################
IN = ctypesForFortran.IN
OUT = ctypesForFortran.OUT
INOUT = ctypesForFortran.INOUT

# shared objects library
########################
shared_objects_library = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      _libs4py)
complete_GRIB_samples_path_from_dynamic_gribapi(shared_objects_library)
ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)


# sub-modules
#############
from . import wfa
from . import wlfi
from . import wlfa
from . import wtransforms
