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

from ctypes import c_longlong, c_char_p, PyDLL, c_bool, c_double
import os
import numpy as np
import subprocess

from . import ctypesForFortran

_libs4py = "libs4py.so"  # local name in the directory


# helpers
#########
def addReturnCode(func):
    def wrapper(*args):
        """
        This decorator adds an integer at the beginning of the "returned"
        signature of the Python function.
        """
        return [(c_longlong(), OUT)] + func(*args)
    wrapper.__name__ = func.__name__
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
    return wrapper


def complete_GRIB_samples_path_from_dynamic_gribapi(lib):
    import re
    ldd_out = subprocess.check_output(['ldd', lib])
    libs_grib_api = {}
    for line in ldd_out.splitlines():
        match = re.match('\s*(libgrib_api.*) => (.*)(/lib/libgrib_api.*\.so.*)\s\(0x.*', line)
        if match:
            libs_grib_api[match.group(1)] = match.group(2)
    gsp = os.getenv('GRIB_SAMPLES_PATH', '.')
    for l in set(libs_grib_api.values()):
        gsp = ':'.join([gsp, '/'.join([l, 'share/grib_api/samples'])])
    os.environ['GRIB_SAMPLES_PATH'] = gsp  # FIXME: seems not to work on Bull: to be exported beforehand ?


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
os.unsetenv('GRIB_DEFINITION_PATH')

so = PyDLL(shared_objects_library)
ctypesFF = ctypesForFortran.ctypesForFortranFactory(so)


# sub-modules
#############
from . import wfa
from . import wlfi
from . import wlfa
from . import wtransforms
