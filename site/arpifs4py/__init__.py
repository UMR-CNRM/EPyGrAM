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
import resource
import numpy as np

import griberies
import ctypesForFortran

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


# common parameters and objects
###############################
IN = ctypesForFortran.IN
OUT = ctypesForFortran.OUT
INOUT = ctypesForFortran.INOUT

# shared objects library
########################
shared_objects_library = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      _libs4py)
ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)


# sub-modules
#############
from . import wfa
from . import wlfi
from . import wlfa
from . import wtransforms


# initialization
################
def init_env(omp_num_threads=None,
             no_mpi=False,
             lfi_C=False,
             mute_FA4py=False,
             ensure_consistent_GRIB_paths=False,
             ignore_gribenv_paths=False):
    """
    Set adequate environment for the inner libraries.

    :param int omp_num_threads: sets OMP_NUM_THREADS
    :param bool no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param bool lfi_C: if True, LFI_HNDL_SPEC set to ':1', to use the C version of LFI
    :param bool mute_FA4py: mute messages from FAIPAR in FA4py library
    :param bool ensure_consistent_GRIB_paths: complete GRIB samples/definition
        paths to be consistent with inner library
    :param bool ignore_gribenv_paths: ignore predefined values of the variables
        GRIB_SAMPLES_PATH and GRIB_DEFINITION_PATH
        (or equivalent ECCODES variables)
    """
    # because arpifs library is compiled with MPI & openMP
    if omp_num_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)
    if no_mpi:
        os.environ['DR_HOOK_NOT_MPI'] = '1'
    # use the C library for LFI
    if lfi_C:
        os.environ['LFI_HNDL_SPEC'] = ':1'
    # option for FA
    if mute_FA4py:
        os.environ['FA4PY_MUTE'] = '1'
    # ensure grib_api/eccodes variables are consistent with inner library
    if ensure_consistent_GRIB_paths:
        if len(griberies.get_samples_paths() + griberies.get_definition_paths()) > 0:
            _GRIB_samples_path_from_dynamic_gribapi(shared_objects_library,
                                                    reset=ignore_gribenv_paths)


def _GRIB_samples_path_from_dynamic_gribapi(lib, reset=False):
    """If needed, set adequate path to the used low level library."""
    libs_grib_api = {}
    for apilib in ('grib_api', 'eccodes'):
        for l, libpath in ctypesForFortran.get_dynamic_libs(lib).items():
            if l.startswith('lib' + apilib):
                griberies.complete_grib_paths(libpath, apilib, reset=reset)
