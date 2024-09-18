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

import ctypesForFortran

_libs4py = "libs4py.so"  # local name in the directory


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


def FALFI_get_dynamic_gribapi_lib_paths():
    """If needed, set adequate path to the used low level library."""
    libs_grib_api = {}
    for apilib in ('grib_api', 'eccodes'):
        for l, libpath in ctypesForFortran.get_dynamic_libs(shared_objects_library).items():
            if l.startswith('lib' + apilib):
                libs_grib_api[apilib] = libpath
    return libs_grib_api


# initialization
################
def FALFI_init_env(omp_num_threads=None,
                   no_mpi=False,
                   lfi_C=False,
                   mute_FA4py=False,
                   unlimited_stack=False):
    """
    Set adequate environment for the inner libraries.

    :param int omp_num_threads: sets OMP_NUM_THREADS
    :param bool no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param bool lfi_C: if True, LFI_HNDL_SPEC set to ':1', to use the C version of LFI
    :param bool mute_FA4py: mute messages from FAIPAR in FA4py library
    :param unlimited_stack: equivalent of 'ulimit -s unlimited'
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
    # unlimited stack
    if unlimited_stack:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY,resource.RLIM_INFINITY))
