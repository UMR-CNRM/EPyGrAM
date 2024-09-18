#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
falfilfa4py:

Contains the interface routines to IO libraries for formats:
- LFI
- FA (overlay of LFI)
- LFA (DDH format)

Note that this is temporary between the former package arpifs4py and a direct python interface to fully externalised
FA_LFI/LFA library.

Actual .so library should be in one of the preinstalled paths or in a directory specified via LD_LIBRARY_PATH
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import os
import resource
import ctypesForFortran


# Shared objects library
########################
shared_objects_library = os.environ.get('FALFILFA4PY_SO', None)
if shared_objects_library is None or not os.path.exists(shared_objects_library):
    # not specified or path does not exist : find in known locations
    so_basename = "falfilfa4py.so.0"  # local name in the directory
    LD_LIBRARY_PATH = [p for p in os.environ.get('LD_LIBRARY_PATH', '').split(':') if p != '']
    potential_locations = LD_LIBRARY_PATH + [
        "/home/common/epygram/public/EPyGrAM/libs4py",  # CNRM
        "/home/gmap/mrpe/mary/public/EPyGrAM/libs4py",  # belenos/taranis
        "/home/acrd/public/EPyGrAM/libs4py",  # ECMWF's Atos aa-ad
        ]
    for _libs4py_dir in potential_locations:
        shared_objects_library = os.path.join(_libs4py_dir, so_basename)
        if os.path.exists(shared_objects_library):
            break
        else:
            shared_objects_library = None
    if shared_objects_library is None:
        msg = ' '.join(["'{}' was not found in any of potential locations: {}.",
                        "You can specify a different location using env var LD_LIBRARY_PATH",
                        "or specify a precise full path using env var FALFILFA4PY_SO."]).format(
                                so_filename, str(potential_locations))
        raise FileNotFoundError(msg)
else:
    so_basename = os.path.basename(shared_objects_library)
ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)

def get_dynamic_eccodes_lib_paths_from_FA():
    """Get paths to the eccodes/grib_api linked for FA purpose in the shared objects library."""
    libs_grib_api = {}
    for apilib in ('grib_api', 'eccodes'):
        for l, libpath in ctypesForFortran.get_dynamic_libs(shared_objects_library).items():
            if l.startswith('lib' + apilib):
                libs_grib_api[apilib] = libpath
    return libs_grib_api

# Initialization
################

def init_env(omp_num_threads=None,
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
    :param unlimited_stack: equivalent to 'ulimit -s unlimited'
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
    # ulimit -s unlimited
    if unlimited_stack:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY,resource.RLIM_INFINITY))


# sub-modules
#############
from . import FA
from . import LFI
from . import LFA
