#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
Enhanced Python for Graphics and Analysis of Meteorological fields
------------------------------------------------------------------

``epygram`` is a package of classes, designed to handle Meteorological Fields,
and various resource Formats from which the Fields can be extracted.

It is distributed along with a series of applicative tools using the package,
that can be used as command line tools, or can be taken as templates for building more complex
applications with the ``epygram`` package.

********************************************************************************

.. _license:

License
-------

Copyright Météo France (2014)

Initial authors:
* A. Mary - Météo France, CNRM/GMAP/COOPE - alexandre.mary@meteo.fr
* S. Riette - Météo France

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

The license text is provided in the LICENSE.txt file of the package.

********************************************************************************
"""

import sys
import io
import os

import footprints

package_rootdir = os.path.dirname(os.path.realpath(__path__[0]))  # realpath to resolve symlinks

__all__ = []

__version__ = "2.0.6"

__license__ = 'CeCILL-C'

__authors__ = ['Alexandre Mary', 'Sébastien Riette']


class epygramError(Exception):
    """Errors class for the package."""
    pass


# : Root log for epygram
epylog = footprints.loggers.getLogger(__name__)


# Check that Python version is compatible
if sys.version_info.major == 3:
    if sys.version_info.minor < 5:
        epylog.warning('*epygram* requires Python3.5 at least. ' +
                       'It may not work properly with older versions.')
else:
    epylog.warning('*epygram* support for Python2 is not maintained !')

# config
from . import config

# COMPONENTS (modules) #
from . import util
from . import base
from . import containers
from . import moves

from . import cli
from . import profiles
from . import spectra

from . import geometries
from . import fields
from . import formats
from . import resources

# Register plugins
from . import _plugins
for p in config.activate_plugins:
    if p in _plugins.available():
        _plugins.activate(p)
    else:
        epylog.info("Plugin '{}' from config.activate_plugins is not available.".format(p))

# Further initializations
if config.hide_footprints_warnings:
    footprints.logger.setLevel(footprints.loggers.logging.ERROR)


# OTHERS #
def showconfig():
    """
    Print current config.
    """
    import copy
    cfg = copy.copy(config.__dict__)
    header = "### " + cfg['__name__'] + " ###"
    print("#" * len(header))
    print(header)
    print("#" * len(header))
    to_be_removed = ['os', 'sys', 'imp', 'platform', 'footprints',
                     '__name__', '__doc__', '__package__', '__file__',
                     '__builtins__']
    for p in to_be_removed:
        cfg.pop(p, None)
    for k in sorted(cfg.keys()):
        print('- ' + k + ' = ' + str(cfg[k]))


def init_env(omp_num_threads=1,
             no_mpi=True,
             lfi_C=True,
             mute_FA4py=None,
             unlimited_stack=True,
             ensure_consistent_GRIB_paths=False,
             ignore_gribenv_paths=False,
             fa_limits=config.FA_limits,
             ):
    """
    A function to modify execution environment (to be called early in
    execution).

    :param omp_num_threads: sets OMP_NUM_THREADS
    :param no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param lfi_C: if True, LFI_HNDL_SPEC set to ':1', to use the C version of LFI
    :param mute_FA4py: mute messages from FAIPAR in FA4py library
    :param unlimited_stack: equivalent to 'ulimit -s unlimited'
    :param bool ensure_consistent_GRIB_paths: complete GRIB samples/definition
        paths to be consistent with inner library
    :param ignore_gribenv_paths: ignore predefined values of the variables
        GRIB_SAMPLES_PATH and GRIB_DEFINITION_PATH
        (or equivalent ECCODES variables)
    :param fa_limits: FA limits
    """
    # 1. falfilfa4py library
    # FA & LFI need some special environment setting
    if any([f in formats.runtime_available_formats for f in ('FA', 'LFI')]):
        import falfilfa4py
        if mute_FA4py is None:
            mute_FA4py = config.FA_mute_FA4py
        falfilfa4py.init_env(omp_num_threads=omp_num_threads, no_mpi=no_mpi,
                             lfi_C=lfi_C, mute_FA4py=mute_FA4py,
                             unlimited_stack=unlimited_stack,
                             fa_limits=fa_limits)
        if ensure_consistent_GRIB_paths: #CLEANME:DEPRECATED:
            from epygram.extra import griberies
            eccodes_libpath = falfilfa4py.get_dynamic_eccodes_lib_path_from_FA()
            griberies.paths.complete_grib_paths(eccodes_libpath, reset=ignore_gribenv_paths)
    # 2. SpectralGeometry inner transformation lib may need some special
    # environment setting, delayed to actual invocation:
    # we simply pass kwargs, initialization is done at first call to the library
    from .geometries.SpectralGeometry import ectrans4py_init_env_kwargs
    ectrans4py_init_env_kwargs.update(omp_num_threads=omp_num_threads,
                                      no_mpi=no_mpi,
                                      unlimited_stack=unlimited_stack,
                                      trigger=True)
    # 3. grib_api or eccodes
    # need some special environment setting
    # ensure eccodes variables are consistent with inner library
    if 'GRIB' in formats.runtime_available_formats and ensure_consistent_GRIB_paths:
        from .formats.GRIB import lowlevelgrib
        lowlevelgrib.init_env(reset=ignore_gribenv_paths)


#: shortcut for epygram.formats.resource()
open = formats.resource

if config.init_at_import:
    init_env()
