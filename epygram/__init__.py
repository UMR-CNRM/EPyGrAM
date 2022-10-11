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
that can be used 'as is', or can be taken as templates for building more complex
applications with the ``epygram`` package.

The package uses extensively the ``footprints`` package (included in epygram
distributions), designed by the VORTEX team (MF) in order (basically) to tackle
once-for-all with Factories for families of similar classes. The basic idea
behind this concept is that similar classes have a similar set of attributes,
but with behavioral differences (values or name/number of attributes), viewed
as their "footprint".

For the needs of the FA, LFI & LFA formats and spectral transforms of fields
from ARPEGE/ALADIN/AROME models, the :mod:`arpifs4py` library is needed. It is also
therein included, already compiled for Mageia4 platforms, and with necessary
stuff for recompiling it on other platforms. Be aware to be recompiled, this
library needs an *arpifs* pack pre-compiled with *gmkpack*, and *gribex*
library.

Complete dependencies to be found in :ref:`Dependancies <dependancies>`.

It is recalled that packages available on Pypi.org can be installed locally
using: ``pip[2|3] install --user <packagename>``

********************************************************************************

.. _license:

License
-------

Copyright Météo France (2014 - 2022), authors:

* A. Mary - Météo France, CNRM/GMAP/COOPE - alexandre.mary@meteo.fr
* S. Riette - Météo France

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.

********************************************************************************
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import sys
import io
import os

import footprints

package_rootdir = os.path.dirname(os.path.realpath(__path__[0]))  # realpath to resolve symlinks

__all__ = []

__version__ = io.open(os.path.join(package_rootdir, 'VERSION'), 'r').read().strip()

__license__ = 'CeCILL-C'

__authors__ = ['Alexandre Mary', 'Sébastien Riette']

__contributors__ = ['Stéphanie Faroux', 'Ghislain Faure']


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
    if sys.version_info.minor < 7:
        epylog.warning('*epygram* requires Python2.7 at least. ' +
                       'It may not work properly with older versions.')

# config
from . import config

# COMPONENTS (modules) #
from . import util
from . import base
from . import containers
from . import moves

from . import args_catalog
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

# User modules  # TODO: ? CLEANME
if len(config.usermodules) > 0:
    footprints.priorities.set_before('debug', 'user')
    for m in config.usermodules:
        if sys.version_info.major == 3 and sys.version_info.minor >= 4:
            import importlib.util as imputil  # @UnresolvedImport
            spec = imputil.spec_from_file_location(m['name'],
                                                   m['abspath'])
            spec.loader.exec_module(imputil.module_from_spec(spec))
            del spec
        else:
            import imp
            imp.load_source(m['name'], m['abspath'])

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
             ensure_consistent_GRIB_paths=True,
             ignore_gribenv_paths=False):
    """
    A function to modify execution environment (to be called early in
    execution).

    :param omp_num_threads: sets OMP_NUM_THREADS
    :param no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param lfi_C: if True, LFI_HNDL_SPEC set to ':1', to use the C version of LFI
    :param mute_FA4py: mute messages from FAIPAR in FA4py library
    :param bool ensure_consistent_GRIB_paths: complete GRIB samples/definition
        paths to be consistent with inner library
    :param ignore_gribenv_paths: ignore predefined values of the variables
        GRIB_SAMPLES_PATH and GRIB_DEFINITION_PATH
        (or equivalent ECCODES variables)
    """
    # 1. arpifs4py library
    # FA & LFI need some special environment setting
    if any([f in formats.runtime_available_formats for f in ('FA', 'LFI')]):
        import arpifs4py
        if mute_FA4py is None:
            mute_FA4py = config.FA_mute_FA4py
        arpifs4py.init_env(omp_num_threads=omp_num_threads, no_mpi=no_mpi,  # common
                           lfi_C=lfi_C, mute_FA4py=mute_FA4py,  # LFI/FA
                           ensure_consistent_GRIB_paths=ensure_consistent_GRIB_paths,  # GRIB paths
                           ignore_gribenv_paths=ignore_gribenv_paths)
    # 2. SpectralGeometry inner transformation lib may need some special
    # environment setting, delayed to actual invocation:
    # we simply pass kwargs, initialization is done at first call to the library
    from .geometries.SpectralGeometry import transforms_lib_init_env_kwargs
    transforms_lib_init_env_kwargs.update(omp_num_threads=omp_num_threads,
                                          no_mpi=no_mpi,
                                          trigger=True)
    # 3. grib_api or eccodes
    # need some special environment setting
    # ensure grib_api/eccodes variables are consistent with inner library
    if 'GRIB' in formats.runtime_available_formats and ensure_consistent_GRIB_paths:
        from .formats.GRIB import lowlevelgrib
        lowlevelgrib.init_env(reset=ignore_gribenv_paths)


if config.init_at_import:
    init_env()
