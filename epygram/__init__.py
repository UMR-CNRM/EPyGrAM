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

Other dependencies,

Mandatory:

- :mod:`numpy`
- :mod:`pyproj`

Functional (needed for specific, more or less usual functionalities):

- :mod:`matplotlib`
- :mod:`mpl_toolkits.basemap`
- :mod:`scipy`
- :mod:`pyresample`
- :mod:`nbsphinx`

it is recalled that packages available on Pypi can be installed locally using:
pip install --user <packagename>

and (if according formats are activated):

- :mod:`dateutil` (only for :class:`epygram.formats.netCDF`)
- :mod:`gribapi` (only for :mod:`epygram.formats.GRIB`)
- :mod:`netCDF4` (only for :class:`epygram.formats.netCDF`)
- :mod:`PIL` (only for :class:`epygram.formats.TIFFMF`)



.. _license:

********************************************************************************

Copyright Météo France (2014 - 2016),
contributors :
A. Mary <Météo France, CNRM/GMAP/COOPE, alexandre.mary@meteo.fr> ||
S. Riette

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

import footprints

__all__ = []

__version__ = '1.2.11'

__license__ = 'CeCILL-C'

__authors__ = ['Alexandre Mary', 'Sébastien Riette']

__contributors__ = ['Stéphanie Faroux', 'Ghislain Faure']


class epygramError(Exception):
    """Errors class for the package."""
    pass

# : Root log for epygram
epylog = footprints.loggers.getLogger(__name__)

# config
from . import config
if config.noninteractive_backend:
    try:
        import matplotlib
        matplotlib.use(config.noninteractive_backend)
    except ImportError:
        pass

# COMPONENTS (modules) #
from . import util
from . import base
from . import containers

from . import args_catalog
from . import profiles
from . import myproj
from . import spectra

from . import geometries
from . import fields
from . import formats
from . import resources

# User modules
if len(config.usermodules) > 0:
    import sys
    footprints.priorities.set_before('debug', 'user')
    for m in config.usermodules:
        if sys.version_info.major == 3 and sys.version_info.minor >= 4:
            import importlib.util as imputil
            spec = imputil.spec_from_file_location(m['name'],
                                                   m['abspath'])
            spec.loader.exec_module(imputil.module_from_spec(spec))
            del spec
        else:
            import imp
            imp.load_source(m['name'], m['abspath'])


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

    toberemoved = ['os', 'sys', 'imp', 'platform', 'footprints', 'copy',
                   '__name__', '__doc__', '__package__', '__file__',
                   '__builtins__']

    for p in toberemoved:
        cfg.pop(p, None)
    for k in sorted(cfg.keys()):
        print('- ' + k + ' = ' + str(cfg[k]))


def init_env(omp_num_threads=1, no_mpi=True, unlimited_stack=True, lfi_C=True):
    """
    A function to modify execution environment (to be called early in
    execution).

    :param no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param omp_num_threads: sets OMP_NUM_THREADS
    :param lfi_C: if True, LFI_HNDL_SPEC set to ':1', to use the C version of LFI
    :param unlimited_stack: stack size unlimited on Bull supercomputers.
    """
    import os
    import resource

    # because arpifs library is compiled with MPI & openMP
    os.putenv('OMP_NUM_THREADS', str(omp_num_threads))
    if no_mpi:
        os.putenv('DR_HOOK_NOT_MPI', '1')
    if lfi_C:
        os.putenv('LFI_HNDL_SPEC', ':1')
    if unlimited_stack and ('beaufix' in os.getenv('HOSTNAME', '') or
                            'prolix' in os.getenv('HOSTNAME', '')):
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY,
                                                   resource.RLIM_INFINITY))

if config.init_at_import:
    init_env()
