#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains equivalences tables for GRIB encoding.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import os

#: Aliases to *productionStatusOfProcessedData* numbers
productionStatusOfProcessedData_dict = {'oper':0,
                                        'dble':1,
                                        'test':2,
                                        'research':2,
                                        'unknown':255}

#: Aliases to *typeOfGeneratingProcess* numbers
typeOfGeneratingProcess_dict = {'Analysis':0,
                                'Initialization':1,
                                'Forecast':2,
                                'Bias corrected forecast':3,
                                'Ensemble forecast':4,
                                'Probability forecast':5,
                                'Forecast error':6,
                                'Analysis error':7,
                                'Observation':8,
                                'Climatological':9,
                                'Probability-weighted forecast':10,
                                'Bias-corrected ensemble forecast':11,
                                'Post-processed analysis':12,
                                'Post-processed forecast':13,
                                'Nowcast':14,
                                'Hindcast':15,
                                'unknown':255}

#: Geoid shapes for pyproj
pyproj_geoid_shapes = {0:{'a':6367470.,
                          'b':6367470.},
                       2:{'a':6378160.,
                          # 'b':6356775.,
                          'rf':297.},
                       4:{'ellps':'GRS80',
                          # 'a':6378137.,
                          # 'b':6356752.314,
                          # 'rf':298.257222101
                          },
                       5:{'ellps':'WGS84',
                          # 'a':6378137.,
                          # 'b':6356752.3142,
                          # 'flattening':298.257223563,
                          },
                       6:{'a':6371229.,
                          'b':6371229.},
                       8:{'a':6371200.,
                          'b':6371200.},
                       9:{'ellps':'airy',  # TOBECHECKED:
                          # 'a':6377563.369,
                          # 'b':6356752.314,
                          # 'rf':299.3249646,
                          # 'lambda0':0,
                          }}

#: Type of statistical process over a duration
statistical_processes = {0:'average',
                         1:'accumulation',
                         2:'maximum',
                         3:'minimum',
                         4:'difference',
                         5:'rms',
                         6:'stdev',
                         7:'covariance',
                         8:'-difference',
                         9:'ratio',
                         10:'summation',
                         11:'standardized anomaly'
                         }


def complete_grib_paths(rootdir, api_name, reset=False):
    """
    Complete GRIB_SAMPLES_PATH and GRIB_DEFINITION_PATH according to **rootdir**
    installation path of GRIB API **api_name**.

    :param rootdir: the directory in which is installed the API
    :param api_name: the name of the GRIB API, among ('eccodes', 'grib_api')

    Reconstructed path are ``$rootdir$/share/$api_name$/samples``
    and ``$rootdir$/share/$api_name$/definitions``
    """
    # FIXME: seems not to work on Bull: to be exported beforehand ?
    if api_name == 'grib_api':
        sp = 'GRIB_SAMPLES_PATH'
        dp = 'GRIB_DEFINITION_PATH'
    elif api_name == 'eccodes':
        sp = 'ECCODES_SAMPLES_PATH'
        dp = 'ECCODES_DEFINITION_PATH'
    loc_samples = [os.path.join(rootdir, 'share', api_name, 'samples')]
    loc_defs = [os.path.join(rootdir, 'share', api_name, 'definitions')]
    if not reset and os.environ.get(sp, False):
        loc_samples.append(os.environ.get(sp))
    if not reset and os.environ.get(dp, False):
        loc_defs.append(os.environ.get(dp))
    os.environ[sp] = os.pathsep.join(loc_samples)
    os.environ[dp] = os.pathsep.join(loc_defs)


def _get_paths(obj):
    paths = os.pathsep.join([os.environ.get('ECCODES_{}_PATH'.format(obj), ''),
                             os.environ.get('GRIB_{}_PATH'.format(obj), '')])
    return [p for p in paths.split(os.pathsep) if p not in ('', '.')]


def get_samples_paths():
    """Get the environment-variable-set path to samples"""
    return _get_paths('SAMPLES')


def get_definition_paths():
    """Get the environment-variable-set path to definitions"""
    return _get_paths('DEFINITION')
    