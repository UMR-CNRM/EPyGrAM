#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities around GRIB paths.
"""

import os
import subprocess


def get_eccodes_from_ldconfig():
    """DEPRECATED.Get eccodes install directory from ldconfig."""
    print("DEPRECATED:get_eccodes_from_ldconfig shouln't be called anymore.")
    out = str(subprocess.check_output(['/sbin/ldconfig', '-p']))
    out_split = out.split(r'\n')
    libs_eccodes = [lib for lib in out_split if str('libeccodes.so') in lib]
    paths = [lib.split('=>')[1].strip() for lib in libs_eccodes]
    dirs = set([os.path.sep.join(lib.split(os.path.sep)[:-2]) for lib in paths])
    assert len(dirs) == 1, "More than one 'libeccodes.so' has been found"
    return dirs.pop()


def complete_grib_paths(rootdir, reset=False):
    """
    Complete ECCODES_SAMPLES_PATH and ECCODES_DEFINITION_PATH
    according to **rootdir** installation path of ECCODES.

    :param rootdir: the directory in which is installed the API
    :param reset: ignore predefined values of the variables

    Reconstructed path are ``$rootdir$/share/eccodes/samples``
    and ``$rootdir$/share/eccodes/definitions``
    """
    complete_grib_samples_paths(rootdir, reset=reset)
    complete_grib_definition_paths(rootdir, reset=reset)


def complete_grib_samples_paths(rootdir, reset=False):
    """
    Complete ECCODES_SAMPLES_PATH according to **rootdir**
    installation path of ECCODES.

    :param rootdir: the directory in which is installed the API
    :param reset: ignore predefined values of the variables

    Reconstructed path is ``$rootdir$/share/eccodes/samples``
    """
    sp = 'ECCODES_SAMPLES_PATH'
    loc_samples = [os.path.join(rootdir, 'share', 'eccodes', 'samples')]
    if not reset and os.environ.get(sp, False):
        loc_samples.append(os.environ.get(sp))
    os.environ[sp] = os.pathsep.join(loc_samples)


def complete_grib_definition_paths(rootdir, reset=False):
    """
    Complete ECCODES_DEFINITION_PATH according to **rootdir**
    installation path of ECCODES.

    :param rootdir: the directory in which is installed the API
    :param reset: ignore predefined values of the variables

    Reconstructed path are ``$rootdir$/share/eccodes/definitions``
    """
    dp = 'ECCODES_DEFINITION_PATH'
    loc_defs = [os.path.join(rootdir, 'share', 'eccodes', 'definitions')]
    if not reset and os.environ.get(dp, False):
        loc_defs.append(os.environ.get(dp))
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

