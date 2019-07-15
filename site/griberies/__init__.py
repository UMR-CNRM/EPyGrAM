#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities around GRIB format.
"""
from __future__ import print_function, absolute_import, unicode_literals, division
import six

import os
import re
import io
import copy

from bronx.syntax.parsing import str2dict

from . import tables


def complete_grib_paths(rootdir, api_name, reset=False):
    """
    Complete GRIB_SAMPLES_PATH and GRIB_DEFINITION_PATH according to **rootdir**
    installation path of GRIB API **api_name**.

    :param rootdir: the directory in which is installed the API
    :param api_name: the name of the GRIB API, among ('eccodes', 'grib_api')
    :param reset: ignore predefined values of the variables

    Reconstructed path are ``$rootdir$/share/$api_name$/samples``
    and ``$rootdir$/share/$api_name$/definitions``
    """
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


def set_definition_path(path, api_name='eccodes', reset=False):
    """
    Set path to GRIB|ECCODES_DEFINITION_PATH.

    :param api_name: the name of the GRIB API, among ('eccodes', 'grib_api')
    :param reset: ignore predefined values of the variables
    """
    paths = [path]
    if api_name == 'grib_api':
        dp = 'GRIB_DEFINITION_PATH'
    elif api_name == 'eccodes':
        dp = 'ECCODES_DEFINITION_PATH'
    if not reset and os.environ.get(dp, False):
        paths.append(os.environ.get(dp))
    os.environ[dp] = os.pathsep.join(paths)


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


def parse_GRIBstr_todict(strfid):
    """Parse and return a dict GRIB fid from a string."""
    fid = str2dict(strfid, try_convert=int)
    return fid


def read_gribdef(filename):
    """Read a grib definition file and return it as a dict."""
    re_name = re.compile('("|\')(?P<name>[\w\.\-\_ ]+)("|\')\s*=')
    re_real = '\+|-?\d*\.\d*e\+|-?\d*'
    re_real_g = '(?P<real>' + re_real + ')'
    re_int = '\+|-?\d+'
    re_int_g = '(?P<int>' + re_int + ')'
    # re_num = '(' + re_real + ')|(' + re_int + ')'
    re_num_g = re_int_g + '|' + re_real_g
    re_keyvalue = re.compile('(?P<key>\w+)\s*=\s*' + re_num_g + '\s*;')
    # read file
    with io.open(filename, 'r') as f:
        lines_unfold = [l.strip() for l in f.readlines()]
        lines = []
        for line in lines_unfold:
            line = line.replace(';',';\n').replace('{ ','{\n')
            line_split = [l.strip() for l in line.split('\n')]
            lines.extend(line_split)
    # find fields declaration
    dico = {}
    indexes = []
    for i, line in enumerate(lines):
        fmatch = re_name.match(line)
        if fmatch:
            field = fmatch.group('name')
            indexes.append((field, i))
            dico[field] = {}
    # loop on fields
    for (j, (field, istart)) in enumerate(indexes):
        if j + 1 == len(indexes):  # last one
            iend = len(lines)
        else:
            iend = indexes[j + 1][1]
        if istart > 0 and lines[istart - 1].startswith("#"):  # this is a comment
            dico[field]['#comment'] = lines[istart - 1][1:].strip().strip('"')
        for i in range(istart, iend):
            kvmatch = re_keyvalue.match(lines[i])
            if kvmatch:
                if kvmatch.groupdict().get('int'):
                    dico[field][kvmatch.group('key')] = int(kvmatch.group('int'))
                elif kvmatch.groupdict().get('real'):
                    dico[field][kvmatch.group('key')] = float(kvmatch.group('real'))
    return dico
