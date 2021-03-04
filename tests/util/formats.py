#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from epygram import config

#: suffixes of files in datadir for each format
suffixes = {'FA': 'fa',
            'GRIB1':'grb1',
            'GRIB2':'grb2',
            'netCDF':'nc',
            'netCDFMNH':'ncMNH',
            'LFI':'lfi'}

#: some formats have more "loss" of precision when r/w/r
tolerances = {'netCDF':2e-12,
              'netCDFMNH':2e-12,
              'GRIB':1e-12}


def activate(*formats):
    """
    Activate format(s) only if it is in epygram.config.implemented_formats.
    """
    fmap = {'GRIB1':'GRIB', 'GRIB2':'GRIB'}
    return [f for f in formats if fmap.get(f, f) in config.implemented_formats]
