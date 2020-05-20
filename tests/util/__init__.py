#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from epygram import config


def _build_test_for_attr(attr):
    def test_(obj):
        obj._test(attr)
    test_.__name__ = str('test_{:s}'.format(attr))
    return test_


def _build_test_for_mtd(mtd, args, kwargs, assertion, expected, suffix=''):
    def test_(obj):
        obj._test_mtd(mtd, args, kwargs, assertion, expected)
    test_.__name__ = str('test_{:s}'.format(mtd))
    if suffix != '':
        test_.__name__ += suffix
    return test_


def add_tests_for_attrs(*attrs):
    def decocls(cls):
        for attr in attrs:
            tst = _build_test_for_attr(attr)
            setattr(cls, tst.__name__, tst)
        return cls
    return decocls


def add_tests_for_mtd(mtd, args, kwargs, assertion, expected, suffix=''):
    def decocls(cls):
        tst = _build_test_for_mtd(mtd, args, kwargs, assertion, expected,
                                  suffix=suffix)
        setattr(cls, tst.__name__, tst)
        return cls
    return decocls


def activate(*formats):
    """
    Activate format(s) only if it is in epygram.config.implemented_formats.
    """
    fmap = {'GRIB1':'GRIB', 'GRIB2':'GRIB'}
    return [f for f in formats if fmap.get(f, f) in config.implemented_formats]


datadir = datadir = './data'
suffixes = {'FA': 'fa',
            'GRIB1':'grb1',
            'GRIB2':'grb2',
            'netCDF':'nc',
            'netCDFMNH':'ncMNH',
            'LFI':'lfi'}

delta_assertAlmostEqual = 1e-12
delta_assertAlmostEqual4pyproj = 1e-8  # because pyproj not reproducible between python2 and python3
