#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals


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


suffixes = {'FA': 'fa',
            'GRIB1':'grb1',
            'GRIB2':'grb2',
            'netCDF':'nc',
            'LFI':'lfi'}

delta_assertAlmostEqual = 1e-12
