#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import main, skipIf

import epygram
from epygram import epylog

from . import abstract_testclasses as abtc

epygram.init_env()
epylog.setLevel('WARNING')
skip_LFI = False


@skipIf('FA' not in epygram.config.implemented_formats, "format not activated")
class TestFA(abtc.TestFMT):
    basename = 'FA'
    len = 134


@skipIf('GRIB' not in epygram.config.implemented_formats, "format not activated")
class TestGRIB(abtc.TestFMT):
    basename = 'GRIB'
    len = 3


@skipIf('netCDF' not in epygram.config.implemented_formats, "format not activated")
class TestnetCDF(abtc.TestFMT):
    basename = 'netCDF'
    len = 9


@skipIf('DDHLFA' not in epygram.config.implemented_formats, "format not activated")
class TestDDHLFA(abtc.TestFMT):
    basename = 'ddh.LFA'
    len = 160


@skipIf('LFI' not in epygram.config.implemented_formats, "format not activated")
@skipIf(skip_LFI, "LFI testing not ready yet (and greedy)")
class TestLFI(abtc.TestFMT):
    basename = 'LFI'
    len = 160


if __name__ == '__main__':
    main(verbosity=2)
