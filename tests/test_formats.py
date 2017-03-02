#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import main
import os

from . import abstract_testclasses

datadir = 'data/formats'


class TestFA(abstract_testclasses.TestFMT):
    filename = os.path.join(datadir, 'FA')
    len = 134


class TestGRIB(abstract_testclasses.TestFMT):
    filename = os.path.join(datadir, 'GRIB')
    len = 3


class TestNetCDF(abstract_testclasses.TestFMT):
    filename = os.path.join(datadir, 'netCDF')
    len = 9


class TestDDHLFA(abstract_testclasses.TestFMT):
    filename = os.path.join(datadir, 'ddh.LFA')
    len = 160


if __name__ == '__main__':
    main(verbosity=2)
