#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import TestCase, main
from numpy import pi
import matplotlib
import datetime
import os

import epygram

datadir = 'data/formats'
fmts = ['FA', 'GRIB', 'netCDF', 'ddh.LFA', 'LFI']
testfiles = {fmt:os.path.join(datadir, fmt)
             for fmt in fmts}


class TestFMT(TestCase):

    def _test_listfields(self, filename, n):
        with epygram.formats.resource(filename, 'r') as r:
            self.assertEqual(len(r.listfields()), n)

    def _test_what(self, filename):
        with epygram.formats.resource(filename, 'r') as r:
            with open('/dev/null', 'ab') as out:
                r.what(out=out)


class TestFA(TestFMT):

    filename = os.path.join(datadir, 'FA')
    len = 134

    def test_listfields(self):
        self._test_listfields(self.filename, self.len)

    def test_what(self):
        self._test_what(self.filename)


class TestGRIB(TestFMT):

    filename = os.path.join(datadir, 'GRIB')
    len = 3

    def test_listfields(self):
        self._test_listfields(self.filename, self.len)

    def test_what(self):
        self._test_what(self.filename)


class TestNetCDF(TestFMT):

    filename = os.path.join(datadir, 'netCDF')
    len = 9

    def test_listfields(self):
        self._test_listfields(self.filename, self.len)

    def test_what(self):
        self._test_what(self.filename)


class TestDDHLFA(TestFMT):

    filename = os.path.join(datadir, 'ddh.LFA')
    len = 160

    def test_listfields(self):
        self._test_listfields(self.filename, self.len)

    def test_what(self):
        self._test_what(self.filename)


if __name__ == '__main__':
    main(verbosity=2)
