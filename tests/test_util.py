#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import TestCase, main
from numpy import pi
import datetime

from epygram import util


class TestAngle(TestCase):

    def test_conversion(self):
        a = util.Angle(0., 'degrees')
        for u in a.units:
            a.get(u)
        a = util.Angle(0., 'radians')
        for u in a.units:
            a.get(u)
        a = util.Angle((0., 0., 0.), 'DMS')
        for u in a.units:
            a.get(u)
        a = util.Angle((1., 0.), 'cos_sin')
        for u in a.units:
            a.get(u)

    def test_equal(self):
        self.assertEqual(util.Angle(0., 'radians'),
                         util.Angle(0., 'degrees'))
        self.assertEqual(util.Angle(0., 'radians'),
                         util.Angle((0., 0., 0.), 'DMS'))
        self.assertEqual(util.Angle(0., 'radians'),
                         util.Angle((1., 0.), 'cos_sin'))


class TestFunctions(TestCase):

    def test_degrees_nearest_mod(self):
        self.assertEqual(util.degrees_nearest_mod(361., 0.),
                         1.)
        self.assertEqual(util.degrees_nearest_mod(-179., 180.),
                         181.)
        self.assertEqual(util.degrees_nearest_mod(361., 180.),
                         1.)
        self.assertEqual(util.degrees_nearest_mod(359., 180.),
                         359.)
        self.assertEqual(util.degrees_nearest_mod(361., -180.),
                         1.)

    def test_positive_longitude(self):
        self.assertEqual(util.positive_longitude(-1., 'degrees'),
                         359.)
        self.assertEqual(util.positive_longitude(-pi, 'radians'),
                         pi)

    def test_datetimes2fieldvaliditylist(self):
        util.datetimes2fieldvaliditylist(datetime.datetime(2000, 1, 1))
        util.datetimes2fieldvaliditylist(datetime.datetime(2000, 1, 1),
                                         basis=datetime.datetime(1999, 12, 31))
        util.datetimes2fieldvaliditylist([datetime.datetime(2000, 1, 1),
                                          datetime.datetime(2000, 1, 2)],
                                         basis=datetime.datetime(1999, 12, 31))
        util.datetimes2fieldvaliditylist([datetime.datetime(2000, 1, 1),
                                          datetime.datetime(2000, 1, 2)],
                                         basis=[datetime.datetime(1999, 12, 31),
                                                datetime.datetime(2000, 1, 1)])


if __name__ == '__main__':
    main(verbosity=2)
