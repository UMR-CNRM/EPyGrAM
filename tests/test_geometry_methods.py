#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import main, skipIf
import numpy

import epygram

from . import abstract_testclasses as abtc
from .util import delta_assertAlmostEqual

epygram.init_env()

epsilon = delta_assertAlmostEqual
fast = False  # skips some slow tests
basemap_ok = True


# H2D
#####
class Test_gauss(abtc.Test_GeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'gauss.fa'
    point = (120, -45)
    transect = ((2, 45), (16, 35))

    def test_name(self):
        self.assertTrue('gauss' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth(*self.transect),
                               128.28873612363486,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance(*self.transect),
                               1626146.444843353,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqualSeq((lons[0, 0], lats[0, 0]),
                                  (2.578310078088704, 46.628102031648154),
                                  delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         181724)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(*self.ij_test),
                                  (-4.1688266426944685, 52.818404733546714),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace(*self.transect, num=3)).flatten(),
            numpy.array([self.transect[0],
                         (9.5164400860259555, 40.210068310928158),
                         self.transect[1]]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertEqual(self.geo.ll2ij(*self.point),
                         (186, 323))

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(*self.point),
                               0.46747883700739834,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (186, 323))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(186, 323), (187, 323),
                                                (198, 322), (199, 322)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':85000}, [(186, 323), (187, 323),
                                                     (198, 322), (199, 322)])

    def test_point_is_inside_domain_ij(self):
        self.assertFalse(self.geo.point_is_inside_domain_ij(100, 0))
        self.assertTrue(self.geo.point_is_inside_domain_ij(0, 100))

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(*self.point),
                               99924.900787803344, delta=epsilon)


class Test_lambert_HS(abtc.Test_RectGeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'lambert_HS.fa'
    point = (-73, -55)
    transect = [(-73, -55), (-72, -54)]

    def test_name(self):
        self.assertTrue('lambert' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth(*self.transect),
                               35.889107796530737,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               6.7197722043960928,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance(*self.transect),
                               128585.10421692148,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqualSeq((lons[0, 0], lats[0, 0]),
                                  (-88.772039618271094, -60.859390190567119),
                                  delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(*self.ij_test),
                                  (-84.95623428711568, -58.475140869628035),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace(*self.transect, num=3)).flatten(),
            numpy.array([self.transect[0],
                         (-72.494243850000004, -54.501150989999999),
                         self.transect[1]]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(*self.point),
                                  (116.36726532276415, 83.253764915097392),
                                  delta=epsilon)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(self.point[1]),
                               1.003964261558473,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (116, 83))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(116, 83), (116, 84),
                                                (117, 83), (117, 84)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':8500}, [(116, 83), (116, 84),
                                                    (117, 83), (117, 84)])

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(*self.point),
                               7967.6821594830908, delta=epsilon)


class Test_mercator(abtc.Test_RectGeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'mercator_HN.fa'
    point = (-61, 16)
    transect = (-61, 16), (-62, 17)

    def test_name(self):
        self.assertTrue('mercator' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth(*self.transect),
                               - 43.795221249130975,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               0.,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance(*self.transect),
                               154053.58670729026,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqualSeq((lons[0, 0], lats[0, 0]),
                                  (-81.51099961956092, -6.323778282067483),
                                  delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(*self.ij_test),
                                  (-78.05772877762075, -1.1566717493466303),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace(*self.transect, num=3)).flatten(),
            numpy.array([self.transect[0],
                         (-61.5, 16.50064626),
                         self.transect[1]]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(*self.point),
                                  (142.55006728429416, 156.71234418713408),
                                  delta=epsilon)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(self.point[1]),
                               1.040299435861602,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (143, 157))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(142, 156), (142, 157),
                                                (143, 156), (143, 157)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':18500}, [(142, 156), (142, 157),
                                                     (143, 156), (143, 157)])

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(*self.point),
                               15371.771945773886, delta=epsilon)


class Test_stereopol(abtc.Test_RectGeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'stereopol_HN.fa'
    point = (40, 80)
    transect = ((-40, 80), (-41, 81))

    def test_name(self):
        self.assertTrue('polar_stereographic' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth(*self.transect),
                               - 3.8723909400079606,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               - 18.423814350311183,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance(*self.transect),
                               112699.14287003415,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqualSeq((lons[0, 0], lats[0, 0]),
                                  (-53.423814350311183, 42.887302601956691),
                                  delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(*self.ij_test),
                                  (-51.28805020370288, 48.00255024785442),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace(*self.transect, num=3)).flatten(),
            numpy.array([self.transect[0],
                         (-40.473562909999998, 80.500177789999995),
                         self.transect[1]]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(*self.point),
                                  (177.03804286989367, 311.38757537582376),
                                  delta=epsilon)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(self.point[1]),
                               1.0076542662455523,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (177, 311))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(177, 311), (177, 312),
                                                (178, 311), (178, 312)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':19500}, [(177, 311), (177, 312),
                                                     (178, 311), (178, 312)])

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(*self.point),
                               15876.558945361001, delta=epsilon)


class Test_regLL(abtc.Test_RectGeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'regLL_small.fa'
    point = (-2, 48)
    transect = ((0, 46), (1, 47))

    def test_name(self):
        self.assertTrue('regular_lonlat' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth(*self.transect),
                               45.0,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance(*self.transect),
                               134992.09125666993,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqualSeq((lons[0, 0], lats[0, 0]),
                                  (-8.0, 38.0),
                                  delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         481401)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(*self.ij_test),
                                  (-7.4, 38.9),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace(*self.transect, num=3)).flatten(),
            numpy.array([self.transect[0],
                         (0.5, 46.5),
                         self.transect[1]]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(*self.point),
                                  (240.0, 400.0),
                                  delta=epsilon)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (240, 400))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(240, 400), (240, 401),
                                                (241, 400), (241, 401)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':3000}, [(239, 400), (240, 399),
                                                    (240, 400), (240, 401),
                                                    (241, 400)])

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(*self.point),
                               1860.0973413526367, delta=epsilon)


if __name__ == '__main__':
    main(verbosity=2)
