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
fast = True  # skips some slow tests
basemap_ok = True


# H2D
#####
class Test_gauss(abtc.Test_GeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'gauss.fa'

    def test_name(self):
        self.assertTrue('gauss' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth((2, 45), (16, 35)),
                               128.28873612363486,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance((2, 45), (16, 35)),
                               1626146.444843353,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqual(lons[0, 0],
                               2.578310078088704,
                               delta=epsilon)
        self.assertAlmostEqual(lats[0, 0],
                               46.628102031648154,
                               delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         181724)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(24, 36),
                                  (-4.1688266426944685, 52.818404733546714),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace((2, 45), (16, 35), 4)).flatten(),
            numpy.array([(2, 45),
                         (7.138076483278202, 41.86030120579446),
                         (11.780188209507074, 38.513496711778785),
                         (16, 35)]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertEqual(self.geo.ll2ij(2, 38),
                         (147, 40))

    @skipIf(fast, 'slow test')
    @skipIf(not basemap_ok, 'need basemap module')
    def _test_make_basemap(self, **kwargs):
        from mpl_toolkits.basemap import Basemap
        self.assertIsInstance(self.geo.make_basemap(**kwargs),
                              Basemap)

    def test_make_basemap(self):
        self._test_make_basemap()

    def test_make_basemap_opts1(self):
        self._test_make_basemap(specificproj=('nsper', {}))

    def test_make_basemap_opts2(self):
        self._test_make_basemap(specificproj=('nsper', {'lon':2.,
                                                        'lat':45.,
                                                        'sat_height':600.}))

    def test_make_basemap_opts3(self):
        self._test_make_basemap(specificproj='ortho')

    def test_make_basemap_opts4(self):
        self._test_make_basemap(zoom={'lonmin':-10, 'lonmax':15,
                                      'latmin':35.5, 'latmax':53})

    def test_make_point_geometry(self):
        self.assertIsInstance(self.geo.make_point_geometry(2, 45),
                              epygram.geometries.PointGeometry)

    def test_make_profile_geometry(self):
        self.assertIsInstance(self.geo.make_profile_geometry(2, 45),
                              epygram.geometries.V1DGeometry)

    def test_make_section_geometry(self):
        self.assertIsInstance(self.geo.make_section_geometry((2, 45), (3, 46)),
                              epygram.geometries.V2DGeometry)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(120, -45),
                               0.46747883700739834,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def _test_nearest_points(self, request, expected):
        self.assertEqual(self.geo.nearest_points(120, -45, request=request),
                         expected)

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
        self.assertAlmostEqual(self.geo.resolution_ll(120, -45),
                               99924.900787803344, delta=epsilon)
        self.assertAlmostEqual(self.geo.resolution_ll(2, 45),
                               14082.872826968507, delta=epsilon)


class Test_lambert_HS(abtc.Test_GeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'lambert_HS.fa'

    def test_name(self):
        self.assertTrue('lambert' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth((-73, -55), (-72, -54)),
                               35.889107796530737,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               6.7197722043960928,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance((-73, -55), (-72, -54)),
                               128585.10421692148,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqual(lons[0, 0],
                               - 88.772039618271094,
                               delta=epsilon)
        self.assertAlmostEqual(lats[0, 0],
                               - 60.859390190567119,
                               delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(24, 36),
                                  (-84.95623428711568, -58.475140869628035),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace((-73, -55), (-72, -54), 4)).flatten(),
            numpy.array([(-73., -55.),
                         (-72.661532609999995, -54.667694949999998),
                         (-72.328234140000006, -54.334351269999999),
                         (-72., -54.)]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(-73, -55),
                                  (116.36726532276415, 83.253764915097392),
                                  delta=epsilon)

    @skipIf(fast, 'slow test')
    @skipIf(not basemap_ok, 'need basemap module')
    def _test_make_basemap(self, **kwargs):
        from mpl_toolkits.basemap import Basemap
        self.assertIsInstance(self.geo.make_basemap(**kwargs),
                              Basemap)

    def test_make_basemap(self):
        self._test_make_basemap()

    def test_make_basemap_opts1(self):
        self._test_make_basemap(specificproj=('nsper', {}))

    def test_make_basemap_opts2(self):
        self._test_make_basemap(specificproj=('nsper', {'lon':-73.,
                                                        'lat':-55.,
                                                        'sat_height':600.}))

    def test_make_basemap_opts3(self):
        self._test_make_basemap(specificproj='ortho')

    def test_make_basemap_opts4(self):
        self._test_make_basemap(zoom={'lonmin':-80, 'lonmax':-65,
                                      'latmin':-65, 'latmax':-42})

    def test_make_point_geometry(self):
        self.assertIsInstance(self.geo.make_point_geometry(-73, -55),
                              epygram.geometries.PointGeometry)

    def test_make_profile_geometry(self):
        self.assertIsInstance(self.geo.make_profile_geometry(-73, -55),
                              epygram.geometries.V1DGeometry)

    def test_make_section_geometry(self):
        self.assertIsInstance(self.geo.make_section_geometry((-73, -55), (-72, -54)),
                              epygram.geometries.V2DGeometry)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(-55),
                               1.003964261558473,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def _test_nearest_points(self, request, expected):
        self.assertEqual(self.geo.nearest_points(-73, -55, request=request),
                         expected)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (116, 83))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(116, 83), (116, 84),
                                                (117, 83), (117, 84)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':8500}, [(116, 83), (116, 84),
                                                    (117, 83), (117, 84)])

    def test_point_is_inside_domain_ll(self):
        self.assertFalse(self.geo.point_is_inside_domain_ll(2, 45))
        self.assertTrue(self.geo.point_is_inside_domain_ll(-73, -55))

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(-73, -55),
                               7967.6821594830908, delta=epsilon)


class Test_mercator(abtc.Test_GeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'mercator_HN.fa'

    def test_name(self):
        self.assertTrue('mercator' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth((-61, 16), (-62, 17)),
                               - 43.795221249130975,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               0.,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance((-61, 16), (-62, 17)),
                               154053.58670729026,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqual(lons[0, 0],
                               - 81.51099961956092,
                               delta=epsilon)
        self.assertAlmostEqual(lats[0, 0],
                               - 6.323778282067483,
                               delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(24, 36),
                                  (-78.05772877762075, -1.1566717493466303),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace((-61, 16), (-62, 17), 4)).flatten(),
            numpy.array([(-61, 16),
                         (-61.333333330000002, 16.333906070000001),
                         (-61.666666669999998, 16.667242829999999),
                         (-62, 17)]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(-61, 16),
                                  (142.55006728429416, 156.71234418713408),
                                  delta=epsilon)

    @skipIf(fast, 'slow test')
    @skipIf(not basemap_ok, 'need basemap module')
    def _test_make_basemap(self, **kwargs):
        from mpl_toolkits.basemap import Basemap
        self.assertIsInstance(self.geo.make_basemap(**kwargs),
                              Basemap)

    def test_make_basemap(self):
        self._test_make_basemap()

    def test_make_basemap_opts1(self):
        self._test_make_basemap(specificproj=('nsper', {}))

    def test_make_basemap_opts2(self):
        self._test_make_basemap(specificproj=('nsper', {'lon':-61.,
                                                        'lat':16.,
                                                        'sat_height':600.}))

    def test_make_basemap_opts3(self):
        self._test_make_basemap(specificproj='ortho')

    def test_make_basemap_opts4(self):
        self._test_make_basemap(zoom={'lonmin':-65, 'lonmax':-55,
                                      'latmin':10, 'latmax':20})

    def test_make_point_geometry(self):
        self.assertIsInstance(self.geo.make_point_geometry(-61, 16),
                              epygram.geometries.PointGeometry)

    def test_make_profile_geometry(self):
        self.assertIsInstance(self.geo.make_profile_geometry(-61, 16),
                              epygram.geometries.V1DGeometry)

    def test_make_section_geometry(self):
        self.assertIsInstance(self.geo.make_section_geometry((-61, 16), (-62, 17)),
                              epygram.geometries.V2DGeometry)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(16),
                               1.040299435861602,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def _test_nearest_points(self, request, expected):
        self.assertEqual(self.geo.nearest_points(-61, 16, request=request),
                         expected)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (143, 157))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(142, 156), (142, 157),
                                                (143, 156), (143, 157)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':18500}, [(142, 156), (142, 157),
                                                     (143, 156), (143, 157)])

    def test_point_is_inside_domain_ll(self):
        self.assertFalse(self.geo.point_is_inside_domain_ll(2, 45))
        self.assertTrue(self.geo.point_is_inside_domain_ll(-61, 16))

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(-61, 16),
                               15371.771945773886, delta=epsilon)


class Test_stereopol(abtc.Test_GeometryMethods):
    fid = 'SURFGEOPOTENTIEL'
    basename = 'stereopol_HN.fa'

    def test_name(self):
        self.assertTrue('polar_stereographic' in self.geo.name)

    def test_azimuth(self):
        self.assertAlmostEqual(self.geo.azimuth((-40, 80), (-41, 81)),
                               - 3.8723909400079606,
                               delta=epsilon)

    def test_compass_grid(self):
        self.assertAlmostEqual(self.geo.compass_grid()[0, 0],
                               - 18.423814350311183,
                               delta=epsilon)

    def test_distance(self):
        self.assertAlmostEqual(self.geo.distance((-40, 80), (-41, 81)),
                               112699.14287003415,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_get_lonlat_grid(self):
        lons, lats = self.geo.get_lonlat_grid()
        self.assertAlmostEqual(lons[0, 0],
                               - 53.423814350311183,
                               delta=epsilon)
        self.assertAlmostEqual(lats[0, 0],
                               - 6.323778282067483,
                               delta=epsilon)

    def test_gridpoints_number(self):
        self.assertEqual(self.geo.gridpoints_number,
                         102400)

    def test_ij2ll(self):
        self.assertAlmostEqualSeq(self.geo.ij2ll(24, 36),
                                  (-51.28805020370288, 48.00255024785442),
                                  delta=epsilon)

    def test_linspace(self):
        self.assertAlmostEqualSeq(
            numpy.array(self.geo.linspace((-40, 80), (-41, 81), 4)).flatten(),
            numpy.array([(-40, 80),
                         (-40.310237149999999, 80.333484839999997),
                         (-40.642749369999997, 80.666831430000002),
                         (-41, 81)]).flatten(),
            delta=epsilon)

    def test_ll2ij(self):
        self.assertAlmostEqualSeq(self.geo.ll2ij(-40, 80),
                                  (103.66322735842627, 260.00997643994754),
                                  delta=epsilon)

    @skipIf(fast, 'slow test')
    @skipIf(not basemap_ok, 'need basemap module')
    def _test_make_basemap(self, **kwargs):
        from mpl_toolkits.basemap import Basemap
        self.assertIsInstance(self.geo.make_basemap(**kwargs),
                              Basemap)

    def test_make_basemap(self):
        self._test_make_basemap()

    def test_make_basemap_opts1(self):
        self._test_make_basemap(specificproj=('nsper', {}))

    def test_make_basemap_opts2(self):
        self._test_make_basemap(specificproj=('nsper', {'lon':-61.,
                                                        'lat':16.,
                                                        'sat_height':600.}))

    def test_make_basemap_opts3(self):
        self._test_make_basemap(specificproj='ortho')

    def test_make_basemap_opts4(self):
        self._test_make_basemap(zoom={'lonmin':-65, 'lonmax':-55,
                                      'latmin':10, 'latmax':20})

    def test_make_point_geometry(self):
        self.assertIsInstance(self.geo.make_point_geometry(-40, 80),
                              epygram.geometries.PointGeometry)

    def test_make_profile_geometry(self):
        self.assertIsInstance(self.geo.make_profile_geometry(-40, 80),
                              epygram.geometries.V1DGeometry)

    def test_make_section_geometry(self):
        self.assertIsInstance(self.geo.make_section_geometry((-40, 80), (-41, 81)),
                              epygram.geometries.V2DGeometry)

    def test_map_factor(self):
        self.assertAlmostEqual(self.geo.map_factor(80),
                               1.0076542662455523,
                               delta=epsilon)

    @skipIf(fast, 'slow test')
    def test_map_factor_field(self):
        self.assertIsInstance(self.geo.map_factor_field(),
                              epygram.fields.H2DField)

    def _test_nearest_points(self, request, expected):
        self.assertEqual(self.geo.nearest_points(-40, 80, request=request),
                         expected)

    def test_nearest_points1(self):
        self._test_nearest_points({'n':'1'}, (104, 260))

    def test_nearest_points2(self):
        self._test_nearest_points({'n':'2*2'}, [(103, 260), (103, 261),
                                                (104, 260), (104, 261)])

    def test_nearest_points3(self):
        self._test_nearest_points({'radius':19500}, [(103, 260), (103, 261),
                                                     (104, 260), (104, 261)])

    def test_point_is_inside_domain_ll(self):
        self.assertFalse(self.geo.point_is_inside_domain_ll(22, 45))
        self.assertTrue(self.geo.point_is_inside_domain_ll(-40, 80))

    def test_resolution_ll(self):
        self.assertAlmostEqual(self.geo.resolution_ll(-40, 80),
                               15876.793533343482, delta=epsilon)


if __name__ == '__main__':
    main(verbosity=2)
