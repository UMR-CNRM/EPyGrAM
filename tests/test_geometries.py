#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import main, skipIf
import os
import tempfile

import epygram

from . import abstract_testclasses as abtc
from . import util
from . import datadir

epygram.init_env()

fast = False  # skips some slow tests


# H2D
#####
fmts = ['FA', 'GRIB2', 'netCDF']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@skipIf(fast, 'slow test')
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_gaussC1(abtc.Test_H2DGeometry):
    fileprefix = 'gaussC1'


fmts = ['FA', 'GRIB2', 'netCDF']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@skipIf(fast, 'slow test')
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_gaussC2p4(abtc.Test_H2DGeometry):
    fileprefix = 'gaussC2.4'


fmts = ['FA', 'GRIB2', 'netCDF', 'LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_lambertHN(abtc.Test_H2DGeometry):
    fileprefix = 'lambert_HN'


fmts = ['FA', 'GRIB2', 'netCDF', 'LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_lambertHS(abtc.Test_H2DGeometry):
    fileprefix = 'lambert_HS'


fmts = ['LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_lambertHS_rot_sec(abtc.Test_H2DGeometry):
    fileprefix = 'lambert_HS_rotated_secant'


fmts = ['FA', 'GRIB1', 'GRIB2', 'netCDF', 'LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_mercatorHN(abtc.Test_H2DGeometry):
    fileprefix = 'mercator_HN'


fmts = ['FA', 'GRIB1', 'GRIB2', 'netCDF']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_mercatorHS(abtc.Test_H2DGeometry):
    fileprefix = 'mercator_HS'


fmts = ['LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_mercatorHN_rot_sec(abtc.Test_H2DGeometry):
    fileprefix = 'mercator_HN_rotated_secant'


fmts = ['FA', 'GRIB2', 'netCDF', 'LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_stereopolHN(abtc.Test_H2DGeometry):
    fileprefix = 'stereopol_HN'


fmts = ['FA', 'GRIB2', 'netCDF', 'LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_stereopolHS(abtc.Test_H2DGeometry):
    fileprefix = 'stereopol_HS'


fmts = ['LFI']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_stereopolHS_rot_sec(abtc.Test_H2DGeometry):
    fileprefix = 'stereopol_HS_rotated_secant'


fmts = ['FA', 'GRIB1', 'GRIB2', 'netCDF']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_regLLsmall(abtc.Test_H2DGeometry):
    fileprefix = 'regLL_small'


fmts = ['FA', 'GRIB1', 'GRIB2', 'netCDF']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_regLLlarge(abtc.Test_H2DGeometry):
    fileprefix = 'regLL_large'


# Special ones
@skipIf(len(util.activate('DDHLFA')) == 0, "format not activated")
@util.add_tests_for_attrs('gridpoint', 'zonalband')
class Test_DDHLFA_PointGeometry(abtc.Test_DDHLFA_Geometry):
    fid_to_test = 'SVGFS01'
    geom = 'point'
    update_pickle = False  # set to True if necessary to make a new pickle


@skipIf(len(util.activate('DDHLFA')) == 0, "format not activated")
@util.add_tests_for_attrs('gridpoint', 'zonalband')
class Test_DDHLFA_V1DGeometry(abtc.Test_DDHLFA_Geometry):
    fid_to_test = 'VCT1'
    geom = 'profile'
    update_pickle = False  # set to True if necessary to make a new pickle


netCDF_Ndimensions = ['0D', 'V1D', 'V2D', 'H2D', '3D', '4D']  # , '5D']


@skipIf(len(util.activate('netCDF')) == 0, "format not activated")
@util.add_tests_for_attrs(*netCDF_Ndimensions)
class Test_netCDF_Ndimensions(abtc.Test_GeometryInterfaces):

    fid_to_test = 'temperature'
    dims = netCDF_Ndimensions
    indexes_to_check = [(i, j, k, l)
                        for i in (0, -1) for j in (0, -1)
                        for k in (0, -1) for l in (0, -1)]
    checks = {'0D':[296.37378837670036] * 16,
              'V1D':([232.92738732086275] * 4 +
                     [297.33992688752369] * 4) * 2,
              'V2D':([232.08349648457917,
                      231.25490097208404] * 2 +
                     [296.55241455590806,
                      299.92105230868816] * 2) * 2,
              'H2D':[296.95550432206772,
                     295.13011625810094,
                     299.80259631674028,
                     299.88585252134163] * 4,
              '3D':[233.38263541213618,
                    233.19968143561593,
                    233.44929161397826,
                    233.18015036341166,
                    297.41617052305361,
                    296.40097719165516,
                    298.95071663196188,
                    299.47477603703089] * 2,
              '4D':[233.50710303258333, 233.65517138068839,
                    232.54884073249218, 232.16818976450023,
                    297.15752229063031, 296.69304743371845,
                    299.11282903070833, 299.46143782920564,
                    233.29663862699829, 233.15782454347234,
                    233.52239611030629, 233.57453276726824,
                    296.89817229190237, 295.62418123637462,
                    298.88509692631214, 299.43758821205216]
              }

    @classmethod
    def setUpClass(cls):
        cls.filename = {d:os.path.join(cls.datadir,
                                       '_'.join(['dims', d])) + '.nc'
                        for d in cls.dims}

    def _test(self, dims):
        filename = os.path.join(self.datadir,
                                '_'.join(['dims', dims])) + '.nc'
        with epygram.formats.resource(filename, 'r') as r:
            f_r = r.readfield(self.fid_to_test)
            checks = [f_r.getdata(d4=True)[i] for i in self.indexes_to_check]
            self.assertEqual(checks,
                             self.checks[dims],
                             str(checks))
        self._test_rwr(filename, self.fid_to_test)


@skipIf(len(util.activate('netCDF')) == 0, "format not activated")
@util.add_tests_for_attrs('A', 'H', 'P', 'hybridP')
class Test_netCDF_VGeometries(abtc.Test_GeometryInterfaces):

    fid_to_test = 'temperature'
    vdict = {'A':102,
             'H':103,
             'P':100,
             'hybridP':119,
             'hybridH':118}
    vgridlen = {'A':90,
                'H':90,
                'P':90,
                'hybridP':91,
                'hybridH':None}

    def _test(self, typ):
        filename = os.path.join(self.datadir,
                                '_'.join(['VGeometry', typ])) + '.nc'
        with epygram.formats.resource(filename, 'r') as r:
            f_r = r.readfield(self.fid_to_test)
            self.assertEqual(f_r.geometry.vcoordinate.typeoffirstfixedsurface,
                             self.vdict[typ])
            self.assertEqual(len(f_r.geometry.vcoordinate.grid['gridlevels']),
                             self.vgridlen[typ])
        self._test_rwr(filename, self.fid_to_test)


fmts = ['FA', 'GRIB2']
@skipIf(len(util.activate(*fmts)) == 0, "format not activated")
@util.add_tests_for_attrs(*util.activate(*fmts))
class Test_VGeometry_hybridP(abtc.Test_GeometryInterfaces):

    fid_to_test = {'FA':'S090TEMPERATURE',
                   'GRIB2':{'shortName':'t',
                            'typeOfFirstFixedSurface':119,
                            'level':90}}

    def _test(self, fmt):
        filename = os.path.join(self.datadir,
                                '.'.join(['VGeometry_hybridP',
                                          abtc.suffixes[fmt]]))
        with epygram.formats.resource(filename, 'r') as r:
            f_r = r.readfield(self.fid_to_test[fmt])
            self.assertEqual(f_r.geometry.vcoordinate.typeoffirstfixedsurface,
                             119)
            self.assertEqual(len(f_r.geometry.vcoordinate.grid['gridlevels']),
                             91)
        self._test_rwr(filename, self.fid_to_test[fmt])


@skipIf(len(util.activate('FA')) == 0, "format not activated")
@util.add_tests_for_attrs('LAM', 'global')
class Test_SpectralGeometry(abtc.Test_GeometryInterfaces):

    fid = 'SPECSURFGEOPOTEN'

    def _test(self, geom):
        filename = os.path.join(self.datadir,
                                '_'.join(['spectralgeometry',
                                          geom])) + '.fa'
        with epygram.formats.resource(filename, 'r') as r:
            f_r = r.readfield(self.fid)
            tmp = tempfile.mktemp()
            out = epygram.formats.resource(tmp, 'w', fmt=r.format)
            out.writefield(f_r, compression=r.fieldscompression[self.fid])
            out.close()
            out.open(openmode='r')
            f_rwr = out.readfield(self.fid)
            out.close()
            os.remove(tmp)
            self.assertEqual(f_r.spectral_geometry, f_rwr.spectral_geometry,
                             '\n<>\n'.join([str(f_r.spectral_geometry),
                                            str(f_rwr.spectral_geometry)]))


@skipIf(len(util.activate('FA')) == 0, "format not activated")
@util.add_tests_for_attrs('gaussC1', 'lambert_HN')
class Test_GeoPointsFull_WR_fromFA(abtc.Test_GeoPoints_WR_fromFA):
    llv = False


@skipIf(len(util.activate('FA')) == 0, "format not activated")
@util.add_tests_for_attrs('gaussC1', 'lambert_HN')
class Test_GeoPointsLLV_WR_fromFA(abtc.Test_GeoPoints_WR_fromFA):
    llv = True


if __name__ == '__main__':
    main(verbosity=2)
