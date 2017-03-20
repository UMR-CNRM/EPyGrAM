#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import TestCase
import tempfile
import os
import cPickle

import epygram

from . import datadir
from .util import suffixes, delta_assertAlmostEqual


# CLASSES
#########
class TestFMT(TestCase):

    datadir = os.path.join(datadir, 'formats')
    basename = ''  # to be redefined in real classes
    len = 0  # to be redefined in real classes

    @property
    def fmt(self):
        return self.__class__.__name__[4:]

    def setUp(self):
        self.filename = os.path.join(self.datadir, self.basename)

    def tearDown(self):
        del self.filename

    def test_listfields(self):
        with epygram.formats.resource(self.filename, 'r', fmt=self.fmt) as r:
            self.assertEqual(len(r.listfields()), self.len)

    def test_what(self):
        with epygram.formats.resource(self.filename, 'r', fmt=self.fmt) as r:
            with open('/dev/null', 'ab') as out:
                r.what(out=out)


class Test_GeometryInterfaces(TestCase):

    datadir = os.path.join(datadir, 'geometries')

    def _test_rwr(self, filename, fid):
        """Generic test Read/Write/Read and check identity."""
        with epygram.formats.resource(filename, 'r') as r:
            f_r = r.readfield(fid)
            tmp = tempfile.mktemp()
            out = epygram.formats.resource(tmp, 'w', fmt=r.format)
            out.writefield(f_r)
            out.close()
            out.open(openmode='r')
            f_rwr = out.readfield(fid)
            out.close()
            os.remove(tmp)
            self.assertEqual(f_r.geometry, f_rwr.geometry,
                             '\n<>\n'.join([str(f_r.geometry),
                                            str(f_rwr.geometry)]))

    def _test_pickled(self, filename, fid, picklename):
        """Generic test whatever the actual domain type."""
        with epygram.formats.resource(filename, 'r') as r:
            fld = r.readfield(fid)
            if isinstance(fld, epygram.base.FieldSet):
                fld = fld[0]
            if not self.update_pickle:
                with open(picklename, 'rb') as pckl:
                    pt_pckld = cPickle.load(pckl)
                self.assertEqual(fld, pt_pckld)
            else:
                with open(picklename, 'wb') as pckl:
                    cPickle.dump(fld, pckl)


class Test_H2DGeometry(Test_GeometryInterfaces):

    fid_to_test = {'FA':'SURFGEOPOTENTIEL',
                   'GRIB1':{'indicatorOfParameter':6},
                   'GRIB2':{'shortName':'orog'},
                   'netCDF':'SURFGEOPOTENTIEL',
                   'LFI':('ZS', 0)}
    fileprefix = ''  # to be defined in real classes

    def _test(self, fmt):
        filename = os.path.join(self.datadir,
                                '.'.join([self.fileprefix,
                                          suffixes[fmt]]))
        self._test_rwr(filename, self.fid_to_test[fmt])


class Test_DDHLFA_Geometry(Test_GeometryInterfaces):

    ddh_domaintype = ['gridpoint', 'zonalband']
    fid_to_test = ''  # will define the actual type of geom to test, Point/V1D
    geom = None  # to be defined in actual classes: point/profile
    update_pickle = False  # set to True if necessary to make a new pickle

    def _test(self, domain):
        filename = os.path.join(self.datadir,
                                '_'.join([self.geom, domain])) + '.ddhlfa'
        picklename = '_'.join([filename, self.fid_to_test]) + '.cPickle'
        self._test_pickled(filename,
                           self.fid_to_test,
                           picklename)


class Test_GeometryMethods(TestCase):

    fid = ''
    basename = ''
    datadir = os.path.join(datadir, 'geometry_methods')

    def setUp(self):
        filename = os.path.join(self.datadir, self.basename)
        r = epygram.formats.resource(filename, 'r')
        self.fld = r.readfield(self.fid)
        self.geo = self.fld.geometry
        del r

    def tearDown(self):
        del self.fld
        del self.geo

    def assertAlmostEqualSeq(self, first, second, delta):
        self.assertEqual(len(first), len(second))
        for i in range(len(first)):
            self.assertAlmostEqual(first[i], second[i], delta=delta)

    def _test_mtd(self, mtd, args, kwargs, assertion, expected):
        out = getattr(self.geo, mtd)(*args, **kwargs)
        if assertion == 'Equal':
            self.assertEqual(out, expected)
        elif assertion == 'AlmostEqual':
            self.assertAlmostEqual(out, expected, delta=delta_assertAlmostEqual)
        elif assertion == 'IsInstance':
            self.assertIsInstance(out, expected)
