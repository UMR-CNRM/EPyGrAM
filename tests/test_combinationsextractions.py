#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

import epygram
from footprints import proxy as fpx
import numpy
import matplotlib.pyplot as plt
import unittest
from unittest import main, TestCase, skipIf
import os
import re

epygram.init_env()

#Test configuration
test_plot = False  #To plot some figures during the test
doTest = range(18)  #list of tests to execute
terms = range(4)  #terms of the run to use
baseDir = "/cnrm/mesonh/data1/riette/Public/epygram/smallRunFa"  #Directory in which ICMSHAROM+* files are stored
coords_profile = (8.77, 41.9)
coords_section = (8.77, 41.9), (8.9, 41.4)
FA_name = 'S*TEMPERATURE'
GRIB_handgrip = {'discipline': 0,
                 'parameterCategory': 0,
                 'parameterNumber': 0,
                 'typeOfFirstFixedSurface': 119}

class TestMisc(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resources = [epygram.formats.resource(os.path.join(baseDir,
                                                               'ICMSHAROM+' + str(term).zfill(4)),
                                                  'r')
                         for term in terms]

    @classmethod
    def tearDownClass(cls):
        for r in cls.resources:
            r.close()

    @property
    def _CL_resources(self):
        return [fpx.resource_modificator(name='CombineLevels',
                                         resource=r,
                                         openmode='r')
                for r in self.resources]

    @property
    def _MV_resource(self):
        return fpx.resource_modificator(name='MultiValidities',
                                        resources=self.resources,
                                        openmode='r')

    @property
    def _CLMV_resource(self):
        # combinelevels(multivalidities)
        return fpx.resource_modificator(name='CombineLevels',
                                        resource=self._MV_resource,
                                        openmode='r')

    @property
    def _MVCL_resource(self):
        # multivalidities(combinelevels)
        return fpx.resource_modificator(name='MultiValidities',
                                        resources=self._CL_resources,
                                        openmode='r')

    def test_V2D_virtualfield(self):
        # Build a virtualField from H1D - one validity
        fields = epygram.base.FieldSet()
        section1 = self.resources[0].extractsection(FA_name,
                                                    *coords_section)
        for k in range(len(section1.geometry.vcoordinate.levels)):
            field = section1.getlevel(k=k)
            fields.append(field)
        section1_bis = fpx.field(fid={'FA':FA_name},
                                 structure='V2D',
                                 fieldset=fields)
        self.assertTrue(numpy.all(section1.getdata() == section1_bis.getdata()))
        self.assertEqual(section1.geometry, section1_bis.geometry)

    def test_V2DT_virtualfield(self):
        # Build a virtualField from H1D - several validities
        fields = epygram.base.FieldSet()
        resource4D = self._MVCL_resource
        sections = resource4D.extractsection(FA_name,
                                             *coords_section)
        for k in range(len(sections.geometry.vcoordinate.levels)):
            field = sections.getlevel(k=k)
            fields.append(field)
        sections_bis = fpx.field(fid={'FA':FA_name},
                                 structure='V2D',
                                 fieldset=fields)
        self.assertTrue(numpy.all(sections.getdata() == sections_bis.getdata()))
        self.assertEqual(sections.geometry, sections_bis.geometry)

    def test_V2DP_virtualfield(self):
        # Build a virtualField from H1D - one validity - sigma levels converted in P
        section1 = self.resources[0].extractsection(FA_name,
                                                    *coords_section,
                                                    vertical_coordinate=100)
        fields = epygram.base.FieldSet()
        for k in range(len(section1.geometry.vcoordinate.levels)):
            field = section1.getlevel(k=k)
            fields.append(field)
        section1_bis = fpx.field(fid={'FA':FA_name},
                                 structure='V2D',
                                 fieldset=fields)
        if test_plot:
            plot = section1_bis.plotfield()
            plt.show()
        self.assertTrue(numpy.all(section1.getdata() == section1_bis.getdata()))
        self.assertEqual(section1.geometry, section1_bis.geometry)

    def test_V2DTP_virtualfield_as_profiles(self):
        # Build a virtualField from H1D - several validities - sigma levels converted in P
        resource4D = self._MVCL_resource
        sections = resource4D.extractsection(FA_name,
                                             *coords_section,
                                             vertical_coordinate=100,
                                             interpolation='nearest')
        if test_plot:
            anim = sections.plotanimation(repeat=True,
                                          interval=500,
                                          zoom={'ymax':900, 'ymin':1020})
            plt.show()  # we can see mountains movement
        profiles = sections.as_profiles()[0]
        lon, lat = sections.geometry.get_lonlat_grid()
        profiles_bis = resource4D.extractprofile(FA_name,
                                                 lon[0], lat[0],
                                                 vertical_coordinate=100,
                                                 interpolation='nearest')
        self.assertTrue(numpy.all(profiles.getdata() == profiles_bis.getdata()))
        self.assertEqual(profiles.geometry, profiles_bis.geometry)
        if test_plot:
            anim = profiles.plotanimation(repeat=True,
                                          interval=500,
                                          zoom={'ymax':900, 'ymin':1020})
            plt.show()

    def test_V2DT_getvalidity(self):
        resource4D = self._MVCL_resource
        sections = resource4D.extractsection(FA_name,
                                             *coords_section,
                                             vertical_coordinate=100,
                                             interpolation='nearest')
        for i in range(len(self.resources)):
            sect1 = self.resources[i].extractsection(FA_name,
                                                     *coords_section,
                                                     vertical_coordinate=100,
                                                     interpolation='nearest')
            sect2 = sections.getvalidity(i)
            self.assertTrue(numpy.all(sect1.getdata() == sect2.getdata()))
            self.assertEqual(sect1.geometry, sect2.geometry)

    def test_V2DTP_virtualfield(self):
        resources4D = self._MVCL_resource
        sections = resources4D.extractsection(FA_name,
                                              *coords_section,
                                              vertical_coordinate=100,
                                              interpolation='nearest')
        fields = epygram.base.FieldSet()
        for k in range(len(sections.geometry.vcoordinate.levels)):
            field = sections.getlevel(k=k)
            fields.append(field)
        sections_bis = fpx.field(fid={'FA':FA_name},
                                 structure='V2D',
                                 fieldset=fields)
        self.assertTrue(numpy.all(sections.getdata() == sections_bis.getdata()))
        self.assertEqual(sections.geometry, sections_bis.geometry)

    @property
    def profile1(self):
        if hasattr(self, '_profile1'):
            profile1 = self._profile1
        else:
            profile1 = self.resources[0].extractprofile(FA_name,
                                                        *coords_profile)
            self._profile1 = profile1
        return profile1

    def test_V1D_SubdomainResource(self):
        # Test Subdomain, CombineLevels and D3VirtualField
        profile_resource = fpx.resource_modificator(name="Subdomain",
                                                    resource=self._CL_resources,
                                                    openmode='r',
                                                    geometry=self.profile1.geometry)
        profile2 = profile_resource.readfield(GRIB_handgrip)
        self.assertTrue(numpy.all(self.profile1.getdata() == profile2.getdata()))
        self.assertEqual(self.profile1.geometry, profile2.geometry)

    def test_V1D_virtualfield_extract_subdomain_resource(self):
        temp_fids = [fid for fid in self.resources[0].listfields()

                     if (fid[0:2] == 'S0' and fid[4:] == 'TEMPERATURE')]
        virtual3D = fpx.field(fid={'FA':FA_name},
                                  structure='3D',
                                  resource=self.resources[0],
                                  resource_fids=temp_fids)
        virtual3D.sp2gp()
        profile3 = virtual3D.extract_subdomain(self.profile1.geometry)
        self.assertTrue(numpy.all(self.profile1.getdata() == profile3.getdata()))
        self.assertEqual(self.profile1.geometry, profile3.geometry)

    @property
    def virtual3D(self):
        if hasattr(self, '_virtual3D'):
            virtual3D = self._virtual3D
        else:
            temp_fids = [fid for fid in self.resources[0].listfields()
                         if (fid[0:2] == 'S0' and fid[4:] == 'TEMPERATURE')]
            fields = epygram.base.FieldSet()
            for fid in temp_fids:
                field = self.resources[0].readfield(fid)
                fields.append(field)
            virtual3D = fpx.field(fid={'FA':FA_name},
                                  structure='3D',
                                  fieldset=fields)
            virtual3D.sp2gp()
            self._virtual3D = virtual3D
        return virtual3D

    def test_V1D_virtualfield_extract_subdomain_fields(self):
        profile3bis = self.virtual3D.extract_subdomain(self.profile1.geometry)
        self.assertTrue(numpy.all(self.profile1.getdata() == profile3bis.getdata()))
        self.assertEqual(self.profile1.geometry, profile3bis.geometry)

    def test_V1D_virtualfield_asreal(self):
        real_virtual3D = self.virtual3D.as_real_field()
        profile4 = real_virtual3D.extract_subdomain(self.profile1.geometry)
        self.assertTrue(numpy.all(self.profile1.getdata() == profile4.getdata()))
        self.assertEqual(self.profile1.geometry, profile4.geometry)

    def test_V1D_CLvirtual(self):
        virtual3D_resource = fpx.resource_modificator(name='CombineLevels',
                                                      resource=self.resources[0],
                                                      openmode='r',
                                                      virtual=True)
        profile5 = virtual3D_resource.readfield(GRIB_handgrip)
        profile5.sp2gp()
        profile5 = profile5.extract_subdomain(self.profile1.geometry)
        if test_plot:
            self.profile1.plotfield()
            plt.show()
        self.assertTrue(numpy.all(self.profile1.getdata() == profile5.getdata()))
        self.assertEqual(self.profile1.geometry, profile5.geometry)



if __name__ == '__main__':
    main(verbosity=2)
