#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy
from unittest import main, TestCase, skipIf

from footprints import proxy as fpx

import epygram
import matplotlib.pyplot as plt

epygram.init_env()

# Test configuration
test_plot = False  # To plot some figures during the test
terms = range(2)  # terms of the run to use
# Distant Directory in which ICMSHAROM+* files are stored:
baseDir = "/cnrm/phynh/data1/riette/Public/epygram/smallRunFa"
coords_profile = (8.77, 41.9)
coords_section = (8.77, 41.9), (8.9, 41.4)
FA_name = 'S*TEMPERATURE'
GRIB_handgrip = {'discipline': 0,
                 'parameterCategory': 0,
                 'parameterNumber': 0,
                 'typeOfFirstFixedSurface': 119}


@skipIf('FA' not in epygram.config.implemented_formats, "format not activated")
class TestMisc(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resources = [epygram.formats.resource('/'.join([baseDir,
                                                            'ICMSHAROM+{:04}'.format(term)]),
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

# Actual tests
##############
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
            section1_bis.plotfield()
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
            sections.plotanimation(repeat=True,
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
            profiles.plotanimation(repeat=True,
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
        if not hasattr(self, '_profile1'):
            self._profile1 = self.resources[0].extractprofile(FA_name,
                                                              *coords_profile)
        return self._profile1

    @property
    def section1(self):
        if not hasattr(self, '_section1'):
            self._section1 = self.resources[0].extractsection(FA_name,
                                                              *coords_section)
        return self._section1

    def test_V1D_SubdomainResource(self):
        # Test Subdomain, CombineLevels and D3VirtualField
        profile_resource = fpx.resource_modificator(name="Subdomain",
                                                    resource=self._CL_resources[0],
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

    def test_H2D_fromvarious_metaresources(self):
        fid = {'discipline': 0,
               'parameterCategory': 1,
               'parameterNumber': 1,
               'scaledValueOfFirstFixedSurface': 2, #level': 2,
               'typeOfFirstFixedSurface': 103}
        virtual3D_resource = fpx.resource_modificator(name='CombineLevels',
                                                      resource=self.resources[0],
                                                      openmode='r',
                                                      virtual=True)
        f1 = virtual3D_resource.readfield(fid)
        if f1.spectral:
            f1.sp2gp()
        f2 = self._CL_resources[0].readfield(fid)
        if f2.spectral:
            f2.sp2gp()
        f3 = self.resources[0].readfield('CLSHUMI.RELATIVE')
        if f3.spectral:
            f3.sp2gp()
        self.assertTrue(numpy.all(f1.data == f2.data))
        self.assertTrue(numpy.all(f1.data == f3.data))
        if test_plot:
            f1.plotfield()
            f2.plotfield()
            f3.plotfield()
            plt.show()

    def _test__variousways(self, geomtype):
        # Test Subdomain, CombineLevels and MultiValidities
        if geomtype == 'V1D':
            geometry = self.profile1.geometry.deepcopy()
        elif geomtype == 'V2D':
            geometry = self.section1.geometry.deepcopy()
        else:
            raise Exception()
        fields = []
        for _ in range(len(geometry.vcoordinate.levels)):
            geometry.vcoordinate.levels.pop()

        # subdo(comb(mult([file])))
        chain1 = fpx.resource_modificator(name="Subdomain",
                                          resource=self._CLMV_resource,
                                          openmode='r',
                                          geometry=geometry.deepcopy(),
                                          interpolation='nearest')
        fields.append(chain1.readfield(GRIB_handgrip))

        if test_plot:
            field = fields[0]
            field.plotanimation(repeat=True, interval=100)
            plt.show()

        # subdo(mult([comb(file)]))
        chain2 = fpx.resource_modificator(name="Subdomain",
                                          resource=self._MVCL_resource,
                                          openmode='r',
                                          geometry=geometry.deepcopy(),
                                          interpolation='nearest')
        fields.append(chain2.readfield(GRIB_handgrip))

        # comb(subdo(mult([file])))
        subdo_mult = fpx.resource_modificator(name="Subdomain",
                                              resource=self._MV_resource,
                                              openmode='r',
                                              geometry=geometry.deepcopy(),
                                              interpolation='nearest')
        chain3 = fpx.resource_modificator(name='CombineLevels',
                                          resource=subdo_mult,
                                          openmode='r')
        fields.append(chain3.readfield(GRIB_handgrip))

        # mult([subdo(comb(file))])
        subdo_comb = [fpx.resource_modificator(name="Subdomain",
                                               resource=r3,
                                               openmode='r',
                                               geometry=geometry.deepcopy(),
                                               interpolation='nearest')
                      for r3 in self._CL_resources]
        chain4 = fpx.resource_modificator(name='MultiValidities',
                                          resources=subdo_comb,
                                          openmode='r')
        fields.append(chain4.readfield(GRIB_handgrip))

        # comb(mult([subdo(file)]))
        subdo = [fpx.resource_modificator(name="Subdomain",
                                          resource=r,
                                          openmode='r',
                                          geometry=geometry.deepcopy(),
                                          interpolation='nearest')
                 for r in self.resources]
        mult_subdo = fpx.resource_modificator(name='MultiValidities',
                                              resources=subdo,
                                              openmode='r')
        chain5 = fpx.resource_modificator(name='CombineLevels',
                                          resource=mult_subdo,
                                          openmode='r')
        fields.append(chain5.readfield(GRIB_handgrip))

        # mult([comb(subdo(file))])
        comb_subdo = [fpx.resource_modificator(name='CombineLevels',
                                               resource=r,
                                               openmode='r')
                      for r in subdo]
        chain6 = fpx.resource_modificator(name='MultiValidities',
                                          resources=comb_subdo,
                                          openmode='r')
        fields.append(chain6.readfield(GRIB_handgrip))

        # comb(mult([file])).extractprofile
        field = self._CLMV_resource.readfield(GRIB_handgrip)
        if field.spectral:
            field.sp2gp()
        fields.append(field.extract_subdomain(geometry.deepcopy(),
                                              interpolation='nearest'))

        # mult(comb([file])).extractprofile
        field = self._MVCL_resource.readfield(GRIB_handgrip)
        if field.spectral:
            field.sp2gp()
        fields.append(field.extract_subdomain(geometry.deepcopy(),
                                              interpolation='nearest'))

        self.assertTrue(all([numpy.all(field.data == fields[0].data)
                             for field in fields]))
        for field in fields[1:]:
            self.assertEqual(field.geometry, fields[0].geometry)

    def test_V1D_variousways(self):
        # Test Subdomain, CombineLevels and MultiValidities for V1D
        self._test__variousways('V1D')

    def test_V2D_variousways(self):
        # Test Subdomain, CombineLevels and MultiValidities for V2D
        self._test__variousways('V2D')

    def test_3D(self):
        # Test CombineLevels
        field3D = self._CL_resources[0].readfield({'discipline': 0,
                                                   'parameterCategory': 193, #'parameterCategory': 2,
                                                   'typeOfFirstFixedSurface': 119,
                                                   'parameterNumber': 28}) #'parameterNumber': 11
        field3D.sp2gp()
        for k in range(60):
            field = self.resources[0].readfield('S{:03}VERTIC.DIVER'.format(k + 1))
            field.sp2gp()
            self.assertTrue(numpy.all(field.getdata() == field3D.getdata()[k, :, :]))

    def test_4D_2ways(self):
        # Test 4D
        fid = {'discipline': 0,
               'parameterCategory': 193, #'parameterCategory': 2,
               'typeOfFirstFixedSurface': 119,
               'parameterNumber': 28} #'parameterNumber': 11
        field_1 = self._MVCL_resource.readfield(fid)
        field_2 = self._CLMV_resource.readfield(fid)
        field_1.fid.pop('MultiValidities', None)
        self.assertEqual(field_1, field_2)

    def test_MV(self):
        # Test multivalidities
        fieldname = 'SPECSURFGEOPOTEN'
        one_field = self.resources[0].readfield(fieldname)
        if one_field.spectral:
            one_field.sp2gp()
        one_data = one_field.getdata(subzone='CI', d4=True)
        mult_field = self._MV_resource.readfield(fieldname)
        if mult_field.spectral:
            mult_field.sp2gp()
        mult_data = mult_field.getdata(subzone='CI')
        self.assertTrue(numpy.all(one_data[0, :, :, :] == mult_data[0, :, :]))

if __name__ == '__main__':
    main(verbosity=2)
