#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains classes for GRIB resource and GRIB individual message,
editions 1 and 2.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

__all__ = ['GRIB']

import datetime
import os
import numpy
import copy
import sys
import six
import tempfile
import uuid

import footprints
from footprints import proxy as fpx, FPDict, FPList
from bronx.meteo.conversion import q2R
from bronx.syntax.parsing import str2dict

from epygram import config, epygramError, util
from epygram.base import FieldSet, FieldValidity
from epygram.resources import FileResource
from epygram.util import (Angle, RecursiveObject,
                          separation_line, write_formatted_dict)
from epygram.fields import H2DField
from epygram.geometries.H2DGeometry import gauss_latitudes
from epygram.geometries.VGeometry import pressure2altitude
from epygram.geometries.SpectralGeometry import (SpectralGeometry,
                                                 gridpoint_dims_from_truncation,
                                                 nearest_greater_FFT992compliant_int)
from . import grib_utilities

import gribapi

epylog = footprints.loggers.getLogger(__name__)


# FIXME: temporary conversion of surface types
onetotwo = {1:1,  # ground or water surface
            20:20,  # isothermal level
            100:100,  # isobaric surface
            102:101,  # mean sea level
            103:102,  # height above mean sea level
            105:103,  # height above ground
            109:119,  # hybrid pressure coord level
            111:106,  # depth below land surface (!!! cm -> m !!!)
            113:107,  # isentropic (theta) level
            115:108,  # level at specified pressure difference from ground to level
            117:109,  # potential vorticity surface
            }
twotoone = {v:k for (k, v) in onetotwo.items()}


def parse_GRIBstr_todict(strfid):
    """Parse and return a dict GRIB fid from a string."""
    fid = str2dict(strfid, try_convert=int)
    return fid


class GRIBmessage(RecursiveObject, dict):
    """
    Class implementing a GRIB message as an object.
    """

    fid_keys = {1:['name',
                   'shortName',
                   'indicatorOfParameter',
                   'paramId',
                   'indicatorOfTypeOfLevel',
                   'typeOfLevel',
                   'level',
                   'topLevel',
                   'bottomLevel',
                   'editionNumber',
                   'table2Version'],
                2:['name',
                   'shortName',
                   'discipline',
                   'parameterCategory',
                   'parameterNumber',
                   'editionNumber',
                   'typeOfFirstFixedSurface',
                   'level',
                   'topLevel',
                   'typeOfSecondFixedSurface',
                   'bottomLevel',
                   'tablesVersion']}

    def __init__(self, source,
                 ordering=config.GRIB_default_ordering,
                 packing=None,
                 sample=None,
                 grib_edition=None,
                 other_GRIB_options=None):
        """
        Initialize a GRIBmessage from either sources.

        :param source: being a tuple of either form:\n
          - ('file', '*filename*' [, *offset_position*])
            *filename* being a relative or absolute path to the file it is read
            in. n = *offset_position*, if given, enables to read the n+1'th GRIB
            message in file. Defaults to 0. Negative value counts from the end.
          - ('field', :class:`epygram.fields.H2DField`)
          - ('gribid', *gribid*)
            *gribid* being an integer, refering to the *gribid* of a GRIB_API
            message in memory.
          - ('sample', '*samplename*')
            *samplename* being the name of the sample from which to be
            generated.

        In case **source** is a field, some options can be forced:

        :param ordering: flattening of 2D data
        :param packing: options of packing and compression in GRIB (dict).
        :param sample: to use a specific sample GRIB.
          Specific syntax 'file:$filename$' takes the first message
          in $filename$ as sample.
        :param grib_edition: to force a GRIB edition number (1, 2).
        :param other_GRIB_options: other options to be specified in GRIB,
          as a dict(GRIBkey=value)
        """
        super(GRIBmessage, self).__init__()
        self.built_from = source[0]
        if self.built_from == 'file':
            self._file = open(source[1], 'r')
            # if position specified, go to position of message
            if len(source) == 3:
                n = source[2]
                if n < 0:
                    N = gribapi.grib_count_in_file(self._file)
                    n = N + n
                for _ in range(n):
                    gid = gribapi.grib_new_from_file(self._file,
                                                     headers_only=True)
                    gribapi.grib_release(gid)
            # load message in memory and save gribid
            self._gid = gribapi.grib_new_from_file(self._file)
        elif self.built_from == 'field':
            self._build_msg_from_field(source[1],
                                       ordering=ordering,
                                       packing=packing,
                                       sample=sample,
                                       grib_edition=grib_edition,
                                       other_GRIB_options=other_GRIB_options)
        elif self.built_from == 'gribid':
            self._gid = source[1]
        elif self.built_from == 'sample':
            self._gid = self._clone_from_sample(source[1])
        else:
            raise NotImplementedError("not yet.")

    def __del__(self):
        try:
            gribapi.grib_release(self._gid)
        except (gribapi.GribInternalError, AttributeError):
            pass
        self.clear()

    def __getitem__(self, item):
        if item not in list(self.keys()):
            self._readattribute(item)
            if item not in list(self.keys()):
                raise KeyError(item + " not in GRIBmessage keys.")

        return super(GRIBmessage, self).__getitem__(item)

    def __setitem__(self, key, value):
        if value is not None:
            if gribapi.__version__ <= '1.10.4':
                if isinstance(value, numpy.float):
                    value = float(value)
                elif isinstance(value, numpy.int):
                    value = int(value)
            if isinstance(value, six.string_types):  # gribapi str/unicode incompatibility
                v = str(value)
            else:
                v = value
            gribapi.grib_set(self._gid, str(key), v)  # gribapi str/unicode incompatibility
        else:
            gribapi.grib_set_missing(self._gid, str(key))  # gribapi str/unicode incompatibility
        super(GRIBmessage, self).__setitem__(key, value)

    def get(self, key, default=None):
        """Same as dict.get(), but try to read attribute first."""
        try:
            value = self.__getitem__(key)
        except (KeyError, gribapi.GribInternalError):
            value = default
        return value

    def _clone_from_sample(self, sample):
        """Clone a sample GRIB message."""
        sample_gid = gribapi.grib_new_from_samples(str(sample))  # gribapi str/unicode incompatibility
        gid = gribapi.grib_clone(sample_gid)
        gribapi.grib_release(sample_gid)
        return gid

    def _clone_from_file(self, filename):
        """Clone first GRIB message from file."""
        f = open(filename, 'r')
        original_gid = gribapi.grib_new_from_file(f)
        gid = gribapi.grib_clone(original_gid)
        gribapi.grib_release(original_gid)
        f.close()
        return gid

    def _readattribute(self, attribute, array=False):
        """Actual access to the attributes."""
        if attribute != 'values':
            if not array:
                try:
                    # force to get as integer
                    # bug in GRIB_API ? 1, 103 & 105 => 'sfc'
                    if attribute in ('typeOfFirstFixedSurface',
                                     'indicatorOfTypeOfLevel',
                                     'typeOfSecondFixedSurface',
                                     # type error when setting key 'centre' if str
                                     'centre', 'originatingCentre'):
                        attr = gribapi.grib_get(self._gid, str(attribute), int)  # gribapi str/unicode incompatibility
                    else:
                        attr = gribapi.grib_get(self._gid, str(attribute))  # gribapi str/unicode incompatibility
                except gribapi.GribInternalError as e:
                    # differenciation not well done... PB in gribapi
                    if str(e) == 'Passed array is too small':
                        attr = gribapi.grib_get_double_array(self._gid, str(attribute))  # gribapi str/unicode incompatibility
                    else:
                        raise type(e)(str(e) + ' : "' + attribute + '"')
            else:
                attr = gribapi.grib_get_double_array(self._gid, str(attribute))  # gribapi str/unicode incompatibility
        else:
            attr = gribapi.grib_get_values(self._gid)
        super(GRIBmessage, self).__setitem__(attribute, attr)

    def _set_array_attribute(self, attribute, value):
        """Setter for array attributes."""
        gribapi.grib_set_array(self._gid, str(attribute), value)  # gribapi str/unicode incompatibility

    def _build_msg_from_field(self, field,
                              ordering=config.GRIB_default_ordering,
                              packing=None,
                              sample=None,
                              grib_edition=None,
                              other_GRIB_options=None):
        """
        Build the GRIB message from the field, using a template.

        :param field: a :class:`epygram.base.Field`.
        :param ordering: flattening of 2D data
        :param packing: options of packing and compression in GRIB (dict).
        :param sample: to use a specific sample GRIB.
          Specific syntax 'file:$filename$' takes the first message
          in $filename$ as sample.
        :param grib_edition: to force a GRIB edition number (1, 2).
        :param other_GRIB_options: other options to be specified in GRIB,
          as a dict(GRIBkey=value)
        """
        if not isinstance(field, H2DField):
            raise NotImplementedError("not yet.")

        # part 0 --- determine grib_edition
        if grib_edition is None:
            try:
                grib_edition = field.fid['GRIB2'].get('editionNumber',
                                                      config.GRIB_default_edition)
            except KeyError:
                try:
                    grib_edition = field.fid['GRIB1'].get('editionNumber',
                                                          config.GRIB_default_edition)
                except KeyError:
                    grib_edition = config.GRIB_default_edition

        # part 1 --- set sample and preset packing
        if packing is None:
            if field.spectral:
                packing = {'packingType':'spectral_simple',
                           'bitsPerValue':24}
            else:
                packing = config.GRIB_default_packing[grib_edition]
        # get packingType...
        required_packingType = packing.get('packingType', None)
        if required_packingType is not None:
            if 'grid' in required_packingType and field.spectral or\
               'spectral' in required_packingType and not field.spectral:
                required_packingType = None
        # ... to determine sample
        if sample is None:
            # try to get an appropriate sample, or a default one
            if required_packingType is None:
                if field.spectral:
                    sample = 'sh_ml_grib' + str(grib_edition)
                else:
                    sample = config.GRIB_default_sample[grib_edition]
            else:
                # if 'gauss' in field.geometry.name:
                #    sample = 'gg_sfc_grib' + str(grib_edition)
                # else:
                    sample = 'GRIB' + str(grib_edition) + required_packingType[4:]  # [4:] to remove leading 'grid_*'
                    if sample + '.tmpl' not in os.listdir(config.GRIB_samples_path):
                        sample = config.GRIB_default_sample[grib_edition]
        # reset packing "on the fly" if field is uniform
        if field.max() - field.min() < config.epsilon:
            packing = {'packingType':'grid_simple'}
            required_packingType = packing['packingType']
            sample = 'GRIB' + str(grib_edition) + required_packingType[4:]
        # clone from sample
        if sample.startswith('file:'):
            sample = sample[5:]
            self._gid = self._clone_from_file(sample)
        else:
            self._gid = self._clone_from_sample(sample)

        # part 2 --- parameter
        if grib_edition == 1:
            param_list = ['table2Version', 'indicatorOfParameter']
            for k in param_list:
                self[k] = field.fid['GRIB1'][k]
        else:
            if 'GRIB2' in field.fid:
                self['tablesVersion'] = field.fid['GRIB2'].get('tablesVersion', config.GRIB_default_tablesVersion)
            else:
                self['tablesVersion'] = field.fid['generic'].get('tablesVersion', config.GRIB_default_tablesVersion)
            param_list = ['discipline', 'parameterCategory', 'parameterNumber']
            for k in param_list:
                if 'GRIB2' in field.fid:
                    self[k] = field.fid['GRIB2'][k]
                else:
                    self[k] = field.fid['generic'][k]

        # part 3 --- context
        if grib_edition == 2:
            for p in config.GRIB_default_production_parameters.keys():
                self[p] = config.GRIB_default_production_parameters[p]
        try:
            process_id = int(field.processtype)
        except ValueError:
            process_id = config.GRIB_default_production_parameters['generatingProcessIdentifier']
        self['generatingProcessIdentifier'] = process_id

        # part 4 --- validity
        if grib_edition == 2:
            self['hoursAfterDataCutoff'] = None
            self['minutesAfterDataCutoff'] = None
        self['dataDate'] = int(field.validity.getbasis(fmt='IntStr')[:8])
        self['dataTime'] = int(field.validity.getbasis(fmt='IntStr')[8:12])
        term_in_seconds = field.validity.term().total_seconds()
        cumulative = field.validity.cumulativeduration()
        if cumulative:
            startStep_in_seconds = (field.validity.term() - cumulative).total_seconds()
            if grib_edition == 1:
                self['timeRangeIndicator'] = 4
            else:
                self['productDefinitionTemplateNumber'] = 8
        else:
            if grib_edition == 1:
                self['timeRangeIndicator'] = 0
            else:
                self['productDefinitionTemplateNumber'] = 0
        if term_in_seconds / 3600. < config.epsilon:
            self['stepUnits'] = 'h'  # hours
            if cumulative:
                self['startStep'] = int(startStep_in_seconds // 3600)
            self['endStep'] = int(term_in_seconds // 3600)
        else:
            self['stepUnits'] = 's'  # seconds
            if cumulative:
                self['startStep'] = int(startStep_in_seconds)
            self['endStep'] = int(term_in_seconds)

        # part 5 --- geometry
        if not field.spectral:
            # 5.1 : earth shape
            if not hasattr(field.geometry, 'geoid'):
                epylog.warning("geometry has no geoid: assume 'shapeOfTheEarth' = 6.")
                self['shapeOfTheEarth'] = 6
            else:
                if (field.geometry.geoid.get('a') == grib_utilities.pyproj_geoid_shapes[6]['a'] and
                    field.geometry.geoid.get('b') == grib_utilities.pyproj_geoid_shapes[6]['b']) or \
                   field.geometry.geoid.get('geoidradius') == grib_utilities.myproj_geoid_shapes[6]['geoidradius']:
                    self['shapeOfTheEarth'] = 6
                else:
                    found = False
                    for s, g in grib_utilities.pyproj_geoid_shapes.items():
                        if s in (0, 2, 4, 5, 6, 8, 9) and field.geometry.geoid == g:
                            self['shapeOfTheEarth'] = s
                            break
                            found = True
                    if not found:
                        radius = field.geometry.geoid.get('geoidradius')
                        if radius is None:
                            radius = field.geometry.geoid.get('a')
                            if radius is None or radius != field.geometry.geoid.get('b'):
                                radius = None
                        if radius is not None:
                            self['shapeOfTheEarth'] = 1
                            self['scaleFactorOfRadiusOfSphericalEarth'] = 1
                            self['scaledValueOfRadiusOfSphericalEarth'] = radius
                        else:
                            a = field.geometry.geoid.get('a', field.geometry.geoid.get('geoidmajor'))
                            b = field.geometry.geoid.get('b', field.geometry.geoid.get('geoidminor'))
                            if a is None or b is None:
                                raise epygramError("unable to encode geoid.")
                            else:
                                self['shapeOfTheEarth'] = 7
                                self['scaleFactorOfMajorAxisOfOblateSpheroidEarth'] = 1
                                self['scaledValueOfMajorAxisOfOblateSpheroidEarth'] = a
                                self['scaleFactorOfMinorAxisOfOblateSpheroidEarth'] = 1
                                self['scaledValueOfMinorAxisOfOblateSpheroidEarth'] = b
            # 5.2: dimensions
            if field.geometry.rectangular_grid:
                self['Ni'] = field.geometry.dimensions['X']
                self['Nj'] = field.geometry.dimensions['Y']
            else:
                pass  # done below
            # 5.3: type of geometry
            if field.geometry.name == 'regular_lonlat':
                self['gridType'] = 'regular_ll'
                self['iDirectionIncrementInDegrees'] = field.geometry.grid['X_resolution'].get('degrees')
                self['jDirectionIncrementInDegrees'] = field.geometry.grid['Y_resolution'].get('degrees')
            elif field.geometry.projected_geometry:
                self['gridType'] = field.geometry.name
                if field.geometry.name == 'lambert':
                    # lambda0 in geoid in case shapeOfTheEarth == 9
                    lambda0 = field.geometry.geoid.get('lambda0', field.geometry.projection['reference_lon'].get('degrees'))
                    self['LoVInDegrees'] = util.positive_longitude(lambda0)
                    lat_0 = field.geometry.projection.get('reference_lat', None)
                    lat_1 = field.geometry.projection.get('secant_lat1', None)
                    lat_2 = field.geometry.projection.get('secant_lat2', None)
                    if lat_0 is None:
                        lat_1 = lat_1.get('degrees')
                        lat_2 = lat_2.get('degrees')
                        lat_0 = (lat_1 + lat_2) / 2.
                    elif lat_1 is None and lat_2 is None:
                        lat_0 = lat_0.get('degrees')
                        lat_1 = lat_2 = lat_0
                    self['LaDInDegrees'] = lat_0
                    self['Latin1InDegrees'] = lat_1
                    self['Latin2InDegrees'] = lat_2
                    if abs(field.geometry.projection['rotation'].get('degrees')) > config.epsilon:
                        raise NotImplementedError("*rotation* attribute of projection != 0.")
                elif field.geometry.name == 'polar_stereographic':
                    if abs(field.geometry.projection['reference_lat'].get('degrees') - 90.) < config.epsilon:
                        self['projectionCentreFlag'] = 0
                    else:
                        self['projectionCentreFlag'] = 1
                    lat_ts = field.geometry.projection.get('secant_lat', field.geometry.projection['reference_lat']).get('degrees')
                    try:
                        self['LaDInDegrees'] = lat_ts
                    except gribapi.GribInternalError:
                        if abs(lat_ts -
                               numpy.copysign(60., field.geometry.projection['reference_lat'].get('degrees'))) > config.epsilon:
                            raise epygramError("unable to write polar stereographic geometry to GRIB1 if *secant_lat* is not +/-60 degrees.")
                    self['orientationOfTheGridInDegrees'] = field.geometry.projection['reference_lon'].get('degrees')
                    if abs(field.geometry.projection['rotation'].get('degrees')) > config.epsilon:
                        raise NotImplementedError("*rotation* attribute of projection != 0.")
                elif field.geometry.name == 'mercator':
                    lat_ts = field.geometry.projection.get('secant_lat', field.geometry.projection['reference_lat']).get('degrees')
                    try:
                        self['LaDInDegrees'] = lat_ts
                    except gribapi.GribInternalError:
                        pass  # TOBECHECKED: in GRIB1, lat_ts=0. ?
                    if abs(field.geometry.projection['rotation'].get('degrees')) > config.epsilon:
                        raise NotImplementedError("*rotation* attribute of projection != 0.")
                if field.geometry.name in ('lambert', 'polar_stereographic'):
                    self['DxInMetres'] = field.geometry.grid['X_resolution']
                    self['DyInMetres'] = field.geometry.grid['Y_resolution']
                else:
                    self['DiInMetres'] = field.geometry.grid['X_resolution']
                    self['DjInMetres'] = field.geometry.grid['Y_resolution']
            elif 'gauss' in field.geometry.name:
                if field.geometry.name == 'rotated_reduced_gauss':
                    self['gridType'] = 'reduced_stretched_rotated_gg'
                    self['longitudeOfStretchingPoleInDegrees'] = field.geometry.grid['pole_lon'].get('degrees')
                    self['latitudeOfStretchingPoleInDegrees'] = field.geometry.grid['pole_lat'].get('degrees')
                    self['stretchingFactor'] = field.geometry.grid['dilatation_coef']
                elif field.geometry.name == 'reduced_gauss':
                    if field.geometry.grid.get('dilatation_coef') != 1.:
                        self['gridType'] = 'reduced_stretched_gg'
                        self['stretchingFactor'] = field.geometry.grid['dilatation_coef']
                    else:
                        self['gridType'] = 'reduced_gg'
                self['global'] = 1
                self['latitudeOfFirstGridPointInDegrees'] = field.geometry.grid['latitudes'][0].get('degrees')
                self['longitudeOfFirstGridPointInDegrees'] = 0.
                self['latitudeOfLastGridPointInDegrees'] = field.geometry.grid['latitudes'][-1].get('degrees')
                self['longitudeOfLastGridPointInDegrees'] = 360. - 360. / field.geometry.dimensions['lon_number_by_lat'][-1]
                self['Nj'] = field.geometry.dimensions['lat_number']
                self['N'] = field.geometry.dimensions['lat_number'] // 2
                self._set_array_attribute('pl', field.geometry.dimensions['lon_number_by_lat'])
            else:
                raise NotImplementedError("not yet.")
        else:
            # spectral case
            if field.spectral_geometry.space == 'bi-fourier':
                # TODO: update when spectral LAM accepted by WMO and implemented in grib_api
                self['gridType'] = 'sh'  # in bi-fourier case, this is a bypass
                self['J'] = max(field.spectral_geometry.truncation['in_X'],
                                field.spectral_geometry.truncation['in_Y'])
                self['K'] = self['J']
                self['M'] = self['J']
            elif field.spectral_geometry.space == 'legendre':
                self['gridType'] = 'sh'
                self['sphericalHarmonics'] = 1
                self['J'] = field.spectral_geometry.truncation['max']
                self['K'] = field.spectral_geometry.truncation['max']
                self['M'] = field.spectral_geometry.truncation['max']
            else:
                raise NotImplementedError('spectral_geometry.name == ' + field.spectral_geometry.name)

        # 5.4: vertical geometry
        if grib_edition == 1:
            self['indicatorOfTypeOfLevel'] = twotoone.get(field.geometry.vcoordinate.typeoffirstfixedsurface, 255)
            if len(field.geometry.vcoordinate.levels) > 1:
                raise epygramError("field has more than one level")
            self['level'] = field.geometry.vcoordinate.levels[0]
            if self['indicatorOfTypeOfLevel'] == 109:
                self['numberOfVerticalCoordinateValues'] = len(field.geometry.vcoordinate.grid['gridlevels']) * 2
                ab = [ab[1]['Ai'] for ab in field.geometry.vcoordinate.grid['gridlevels']] + \
                     [ab[1]['Bi'] for ab in field.geometry.vcoordinate.grid['gridlevels']]
                self._set_array_attribute('pv', ab)
            else:
                pass
                # !!! this should be done but it changes gridType !!!
                # self['numberOfVerticalCoordinateValues'] = 0
        elif grib_edition == 2:
            if len(field.geometry.vcoordinate.levels) > 1:
                raise epygramError("field has more than one level")
            self['scaleFactorOfFirstFixedSurface'] = 0
            self['scaleFactorOfSecondFixedSurface'] = 0
            self['level'] = field.geometry.vcoordinate.levels[0]
            self['typeOfFirstFixedSurface'] = field.geometry.vcoordinate.typeoffirstfixedsurface
            self['level'] = field.geometry.vcoordinate.levels[0]
            if self['typeOfFirstFixedSurface'] == 119:
                self['numberOfVerticalCoordinateValues'] = len(field.geometry.vcoordinate.grid['gridlevels'])
                ab = [ab[1]['Ai'] for ab in field.geometry.vcoordinate.grid['gridlevels']] + \
                     [ab[1]['Bi'] for ab in field.geometry.vcoordinate.grid['gridlevels']]
                self._set_array_attribute('pv', ab)
            if hasattr(field.geometry.vcoordinate, 'typeofsecondfixedsurface'):
                self['typeOfSecondFixedSurface'] = field.geometry.vcoordinate.typeofsecondfixedsurface
            else:
                self['typeOfSecondFixedSurface'] = 255
            if hasattr(field.geometry.vcoordinate, 'toplevel'):
                self['topLevel'] = field.geometry.vcoordinate.toplevel
            else:
                if 'GRIB2' in field.fid:
                    self['topLevel'] = field.fid['GRIB2'].get('topLevel', field.geometry.vcoordinate.levels[0])
                else:
                    self['topLevel'] = field.fid['generic'].get('topLevel', field.geometry.vcoordinate.levels[0])
            if hasattr(field.geometry.vcoordinate, 'bottomlevel'):
                self['bottomLevel'] = field.geometry.vcoordinate.bottomlevel
            else:
                if 'GRIB2' in field.fid:
                    self['bottomLevel'] = field.fid['GRIB2'].get('bottomLevel', field.geometry.vcoordinate.levels[0])
                else:
                    self['bottomLevel'] = field.fid['generic'].get('bottomLevel', field.geometry.vcoordinate.levels[0])

        # part 6 --- ordering
        if not field.spectral:
            self.update(ordering)

        # part 7 --- other options
        if other_GRIB_options is not None:
            try:
                for k, v in other_GRIB_options.items():
                    self[k] = v
            except gribapi.GribInternalError as e:
                epylog.warning('set ' + k + ' failed: ' + str(e))

        # part 8 --- set first gridpoint according to ordering
        if not field.spectral:
            if field.geometry.rectangular_grid:
                corners = field.geometry.gimme_corners_ll()
                if self['iScansNegatively'] == 0 and \
                   self['jScansPositively'] == 0:
                    self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                    self['latitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ul'][1])
                    if self['gridType'] == 'regular_ll':
                        self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                        self['latitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['lr'][1])
                elif self['iScansNegatively'] == 0 and \
                     self['jScansPositively'] == 1:
                    self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                    self['latitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ll'][1])
                    if self['gridType'] == 'regular_ll':
                        self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                        self['latitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ur'][1])
                elif self['iScansNegatively'] == 1 and \
                     self['jScansPositively'] == 0:
                    self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                    self['latitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ur'][1])
                    if self['gridType'] == 'regular_ll':
                        self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                        self['latitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ll'][1])
                elif self['iScansNegatively'] == 1 and \
                     self['jScansPositively'] == 1:
                    self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                    self['latitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['lr'][1])
                    if self['gridType'] == 'regular_ll':
                        self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                        self['latitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ul'][1])
                else:
                    raise NotImplementedError('this ordering: not yet.')
            else:
                if not(self['iScansNegatively'] == 0 and
                       self['jScansPositively'] == 0):
                    raise NotImplementedError('this ordering: not yet.')

        # part 9 --- values
        values = field.getdata(d4=True).copy()
        if isinstance(values, numpy.ma.masked_array) and values.mask.any():
            if self['editionNumber'] == 2:
                self['bitMapIndicator'] = 0
                self['bitmapPresent'] = 1
                self['missingValue'] = values.fill_value
                values = field.geometry.fill_maskedvalues(values)
            else:
                # bitmap in GRIB1 ?
                raise NotImplementedError("didn't succeed to make this work")
        values = values.squeeze()
        if not field.spectral:
            # is it necessary to pre-write values ? (packingType != from sample)
            # Yes it is (don't really know why...)
            """if required_packingType is None:  # unable to guess
                pre_write = True
            else:
                try:
                    pre_write = required_packingType != self['packingType']
                except gribapi.GribInternalError:
                    pre_write = True"""
            pre_write = True
            if pre_write:
                self.set_values(values)
        self.set_packing(packing)
        self.set_values(values)

    def set_packing(self, packing):
        """
        Specific method to set **packing** because the order of the elements is
        important.

        :param dict packing: GRIB keys from packing.
        """
        packing = copy.copy(packing)
        if packing.get('bitsPerValue') is not None and \
           packing.get('bitsPerValue') > config.GRIB_max_bitspervalue:
            # problem with bitsPerValue = 30 at least
            epylog.warning(('GRIB encoding higher than {} ' +
                            '(bits per value): {} : may be untrustful.').
                           format(config.GRIB_max_bitspervalue,
                                  packing.get('bitsPerValue')))
            if config.GRIB_force_bitspervalue:
                packing['bitsPerValue'] = config.GRIB_max_bitspervalue
        order = ['packingType', 'complexPacking', 'boustrophedonicOrdering',
                 'bitsPerValue']
        for k in order:
            if k in packing:
                try:
                    self[k] = packing.pop(k)
                except gribapi.GribInternalError:
                    if config.GRIB_packing_fatal:
                        raise
        for k, v in packing.items():  # remaining items
            self[k] = v

    @property
    def grib_edition(self):
        return self['editionNumber']

    def set_values(self, values):
        """
        Wrapper to set **values** as a 2D array if gridpoint or 1D if spectral.
        """
        if len(values.shape) == 1:
            gribapi.grib_set_values(self._gid, values)
        elif len(values.shape) == 2:
            self.set_2Dvalues(values)

    def set_2Dvalues(self, values):
        """
        Wrapper to set **values** as a 2D array, in coherence with ordering
        parameters already set beforehand.
        """
        if self['iScansNegatively'] == 0 and \
           self['jScansPositively'] == 0 and \
           self['jPointsAreConsecutive'] == 0:
            if isinstance(values, numpy.ma.MaskedArray):
                data1d = values[:, :].compressed()
            else:
                data1d = values[::-1, :].flatten(order='C')
        elif self['iScansNegatively'] == 0 and \
             self['jScansPositively'] == 1 and \
             self['jPointsAreConsecutive'] == 0:
            if isinstance(values, numpy.ma.MaskedArray):
                data1d = values[::-1, :].compressed()
            else:
                data1d = values[:, :].flatten(order='C')
        elif self['iScansNegatively'] == 1 and \
             self['jScansPositively'] == 0 and \
             self['jPointsAreConsecutive'] == 0:
            if isinstance(values, numpy.ma.MaskedArray):
                data1d = values[:, ::-1].compressed()
            else:
                data1d = values[::-1, ::-1].flatten(order='C')
        elif self['iScansNegatively'] == 1 and \
             self['jScansPositively'] == 1 and \
             self['jPointsAreConsecutive'] == 0:
            if isinstance(values, numpy.ma.MaskedArray):
                data1d = values[::-1, ::-1].compressed()
            else:
                data1d = values[:, ::-1].flatten(order='C')
        else:
            raise NotImplementedError('this ordering: not yet.')
        gribapi.grib_set_values(self._gid, data1d)

    def _read_spectralgeometry(self):
        """
        Returns a SpectralGeometry object containing
        the spectral geometry information of the GRIB message.
        """
        if self['gridType'] == 'sh':
            if not self['J'] == self['K'] == self['M']:
                raise NotImplementedError("case: not J==K==M")
            spgeom = SpectralGeometry(space='legendre',
                                      truncation={'max':self['J']})
        else:
            raise NotImplementedError('gridType==' + self['gridType'])
        return spgeom

    def _read_geometry(self):
        """
        Returns a geometry object containing
        the geometry information of the GRIB message.
        """
        geoid = config.default_geoid
        try:
            geoid = grib_utilities.pyproj_geoid_shapes[self['shapeOfTheEarth']]
        except KeyError:
            if self['shapeOfTheEarth'] == 1:
                radius = (float(self['scaledValueOfRadiusOfSphericalEarth']) /
                          10 ** self['scaleFactorOfRadiusOfSphericalEarth'])
                geoid = {'a':radius, 'b':radius}
            elif self['shapeOfTheEarth'] in (3, 7):
                a = (float(self['scaledValueOfMajorAxisOfOblateSpheroidEarth']) /
                     10 ** self['scaleFactorOfMajorAxisOfOblateSpheroidEarth'])
                b = (float(self['scaledValueOfMinorAxisOfOblateSpheroidEarth']) /
                     10 ** self['scaleFactorOfMinorAxisOfOblateSpheroidEarth'])
                if self['shapeOfTheEarth'] == 3:
                    a *= 1000
                    b *= 1000
                geoid = {'a':a, 'b':b}
            else:
                raise
        except gribapi.GribInternalError:
            pass

        kwargs_vcoord = {'structure': 'V'}
        kwargs_vcoord['position_on_grid'] = 'mass'
        if self.grib_edition == 1:
            kwargs_vcoord['typeoffirstfixedsurface'] = onetotwo.get(self['indicatorOfTypeOfLevel'], 255)
            kwargs_vcoord['levels'] = [self['level']]
            if self['indicatorOfTypeOfLevel'] in (112,):
                kwargs_vcoord['toplevel'] = self['topLevel']
                kwargs_vcoord['bottomlevel'] = self['bottomLevel']
        elif self.grib_edition == 2:
            kwargs_vcoord['typeoffirstfixedsurface'] = self['typeOfFirstFixedSurface']
            kwargs_vcoord['levels'] = [self['level']]
            if self['typeOfSecondFixedSurface'] != 255:
                kwargs_vcoord['typeofsecondfixedsurface'] = self['typeOfSecondFixedSurface']
                kwargs_vcoord['toplevel'] = self['topLevel']
                kwargs_vcoord['bottomlevel'] = self['bottomLevel']
        if kwargs_vcoord['typeoffirstfixedsurface'] == 119:
            try:
                self._readattribute('pv', array=True)
                A_and_B = self['pv']
            except gribapi.GribInternalError:
                epylog.warning('Error while reading A/B vertical levels coefficients ! Ignore.')
                A_and_B = []
            Ai = A_and_B[:len(A_and_B) // 2]
            Bi = A_and_B[len(A_and_B) // 2:]
            kwargs_vcoord['grid'] = {'gridlevels': tuple([(i + 1, FPDict({'Ai':Ai[i], 'Bi':Bi[i]})) for
                                                          i in range(len(Ai))]),
                                     'ABgrid_position':'flux'}
        vcoordinate = fpx.geometry(**kwargs_vcoord)

        if self['gridType'] != 'sh' and 'gg' not in self['gridType']:
            dimensions = {'X':self.get('Nx', self['Ni']),
                          'Y':self.get('Ny', self['Nj'])}
            if self['iScansNegatively'] == 0 and self['jScansPositively'] == 0:
                input_position = (0, dimensions['Y'] - 1)
            elif self['iScansNegatively'] == 0 and self['jScansPositively'] == 1:
                input_position = (0, 0)
            elif self['iScansNegatively'] == 1 and self['jScansPositively'] == 0:
                input_position = (dimensions['X'] - 1, dimensions['Y'] - 1)
            elif self['iScansNegatively'] == 1 and self['jScansPositively'] == 1:
                input_position = (dimensions['X'] - 1, 0)

        if self['gridType'] == 'regular_ll':
            geometryname = 'regular_lonlat'
            grid = {'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'X_resolution':Angle(self['iDirectionIncrementInDegrees'],
                                         'degrees'),
                    'Y_resolution':Angle(self['jDirectionIncrementInDegrees'],
                                         'degrees')
                    }
            projection = None
        elif self['gridType'] in ('polar_stereographic',):
            geometryname = self['gridType']
            if self['projectionCentreFlag'] == 0:
                lat_0 = 90.
            else:
                lat_0 = -90.
            # In GRIB1, LaD is not present but resolution is supposed to
            # be given at 60deg
            lat_ts = self.get('LaDInDegrees', numpy.copysign(60, lat_0))
            # !!! dirty bypass of erroneous GRIBs from OSI-SAF !..
            if self['centreDescription'] == 'Oslo':
                epylog.warning(' '.join(['for centre:',
                                         self['centreDescription'],
                                         "(OSI-SAF)"
                                         "and polar stereographic" +
                                         "projection, lat_ts is taken in" +
                                         "3rd position in attribute 'pv[]'" +
                                         "of GRIB message !"]))
                (a, b, lat_ts) = self['pv']
                geoid = {'a':a, 'b':b}
            projection = {'reference_lon':Angle(float(self['orientationOfTheGridInDegrees']), 'degrees'),
                          'reference_lat':Angle(lat_0, 'degrees'),
                          'secant_lat':Angle(lat_ts, 'degrees'),
                          'rotation':Angle(0., 'degrees')
                          }
            # TOBEDELETED: in GRIB stereo-polar projection, the resolution is supposed
            # to be given at 60deg (GRIB1) and 'LaD'deg (GRIB2)
            # m = (1. + numpy.copysign(1., lat_0) * numpy.sin(numpy.radians(lat_0))) / \
            #     (1. + numpy.copysign(1., lat_0) * numpy.sin(numpy.radians(lat_ts)))
            m = 1
            grid = {'X_resolution':float(self['DxInMetres']) * m,
                    'Y_resolution':float(self['DyInMetres']) * m,
                    'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'LAMzone':None}
        elif self['gridType'] in ('lambert',):
            geometryname = self['gridType']
            lon_0 = geoid.get('lambda0', self['LoVInDegrees'])  # lambda0 in geoid in case shapeOfTheEarth == 9
            lat_0 = self['LaDInDegrees']
            lat_1 = self['Latin1InDegrees']
            lat_2 = self['Latin2InDegrees']
            projection = {'reference_lon':Angle(lon_0, 'degrees'),
                          'rotation':Angle(0., 'degrees')}
            if abs(lat_1 - lat_2) < config.epsilon:
                projection['reference_lat'] = Angle(lat_0, 'degrees')
            else:
                projection['secant_lat1'] = Angle(lat_1, 'degrees')
                projection['secant_lat2'] = Angle(lat_2, 'degrees')
            grid = {'X_resolution':float(self['DxInMetres']),
                    'Y_resolution':float(self['DyInMetres']),
                    'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'LAMzone':None}
        elif self['gridType'] in ('mercator',):
            geometryname = self['gridType']
            lat_ts = self['LaDInDegrees']
            projection = {'reference_lon':Angle(0., 'degrees'),
                          'reference_lat':Angle(0., 'degrees'),
                          'rotation':Angle(0., 'degrees')}
            if abs(lat_ts) > config.epsilon:
                projection['secant_lat'] = lat_ts
            grid = {'X_resolution':float(self['DiInMetres']),
                    'Y_resolution':float(self['DjInMetres']),
                    'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'LAMzone':None}
        elif 'gg' in self['gridType']:
            projection = None
            latitudes = gauss_latitudes(self['Nj'])
            grid = {'latitudes':FPList([Angle(l, 'degrees') for l in latitudes])}
            if self['gridType'] == 'reduced_gg':
                geometryname = 'reduced_gauss'
                grid['dilatation_coef'] = 1.
            elif self['gridType'] == 'regular_gg':
                geometryname = 'regular_gauss'
                grid['dilatation_coef'] = 1.
            else:
                if 'stretched' in self['gridType']:
                    grid['dilatation_coef'] = self['stretchingFactor']
                else:
                    grid['dilatation_coef'] = 1.
                if 'rotated' in self['gridType']:
                    grid['pole_lon'] = Angle(self['longitudeOfStretchingPoleInDegrees'], 'degrees')
                    grid['pole_lat'] = Angle(self['latitudeOfStretchingPoleInDegrees'], 'degrees')
                    geometryname = 'rotated_reduced_gauss'

            if 'reduced' in self['gridType']:
                self._readattribute('pl', array=True)  # pre-load with array=True to bypass gribapi error
                lon_number_by_lat = self['pl']
            elif 'regular' in self['gridType']:
                lon_number_by_lat = [self['Ni'] for _ in range(self['Nj'])]
            else:
                raise NotImplementedError('gauss grid of that type' + self['gridType'])
            dimensions = {'max_lon_number':int(max(lon_number_by_lat)),
                          'lat_number':len(latitudes),
                          'lon_number_by_lat':FPList([int(n) for n in
                                                      lon_number_by_lat])
                          }
        elif self['gridType'] == 'sh':
            # spherical harmonics: => forced to a linear gauss grid
            projection = None
            spgeom = self._read_spectralgeometry()
            gpdims = gridpoint_dims_from_truncation(spgeom.truncation,
                                                    grid='linear')
            latitudes = gauss_latitudes(gpdims['lat_number'])
            grid = {'latitudes':FPList([Angle(l, 'degrees')
                                        for l in latitudes]),
                    'dilatation_coef':1.}
            dimensions = gpdims
            geometryname = 'reduced_gauss'
            # try to have roughly the same zonal resolution as on equator
            lon_number_by_lat = 2 * gpdims['lat_number'] * numpy.cos(numpy.radians(latitudes))
            lon_number_by_lat = [min(nearest_greater_FFT992compliant_int(n),
                                     dimensions['max_lon_number'])
                                 for n in lon_number_by_lat]
            dimensions['lon_number_by_lat'] = FPList(lon_number_by_lat)
        else:
            raise NotImplementedError("not yet !")
        # Make geometry object
        kwargs_geom = dict(structure='H2D',
                           name=geometryname,
                           grid=grid,
                           dimensions=dimensions,
                           vcoordinate=vcoordinate,
                           projection=projection,
                           position_on_horizontal_grid='center',
                           geoid=geoid)
        geometry = fpx.geometry(**kwargs_geom)

        return geometry

    def _read_validity(self):
        """
        Returns a :class:`epygram.base.FieldValidity` object containing the
        validity of the GRIB message.
        """
        accepted_tRI = (0, 1, 2, 4, 10, 113, 123)
        year = int(str(self['dataDate'])[0:4])
        month = int(str(self['dataDate'])[4:6])
        day = int(str(self['dataDate'])[6:8])
        dataTime = '{:0>{width}}'.format(self['dataTime'], width=4)
        hour = int(dataTime[0:2])
        minutes = int(dataTime[2:4])
        basis = datetime.datetime(year, month, day, hour, minutes)
        cum = None
        if self['stepUnits'] == 1:
            timeunitfactor = 3600
        if self['timeRangeIndicator'] in accepted_tRI:
            term = self['endStep']
            term = datetime.timedelta(seconds=term * timeunitfactor)
            if self['timeRangeIndicator'] in (2, 4) or self['productDefinitionTemplateNumber'] == 8:
                cum = self['endStep'] - self['startStep']
                cum = datetime.timedelta(seconds=cum * timeunitfactor)
            elif self['timeRangeIndicator'] in (113, 123,):
                epylog.warning('not able to interpret timeRangeIndicator={}'.format(self['timeRangeIndicator']))
        else:
            raise NotImplementedError("'timeRangeIndicator' not in {}.".format(accepted_tRI))
        validity = FieldValidity(basis=basis, term=term, cumulativeduration=cum)

        return validity

    def update(self, E, **F):
        """
        M.update([E, ]**F) -> None. Update D from dict/iterable E and F.
        If E present and has a .keys() method, does: for k in E: D[k] = E[k]
        If E present and lacks .keys() method, does: for (k, v) in E: D[k] = v
        In either case, this is followed by: for k in F: D[k] = F[k]
        """
        if isinstance(E, dict):
            items = list(E.items())
        elif '__iter__' in dir(E):
            items = E
        for (k, v) in items:
            self[k] = v
        for k in F:
            self[k] = F[k]

    def genfid(self):
        """Generates and returns a GRIB-type epygram fid from the message."""
        fid_keys = copy.copy(self.fid_keys[self.grib_edition])
        if self['topLevel'] == self['bottomLevel']:
            fid_keys.pop(fid_keys.index('topLevel'))
            fid_keys.pop(fid_keys.index('bottomLevel'))
        fid = {'GRIB' + str(self.grib_edition):FPDict({k:self[k] for k in fid_keys})}
        if self.grib_edition == 1:
            # Here we should complete the generic part of fid
            pass
        elif self.grib_edition == 2:
            fid['generic'] = copy.copy(fid['GRIB2'])
        return fid

    def readkeys(self, namespace=None):
        """
        Reads and returns the available keys of the message.

        :param namespace: the namespace of keys to be read, among:\n
          - **None**: to get all keys present in message,
          - 'ls': to get the same default keys as the grib_ls,
          - 'mars': to get the keys used by MARS.
        """
        if isinstance(namespace, six.string_types):
            namespace = str(namespace)  # gribapi str/unicode incompatibility
        key_iter = gribapi.grib_keys_iterator_new(self._gid,
                                                  namespace=namespace)
        namespace = []
        while gribapi.grib_keys_iterator_next(key_iter):
            namespace.append(gribapi.grib_keys_iterator_get_name(key_iter))
        gribapi.grib_keys_iterator_delete(key_iter)

        return namespace

    def readmessage(self, namespace=None):
        """
        Reads the meta-data of the message.

        :param namespace: the namespace of keys to be read, among:\n
          - **None**: to get all keys present in message,
          - ['myKey1', 'myKey2', ...] for any custom namespace,
          - 'ls': to get the same default keys as the grib_ls,
          - 'mars': to get the keys used by MARS.
        """
        self.clear()
        if namespace in (None, 'ls', 'mars'):
            namespace = self.readkeys(namespace=namespace)
        for k in namespace:
            self._readattribute(k)

    def asfield(self,
                getdata=True,
                footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                get_info_as_json=None):
        """
        Returns an :class:`epygram.base.Field` made out from the GRIB message.

        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        :param footprints_proxy_as_builder: if **True**, uses footprints.proxy
          to build fields.
        :param get_info_as_json: if not **None**, writes the keys given in
          *get_info_as_json* as json in field.comment.
        """
        if footprints_proxy_as_builder:
            builder = fpx.field
        else:
            builder = H2DField

        field_kwargs = {}
        field_kwargs['fid'] = self.genfid()
        try:
            field_kwargs['validity'] = self._read_validity()
        except (epygramError, NotImplementedError):
            if config.GRIB_ignore_validity_decoding_errors:
                field_kwargs['validity'] = FieldValidity()
            else:
                raise
        field_kwargs['spectral_geometry'] = None
        field_kwargs['processtype'] = self['generatingProcessIdentifier']
        geometry = self._read_geometry()
        if self['gridType'] == 'sh':
            field_kwargs['spectral_geometry'] = self._read_spectralgeometry()
        field_kwargs['geometry'] = geometry
        field_kwargs['structure'] = geometry.structure
        try:
            field_kwargs['units'] = self['units']
        except gribapi.GribInternalError:
            try:
                field_kwargs['units'] = self['parameterUnits']
            except gribapi.GribInternalError:
                pass
        if get_info_as_json is not None:
            import json
            comment = {}
            for k in get_info_as_json:
                try:
                    comment[k] = self[k]
                except gribapi.GribInternalError:
                    pass
            field_kwargs['comment'] = json.dumps(comment)
        field = builder(**field_kwargs)
        if getdata:
            data1d = self['values']
            if 'gauss' in geometry.name:
                if self['gridType'] == 'sh':
                    data2d = data1d
                else:
                    if self['iScansNegatively'] == 0 and \
                       self['jScansPositively'] == 0 and \
                       self['jPointsAreConsecutive'] == 0:
                        data2d = geometry.reshape_data(data1d[:])
                    else:
                        raise NotImplementedError("not yet !")
            else:
                if self['iScansNegatively'] == 0 and \
                   self['jScansPositively'] == 0 and \
                   self['jPointsAreConsecutive'] == 0:
                    data2d = data1d.reshape((self['Nj'], self['Ni']),
                                            order='C')[::-1, :]
                elif self['iScansNegatively'] == 0 and \
                     self['jScansPositively'] == 1 and \
                     self['jPointsAreConsecutive'] == 0:
                    data2d = data1d.reshape((self['Nj'], self['Ni']),
                                            order='C')[:, :]
                elif self['iScansNegatively'] == 1 and \
                     self['jScansPositively'] == 0 and \
                     self['jPointsAreConsecutive'] == 0:
                    data2d = data1d.reshape((self['Nj'], self['Ni']),
                                            order='C')[::-1, ::-1]
                elif self['iScansNegatively'] == 1 and \
                     self['jScansPositively'] == 1 and \
                     self['jPointsAreConsecutive'] == 0:
                    data2d = data1d.reshape((self['Nj'], self['Ni']),
                                            order='C')[:, ::-1]
                else:
                    raise NotImplementedError("not yet !")

            if ((self['editionNumber'] == 2 and
                 self['bitMapIndicator'] == 0) or
                (self['editionNumber'] == 1 and
                 self['bitmapPresent'] == 1)):
                data2d = numpy.ma.masked_equal(data2d, self['missingValue'])
            field.setdata(data2d)

        return field

    def write_to_file(self, ofile):
        """
        *ofile* being an open file-like object designing the physical GRIB
        file to be written to.
        """
        gribapi.grib_write(self._gid, ofile)


class GRIB(FileResource):
    """Class implementing all specificities for GRIB resource format."""

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['GRIB']),
                default='GRIB')
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        super(GRIB, self).__init__(*args, **kwargs)
        if self.openmode in ('r', 'a'):
            _file = open(self.container.abspath, 'r')
            isgrib = _file.readline()
            _file.close()
            if isgrib[0:4] != 'GRIB':
                raise IOError("this resource is not a GRIB one.")
        if not self.fmtdelayedopen:
            self.open()

    def open(self, openmode=None):
        """
        Opens a GRIB and initializes some attributes.

        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(GRIB, self).open(openmode=openmode)
        self._file = open(self.container.abspath, self.openmode)
        self.isopen = True

    def close(self):
        """
        Closes a GRIB.
        """
        if hasattr(self, '_file'):
            self._file.close()
        self.isopen = False

    @property
    @FileResource._openbeforedelayed
    def messages_number(self):
        """Counts the number of messages in file."""
        return gribapi.grib_count_in_file(self._file)

    def listfields(self, onlykey=None, select=None, complete=False):
        """
        Returns a list containing the GRIB identifiers of all the fields of the
        resource.

        :param onlykey: can be specified as a string or a tuple of strings,
          so that only specified keys of the fid will returned.
        :param select: can be specified as a dict(key=value) to restrain
          the list of fields to those that match the key:value pairs.
        :param complete: list fields with their natural fid + generic one.
        """
        if select is not None:
            additional_keys = [k for k in select.keys() if k not in
                               GRIBmessage.fid_keys[1] + GRIBmessage.fid_keys[2]]
        else:
            additional_keys = []
        fidlist = super(GRIB, self).listfields(additional_keys=additional_keys)
        if select is not None:
            fidlist = [f for f in fidlist if all([f[k] == select[k] for k in select.keys()])]
        if onlykey is not None:
            if isinstance(onlykey, six.string_types):
                fidlist = [f[onlykey] for f in fidlist]
            elif isinstance(onlykey, tuple):
                fidlist = [{k:f[k] for k in onlykey} for f in fidlist]
        if complete:
            fidlist = [{'GRIB' + str(f['editionNumber']):f} for f in fidlist]
            for f in fidlist:
                if 'GRIB2' in f:
                    f['generic'] = f['GRIB2']

        return fidlist

    def _listfields(self, additional_keys=[]):
        """Returns a list of GRIB-type fid of the fields inside the resource."""
        fidlist = []
        _file = open(self.container.abspath, 'r')
        while True:
            fid = {}
            gid = gribapi.grib_new_from_file(_file, headers_only=True)
            if gid is None:
                break
            n = gribapi.grib_get(gid, b'editionNumber')  # gribapi str/unicode incompatibility
            for k in GRIBmessage.fid_keys[n] + additional_keys:
                # bug in GRIB_API ? 1, 103 & 105 => 'sfc'
                if k in ('typeOfFirstFixedSurface',
                         'indicatorOfTypeOfLevel',
                         'typeOfSecondFixedSurface',
                         # type error when setting key 'centre' if str
                         'centre', 'originatingCentre'):
                    fid[k] = gribapi.grib_get(gid, str(k), int)  # gribapi str/unicode incompatibility
                elif k in ('topLevel', 'bottomLevel'):
                    if gribapi.grib_get(gid, b'topLevel') != gribapi.grib_get(gid, b'bottomLevel'):  # gribapi str/unicode incompatibility
                        fid[k] = gribapi.grib_get(gid, str(k))  # gribapi str/unicode incompatibility
                else:
                    fid[k] = gribapi.grib_get(gid, str(k))  # gribapi str/unicode incompatibility
            gribapi.grib_release(gid)
            fidlist.append(fid)
        _file.close()

        return fidlist

    def split_UV(self, fieldseed):
        """
        Return two lists of fids corresponding respectively to U and V
        components of wind, given a *fieldseed*.
        Syntax example: 'shortName':'u+v', or 'indicatorOfParameter':'33+34'
        """
        if isinstance(fieldseed, six.string_types):
            fieldseed = parse_GRIBstr_todict(fieldseed)

        seeds = [fieldseed.copy(), fieldseed.copy()]
        for i in (0, 1):
            for k, v in seeds[i].items():
                if isinstance(v, six.string_types) and '+' in v:
                    v = v.split('+')[i]
                    try:
                        seeds[i][k] = int(v)
                    except ValueError:
                        seeds[i][k] = v
        Ufid = self.find_fields_in_resource(seeds[0])
        Vfid = self.find_fields_in_resource(seeds[1])

        return (sorted(Ufid), sorted(Vfid))

    def get_message_at_position(self, position):
        """
        Returns the message at position *position*, from 0 (first message)
        to messages_number-1 (last message).
        Negative position starts from the end: -1 is last message.

        Should not be used sequentially, probably very inefficient.
        """
        return GRIBmessage(('file', self.container.abspath, position))

    def iter_messages(self, headers_only=True):
        """
        Iterates sequentially on messages, returning messages.

        :param headers_only: if False, read data from messages.
        """
        try:
            ok = not self._sequential_file.closed
            if not ok:
                self._sequential_file.open()
        except AttributeError:
            self._sequential_file = open(self.container.abspath, 'r')

        gid = gribapi.grib_new_from_file(self._sequential_file,
                                         headers_only=headers_only)
        if gid is not None:
            msg = GRIBmessage(('gribid', gid))
        else:
            msg = None
            self._sequential_file.close()
        return msg

    def iter_fields(self, getdata=True, **kwargs):
        """
        Iterates sequentially on messages, returning fields.

        :param getdata: if False, do not read data from the messages.
        """
        fld = self.iter_messages(headers_only=False)
        if fld is not None:
            fld = fld.asfield(getdata=getdata, **kwargs)
        return fld

    def sortfields(self, sortingkey, onlykey=None):
        """
        Returns a sorted list of fields with regards to the given *sortingkey*
        of their fid, as a dict of lists.

        :param sortingkey: key on which to sort out fields
        :param onlykey: specified as a string or a tuple of strings,
          so that only specified keys of the fid will returned.
        """
        sortedfields = {}
        listoffields = self.listfields()
        onlykeylistoffields = self.listfields(onlykey=onlykey)
        for f in range(len(listoffields)):
            try:
                category = listoffields[f][sortingkey]
            except KeyError:
                category = 'None'
            field = onlykeylistoffields[f]
            if category in sortedfields:
                sortedfields[category].append(field)
            else:
                sortedfields[category] = [field]

        return sortedfields

    def find_fields_in_resource(self, seed=None, generic=False, **_):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        :param seed: might be:\n
          - a 'handgrip', i.e. a dict where you can store all requested GRIB
            keys, e.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl',
                        'level':850},
          - a list of handgrips
          - a string like "{'shortName':'t', 'level':850}", that would be
            converted to a dict handgrip
          - *None*. If *None* (default), returns the list of all fields in
            resource.
        :param generic: if True, returns complete fid's,
          union of {'FORMATname':fieldname} and the according generic fid of
          the fields.
        """
        if seed is None or isinstance(seed, dict):
            fieldslist = self.listfields(select=seed)
        elif isinstance(seed, list):
            fieldslist = []
            for s in seed:
                if isinstance(s, six.string_types):
                    s = parse_GRIBstr_todict(s)
                fieldslist.extend(self.listfields(select=s))
        elif isinstance(seed, six.string_types):
            fieldslist = self.listfields(select=parse_GRIBstr_todict(seed))
        else:
            raise epygramError("unknown type for seed: " + str(type(seed)))
        if fieldslist == []:
            raise epygramError("no field matching '" + str(seed) +
                               "' was found in resource " +
                               self.container.abspath)
        if generic:
            fieldslist2 = [copy.copy(f) for f in fieldslist]
            for f in fieldslist2:
                if f['editionNumber'] == 1:
                    f['typeOfFirstFixedSurface'] = onetotwo.get(f['indicatorOfTypeOfLevel'], 255)
            fieldslist = [(fieldslist[i], fieldslist2[i]) for i in range(len(fieldslist))]

        return fieldslist

    def readfield(self, handgrip,
                  getdata=True,
                  footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                  get_info_as_json=None):
        """
        Finds in GRIB the message that correspond to the *handgrip*,
        and returns it as a :class:`epygram.base.Field`.
        If several messages meet the requirements, raises error (use
        readfields() method instead).

        :param dict handgrip: a dict where you can store all requested GRIB
          keys for discrimination...
          E.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl', 'level':850}
          will return the Temperature at 850hPa field.
        :param getdata: if False, the data is not read, the field consist
          in the meta-data only.
        :param footprints_proxy_as_builder: if True, uses footprints.proxy
          to build fields. True decreases performance.
        :param get_info_as_json: if not None, writes the keys given in
          *get_info_as_json* as json in field.comment.
        """
        if isinstance(handgrip, six.string_types):
            handgrip = parse_GRIBstr_todict(handgrip)

        matchingfields = self.readfields(handgrip,
                                         getdata=getdata,
                                         footprints_proxy_as_builder=footprints_proxy_as_builder,
                                         get_info_as_json=get_info_as_json)
        # filter out those unfiltered by the below GRIBAPI bug
        filtered_matchingfields = FieldSet()
        for field in matchingfields:
            fid = field.fid.get('GRIB1', field.fid.get('GRIB2'))
            if all([fid.get(k, None) in (v, None) for k, v in handgrip.items()]):
                filtered_matchingfields.append(field)
        if len(filtered_matchingfields) > 1:
            raise epygramError("several fields found for that *handgrip*;" +
                               " please refine.")
        elif len(filtered_matchingfields) == 0:
            raise epygramError("inconsistency in *handgrip*; check again" +
                               " values and types of values")

        return filtered_matchingfields[0]

    def readfields(self, handgrip,
                   getdata=True,
                   footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                   get_info_as_json=None):
        """
        Finds in GRIB the message(s) that correspond to the *handgrip*,
        and returns it as a :class:`epygram.base.FieldSet` of
        :class:`epygram.base.Field`.

        :param dict handgrip: a dict where you can store all requested GRIB
          keys for discrimination...
          E.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl'}
          will return all the Temperature fields on Pressure levels.
        :param getdata: if False, the data is not read, the field consist
          in the meta-data only.
        :param footprints_proxy_as_builder: if True, uses footprints.proxy
          to build fields. True decreases performance.
        :param get_info_as_json: if not None, writes the keys given in
          *get_info_as_json* as json in field.comment.
        *handgrip* is a dict where you can store all requested GRIB keys...
        """
        if config.GRIB_safe_indexes:  # FIXME: well not me, gribapi: grib index workaround
            # find an available AND unique filename
            self._index_alias = str(tempfile.mkstemp(dir=config.GRIB_safe_indexes,
                                                     suffix=str(uuid.uuid4()))[1])
            os.remove(self._index_alias)
            os.symlink(self.container.abspath, self._index_alias)
            print(self._index_alias)
        else:
            self._index_alias = self.container.abspath
        # try:finally: ensure the removal of the temporary link
        try:
            matchingfields = self._readfields(handgrip,
                                              getdata=getdata,
                                              footprints_proxy_as_builder=footprints_proxy_as_builder,
                                              get_info_as_json=get_info_as_json)
        finally:
            if config.GRIB_safe_indexes:
                os.unlink(self._index_alias)
            del self._index_alias
        return matchingfields

    def _readfields(self, handgrip,
                    getdata=True,
                    footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                    get_info_as_json=None):
        """Actual method."""
        matchingfields = FieldSet()
        idx = gribapi.grib_index_new_from_file(self._index_alias,
                                               [str(k) for k in handgrip.keys()])  # gribapi str/unicode incompatibility
        # filter
        for k, v in handgrip.items():
            # BUG in gribapi ? type conversion seems not to work for index
            if k == 'indicatorOfTypeOfLevel' and isinstance(v, int):  # GRIB1
                type_conv_GRIB1 = {1:'sfc',
                                   8:'sfc',
                                   # 20:'20',
                                   100:'pl',
                                   102:'sfc',
                                   103:'sfc',
                                   105:'sfc',
                                   109:'ml',
                                   111:'sfc',
                                   112:'sfc',
                                   # 113:'sfc',
                                   # 115:'115',
                                   117:'pv',
                                   200:'sfc',
                                   }
                v = type_conv_GRIB1.get(v, str(v))
            elif k == 'typeOfFirstFixedSurface' and isinstance(v, int):  # GRIB2
                type_conv_GRIB2 = {1:'sfc',
                                   100:'pl',
                                   102:'sfc',
                                   103:'sfc',
                                   109:'pv',
                                   119:'hpl', }
                v = type_conv_GRIB2.get(v, str(v))
            if isinstance(v, six.string_types):  # gribapi str/unicode incompatibility
                v = str(v)
            gribapi.grib_index_select(idx, str(k), v)  # gribapi str/unicode incompatibility
        # load messages
        while True:
            gid = gribapi.grib_new_from_index(idx)
            if gid is None:
                break
            msg = GRIBmessage(('gribid', gid))
            matchingfields.append(msg.asfield(getdata=getdata,
                                              footprints_proxy_as_builder=footprints_proxy_as_builder,
                                              get_info_as_json=get_info_as_json))
            del msg
        gribapi.grib_index_release(idx)
        if len(matchingfields) == 0:
            raise epygramError("no field matching *handgrip* was found.")
        return matchingfields

    @FileResource._openbeforedelayed
    def writefield(self, field,
                   ordering=config.GRIB_default_ordering,
                   packing=None,
                   sample=None,
                   grib_edition=None,
                   other_GRIB_options=None):
        """
        Writes a Field as a GRIBmessage into the GRIB resource.

        :param field: a :class:`epygram.base.Field` instance
        :param ordering: way of ordering data in GRIB, dict of GRIB keys.
        :param packing: options of packing and compression in GRIB (dict).
        :param sample: to use a specific sample GRIB
        :param grib_edition: to force a GRIB edition number (1, 2).
        :param other_GRIB_options: other options to be specified in GRIB,
          as a dict(GRIBkey=value)
        """
        if not isinstance(field, H2DField):
            raise NotImplementedError("'field' argument other than a H2DField.")
        m = GRIBmessage(('field', field),
                        ordering=ordering,
                        packing=packing,
                        sample=sample,
                        grib_edition=grib_edition,
                        other_GRIB_options=other_GRIB_options)
        m.write_to_file(self._file)

    def extractprofile(self, handgrip, lon=None, lat=None,
                       geometry=None,
                       vertical_coordinate=None,
                       interpolation='nearest',
                       cheap_height=True,
                       external_distance=None):
        """
        Extracts a vertical profile from the GRIB resource, given a handgrip
        and the geographic location (*lon*/*lat*) of the profile.

        :param handgrip: MUST define the parameter and the type of levels
        :param lon: the longitude of the desired point.
        :param lat: the latitude of the desired point.
        :param geometry: the geometry on which extract data. If None, it is built from
          lon/lat.
        :param vertical_coordinate defines the requested vertical coordinate of the
          V1DField (cf. `epygram.geometries.vertical_coordinates` possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile at requested lon/lat from the fields grid:\n
          - if 'nearest' (default), extracts profile at the horizontal nearest neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear spline interpolation;
          - if 'cubic', computes profile with horizontal cubic spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = |target_value - external_field.data|
        """
        field3d = fpx.field(fid={'GRIB':handgrip},
                            structure='3D',
                            resource=self, resource_fids=[handgrip])

        if geometry is None:
            if None in [lon, lat]:
                raise epygramError("You must give a geometry or lon *and* lat")
            pointG = field3d.geometry.make_profile_geometry(lon, lat)
        else:
            if lon is not None or lat is not None:
                raise epygramError("You cannot provide lon or lat when geometry is given")
            if geometry.structure != "V1D":
                raise epygramError("geometry must be a V1D")
            pointG = geometry

        profile = self.extract_subdomain(handgrip, pointG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         external_distance=external_distance,
                                         cheap_height=cheap_height,
                                         field3d=field3d)

        return profile

    def extractsection(self, handgrip, end1=None, end2=None,
                       geometry=None,
                       points_number=None,
                       resolution=None,
                       vertical_coordinate=None,
                       interpolation='linear',
                       cheap_height=True):
        """
        Extracts a vertical section from the GRIB resource, given its handgrip
        and the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param handgrip: MUST define the parameter and the type of levels
        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param geometry: is the geometry on which extract data. If None, defaults to
          linearily spaced positions computed from  *points_number*.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. `epygram.geometries.vertical_coordinates` possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
        """

        field3d = fpx.field(fid={'GRIB':handgrip},
                            structure='3D',
                            resource=self, resource_fids=[handgrip])

        if geometry is None:
            if None in [end1, end2]:
                raise epygramError("You must give a geometry or end1 *and* end2")
            sectionG = field3d.geometry.make_section_geometry(end1, end2,
                                                              points_number=points_number,
                                                              resolution=resolution)
        else:
            if end1 is not None or end2 is not None:
                raise epygramError("You cannot provide end1 or end2 when geometry is given")
            if geometry.structure != "V2D":
                raise epygramError("geometry must be a V2D")
            sectionG = geometry

        section = self.extract_subdomain(handgrip, sectionG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         cheap_height=cheap_height,
                                         field3d=field3d)

        return section

    def extract_subdomain(self, handgrip, geometry,
                          vertical_coordinate=None,
                          interpolation='linear',
                          cheap_height=True,
                          external_distance=None,
                          field3d=None):
        """
        Extracts a subdomain from the GRIB resource, given its handgrip
        and the geometry to use.

        :param handgrip: MUST define the parameter and the type of levels
        :param geometry: is the geometry on which extract data.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. `epygram.geometries.vertical_coordinates` possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
          # TODO: not implemented yet
        """
        if field3d is None:
            field3d = fpx.field(fid={'GRIB':handgrip},
                                structure='3D',
                                resource=self, resource_fids=[handgrip])

        subdomain = field3d.extract_subdomain(geometry,
                                              interpolation=interpolation,
                                              exclude_extralevels=True,
                                              external_distance=external_distance)

        # preparation for vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):
            # P => H necessary profiles
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 100 and \
               vertical_coordinate in (102, 103):
                side_profiles = {'t':'t',
                                 'q':'q'}
                for p in sorted(side_profiles.keys(), reverse=True):  # reverse to begin by t
                    try:
                        # try to extract profiles for each lon/lat and each parameter
                        if handgrip.get('shortName') == side_profiles[p]:
                            # already extracted as requested profile
                            side_profiles[p] = subdomain
                        else:
                            if 'typeOfFirstFixedSurface' in handgrip:
                                side_handgrip = {'typeOfFirstFixedSurface':handgrip.get('typeOfFirstFixedSurface')}
                            elif 'indicatorOfTypeOfLevel' in handgrip:
                                side_handgrip = {'indicatorOfTypeOfLevel':handgrip.get('indicatorOfTypeOfLevel')}
                            side_handgrip['shortName'] = p
                            side_profiles[p] = self.extract_subdomain(side_handgrip, geometry,
                                                                      interpolation=interpolation,
                                                                      external_distance=external_distance)
                        side_profiles[p] = side_profiles[p].getdata()
                    except epygramError:
                        # fields not present in file
                        if p in ('t', 'q'):
                            raise epygramError("Temperature and Specific" +
                                               " Humidity must be in" +
                                               " resource.")
                        else:
                            side_profiles[p] = numpy.zeros(side_profiles['t'].shape)
                R = q2R(*[side_profiles[p] for p in
                          ['q']])
                if vertical_coordinate == 102:
                    raise NotImplementedError("vertical_coordinate=={}.".format(vertical_coordinate))
                    # surface_geopotential = geopotential.getvalue_ll(*geometry.get_lonlat_grid(),
                    #                                                 interpolation=interpolation,
                    #                                                 one=False,
                    #                                                 external_distance=external_distance)
                    # del geopotential
                else:
                    surface_geopotential = numpy.zeros(geometry.get_lonlat_grid()[0].size)

            # effective vertical coords conversion
            vertical_mean = 'geometric'
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 100 and \
                 vertical_coordinate in (102, 103):
                subdomain.geometry.vcoordinate = pressure2altitude(subdomain.geometry.vcoordinate,
                                                                   R,
                                                                   side_profiles['t'],
                                                                   vertical_mean,
                                                                   Pdep=side_profiles['pdep'],
                                                                   Phi_surf=surface_geopotential)
            else:
                raise NotImplementedError("this vertical coordinate" +
                                          " conversion.")

        return subdomain

    def what(self, out=sys.stdout,
             mode='one+list',
             sortfields=None,
             details=None,
             **_):
        """
        Writes in file a summary of the contents of the GRIB.

        :param out: the output open file-like object
        :param mode: among ('one+list', 'fid_list', 'what', 'ls', 'mars'), \n
          - 'one+list' = gives the validity/geometry of the first field in
            GRIB, plus the list of fid.
          - 'fid_list' = gives only the fid of each field in GRIB.
          - 'what' = gives the values of the keys from each GRIB message that
            are used to generate an **epygram** field from the message (slower).
          - 'ls' = gives the values of the 'ls' keys from each GRIB message.
          - 'mars' = gives the values of the 'mars' keys from each GRIB message.
        :param sortfields: name of the fid key used to sort fields; e.g. 'typeOfLevel';
          only for *mode* = 'one+list' or 'fid_list'.
        :param details: if 'compression', gives the 'packingType' and 'bitsPerValue'
          parameters of field packing. Only with 'what' mode.
        """
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")

        if mode == 'one+list':
            onefield = self.get_message_at_position(0).asfield(getdata=False)
            g0 = onefield.geometry.footprint_as_dict()
            g0.pop('vcoordinate')
            while True:
                m = self.iter_messages()
                if m is None:
                    break
                try:
                    v = m._read_validity()
                except (epygramError, NotImplementedError):
                    if config.GRIB_ignore_validity_decoding_errors:
                        v = FieldValidity()
                    else:
                        raise
                if v.get() != onefield.validity.get():
                    epylog.error(str(m._read_validity()))
                    epylog.error(str(onefield.validity))
                    raise epygramError("several validities found in file: " +
                                       "'one+list' mode disabled; " +
                                       "try mode 'ls' or 'mars'.")
            onefield.what(out, vertical_geometry=False,
                          cumulativeduration=False,
                          fid=True)

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        listoffields = self.listfields()
        out.write("Number: " + str(len(listoffields)) + "\n")
        if mode in ('what', 'ls', 'mars'):
            out.write(separation_line)
            while True:
                m = self.iter_messages()
                if m is None:
                    break  # end of file
                if mode == 'what':
                    _ = m.asfield(getdata=False)
                    if details == 'compression':
                        m._readattribute('packingType')
                        m._readattribute('bitsPerValue')
                elif mode in ('ls', 'mars'):
                    m.readmessage(mode)
                    m._readattribute('name')
                m_dict = dict(m)
                write_formatted_dict(out, m_dict)
        elif sortfields:
            out.write("sorted by: " + sortfields + "\n")
            out.write(separation_line)
            sortedfields = self.sortfields(sortfields)
            for category in sorted(sortedfields.keys()):
                out.write(sortfields + ": " + str(category) + '\n')
                out.write('--------------------' + '\n')
                for f in sortedfields[category]:
                    f.pop(sortfields)
                    write_formatted_dict(out, f)
                out.write('--------------------' + '\n')
        else:
            out.write(separation_line)
            for f in listoffields:
                write_formatted_dict(out, f)
        out.write(separation_line)
