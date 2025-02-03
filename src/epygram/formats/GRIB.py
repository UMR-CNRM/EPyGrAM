#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains classes for GRIB resource and GRIB individual message,
editions 1 and 2.
"""

__all__ = ['GRIB']

import datetime
import os
import numpy
import copy
import sys
import tempfile
import uuid
import io
import json

import footprints
from footprints import proxy as fpx, FPDict, FPList
from bronx.meteo.conversion import q2R
from bronx.meteo import constants
from epygram.extra import griberies

from epygram import config, epygramError, util
from epygram.base import FieldSet, FieldValidity
from epygram.resources import FileResource
from epygram.util import (Angle, RecursiveObject,
                          separation_line, write_formatted_dict)
from epygram.fields import H2DField
from epygram.geometries import (VGeometry, gauss_latitudes,
                                RegLLGeometry, RotLLGeometry,
                                ProjectedGeometry, GaussGeometry)
from epygram.geometries.VGeometry import pressure2altitude, altitude2height, height2altitude
from epygram.geometries.SpectralGeometry import (SpectralGeometry,
                                                 gridpoint_dims_from_truncation,
                                                 nearest_greater_FFT992compliant_int)

epylog = footprints.loggers.getLogger(__name__)


class LowLevelGRIB(object):
    def __init__(self, GRIB_lowlevel_api):
        self.api_name = GRIB_lowlevel_api.lower()
        if self.api_name == 'eccodes':
            self.samples_path_var = 'ECCODES_SAMPLES_PATH'
            import eccodes  # @UnresolvedImport
            self.api = eccodes
            self.count_in_file = eccodes.codes_count_in_file
            self.any_new_from_file = eccodes.codes_any_new_from_file
            self.release = eccodes.codes_release
            self._set = eccodes.codes_set
            self._set_missing = eccodes.codes_set_missing
            self._new_from_samples = eccodes.codes_grib_new_from_samples
            self.clone = eccodes.codes_clone
            self._get = eccodes.codes_get
            self._get_double_array = eccodes.codes_get_double_array
            self.get_values = eccodes.codes_get_values
            self._set_array = eccodes.codes_set_array
            self.set_values = eccodes.codes_set_values
            self.keys_iterator_new = eccodes.codes_keys_iterator_new
            self.keys_iterator_next = eccodes.codes_keys_iterator_next
            self.keys_iterator_get_name = eccodes.codes_keys_iterator_get_name
            self.keys_iterator_delete = eccodes.codes_keys_iterator_delete
            self.write = eccodes.codes_write
            self._index_new_from_file = eccodes.codes_index_new_from_file
            self._index_select = eccodes.codes_index_select
            self.new_from_index = eccodes.codes_new_from_index
            self.index_release = eccodes.codes_index_release
            self.InternalError = eccodes.CodesInternalError
            self.version = eccodes.__version__
        else:
            raise NotImplementedError("GRIB_lowlevel_api:", GRIB_lowlevel_api)

    @property
    def install_dir(self):
        if 'ECCODES_DIR' in os.environ:
            install_dir = os.environ['ECCODES_DIR']
        else:
            install_dir = griberies.paths.get_eccodes_from_ldconfig()
        return install_dir

    def init_env(self, reset=False):
        """Ensure grib_api/eccodes variables are consistent with inner library."""
        # from python package eccodes 2.38.0, samples and packages are provided as eccodes_memfs
        if self.api.codes_get_api_version() < '2.38.0':
            griberies.paths.complete_grib_samples_paths(self.install_dir,
                                                        reset=reset)
            if len(griberies.paths.get_definition_paths()) > 0:
                griberies.paths.complete_grib_definition_paths(self.install_dir,
                                                               reset=reset)

    #  BELOW: gribapi str/unicode incompatibility
    def set(self, gid, key, value):
        if isinstance(value, (numpy.int64, numpy.int32)):
            value = int(value)  # eccodes does not accept numpy.int...
        return self._set(gid, str(key), value)

    def set_missing(self, gid, key):
        return self._set_missing(gid, str(key))

    def get(self, gid, key, *args):
        return self._get(gid, str(key), *args)

    def get_double_array(self, gid, key):
        return self._get_double_array(gid, str(key))

    def set_array(self, gid, key, *args):
        return self._set_array(gid, str(key), *args)

    def index_new_from_file(self, filename, *args):
        return self._index_new_from_file(str(filename), *args)

    def index_select(self, idx, key, value):
        return self._index_select(idx, str(key), value)

    def new_from_samples(self, sample):
        return self._new_from_samples(str(sample))


lowlevelgrib = LowLevelGRIB(config.GRIB_lowlevel_api)

# conversion of surface types for geometry purposes only
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


def sorted_GRIB2_fid(fid):
    allsorted = copy.copy(GRIBmessage._fid_keys[2])
    allsorted.insert(allsorted.index('topLevel') + 1,
                     'scaledValueOfFirstFixedSurface')
    allsorted.insert(allsorted.index('topLevel') + 2,
                     'scaleFactorOfFirstFixedSurface')
    allsorted.insert(allsorted.index('bottomLevel') + 1,
                     'scaledValueOfSecondFixedSurface')
    allsorted.insert(allsorted.index('topLevel') + 2,
                     'scaleFactorOfSecondFixedSurface')
    allsorted.append('lengthOfTimeRange')
    allsorted.append('indicatorOfUnitForTimeRange')
    allsorted.append('typeOfStatisticalProcessing')
    allsorted.extend(GRIBmessage._satellite_imagery_keys)
    _sorted = []
    for k in allsorted:
        if k in fid:
            _sorted.append(k)
    return _sorted


class GRIBmessage(RecursiveObject, dict):
    """
    Class implementing a GRIB message as an object.
    """

    _fid_keys = {1:['name',
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
                 2:['editionNumber',
                    'name',
                    'shortName',
                    'discipline',
                    'parameterCategory',
                    'parameterNumber',
                    'typeOfFirstFixedSurface',
                    'level',
                    'topLevel',
                    'typeOfSecondFixedSurface',
                    'bottomLevel',
                    'tablesVersion',
                    'productDefinitionTemplateNumber']}
    _satellite_imagery_keys = ['NB',
                               'satelliteSeries',
                               'satelliteNumber',
                               'instrumentType',
                               'scaleFactorOfCentralWaveNumber',
                               'scaledValueOfCentralWaveNumber']
    # Class methods ------------------------------------------------------------

    @classmethod
    def specific_fid_keys_for(cls,
                              productDefinitionTemplateNumber=0,
                              scaleFactorOfFirstFixedSurface=0,
                              scaleFactorOfSecondFixedSurface=0):
        """
        Get specific fid keys according to **productDefinitionTemplateNumber**
        and **scaleFactorOfFirstFixedSurface**
        and **scaleFactorOfSecondFixedSurface**.

        (GRIB2 only).
        """
        if productDefinitionTemplateNumber in (32, 33):
            specific_keys = cls._satellite_imagery_keys
        elif productDefinitionTemplateNumber in (8, 11):
            specific_keys = ['lengthOfTimeRange',
                             # 'indicatorOfUnitForTimeRange',  # FIXME: pb with eccodes index
                             'typeOfStatisticalProcessing',  # FIXME: pb with eccodes index
                             ]
        else:
            specific_keys = []
        if scaleFactorOfFirstFixedSurface != 0:
            specific_keys.append('scaleFactorOfFirstFixedSurface')
            specific_keys.append('scaledValueOfFirstFixedSurface')
        if scaleFactorOfSecondFixedSurface != 0:
            specific_keys.append('scaleFactorOfSecondFixedSurface')
            specific_keys.append('scaledValueOfSecondFixedSurface')
        return specific_keys

    @classmethod
    def fid_keys_for(cls, editionNumber,
                     productDefinitionTemplateNumber=0,
                     scaleFactorOfFirstFixedSurface=0,
                     scaleFactorOfSecondFixedSurface=0):
        """
        Get fid keys according to
        **editionNumber** and **productDefinitionTemplateNumber**,
        **scaleFactorOfFirstFixedSurface**
        and **scaleFactorOfSecondFixedSurface**.
        """
        fid_keys = copy.copy(cls._fid_keys[editionNumber])
        add_keys = []
        remove_keys = []
        if productDefinitionTemplateNumber in (32, 33):
            remove_keys = ['typeOfFirstFixedSurface',
                           'level',
                           'topLevel',
                           'typeOfSecondFixedSurface',
                           'bottomLevel']
            add_keys = cls.specific_fid_keys_for(productDefinitionTemplateNumber=productDefinitionTemplateNumber)
        elif productDefinitionTemplateNumber in (8, 11):
            add_keys = cls.specific_fid_keys_for(productDefinitionTemplateNumber=productDefinitionTemplateNumber)
        if scaleFactorOfFirstFixedSurface != 0 or scaleFactorOfSecondFixedSurface != 0:
            add_keys = cls.specific_fid_keys_for(scaleFactorOfFirstFixedSurface=scaleFactorOfFirstFixedSurface,
                                                 scaleFactorOfSecondFixedSurface=scaleFactorOfSecondFixedSurface)
            remove_keys = ['level']
        for k in remove_keys:
            fid_keys.remove(k)
        fid_keys.extend(add_keys)
        return fid_keys

    @classmethod
    def _GRIB1_default_sample(cls, field, required_packingType=None):
        """
        Get default sample, according in order to:
        1/ spectralness
        2/ model levels
        3/ packing
        """
        sample = griberies.defaults.GRIB1_sample
        if field.spectral:
            sample = 'sh_ml_grib1'
        elif field.geometry.vcoordinate.typeoffirstfixedsurface == 119:
            sample = 'reduced_rotated_gg_ml_grib1'
        elif required_packingType is not None:
            guess_sample = 'GRIB1_{}'.format(required_packingType)  # an epygram sample with this packing
            if guess_sample + '.tmpl' in os.listdir(config.GRIB_epygram_samples_path):  # sample is available with this packing
                sample = guess_sample
        return sample

    # Special methods ----------------------------------------------------------

    def __init__(self, source,
                 ordering=griberies.defaults.GRIB1_ordering,
                 packing=None,
                 sample=None,
                 grib_edition=None,
                 other_GRIB_options={},
                 interpret_comment=False,
                 set_misc_metadata=True):
        """
        Initialize a GRIBmessage from either sources.

        :param source: being a tuple of either form:\n
          - ('file', '*filename*' [, *offset_position*])
            *filename* being a relative or absolute path to the file it is read
            in. n = *offset_position*, if given, enables to read the n+1'th GRIB
            message in file. Defaults to 0. Negative value counts from the end.
          - ('field', :class:`epygram.fields.H2DField`)
          - ('gribid', *gribid*)
            *gribid* being an integer, refering to the *gribid* of an eccodes
            message in memory.
          - ('sample', '*samplename*')
            *samplename* being the name of the sample from which to be
            generated.

        In case **source** is a field, some options can be forced:

        :param sample: to use a specific sample GRIB.
          Specific syntax 'file:$filename$' takes the first message
          in $filename$ as sample.
        :param grib_edition: to force a GRIB edition number (1, 2).
        :param other_GRIB_options: other options to be specified in GRIB,
            as a dict(GRIBkey=value).
            From v1.3.9, any GRIB2 key should be given through this rather than
            packing/ordering.

        GRIB1

        :param ordering: flattening of 2D data
        :param packing: options of packing and compression in GRIB (dict).

        GRIB2 (new way)

        :param interpret_comment: set additional key/values taken from
            field.comment (interpreted as json)
        :param set_misc_metadata: set additional key/values taken from
            field.misc_metadata
        """
        super(GRIBmessage, self).__init__()
        built_from = source[0]
        if built_from == 'file':
            filename = source[1]
            msg_position = source[2] if len(source) == 3 else None
            self._init_from_file(filename, msg_position=msg_position)
        elif built_from == 'field':
            field = source[1]
            self._init_from_field(field,
                                  grib_edition,
                                  sample,
                                  other_GRIB_options,
                                  packing,
                                  ordering,
                                  interpret_comment,
                                  set_misc_metadata)
        elif built_from == 'gribid':
            self._gid = source[1]
        elif built_from == 'sample':
            self._gid = self._clone_from_sample(source[1])
        else:
            raise NotImplementedError("source={}".format(str(source)))

    def __del__(self):
        try:
            lowlevelgrib.release(self._gid)
        except (lowlevelgrib.InternalError, AttributeError):
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
            if lowlevelgrib.version <= '1.10.4':
                if any([isinstance(value, t) for t in [numpy.float32, numpy.float64]]):
                    value = float(value)
                elif any([isinstance(value, t) for t in [numpy.int8, numpy.int16, numpy.int32, numpy.int64]]):
                    value = int(value)
            if isinstance(value, str):  # gribapi str/unicode incompatibility
                v = str(value)
            else:
                v = value
            lowlevelgrib.set(self._gid, key, v)
        else:
            lowlevelgrib.set_missing(self._gid, key)
        super(GRIBmessage, self).__setitem__(key, value)
        epylog.debug("DEBUG: set {} = {}".format(key, value))

    # Low level methods --------------------------------------------------------

    def _init_from_file(self, filename, msg_position=None):
        """Initialize message from message in file, optionally with given position in file."""
        self._file = open(filename, 'rb')
        # if position specified, go to position of message
        if msg_position is not None:
            n = msg_position
            if n < 0:
                N = lowlevelgrib.count_in_file(self._file)
                n = N + n
            for _ in range(n):
                gid = lowlevelgrib.any_new_from_file(self._file,
                                                     headers_only=True)
                lowlevelgrib.release(gid)
        # load message in memory and save gribid
        self._gid = lowlevelgrib.any_new_from_file(self._file)
        self._file.close()

    def _init_from_field(self, field,
                         grib_edition=None,
                         sample=None,
                         other_GRIB_options={},
                         # GRIB1/old way
                         packing=None,
                         ordering=None,
                         # GRIB2/new way
                         interpret_comment=False,
                         set_misc_metadata=True):
        """Initialize message from epygram field (and sample)."""
        if grib_edition == 2 or (grib_edition is None and
                                     'GRIB2' in field.fid):
            # new way for GRIB2
            # first, transfer options
            other_GRIB_options = copy.copy(other_GRIB_options)
            if packing is not None:
                for k,v in packing.items():
                    other_GRIB_options.setdefault(k,v)
            if ordering is not None:
                for k,v in ordering.items():
                    other_GRIB_options.setdefault(k,v)
            self._GRIB2_set(field, sample,
                            interpret_comment=interpret_comment,
                            set_misc_metadata=set_misc_metadata,
                            **other_GRIB_options)
        elif grib_edition == 1 or (grib_edition is None and
                                   'GRIB1' in field.fid):
            # old way for GRIB1
            self._GRIB1_set(field,
                            ordering=ordering,
                            packing=packing,
                            sample=sample,
                            other_GRIB_options=other_GRIB_options)
        else:
            raise epygramError("field.fid must have a 'GRIB1' or 'GRIB2' identifier")

    def _clone_from_sample(self, sample):
        """Clone a sample GRIB message."""
        sample_gid = lowlevelgrib.new_from_samples(sample)
        gid = lowlevelgrib.clone(sample_gid)
        lowlevelgrib.release(sample_gid)
        return gid

    def _clone_from_file(self, filename):
        """Clone first GRIB message from file."""
        with open(filename, 'rb') as f:
            original_gid = lowlevelgrib.any_new_from_file(f)
            gid = lowlevelgrib.clone(original_gid)
            lowlevelgrib.release(original_gid)
        return gid

    def write_to_file(self, ofile):
        """
        *ofile* being an open file-like object designing the physical GRIB
        file to be written to.
        """
        lowlevelgrib.write(self._gid, ofile)

    def _readattribute(self, attribute, array=False, fatal=True):
        """
        Actual access to the attributes.
        :param attribute: key to read
        :param array: attribute is an array
        :param fatal: raise errors or ignore failed keys
        """
        attr = None
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
                        attr = lowlevelgrib.get(self._gid, attribute, int)
                    else:
                        attr = lowlevelgrib.get(self._gid, attribute)
                except lowlevelgrib.InternalError as e:
                    # differenciation not well done... PB in gribapi
                    if str(e) == 'Passed array is too small':
                        try:
                            attr = lowlevelgrib.get_double_array(self._gid, attribute)
                        except lowlevelgrib.InternalError as e:
                            if fatal:
                                raise e
                            else:
                                epylog.warning('GRIB error ignored:' + str(e))
                    elif fatal:
                        raise type(e)(str(e) + ' : "' + attribute + '"')
            else:
                try:
                    attr = lowlevelgrib.get_double_array(self._gid, attribute)
                except lowlevelgrib.InternalError as e:
                    if fatal:
                        raise e
                    else:
                        epylog.warning('GRIB error ignored:' + str(e))
        else:
            try:
                attr = lowlevelgrib.get_values(self._gid)
            except lowlevelgrib.InternalError as e:
                if fatal:
                    raise e
                else:
                    epylog.warning('GRIB error ignored:' + str(e))
        if attr is not None:
            super(GRIBmessage, self).__setitem__(attribute, attr)

    def _set_array_attribute(self, attribute, value):
        """Setter for array attributes."""
        lowlevelgrib.set_array(self._gid, attribute, value)

    def get(self, key, default=None):
        """Same as dict.get(), but try to read attribute first."""
        try:
            value = self.__getitem__(key)
        except (KeyError, lowlevelgrib.InternalError):
            value = default
        return value

    # Write GRIB2 (new way) ----------------------------------------------------

    def _GRIB2_set(self, field,
                   sample=None,
                   interpret_comment=False,
                   set_misc_metadata=True,
                   **other_GRIB_options):
        """Set up a GRIB2 message from a field."""
        # clone from sample or file
        if sample is None:
            sample = field.geometry.suggested_GRIB2_sample(field.spectral)
        if sample.startswith('file:'):
            # clone from file
            self._gid = self._clone_from_file(sample.replace('file:', ''))
        else:
            try:
                self._gid = self._clone_from_sample(sample)
            except Exception:
                # sample not genuinely found by eccodes:
                if sample + '.tmpl' in os.listdir(config.GRIB_epygram_samples_path):
                    # try in epygram samples
                    fsample = os.path.join(config.GRIB_epygram_samples_path, sample) + '.tmpl'
                    self._gid = self._clone_from_file(fsample)  # clone it from file
                else:
                    # or try manually exploring other paths from SAMPLES_PATH
                    for d in griberies.paths.get_samples_paths():
                        fsample = os.path.join(d, sample) + '.tmpl'
                        try:
                            self._gid = self._clone_from_file(fsample)
                        except Exception:
                            raise epygramError('sample {} not found; check {} environment variable'.
                                               format(sample, lowlevelgrib.samples_path_var))
        # interpret comment
        if interpret_comment and field.comment is not None:
            comment_options = json.loads(field.comment)
            for k, v in comment_options.items():
                other_GRIB_options.setdefault(k, v)
        if set_misc_metadata:
            for k, v in field.misc_metadata.items():
                other_GRIB_options.setdefault(k, v)
        # then set keys/values
        self._untouch = set()
        self._GRIB2_set_sections(field, **other_GRIB_options)
        # set other_GRIB_options which have not been set
        for k, v in other_GRIB_options.items():
            if k not in self._untouch:
                try:
                    self[k] = v
                except Exception as e:
                    epylog.info("failed to set {}={}: {}".format(k, v, str(e)))
                    pass

    def _GRIB2_set_sections(self, field, **other_GRIB_options):
        """Set k/v for sections of GRIB2, one by one."""
        if not isinstance(field, H2DField):
            raise NotImplementedError("writing of non-H2DField.")
        self._GRIB2_set_section0(field)
        self._GRIB2_set_section1(field, **other_GRIB_options)
        self._GRIB2_set_section3(field, **other_GRIB_options)
        self._GRIB2_set_section4(field, **other_GRIB_options)
        self._GRIB2_set_section5(field, **other_GRIB_options)
        self._GRIB2_set_section6_7(field, **other_GRIB_options)

    def _GRIB2_set_section0(self, field):
        self['editionNumber'] = 2
        self['discipline'] = field.fid['GRIB2']['discipline']

    def _GRIB2_set_section1(self, field, **other_GRIB_options):
        """Set k/v for section 1."""
        for k in ['centre', 'subCentre',
                  'tablesVersion', 'localTablesVersion']:
            self.set_from(k, [other_GRIB_options,  # forcing k/v argument
                              field.fid['GRIB2'],  # in fid, if present/relevant
                              {'centre':config.GRIB_default_centre},
                              griberies.defaults.GRIB2_keyvalue[1]])  # if present in default
        # basis cutoff
        assert len(field.validity) == 1
        if other_GRIB_options.get('significanceOfReferenceTime', 1) == 1:
            basis = field.validity.getbasis(fmt='IntStr')
            self['year'] = int(basis[:4])
            self['month'] = int(basis[4:6])
            self['day'] = int(basis[6:8])
            self['hour'] = int(basis[8:10])
            self['minute'] = int(basis[10:12])
            self['second'] = int(basis[12:14])
        else:
            raise NotImplementedError("significanceOfReferenceTime != 1")
        # others
        for k in ['productionStatusOfProcessedData',
                  'typeOfProcessedData']:
            self.set_from(k, [other_GRIB_options,  # forcing k/v argument
                              griberies.defaults.GRIB2_keyvalue[1]])  # if present in default

    def _GRIB2_set_section3(self, field, **other_GRIB_options):
        """Set k/v for section 3."""
        self._GRIB2_set_earth_geometry(field.geometry)
        if field.spectral:
            self._GRIB2_set_spectral_geometry(field.geometry)
        else:
            self._GRIB2_set_horizontal_geometry(field.geometry,
                                                **other_GRIB_options)
            self._GRIB2_set_data_ordering(field.geometry, **other_GRIB_options)

    def _GRIB2_set_earth_geometry(self, geometry):
        """Set earth geometry."""
        if not hasattr(geometry, 'geoid'):
            shapeOfTheEarth = griberies.defaults.GRIB2_keyvalue[3].get('shapeOfTheEarth')
            epylog.warning("geometry has no geoid: assume 'shapeOfTheEarth' = {}.".format(shapeOfTheEarth))
            self['shapeOfTheEarth'] = shapeOfTheEarth
        else:
            if all([geometry.geoid.get(axis) == griberies.tables.pyproj_geoid_shapes[6][axis]
                    for axis in ('a', 'b')]):
                self['shapeOfTheEarth'] = 6
            else:
                found = False
                for s, g in griberies.tables.pyproj_geoid_shapes.items():
                    if s in (0, 2, 4, 5, 6, 8, 9) and geometry.geoid == g:
                        self['shapeOfTheEarth'] = s
                        found = True
                        break
                if not found:
                    radius = geometry.geoid.get('geoidradius')
                    if radius is None:
                        radius = geometry.geoid.get('a')
                        if radius is None or radius != geometry.geoid.get('b'):
                            radius = None
                    if radius is not None:
                        self['shapeOfTheEarth'] = 1
                        self['scaleFactorOfRadiusOfSphericalEarth'] = 1
                        self['scaledValueOfRadiusOfSphericalEarth'] = radius
                    else:
                        a = geometry.geoid.get('a', geometry.geoid.get('geoidmajor'))
                        b = geometry.geoid.get('b', geometry.geoid.get('geoidminor'))
                        if a is None or b is None:
                            raise epygramError("unable to encode geoid.")
                        else:
                            self['shapeOfTheEarth'] = 7
                            self['scaleFactorOfMajorAxisOfOblateSpheroidEarth'] = 1
                            self['scaledValueOfMajorAxisOfOblateSpheroidEarth'] = a
                            self['scaleFactorOfMinorAxisOfOblateSpheroidEarth'] = 1
                            self['scaledValueOfMinorAxisOfOblateSpheroidEarth'] = b

    def _GRIB2_set_horizontal_geometry(self, geometry, **other_GRIB_options):
        """Set H2D horizontal geometry."""
        for k in ['gridDefinitionTemplateNumber',]:
            self.set_from(k, [other_GRIB_options])
        # dimensions -----------------------------------------------------------
        if geometry.rectangular_grid:
            self['Ni'] = geometry.dimensions['X']
            self['Nj'] = geometry.dimensions['Y']
            # if lower than sample, it is not modified by ecCodes !
            self['numberOfDataPoints'] = geometry.dimensions['X'] * geometry.dimensions['Y']
        else:
            pass  # done below
        # regular_lonlat -------------------------------------------------------
        if geometry.name == 'regular_lonlat':
            self['gridType'] = 'regular_ll'
            self['iDirectionIncrementInDegrees'] = geometry.grid['X_resolution'].get('degrees')
            self['jDirectionIncrementInDegrees'] = geometry.grid['Y_resolution'].get('degrees')
        # rotated_lonlat -------------------------------------------------------
        elif geometry.name == 'rotated_lonlat':
            self['gridType'] = 'rotated_ll'
            self['iDirectionIncrementInDegrees'] = geometry.grid['X_resolution'].get('degrees')
            self['jDirectionIncrementInDegrees'] = geometry.grid['Y_resolution'].get('degrees')
            self['longitudeOfSouthernPoleInDegrees'] = geometry.grid['southern_pole_lon'].get('degrees')
            self['latitudeOfSouthernPoleInDegrees'] = geometry.grid['southern_pole_lat'].get('degrees')
            self['angleOfRotation'] = geometry.grid['rotation'].get('degrees')
        # projected ------------------------------------------------------------
        elif geometry.projected_geometry:
            self['gridType'] = geometry.name
            # lambert ----------------------------------------------------------
            if geometry.name == 'lambert':
                # lambda0 in geoid in case shapeOfTheEarth == 9
                lambda0 = geometry.geoid.get('lambda0', geometry.projection['reference_lon'].get('degrees'))
                self['LoVInDegrees'] = util.positive_longitude(lambda0)
                lat_0 = geometry.projection.get('reference_lat', None)
                lat_1 = geometry.projection.get('secant_lat1', None)
                lat_2 = geometry.projection.get('secant_lat2', None)
                if lat_0 is None:
                    lat_1 = lat_1.get('degrees')
                    lat_2 = lat_2.get('degrees')
                    lat_0 = (lat_1 + lat_2) / 2.
                elif lat_1 is None and lat_2 is None:
                    lat_0 = lat_0.get('degrees')
                    lat_1 = lat_2 = lat_0
                self['LaD'] = int(lat_0 * 1e6)  # LaDInDegrees sometimes caused errors
                self['Latin1InDegrees'] = lat_1
                self['Latin2InDegrees'] = lat_2
                if abs(geometry.projection['rotation'].get('degrees')) > config.epsilon:
                    raise NotImplementedError("*rotation* attribute of projection != 0.")
            # polar_stereographic ----------------------------------------------
            elif geometry.name == 'polar_stereographic':
                if abs(geometry.projection['reference_lat'].get('degrees') - 90.) < config.epsilon:
                    self['projectionCentreFlag'] = 0
                else:
                    self['projectionCentreFlag'] = 1
                lat_ts = geometry.projection.get('secant_lat', geometry.projection['reference_lat']).get('degrees')
                try:
                    self['LaDInDegrees'] = lat_ts
                except lowlevelgrib.InternalError:
                    if abs(lat_ts -
                           numpy.copysign(60., geometry.projection['reference_lat'].get('degrees'))) > config.epsilon:
                        raise epygramError("unable to write polar stereographic geometry to GRIB1 if *secant_lat* is not +/-60 degrees.")
                self['orientationOfTheGridInDegrees'] = geometry.projection['reference_lon'].get('degrees')
                if abs(geometry.projection['rotation'].get('degrees')) > config.epsilon:
                    raise NotImplementedError("*rotation* attribute of projection != 0.")
            # mercator ---------------------------------------------------------
            elif geometry.name == 'mercator':
                lat_ts = geometry.projection.get('secant_lat', geometry.projection['reference_lat']).get('degrees')
                try:
                    self['LaDInDegrees'] = lat_ts
                except lowlevelgrib.InternalError:
                    pass  # TOBECHECKED: in GRIB1, lat_ts=0. ?
                if abs(geometry.projection['rotation'].get('degrees')) > config.epsilon:
                    raise NotImplementedError("*rotation* attribute of projection != 0.")
            if geometry.name in ('lambert', 'polar_stereographic'):
                self['DxInMetres'] = geometry.grid['X_resolution']
                self['DyInMetres'] = geometry.grid['Y_resolution']
            else:
                self['DiInMetres'] = geometry.grid['X_resolution']
                self['DjInMetres'] = geometry.grid['Y_resolution']
        # gauss ----------------------------------------------------------------
        elif 'gauss' in geometry.name:
            if geometry.name == 'rotated_reduced_gauss':
                self['gridType'] = 'reduced_stretched_rotated_gg'
                self['longitudeOfStretchingPoleInDegrees'] = geometry.grid['pole_lon'].get('degrees')
                self['latitudeOfStretchingPoleInDegrees'] = geometry.grid['pole_lat'].get('degrees')
                self['stretchingFactor'] = geometry.grid['dilatation_coef']
            elif geometry.name == 'reduced_gauss':
                if geometry.grid.get('dilatation_coef') != 1.:
                    self['gridType'] = 'reduced_stretched_gg'
                    self['stretchingFactor'] = geometry.grid['dilatation_coef']
                else:
                    self['gridType'] = 'reduced_gg'
            self['global'] = 1
            self['Nj'] = geometry.dimensions['lat_number']
            self['N'] = geometry.dimensions['lat_number'] // 2
            self._set_array_attribute('pl', geometry.dimensions['lon_number_by_lat'])
        else:
            raise NotImplementedError("not yet.")

    def _GRIB2_set_spectral_geometry(self, spectral_geometry):
        """Set spectral geometry."""
        if spectral_geometry.space == 'bi-fourier':
            # FIXME: update when spectral LAM accepted by WMO and implemented in grib_api
            self['gridType'] = 'sh'  # in bi-fourier case, this is a bypass
            self['J'] = max(spectral_geometry.truncation['in_X'],
                            spectral_geometry.truncation['in_Y'])
            self['K'] = self['J']
            self['M'] = self['J']
        elif spectral_geometry.space == 'legendre':
            self['gridType'] = 'sh'
            self['sphericalHarmonics'] = 1
            self['J'] = spectral_geometry.truncation['max']
            self['K'] = spectral_geometry.truncation['max']
            self['M'] = spectral_geometry.truncation['max']
        else:
            raise NotImplementedError('spectral_geometry.name == ' + spectral_geometry.name)

    def _GRIB2_set_data_ordering(self, geometry, **other_GRIB_options):
        """Set data ordering and coordinates of first/last gridpoints."""
        for k in ['iScansNegatively',
                  'jScansPositively',
                  'jPointsAreConsecutive',
                  'scanningMode']:
            self.set_from(k, [other_GRIB_options,
                              griberies.defaults.GRIB2_keyvalue[3]])
        ordering_str = 'ordering: {}={}, {}={}'.format('iScansNegatively',
                                                       self['iScansNegatively'],
                                                       'jScansPositively',
                                                       self['jScansPositively'])
        if geometry.rectangular_grid:
            corners = geometry.gimme_corners_ll()
            if geometry.name == 'rotated_lonlat':  # special case: coordinates to be given in the rotated referential
                for k, v in corners.items():
                    corners[k] = geometry.ll2xy(*v)
            if self['iScansNegatively'] == 0 and \
               self['jScansPositively'] == 0:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ul'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['lr'][1]
            elif self['iScansNegatively'] == 0 and \
                 self['jScansPositively'] == 1:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ll'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ur'][1]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 0:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ur'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ll'][1]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 1:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['lr'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ul'][1]
            else:
                raise NotImplementedError(ordering_str)
        else:
            if self['iScansNegatively'] == 0 and \
               self['jScansPositively'] == 0:
                self['latitudeOfFirstGridPointInDegrees'] = geometry.grid['latitudes'][0].get('degrees')
                self['longitudeOfFirstGridPointInDegrees'] = 0.
                self['latitudeOfLastGridPointInDegrees'] = geometry.grid['latitudes'][-1].get('degrees')
                self['longitudeOfLastGridPointInDegrees'] = 360. - 360. / geometry.dimensions['lon_number_by_lat'][-1]
            else:
                raise NotImplementedError(ordering_str)

    def _GRIB2_set_section4(self, field, **other_GRIB_options):
        """Set k/v for section 4."""
        self._GRIB2_set_product(field, **other_GRIB_options)
        template = field.fid['GRIB2'].get('productDefinitionTemplateNumber', 0)
        if template in (32, 33):  # 32,33 == simulated satellite imagery
            self._GRIB2_set_simulated_satellite_imagery(field.fid['GRIB2'],
                                                        template,
                                                        **other_GRIB_options)
        else:
            self._GRIB2_set_vertical_geometry(field.geometry, field.fid['GRIB2'],
                                              **other_GRIB_options)
        self._GRIB2_set_validity(field, **other_GRIB_options)

    def _GRIB2_set_product(self, field, **other_GRIB_options):
        """Set product and process."""
        if field.validity.cumulativeduration():
            tpln = field.fid['GRIB2'].get('productDefinitionTemplateNumber',
                                          other_GRIB_options.get('productDefinitionTemplateNumber', 8))
        else:
            tpln = field.fid['GRIB2'].get('productDefinitionTemplateNumber',
                                          other_GRIB_options.get('productDefinitionTemplateNumber', 0))
        self['productDefinitionTemplateNumber'] = tpln
        self['parameterCategory'] = field.fid['GRIB2']['parameterCategory']
        self['parameterNumber'] = field.fid['GRIB2']['parameterNumber']
        for k in ['typeOfGeneratingProcess',
                  'backgroundProcess',
                  'generatingProcessIdentifier']:
            self.set_from(k, [other_GRIB_options,
                              griberies.defaults.GRIB2_keyvalue[4]])
        self._untouch.add('productDefinitionTemplateNumber')

    def _GRIB2_set_simulated_satellite_imagery(self, fid, template_number,
                                               **other_GRIB_options):
        """Set specific keys for simulated satellite imagery."""
        for k in self.specific_fid_keys_for(template_number):
            self.set_from(k, [other_GRIB_options,
                              fid,
                              griberies.defaults.GRIB2_keyvalue[4]],
                          fatal=False)  # FIXME: known issue...

    def _GRIB2_set_vertical_geometry(self, geometry, fid, **other_GRIB_options):
        """Set vertical geometry."""
        if len(geometry.vcoordinate.levels) > 1:
            raise epygramError("field has more than one level")
        for k in ['scaleFactorOfFirstFixedSurface', 'scaleFactorOfSecondFixedSurface']:
            self.set_from(k, [other_GRIB_options,
                              griberies.defaults.GRIB2_keyvalue[4]])
        # CLEANME: if not needed for some obscure GRIB reasons before setting typeOfFirstFixedSurface
        # self['level'] = geometry.vcoordinate.levels[0]
        self['typeOfFirstFixedSurface'] = geometry.vcoordinate.typeoffirstfixedsurface
        if self['typeOfFirstFixedSurface'] == 119:  # HybridPressure
            self['numberOfVerticalCoordinateValues'] = len(geometry.vcoordinate.grid['gridlevels'])
            ab = [ab[1]['Ai'] for ab in geometry.vcoordinate.grid['gridlevels']] + \
                 [ab[1]['Bi'] for ab in geometry.vcoordinate.grid['gridlevels']]
            self._set_array_attribute('pv', ab)
        if hasattr(geometry.vcoordinate, 'typeofsecondfixedsurface'):
            self['typeOfSecondFixedSurface'] = geometry.vcoordinate.typeofsecondfixedsurface
        else:
            self['typeOfSecondFixedSurface'] = griberies.defaults.GRIB2_keyvalue[4]['typeOfSecondFixedSurface']
        # set levels
        if hasattr(geometry.vcoordinate, 'toplevel'):
            level = geometry.vcoordinate.toplevel
        else:
            level = fid.get('topLevel', geometry.vcoordinate.levels[0])
        self['scaledValueOfFirstFixedSurface'] = level * 10**self['scaleFactorOfFirstFixedSurface']
        if self['typeOfSecondFixedSurface'] != 255:
            if hasattr(geometry.vcoordinate, 'bottomlevel'):
                level = geometry.vcoordinate.bottomlevel
            else:
                level = fid.get('bottomLevel', geometry.vcoordinate.levels[0])
            self['scaledValueOfSecondFixedSurface'] = level * 10**self['scaleFactorOfSecondFixedSurface']

    def _GRIB2_set_validity(self, field, **other_GRIB_options):
        """Set temporal validity."""
        for k in ['hoursAfterDataCutoff',
                  'minutesAfterDataCutoff']:
            self.set_from(k, [other_GRIB_options,
                              griberies.defaults.GRIB2_keyvalue[4]])
        term = field.validity.term()
        unit = other_GRIB_options.get('indicatorOfUnitOfTimeRange',
                                      griberies.defaults.GRIB2_keyvalue[4]['indicatorOfUnitOfTimeRange'])
        cumulative = field.validity.cumulativeduration()
        if cumulative:
            cumul_unit = other_GRIB_options.get('indicatorOfUnitForTimeRange',
                                                griberies.defaults.GRIB2_keyvalue[4]['indicatorOfUnitForTimeRange'])
            start = (term - cumulative).total_seconds()
            length = cumulative.total_seconds()
            if cumul_unit == 0:
                start = start // 60
                length = length // 60
            elif cumul_unit == 1:
                start = start // 3600
                length = length // 3600
            elif cumul_unit != 13:
                raise NotImplementedError("indicatorOfUnitForTimeRange=={}".format(cumul_unit))
            self['lengthOfTimeRange'] = length
            self['indicatorOfUnitForTimeRange'] = cumul_unit
            self['forecastTime'] = start
            end = field.validity.get()
            self['yearOfEndOfOverallTimeInterval'] = end.year
            self['monthOfEndOfOverallTimeInterval'] = end.month
            self['dayOfEndOfOverallTimeInterval'] = end.day
            self['hourOfEndOfOverallTimeInterval'] = end.hour
            self['minuteOfEndOfOverallTimeInterval'] = end.minute
            self['secondOfEndOfOverallTimeInterval'] = end.second
            statistical_process = field.validity[0].statistical_process_on_duration(asGRIB2code=True)
            if statistical_process is not None:
                self['typeOfStatisticalProcessing'] = statistical_process
            time_increment = field.validity.statistical_time_increment()
            time_increment_unit = other_GRIB_options.get('indicatorOfUnitForTimeIncrement',
                                                         griberies.defaults.GRIB2_keyvalue[4]['indicatorOfUnitForTimeIncrement'])
            if time_increment:
                time_increment = time_increment.total_seconds()
                if time_increment_unit == 0:
                    time_increment = time_increment // 60
                if time_increment_unit == 1:
                    time_increment = time_increment // 3600
                elif time_increment_unit != 13:
                    raise NotImplementedError("indicatorOfUnitForTimeIncrement=={}".format(time_increment_unit))
                self['timeIncrement'] = time_increment
                self['indicatorOfUnitForTimeIncrement'] = time_increment_unit
        else:
            term = term.total_seconds()
            if unit == 0:
                term = term // 60
            elif unit == 1:
                term = term // 3600
            elif unit != 13:
                raise NotImplementedError("indicatorOfUnitOfTimeRange=={}".format(unit))
            self['forecastTime'] = term
        self['indicatorOfUnitOfTimeRange'] = unit
        for k in ['indicatorOfUnitOfTimeRange',
                  'indicatorOfUnitForTimeRange'
                  'indicatorOfUnitForTimeIncrement']:
            self._untouch.add(k)

    def _GRIB2_set_section5(self, field, **other_GRIB_options):
        """Set k/v for section 5."""
        order = ['packingType', 'complexPacking', 'boustrophedonicOrdering',
                 'bitsPerValue']
        packing = {k:other_GRIB_options.get(k, griberies.defaults.GRIB2_keyvalue[5].get(k, None))
                   for k in order}
        packing = {k:v for k,v in packing.items() if v is not None}
        if field.spectral:
            packing = {'packingType':'spectral_simple',
                       'bitsPerValue':24}
            epylog.info("field is spectral: shifting to {}".format(str(packing)))
        # reset packing "on the fly" if field is uniform
        elif field.max() - field.min() < config.epsilon:
            packing = {'packingType':'grid_simple'}
            epylog.info("field is uniform: shifting to {}".format(str(packing)))
        else:
            bitsPerValue = packing.get('bitsPerValue')
            # problem with bitsPerValue = 30 at least
            if bitsPerValue is not None and bitsPerValue > config.GRIB_max_bitspervalue:
                msg = "requested bitsPerValue={} higher than {} : may be untrustful.".format(
                    bitsPerValue, config.GRIB_max_bitspervalue)
                epylog.warning(msg)
                if config.GRIB_force_bitspervalue:
                    msg = "config.GRIB_force_bitspervalue is True => force bitsPerValue={}".format(config.GRIB_max_bitspervalue)
                    epylog.warning(msg)
                    packing['bitsPerValue'] = config.GRIB_max_bitspervalue
        # FIXME: ? pre-set values before packing seems necessary if packing different from sample
        if any([self[k] != v for k, v in packing.items()]):
            self._GRIB2_set_section6_7(field, **other_GRIB_options)
        if packing.get('packingType') == 'grid_simple':
            self._untouch.add('dataRepresentationTemplateNumber')
            self._untouch.add('bitsPerValue')
        for k in order:
            if k in packing:
                try:
                    self[k] = packing[k]
                except lowlevelgrib.InternalError:
                    if config.GRIB_packing_fatal:
                        raise
                else:
                    self._untouch.add(k)

    def _GRIB2_set_section6_7(self, field, **other_GRIB_options):
        """Set k/v for sections 6 and 7."""
        # bitmap
        missingValue = other_GRIB_options.get('missingValue', None)
        values = field.getdata(d4=True).copy()
        if isinstance(values, numpy.ma.masked_array) and values.mask.any():
            self['bitMapIndicator'] = 0
            self['bitmapPresent'] = 1
            self['missingValue'] = values.fill_value
            values = field.geometry.fill_maskedvalues(values, missingValue)
        values = values.squeeze()
        # values
        self.set_values(values)
        self._untouch.add('missingValue')

    # Write GRIB1 (old way) ----------------------------------------------------

    def _GRIB1_set(self, field,
                   ordering=griberies.defaults.GRIB1_ordering,
                   packing=None,
                   sample=None,
                   other_GRIB_options={}):
        """
        Build the GRIB message from the field, using a template.

        :param field: a :class:`epygram.base.Field`.
        :param ordering: flattening of 2D data
        :param packing: options of packing and compression in GRIB (dict).
        :param sample: to use a specific sample GRIB.
          Specific syntax 'file:$filename$' takes the first message
          in $filename$ as sample.
        :param other_GRIB_options: other options to be specified in GRIB,
          as a dict(GRIBkey=value)
        """
        if not isinstance(field, H2DField):
            raise NotImplementedError("not yet.")
        if other_GRIB_options is None:
            other_GRIB_options = {}

        # part 1 --- set sample and preset packing
        if packing is None:
            if field.spectral:
                packing = {'packingType':'spectral_simple',
                           'bitsPerValue':24}
            else:
                packing = griberies.defaults.GRIB1_packing
        # get packingType...
        required_packingType = packing.get('packingType', None)
        if required_packingType is not None:
            if 'grid' in required_packingType and field.spectral or\
               'spectral' in required_packingType and not field.spectral:
                required_packingType = None
        # ... to determine sample
        if sample is None:
            # try to get an appropriate sample, or a default one
            sample = self._GRIB1_default_sample(field, required_packingType)
        # reset packing "on the fly" if field is uniform
        if field.max() - field.min() < config.epsilon:
            packing = {'packingType':'grid_simple'}
            required_packingType = packing['packingType']
            sample = self._GRIB1_default_sample(field, required_packingType)
        # clone from sample
        if sample.startswith('file:'):
            self._gid = self._clone_from_file(sample.replace('file:', ''))
        else:
            if all([sample + '.tmpl' not in os.listdir(pth)
                    for pth in [p for p in griberies.paths.get_samples_paths()
                                if (os.path.exists(p) and os.path.isdir(p))]]):  # this sample is not found in samples path
                if sample + '.tmpl' in os.listdir(config.GRIB_epygram_samples_path):  # but it is in epygram samples
                    fsample = os.path.join(config.GRIB_epygram_samples_path, sample) + '.tmpl'
                    self._gid = self._clone_from_file(fsample)  # clone it from file
                else:
                    raise epygramError('sample {} not found; check {} environment variable'.
                                       format(sample, lowlevelgrib.samples_path_var))
            else:
                self._gid = self._clone_from_sample(sample)

        # part 2 --- parameter
        self._GRIB1_set_fid(field)

        # part 3 --- context
        try:
            process_id = int(field.processtype)
        except ValueError:
            process_id = griberies.defaults.GRIB1_keyvalue['generatingProcessIdentifier']
        self['generatingProcessIdentifier'] = process_id

        # part 4 --- validity
        self._GRIB1_set_validity(field)

        # part 5 --- geometry
        self._GRIB1_set_geometry(field)

        # part 6 --- ordering
        if not field.spectral:
            self.update(ordering)

        # part 7 --- other options
        if len(other_GRIB_options) > 0:
            try:
                for k, v in other_GRIB_options.items():
                    self[k] = v
            except lowlevelgrib.InternalError as e:
                epylog.warning('set ' + k + ' failed: ' + str(e))

        # part 8 --- set first gridpoint according to ordering
        if not field.spectral:
            self._GRIB1_set_data_ordering(field)

        # part 9 --- values
        values = field.getdata(d4=True).copy()
        if isinstance(values, numpy.ma.masked_array) and values.mask.any():
            # bitmap in GRIB1 ?
            self['bitmapPresent'] = 1
            values = field.geometry.fill_maskedvalues(values)
            # TODO: ok now ? raise NotImplementedError("didn't succeed to make this work")
        values = values.squeeze()
        if not field.spectral:
            self.set_values(values)
        self._GRIB1_set_packing(packing)
        self.set_values(values)

    def _GRIB1_set_fid(self, field):
        """Set GRIB1 fid k/v."""
        param_list = ['table2Version', 'indicatorOfParameter']
        for k in param_list:
            self[k] = field.fid['GRIB1'][k]

    def _GRIB1_set_validity(self, field):
        """Set validity k/v."""
        # cutoff
        self['dataDate'] = int(field.validity.getbasis(fmt='IntStr')[:8])
        self['dataTime'] = int(field.validity.getbasis(fmt='IntStr')[8:12])
        # term and cumulative duration
        term_in_seconds = field.validity.term().total_seconds()
        cumulative = field.validity.cumulativeduration()
        if cumulative:
            startStep_in_seconds = (field.validity.term() - cumulative).total_seconds()
            self['timeRangeIndicator'] = 4
        else:
            self['timeRangeIndicator'] = 0
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

    def _GRIB1_set_geometry(self, field):
        """Set geometry k/v."""
        if not field.spectral:
            # earth shape ------------------------------------------------------
            if not hasattr(field.geometry, 'geoid'):
                epylog.warning("geometry has no geoid: assume 'shapeOfTheEarth' = 6.")
                self['shapeOfTheEarth'] = 6
            else:
                if all([field.geometry.geoid.get(axis) == griberies.tables.pyproj_geoid_shapes[6][axis]
                        for axis in ('a', 'b')]):
                    self['shapeOfTheEarth'] = 6
                else:
                    found = False
                    for s, g in griberies.tables.pyproj_geoid_shapes.items():
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
            # dimensions -------------------------------------------------------
            if field.geometry.rectangular_grid:
                self['Ni'] = field.geometry.dimensions['X']
                self['Nj'] = field.geometry.dimensions['Y']
            else:
                pass  # done below
            # type of geometry  ------------------------------------------------
            if field.geometry.name == 'regular_lonlat':
                self['gridType'] = 'regular_ll'
                self['iDirectionIncrementInDegrees'] = field.geometry.grid['X_resolution'].get('degrees')
                self['jDirectionIncrementInDegrees'] = field.geometry.grid['Y_resolution'].get('degrees')
            elif field.geometry.name == 'rotated_lonlat':
                self['gridType'] = 'rotated_ll'
                self['iDirectionIncrementInDegrees'] = field.geometry.grid['X_resolution'].get('degrees')
                self['jDirectionIncrementInDegrees'] = field.geometry.grid['Y_resolution'].get('degrees')
                self['longitudeOfSouthernPoleInDegrees'] = field.geometry.grid['southern_pole_lon'].get('degrees')
                self['latitudeOfSouthernPoleInDegrees'] = field.geometry.grid['southern_pole_lat'].get('degrees')
                self['angleOfRotation'] = field.geometry.grid['rotation'].get('degrees')
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
                    self['LaD'] = int(lat_0 * 1e6)  # LaDInDegrees sometimes caused errors
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
                    except lowlevelgrib.InternalError:
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
                    except lowlevelgrib.InternalError:
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
            # spectral case  ---------------------------------------------------
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

        # vertical geometry ----------------------------------------------------
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

    def _GRIB1_set_data_ordering(self, field):
        if field.geometry.rectangular_grid:
            corners = field.geometry.gimme_corners_ll()
            if field.geometry.name == 'rotated_lonlat':  # special case: coordinates to be given in the rotated referential
                for k, v in corners.items():
                    corners[k] = field.geometry.ll2xy(*v)
            if self['iScansNegatively'] == 0 and \
               self['jScansPositively'] == 0:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ul'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['lr'][1]
            elif self['iScansNegatively'] == 0 and \
                 self['jScansPositively'] == 1:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ll'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ur'][1]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 0:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['ur'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['ur'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ll'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ll'][1]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 1:
                self['longitudeOfFirstGridPointInDegrees'] = util.positive_longitude(corners['lr'][0])
                self['latitudeOfFirstGridPointInDegrees'] = corners['lr'][1]
                if self['gridType'] == 'regular_ll':
                    self['longitudeOfLastGridPointInDegrees'] = util.positive_longitude(corners['ul'][0])
                    self['latitudeOfLastGridPointInDegrees'] = corners['ul'][1]
            else:
                raise NotImplementedError('this ordering: not yet.')
        else:
            if not(self['iScansNegatively'] == 0 and
                   self['jScansPositively'] == 0):
                raise NotImplementedError('this ordering: not yet.')

    def _GRIB1_set_packing(self, packing):
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
                except lowlevelgrib.InternalError:
                    if config.GRIB_packing_fatal:
                        raise
        for k, v in packing.items():  # remaining items
            self[k] = v

    # Read ---------------------------------------------------------------------

    def as_field(self,
                 getdata=True,
                 read_misc_metadata=griberies.defaults.GRIB2_metadata_to_embark):
        """
        Make an epygram H2DField from message.

        :param getdata: if *False*, only metadata are read, the field do not
            contain data.
        :param read_misc_metadata: read the specified keys, and store it in
            field.misc_metadata
        """
        field_kwargs = dict(
            fid=self.genfid(),
            validity=self._read_validity(),
            geometry=self._read_geometry(),
            misc_metadata=self._read_misc(read_misc_metadata))
        # additional
        field_kwargs['structure'] = field_kwargs['geometry'].structure
        field_kwargs['processtype'] = self['generatingProcessIdentifier']
        if self['gridType'] == 'sh':
            field_kwargs['spectral_geometry'] = self._read_spectralgeometry()
        # build field
        field = H2DField(**field_kwargs)
        # getdata
        if getdata:
            field.setdata(self._read_data(field.geometry))
        return field

    def asfield(self,
                getdata=True,
                footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                get_info_as_json=None):
        """
        .. deprecated:: 1.3.9

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
        except lowlevelgrib.InternalError:
            try:
                field_kwargs['units'] = self['parameterUnits']
            except lowlevelgrib.InternalError:
                pass
        if get_info_as_json is not None:
            comment = {}
            for k in get_info_as_json:
                try:
                    comment[k] = self[k]
                except lowlevelgrib.InternalError:
                    pass
            field_kwargs['comment'] = json.dumps(comment)
        field = builder(**field_kwargs)
        if getdata:
            field.setdata(self._read_data(field.geometry))
        return field

    def _read_earth_geometry(self):
        """Read info about geoid."""
        geoid = config.default_geoid
        try:
            geoid = griberies.tables.pyproj_geoid_shapes[self['shapeOfTheEarth']]
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
        except lowlevelgrib.InternalError:
            pass
        finally:
            if self['centreDescription'] == 'Oslo':
                try:
                    (a, b, _) = self['pv']
                    geoid = {'a':a, 'b':b}
                except KeyError:
                    pass
        return geoid

    def _read_misc(self, read_misc_metadata):
        """Read specified additional keys."""
        misc_metadata = FPDict({})
        for k in read_misc_metadata:
            v = self.get(k)
            if v is not None:
                misc_metadata[k] = v
        return misc_metadata

    def _read_validity(self):
        """
        Returns a :class:`epygram.base.FieldValidity` object containing the
        validity of the GRIB message.
        """
        if self.grib_edition == 2:
            return self._GRIB2_read_validity()
        else:
            return self._GRIB1_read_validity()

    def _read_geometry(self):
        """
        Returns a geometry object containing
        the geometry information of the GRIB message.
        """
        geoid = self._read_earth_geometry()
        vcoordinate = self._read_vertical_geometry()
        geometryname, geometryclass, dimensions, grid, projection = self._read_horizontal_geometry()
        # Make geometry object
        kwargs_geom = dict(name=geometryname,
                           grid=grid,
                           dimensions=dimensions,
                           vcoordinate=vcoordinate,
                           position_on_horizontal_grid='center',
                           geoid=geoid)
        if projection is not None:
            kwargs_geom['projection'] = projection
        geometry = geometryclass(**kwargs_geom)
        return geometry

    def _read_vertical_geometry(self):
        """Read vertical part of geometry."""
        if self.grib_edition == 2:
            return self._GRIB2_read_vertical_geometry()
        else:
            return self._GRIB1_read_vertical_geometry()

    def _read_hybridP_levels(self):
        """Read A and B coefficients of HybridPressure coordinate."""
        try:
            self._readattribute('pv', array=True)
            A_and_B = self['pv']
        except lowlevelgrib.InternalError:
            epylog.warning('Error while reading A/B vertical levels coefficients ! Ignore.')
            A_and_B = []
        Ai = A_and_B[:len(A_and_B) // 2]
        Bi = A_and_B[len(A_and_B) // 2:]
        return {'gridlevels': tuple([(i + 1, FPDict({'Ai':Ai[i], 'Bi':Bi[i]}))
                                     for i in range(len(Ai))]),
                'ABgrid_position':'flux'}

    def _read_horizontal_geometry(self):
        """Read horizontal geometry."""
        # rectangular grid dimensions
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
        # ----------------------------------------------------------------------
        if self['gridType'] == 'regular_ll':
            geometryclass = RegLLGeometry
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
        # ----------------------------------------------------------------------
        elif self['gridType'] == 'rotated_ll':
            geometryclass = RotLLGeometry
            geometryname = 'rotated_lonlat'
            grid = {'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'X_resolution':Angle(self['iDirectionIncrementInDegrees'], 'degrees'),
                    'Y_resolution':Angle(self['jDirectionIncrementInDegrees'], 'degrees'),
                    'southern_pole_lon':Angle(self['longitudeOfSouthernPoleInDegrees'], 'degrees'),
                    'southern_pole_lat':Angle(self['latitudeOfSouthernPoleInDegrees'], 'degrees'),
                    'rotation':Angle(self['angleOfRotation'], 'degrees')
                    }
            projection = None
        # ----------------------------------------------------------------------
        elif self['gridType'] in ('polar_stereographic',):
            geometryclass = ProjectedGeometry
            geometryname = self['gridType']
            if self['projectionCentreFlag'] == 0:
                lat_0 = 90.
            else:
                lat_0 = -90.
            # In GRIB1, LaD is not present but resolution is supposed to
            # be given at +/-60deg
            lat_ts = numpy.copysign(self.get('LaDInDegrees', 60.), lat_0)
            # !!! dirty bypass of erroneous GRIBs from OSI-SAF !..
            if self['centreDescription'] == 'Oslo':
                epylog.warning(' '.join(['for centre:',
                                         self['centreDescription'],
                                         "(OSI-SAF)"
                                         "and polar stereographic" +
                                         "projection, lat_ts is taken in" +
                                         "3rd position in attribute 'pv[]'" +
                                         "of GRIB message !"]))
                lat_ts = self['pv'][3]
            projection = {'reference_lon':Angle(float(self['orientationOfTheGridInDegrees']), 'degrees'),
                          'reference_lat':Angle(lat_0, 'degrees'),
                          'secant_lat':Angle(lat_ts, 'degrees'),
                          'rotation':Angle(0., 'degrees')
                          }
            m = 1
            grid = {'X_resolution':float(self['DxInMetres']) * m,
                    'Y_resolution':float(self['DyInMetres']) * m,
                    'input_lon':Angle(self['longitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_lat':Angle(self['latitudeOfFirstGridPointInDegrees'], 'degrees'),
                    'input_position':input_position,
                    'LAMzone':None}
        # ----------------------------------------------------------------------
        elif self['gridType'] in ('lambert', 'lambert_lam'):
            geometryname = 'lambert'  # self['gridType'] # account for lambert_lam
            geometryclass = ProjectedGeometry
            lon_0 = self['LoVInDegrees']
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
        # ----------------------------------------------------------------------
        elif self['gridType'] in ('mercator',):
            geometryclass = ProjectedGeometry
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
        # ----------------------------------------------------------------------
        elif 'gg' in self['gridType']:
            geometryclass = GaussGeometry
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
        # ----------------------------------------------------------------------
        elif self['gridType'] == 'sh':
            # spherical harmonics: => forced to a linear gauss grid
            projection = None
            spgeom = self._read_spectralgeometry()
            gpdims = gridpoint_dims_from_truncation(spgeom.truncation,
                                                    grid='linear',
                                                    stretching_coef=1.)
            latitudes = gauss_latitudes(gpdims['lat_number'])
            grid = {'latitudes':FPList([Angle(l, 'degrees')
                                        for l in latitudes]),
                    'dilatation_coef':1.}
            dimensions = gpdims
            geometryclass = GaussGeometry
            geometryname = 'reduced_gauss'
            # try to have roughly the same zonal resolution as on equator
            lon_number_by_lat = 2 * gpdims['lat_number'] * numpy.cos(numpy.radians(latitudes))
            lon_number_by_lat = [min(nearest_greater_FFT992compliant_int(n),
                                     dimensions['max_lon_number'])
                                 for n in lon_number_by_lat]
            dimensions['lon_number_by_lat'] = FPList(lon_number_by_lat)
        else:
            raise NotImplementedError("gridType == {} : not yet !".
                                      format(self['gridType']))
        return geometryname, geometryclass, dimensions, grid, projection

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

    def _read_data(self, geometry):
        """Read data then reorder/reshape."""
        data1d = self['values']
        if 'gauss' in geometry.name:
            if self['gridType'] == 'sh':
                data = data1d
            else:
                if self['iScansNegatively'] == 0 and \
                   self['jScansPositively'] == 0 and \
                   self['jPointsAreConsecutive'] == 0:
                    data = geometry.reshape_data(data1d[:])
                else:
                    raise NotImplementedError("not yet !")
        else:
            if self['iScansNegatively'] == 0 and \
               self['jScansPositively'] == 0 and \
               self['jPointsAreConsecutive'] == 0:
                data = data1d.reshape((self['Nj'], self['Ni']),
                                        order='C')[::-1, :]
            elif self['iScansNegatively'] == 0 and \
                 self['jScansPositively'] == 1 and \
                 self['jPointsAreConsecutive'] == 0:
                data = data1d.reshape((self['Nj'], self['Ni']),
                                        order='C')[:, :]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 0 and \
                 self['jPointsAreConsecutive'] == 0:
                data = data1d.reshape((self['Nj'], self['Ni']),
                                        order='C')[::-1, ::-1]
            elif self['iScansNegatively'] == 1 and \
                 self['jScansPositively'] == 1 and \
                 self['jPointsAreConsecutive'] == 0:
                data = data1d.reshape((self['Nj'], self['Ni']),
                                        order='C')[:, ::-1]
            else:
                raise NotImplementedError("not yet !")

        if ((self['editionNumber'] == 2 and
             self['bitMapIndicator'] == 0) or
            (self['editionNumber'] == 1 and
             self['bitmapPresent'] == 1)):
            data = numpy.ma.masked_equal(data, self['missingValue'])
        return data

    def _GRIB2_read_validity(self):
        """Read validity."""
        try:
            # defaults
            cum = None
            time_increment = None
            statistical_process = None
            # basis
            basis = datetime.datetime(
                self['year'], self['month'], self['day'],
                self['hour'], self['minute'], self['second'])
            # term and accumulation
            units = {0:dict(seconds=60),  # minute
                     1:dict(seconds=3600),  # hour
                     2:dict(days=1),  # day
                     13:dict(seconds=1),  # seconds in GRIB2
                     }
            unit = units[self['indicatorOfUnitOfTimeRange']]
            forecastTime_unit = {k:v * self['forecastTime'] for k,v in unit.items()}
            forecastTime = datetime.timedelta(**forecastTime_unit)
            if self['productDefinitionTemplateNumber'] in (8, 11):
                # case of an accumulated field
                term = datetime.datetime(
                    self['yearOfEndOfOverallTimeInterval'],
                    self['monthOfEndOfOverallTimeInterval'],
                    self['dayOfEndOfOverallTimeInterval'],
                    self['hourOfEndOfOverallTimeInterval'],
                    self['minuteOfEndOfOverallTimeInterval'],
                    self['secondOfEndOfOverallTimeInterval'])
                beginning = basis + forecastTime
                cum = term - beginning
                term = term - basis
                statistical_process = self['typeOfStatisticalProcessing']
                time_increment_unit = {k:v * self['timeIncrement'] for k,v in unit.items()}
                time_increment = datetime.timedelta(**time_increment_unit)
            else:
                term = forecastTime
        except (epygramError, NotImplementedError):
            # failed
            if config.GRIB_ignore_validity_decoding_errors:
                validity = FieldValidity()
            else:
                raise
        else:
            # success
            validity = FieldValidity(basis=basis,
                                     term=term,
                                     cumulativeduration=cum,
                                     statistical_process_on_duration=statistical_process,
                                     statistical_time_increment=time_increment)
        finally:
            return validity

    def _GRIB1_read_validity(self):
        """
        Returns a :class:`epygram.base.FieldValidity` object containing the
        validity of the GRIB message.
        """
        accepted_tRI = (0, 1, 2, 3, 4, 10, 113, 123)
        basis = datetime.datetime(
            self['year'], self['month'], self['day'],
            self['hour'], self['minute'], self['second'])
        cum = None
        units = {0:dict(seconds=60),  # minute
                 1:dict(seconds=3600),  # hour
                 2:dict(days=1),  # day
                 13:dict(seconds=1 if self.grib_edition == 2 else 15 * 60),  # seconds in GRIB2, 1/4h in GRIB1
                 }
        unit = units[self['stepUnits']]
        if self['timeRangeIndicator'] in accepted_tRI:
            term = self['endStep']
            termunit = {k:v * term for k,v in unit.items()}
            term = datetime.timedelta(**termunit)
            if self['timeRangeIndicator'] in (2, 3, 4) or self['productDefinitionTemplateNumber'] in (8, 11):
                cum = self['endStep'] - self['startStep']
                cumunit = {k:v * cum for k,v in unit.items()}
                cum = datetime.timedelta(**cumunit)
            elif self['timeRangeIndicator'] in (113, 123,):
                epylog.warning('not able to interpret timeRangeIndicator={}'.format(self['timeRangeIndicator']))
        else:
            if config.GRIB_ignore_validity_decoding_errors:
                try:
                    term = self['endStep']
                    termunit = {k:v * term for k,v in unit.items()}
                    term = datetime.timedelta(**termunit)
                except Exception:
                    term = None
                epylog.warning('not able to interpret timeRangeIndicator={}'.format(self['timeRangeIndicator']))
            else:
                raise NotImplementedError("'timeRangeIndicator' not in {}.".format(accepted_tRI))
        validity = FieldValidity(basis=basis, term=term, cumulativeduration=cum)
        return validity

    def _GRIB2_read_vertical_geometry(self):
        """Read vertical part of geometry in GRIB2."""
        kwargs_vcoord = dict(position_on_grid='mass')
        kwargs_vcoord['typeoffirstfixedsurface'] = self.get('typeOfFirstFixedSurface', 255)
        if kwargs_vcoord['typeoffirstfixedsurface'] != 255:
            first_level = self['scaledValueOfFirstFixedSurface'] * 10**(-self['scaleFactorOfFirstFixedSurface'])
            kwargs_vcoord['levels'] = [first_level]
        else:
            kwargs_vcoord['levels'] = [0]
        if self.get('typeOfSecondFixedSurface', 255) != 255:
            kwargs_vcoord['typeofsecondfixedsurface'] = self['typeOfSecondFixedSurface']
            second_level = self['scaledValueOfSecondFixedSurface'] * 10**(-self['scaleFactorOfSecondFixedSurface'])
            kwargs_vcoord['toplevel'] = first_level
            kwargs_vcoord['bottomlevel'] = second_level
        if kwargs_vcoord['typeoffirstfixedsurface'] == 119:
            kwargs_vcoord['grid'] = self._read_hybridP_levels()
        return VGeometry(**kwargs_vcoord)

    def _GRIB1_read_vertical_geometry(self):
        """Read vertical part of geometry in GRIB1."""
        kwargs_vcoord = {}
        kwargs_vcoord['position_on_grid'] = 'mass'
        kwargs_vcoord['typeoffirstfixedsurface'] = onetotwo.get(self['indicatorOfTypeOfLevel'], 255)
        kwargs_vcoord['levels'] = [self['level']]
        if self['indicatorOfTypeOfLevel'] in (112,):
            kwargs_vcoord['toplevel'] = self['topLevel']
            kwargs_vcoord['bottomlevel'] = self['bottomLevel']
        if kwargs_vcoord['typeoffirstfixedsurface'] == 119:
            kwargs_vcoord['grid'] = self._read_hybridP_levels()
        return VGeometry(**kwargs_vcoord)

    # Utilities methods --------------------------------------------------------

    def set_from(self, k, series_of_dicts, fatal=True):
        """Set key k from its value first found in one of **series_of_dicts**."""
        for d in series_of_dicts:
            if k in d:
                try:
                    self[k] = d[k]
                    if hasattr(self, '_untouch'):
                        self._untouch.add(k)
                except Exception:
                    if fatal:
                        raise
                    else:
                        epylog.warning('Failed to set {}={}'.format(k, d[k]))
                break

    def set_values(self, values):
        """
        Wrapper to set **values** as a 1D if spectral, or 2D array if gridpoint,
        with consistency with ordering parameters already set beforehand.
        """
        if len(values.shape) == 1:
            lowlevelgrib.set_values(self._gid, values)
        elif len(values.shape) == 2:
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
            lowlevelgrib.set_values(self._gid, data1d)

    def actual_fid_keys(self):
        """Adjust fid keys to specific grib edition and templates."""
        template = None
        scaleFactorOfFirstFixedSurface = 0
        scaleFactorOfSecondFixedSurface = 0
        if self.grib_edition == 2:
            template = self['productDefinitionTemplateNumber']
            if template not in (32, 33):
                scaleFactorOfFirstFixedSurface = self['scaleFactorOfFirstFixedSurface']
                if self['typeOfSecondFixedSurface'] != 255:
                    scaleFactorOfSecondFixedSurface = self['scaleFactorOfSecondFixedSurface']
        return self.fid_keys_for(self.grib_edition,
                                 productDefinitionTemplateNumber=template,
                                 scaleFactorOfFirstFixedSurface=scaleFactorOfFirstFixedSurface,
                                 scaleFactorOfSecondFixedSurface=scaleFactorOfSecondFixedSurface)

    def update(self, E, **F):
        """
        D.update([E, ]**F) -> None. Update D from dict/iterable E and F.
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
        fid_keys = copy.copy(self.actual_fid_keys())
        try:
            onelevel = self['topLevel'] == self['bottomLevel']
        except lowlevelgrib.InternalError:
            onelevel = True
        if onelevel:
            if 'topLevel' in fid_keys:
                fid_keys.remove('topLevel')
            if 'bottomLevel' in fid_keys:
                fid_keys.remove('bottomLevel')
        fid = {self.fmt:FPDict({k:self[str(k)] for k in fid_keys})}
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
        if isinstance(namespace, str):
            namespace = str(namespace)  # gribapi str/unicode incompatibility
        key_iter = lowlevelgrib.keys_iterator_new(self._gid,
                                                  namespace=namespace)
        namespace = []
        while lowlevelgrib.keys_iterator_next(key_iter):
            namespace.append(lowlevelgrib.keys_iterator_get_name(key_iter))
        lowlevelgrib.keys_iterator_delete(key_iter)
        return namespace

    def readmessage(self, namespace=None, fatal=True):
        """
        Reads the meta-data of the message.

        :param namespace: the namespace of keys to be read, among:\n
          - **None**: to get all keys present in message,
          - ['myKey1', 'myKey2', ...] for any custom namespace,
          - 'ls': to get the same default keys as the grib_ls,
          - 'mars': to get the keys used by MARS.
        :param fatal: raise errors or ignore failed keys
        """
        self.clear()
        if namespace in (None, 'ls', 'mars'):
            namespace = self.readkeys(namespace=namespace)
        for k in namespace:
            self._readattribute(k, fatal=fatal)

    @property
    def grib_edition(self):
        return self['editionNumber']

    @property
    def fmt(self):
        return 'GRIB' + str(self.grib_edition)


# File containing a Collection of GRIB messages --------------------------------

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
        self._open_through = str(self.container.abspath)  # gribapi str/unicode incompatibility
        if self.openmode in ('r', 'a'):
            _file = io.open(self.container.abspath, 'rb')
            isgrib = _file.readline()[:4]
            try:
                isgrib = isgrib.decode('utf-8')
            except UnicodeDecodeError:
                raise IOError("cannot decode this resource as GRIB.")
            _file.close()
            if isgrib != 'GRIB':
                raise IOError("this resource is not a GRIB one.")
        if not self.fmtdelayedopen:
            self.open()

    def __iter__(self):
        self.close()
        self.open()
        for _ in range(len(self)):
            yield self.iter_fields(getdata=True)

    def open(self, openmode=None):
        """
        Opens a GRIB and initializes some attributes.

        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(GRIB, self).open(openmode=openmode)
        if self.openmode != 'w' and config.GRIB_safe_indexes:
            # grib index bug workaround # FIXME: well not me, gribapi:
            # find an available AND unique filename
            fd, fn = tempfile.mkstemp(dir=config.GRIB_safe_indexes,
                                      suffix=str(uuid.uuid4()))
            os.close(fd)
            self._open_through = str(fn)
            os.remove(self._open_through)
            os.symlink(self.container.abspath, self._open_through)
        self._file = open(self._open_through, self.openmode + 'b')
        self.isopen = True

    def close(self):
        """
        Closes a GRIB.
        """
        if hasattr(self, '_file'):
            self._file.close()
        if config.GRIB_safe_indexes and \
           hasattr(self, '_open_through') and \
           os.path.exists(self._open_through) and \
           self._open_through != self.container.abspath:  # ceinture et bretelles
                os.unlink(self._open_through)
                self._open_through = str(self.container.abspath)  # keep a valid filename
        if hasattr(self, '_sequential_file'):
            if hasattr(self._sequential_file, 'closed') and \
               not self._sequential_file.closed:
                self._sequential_file.close()
        self.isopen = False

    @property
    @FileResource._openbeforedelayed
    def messages_number(self):
        """Counts the number of messages in file."""
        return lowlevelgrib.count_in_file(self._file)

    def listfields(self, onlykey=None, select=None, complete=False,
                   additional_keys=[]):
        """
        Returns a list containing the GRIB identifiers of all the fields of the
        resource.

        :param onlykey: can be specified as a string or a tuple of strings,
          so that only specified keys of the fid will returned.
        :param select: can be specified as a dict(key=value) to restrain
          the list of fields to those that match the key:value pairs.
        :param complete: list fields with their natural fid + generic one.
        :param additional_keys: add given keys in fids
        """
        if select is not None:
            select_keys = list(select.keys())
        else:
            select_keys = []
        fidlist = super(GRIB, self).listfields(additional_keys=additional_keys,
                                               select_keys=select_keys)
        if select is not None:
            fidlist = [f for f in fidlist if all([f[k] == select[k] for k in select.keys()])]
        if onlykey is not None:
            if isinstance(onlykey, str):
                fidlist = [f[onlykey] for f in fidlist]
            elif isinstance(onlykey, tuple):
                fidlist = [{k:f[k] for k in onlykey} for f in fidlist]
        if complete:
            fidlist = [{'GRIB' + str(f['editionNumber']):f} for f in fidlist]
            for f in fidlist:
                if 'GRIB2' in f:
                    f['generic'] = f['GRIB2']
        return fidlist

    @FileResource._openbeforedelayed
    def _listfields(self, additional_keys=[], select_keys=[]):
        """Returns a list of GRIB-type fid of the fields inside the resource."""
        def gid_key_to_fid(gid, k, fid):
            # bug in GRIB_API ? 1, 103 & 105 => 'sfc'
            if k in ('typeOfFirstFixedSurface',
                     'indicatorOfTypeOfLevel',
                     'typeOfSecondFixedSurface',
                     # type error when setting key 'centre' if str
                     'centre', 'originatingCentre'):
                fid[k] = lowlevelgrib.get(gid, k, int)
            elif k in ('topLevel', 'bottomLevel'):
                if lowlevelgrib.get(gid, 'topLevel') != lowlevelgrib.get(gid, 'bottomLevel'):
                    fid[k] = lowlevelgrib.get(gid, k)
            else:
                fid[k] = lowlevelgrib.get(gid, k)

        fidlist = []
        _file = open(self._open_through, 'rb')
        while True:
            filter_out = False
            fid = {}
            gid = lowlevelgrib.any_new_from_file(_file, headers_only=True)
            if gid is None:
                break
            n = lowlevelgrib.get(gid, 'editionNumber')
            particularities = {}
            if n == 2:
                t = lowlevelgrib.get(gid, 'productDefinitionTemplateNumber')
                if t not in (32, 33):
                    particularities['scaleFactorOfFirstFixedSurface'] = lowlevelgrib.get(gid, 'scaleFactorOfFirstFixedSurface')
                    if lowlevelgrib.get(gid, 'typeOfSecondFixedSurface') != 255:
                        particularities['scaleFactorOfSecondFixedSurface'] = lowlevelgrib.get(gid, 'scaleFactorOfSecondFixedSurface')
            else:
                t = None
            for k in GRIBmessage.fid_keys_for(n, t, **particularities):
                gid_key_to_fid(gid, k, fid)
            for k in additional_keys:  # FIXME: here we mix additional, if present, and select/filter
                try:
                    gid_key_to_fid(gid, k, fid)
                except lowlevelgrib.InternalError:
                    pass
            for k in select_keys:  # FIXME: here we mix additional, if present, and select/filter
                try:
                    gid_key_to_fid(gid, k, fid)
                except lowlevelgrib.InternalError:
                    filter_out = True
                    break
            lowlevelgrib.release(gid)
            if not filter_out:
                fidlist.append(fid)
        _file.close()
        return fidlist

    def split_UV(self, fieldseed):
        """
        Return two lists of fids corresponding respectively to U and V
        components of wind, given a *fieldseed*.
        Syntax example: 'shortName':'u+v', or 'indicatorOfParameter':'33+34'
        """
        if isinstance(fieldseed, str):
            fieldseed = griberies.parse_GRIBstr_todict(fieldseed)
        seeds = [fieldseed.copy(), fieldseed.copy()]
        for i in (0, 1):
            for k, v in seeds[i].items():
                if isinstance(v, str) and '+' in v:
                    v = v.split('+')[i]
                    try:
                        seeds[i][k] = int(v)
                    except ValueError:
                        seeds[i][k] = v
        Ufid = self.find_fields_in_resource(seeds[0])
        Vfid = self.find_fields_in_resource(seeds[1])
        return (sorted(Ufid), sorted(Vfid))

    @FileResource._openbeforedelayed
    def get_message_at_position(self, position):
        """
        Returns the message at position *position*, from 0 (first message)
        to messages_number-1 (last message).
        Negative position starts from the end: -1 is last message.

        Should not be used sequentially, probably very inefficient.
        """
        return GRIBmessage(('file', self._open_through, position))

    @FileResource._openbeforedelayed
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
            self._sequential_file = open(self._open_through, 'rb')

        gid = lowlevelgrib.any_new_from_file(self._sequential_file,
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
        msg = self.iter_messages(headers_only=False)
        if msg is not None:
            if 'footprints_proxy_as_builder' in kwargs or 'get_info_as_json' in kwargs:
                fld = msg.asfield(getdata=getdata, **kwargs)  # CLEANME: deprecated
            else:
                fld = msg.as_field(getdata=getdata, **kwargs)
        else:
            fld = None
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
                if isinstance(s, str):
                    s = griberies.parse_GRIBstr_todict(s)
                fieldslist.extend(self.listfields(select=s))
        elif isinstance(seed, str):
            fieldslist = self.listfields(select=griberies.parse_GRIBstr_todict(seed))
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
                  get_info_as_json=None,
                  read_misc_metadata=griberies.defaults.GRIB2_metadata_to_embark):
        """
        Finds in GRIB the message that correspond to the *handgrip*,
        and returns it as a :class:`epygram.base.Field`.
        If several messages meet the requirements, raises error (use
        readfields() method instead).

        :param dict handgrip: a dict where you can store all requested GRIB
            keys for discrimination...
            E.g. {'shortName':'t', 'typeOfFirstFixedSurface':100, 'level':850}
            will return the Temperature at 850hPa field.
        :param getdata: if False, the data is not read, the field consist
            in the meta-data only.
        :param footprints_proxy_as_builder: if True, uses footprints.proxy
            to build fields. True decreases performance.
            .. deprecated:: 1.3.9
        :param get_info_as_json: if not None, writes the keys given in
            *get_info_as_json* as json in field.comment.
            .. deprecated:: 1.3.9
        :param read_misc_metadata: read the specified keys, and store it in
            field.misc_metadata
        """
        if isinstance(handgrip, str):
            handgrip = griberies.parse_GRIBstr_todict(handgrip)
        matchingfields = self.readfields(handgrip,
                                         getdata=getdata,
                                         footprints_proxy_as_builder=footprints_proxy_as_builder,
                                         get_info_as_json=get_info_as_json,
                                         read_misc_metadata=read_misc_metadata)
        # filter out those unfiltered by the below GRIBAPI bug
        filtered_matchingfields = FieldSet()
        for field in matchingfields:
            fid = field.fid.get('GRIB1', field.fid.get('GRIB2'))
            if all([fid.get(k, None) in (v, None) for k, v in handgrip.items()]):
                filtered_matchingfields.append(field)
        if len(filtered_matchingfields) > 1:
            raise epygramError("several fields found for that *handgrip*;" +
                               " please refine:" + str(handgrip))
        elif len(filtered_matchingfields) == 0:
            raise epygramError("inconsistency in *handgrip*; check again" +
                               " values and types of values")
        else:
            filtered_matchingfields[0].fid['short'] = FPDict(handgrip)
        return filtered_matchingfields[0]

    @FileResource._openbeforedelayed
    def readfields(self, handgrip,
                   getdata=True,
                   footprints_proxy_as_builder=config.footprints_proxy_as_builder,
                   get_info_as_json=None,
                   read_misc_metadata=griberies.defaults.GRIB2_metadata_to_embark):
        """
        Finds in GRIB the message(s) that correspond to the *handgrip*,
        and returns it as a :class:`epygram.base.FieldSet` of
        :class:`epygram.base.Field`.

        :param dict handgrip: a dict where you can store all requested GRIB
            keys for discrimination...
            E.g. {'shortName':'t', 'typeOfFirstFixedSurface':100}
            will return all the Temperature fields on Pressure levels.
        :param getdata: if False, the data is not read, the field consist
            in the meta-data only.
        :param footprints_proxy_as_builder: if True, uses footprints.proxy
            to build fields. True decreases performance.
            .. deprecated:: 1.3.9
        :param get_info_as_json: if not None, writes the keys given in
            *get_info_as_json* as json in field.comment.
            .. deprecated:: 1.3.9
        :param read_misc_metadata: read the specified keys, and store it in
            field.misc_metadata
        """
        if isinstance(handgrip, str):
            handgrip = griberies.parse_GRIBstr_todict(handgrip)
        matchingfields = FieldSet()
        filtering_keys = [str(k) for k in handgrip.keys()]  # gribapi str/unicode incompatibility
        for i, k in enumerate(filtering_keys):
            if isinstance(handgrip[k], int):
                filtering_keys[i] = str(k + ':l')  # force key to be selected as int
        idx = lowlevelgrib.index_new_from_file(self._open_through,
                                               filtering_keys)
        # filter
        for k, v in handgrip.items():
            if isinstance(v, str):  # gribapi str/unicode incompatibility
                v = str(v)
            lowlevelgrib.index_select(idx, k, v)
        # load messages
        while True:
            gid = lowlevelgrib.new_from_index(idx)
            if gid is None:
                break
            msg = GRIBmessage(('gribid', gid))
            if get_info_as_json or footprints_proxy_as_builder:
                # CLEANME: deprecated
                matchingfields.append(msg.asfield(getdata=getdata,
                                                  footprints_proxy_as_builder=footprints_proxy_as_builder,
                                                  get_info_as_json=get_info_as_json))
                if get_info_as_json:
                    epylog.warning("Deprecated argument: get_info_as_json; use read_misc_metadata instead.")
                if footprints_proxy_as_builder:
                    epylog.warning("Deprecated argument: footprints_proxy_as_builder.")
            else:
                matchingfields.append(msg.as_field(getdata=getdata,
                                                   read_misc_metadata=read_misc_metadata))
            del msg
        lowlevelgrib.index_release(idx)
        if len(matchingfields) == 0:
            raise epygramError("no field matching *handgrip* was found ({}).".
                               format(handgrip))
        return matchingfields

    @FileResource._openbeforedelayed
    def writefield(self, field,
                   ordering=griberies.defaults.GRIB1_ordering,
                   packing=None,
                   sample=None,
                   grib_edition=None,
                   other_GRIB_options={},
                   interpret_comment=False):
        """
        Writes a Field as a GRIBmessage into the GRIB resource.

        :param field: a :class:`epygram.base.Field` instance
        :param ordering: way of ordering data in GRIB, dict of GRIB keys.
        :param packing: options of packing and compression in GRIB (dict).
        :param sample: to use a specific sample GRIB
        :param grib_edition: to force a GRIB edition number (1, 2).
        :param other_GRIB_options: other options to be specified in GRIB,
            as a dict(GRIBkey=value).
            From v1.3.9, any GRIB2 key should be given through this rather than
            packing/ordering.
        :param interpret_comment: set additional key/values taken from
            field.comment (interpreted as json)
        """
        if not isinstance(field, H2DField):
            raise NotImplementedError("'field' argument other than a H2DField.")
        if 'generic' in field.fid and 'GRIB2' not in field.fid and 'GRIB1' not in field.fid:
            field.fid['GRIB2'] = field.fid['generic']
            remove_fid = True
        else:
            remove_fid = False
        m = GRIBmessage(('field', field),
                        ordering=ordering,
                        packing=packing,
                        sample=sample,
                        grib_edition=grib_edition,
                        other_GRIB_options=other_GRIB_options,
                        interpret_comment=interpret_comment)
        m.write_to_file(self._file)
        if remove_fid:
            del field.fid['GRIB2']

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
          distance = abs(target_value - external_field.data)
        """
        if isinstance(handgrip, str):
            handgrip = griberies.parse_GRIBstr_todict(handgrip)

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
                       cheap_height=True,
                       global_shift_center=None):
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
        :param global_shift_center: for global lon/lat grids, shift the center by the
            requested angle (in degrees). Enables a [0,360] grid
            to be shifted to a [-180,180] grid, for instance (with -180 argument).
        """
        if isinstance(handgrip, str):
            handgrip = griberies.parse_GRIBstr_todict(handgrip)

        field3d = fpx.field(fid={'GRIB':handgrip},
                            structure='3D',
                            resource=self, resource_fids=[handgrip])
        if global_shift_center is not None:
            field3d.global_shift_center(global_shift_center)
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
                          external_distance=None,
                          field3d=None,
                          **_):
        """
        Extracts a subdomain from the GRIB resource, given its handgrip
        and the geometry to use.

        :param handgrip: MUST define the parameter and the type of levels
        :param geometry: is the geometry on which extract data.
                         None to keep the geometry untouched.
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
        """
        if isinstance(handgrip, str):
            handgrip = griberies.parse_GRIBstr_todict(handgrip)

        if field3d is None:
            field3d = fpx.field(fid={'GRIB':handgrip},
                                structure='3D',
                                resource=self, resource_fids=[handgrip])

        if geometry is None or geometry == field3d.geometry:
            subdomain = field3d
            geometry = field3d.geometry
        else:
            subdomain = field3d.extract_subdomain(geometry,
                                                  interpolation=interpolation,
                                                  external_distance=external_distance)

        # preparation for vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):
            # surface height
            if (subdomain.geometry.vcoordinate.typeoffirstfixedsurface, vertical_coordinate) in ((102, 103), (103, 102), (100, 103)):
                fids = self.listfields(complete=True)
                zs = None
                for fid in fids:
                    hg = (fid['generic'].get('discipline', 255),
                          fid['generic'].get('parameterCategory', 255),
                          fid['generic'].get('parameterNumber', 255),
                          fid['generic'].get('typeOfFirstFixedSurface', 255))
                    if hg in ((0, 3, 4, 1), (0, 193, 5, 1), (0, 3, 5, 1), (0, 3, 6, 1), (2, 0, 7, 1)):
                        #(0, 193, 5) is for SPECSURFGEOPOTEN
                        fmt = [fk for fk in fid.keys() if fk != 'generic'][0]
                        zs = self.extract_subdomain(fid[fmt], geometry,
                                                    interpolation=interpolation,
                                                    external_distance=external_distance)
                        if zs.spectral:
                            zs.sp2gp()
                        if hg in ((0, 3, 4, 1), (0, 193, 5, 1)):
                            zs.setdata(zs.getdata() / constants.g0)
                        zs = zs.getdata()
                        break
                if zs is None:
                    raise epygramError("No terrain height field found, conversion height/altitude cannot be done")

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
                    surface_geopotential = zs * constants.g0
                else:
                    surface_geopotential = None

            # effective vertical coords conversion
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 100 and \
                 vertical_coordinate in (102, 103):
                vertical_mean = 'geometric'
                subdomain.geometry.vcoordinate = pressure2altitude(subdomain.geometry.vcoordinate,
                                                                   R,
                                                                   side_profiles['t'],
                                                                   vertical_mean,
                                                                   Pdep=side_profiles['pdep'],
                                                                   Phi_surf=surface_geopotential)
            elif (subdomain.geometry.vcoordinate.typeoffirstfixedsurface, vertical_coordinate) in ((102, 103), (103, 102)):
                if vertical_coordinate == 102:
                    subdomain.geometry.vcoordinate = height2altitude(subdomain.geometry.vcoordinate, zs)
                else:
                    subdomain.geometry.vcoordinate = altitude2height(subdomain.geometry.vcoordinate, zs)
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
        args = locals().copy()
        args.pop('self')
        try:
            self._what(**args)
        except Exception as e:
            print('!!! Exploration of the resource failed: use "grib_dump -O" for a raw decoding of a GRIB file !!!')
            raise e

    def _what(self, out=sys.stdout,
              mode='one+list',
              sortfields=None,
              details=None,
              **_):
        """Actual what method."""
        out.write("### FORMAT: " + self.format + "\n")
        out.write("(For a more thorough insight into GRIB files, " +
                  "use 'grib_dump' or (even better) 'grib_dump -O')")
        out.write("\n")

        if mode == 'one+list':
            onefield = self.get_message_at_position(0).as_field(getdata=False)
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
        out.write("There are: {} fields in this file.\n".format(len(listoffields)))
        if mode in ('what', 'ls', 'mars'):
            out.write(separation_line)
            while True:
                m = self.iter_messages()
                if m is None:
                    break  # end of file
                if mode == 'what':
                    _ = m.as_field(getdata=False)
                    if details == 'compression':
                        m._readattribute('packingType')
                        m._readattribute('bitsPerValue')
                elif mode in ('ls', 'mars'):
                    m.readmessage(mode, fatal=False)
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
            n = 0
            for f in listoffields:
                n += 1
                out.write('{:-^50}\n'.format(' Message: {:<4d} '.format(n)))
                for k in sorted_GRIB2_fid(f):
                    v = f[k]
                    if v != 'unknown':
                        if isinstance(v, str):
                            v = "'{}'".format(v)
                        out.write('{}: {},\n'.format(k, v))
        out.write(separation_line)
