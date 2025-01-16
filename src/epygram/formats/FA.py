#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for FA format.
"""

import datetime
import os
import copy
import numpy
import math
import re
import sys
import io

import footprints
from footprints import FPDict, FPList, proxy as fpx
from bronx.meteo.conversion import q2R
from bronx.syntax.arrays import stretch_array

from falfilfa4py import FA as FA4py
from falfilfa4py import LFI as LFI4py

from epygram import config, epygramError, util
from epygram.util import Angle, separation_line, write_formatted_fields
from epygram.base import FieldSet, FieldValidity, FieldValidityList
from epygram.resources import FileResource
from epygram.geometries import (Geometry, SpectralGeometry, ProjectedGeometry, GaussGeometry,
                                RegLLGeometry, VGeometry, AcademicGeometry, truncation_from_gridpoint_dims)
from epygram.geometries.VGeometry import (hybridP2pressure, hybridP2altitude,
                                          pressure2altitude)
from epygram.fields import MiscField, H2DField
from epygram.formats.fafields import SfxFldDesc_Mod, get_generic_fid
from epygram.extra.griberies.definitions.fa import fagribdef

__all__ = []

epylog = footprints.loggers.getLogger(__name__)


def _gen_headername():
    """Generates a random headername for the FA software."""
    import uuid
    return str(uuid.uuid4()).replace('-', '')[0:16]


def _create_header_from_geometry(geometry, spectral_geometry=None):
    """
    Create a header and returns its name, from a geometry and (preferably)
    a SpectralGeometry.

    :param geometry: a Geometry (or heirs) instance, from which the header is set.
    :param spectral_geometry: optional, a SpectralGeometry instance, from which
      truncation is set in header. If not provided, in LAM case the X/Y
      truncations are computed from field dimension, as linear grid.
      In global case, an error is raised.

    Contains a call to FA4py.wfacade (wrapper for FACADE routine).
    """
    assert isinstance(geometry, Geometry), \
           "geometry must be a Geometry (or heirs) instance."
    if spectral_geometry is not None and\
       not isinstance(spectral_geometry, SpectralGeometry):
        raise epygramError("spectral_geometry must be a SpectralGeometry" +
                           " instance.")
    if geometry.projected_geometry:
        assert geometry.projection['rotation'] == Angle(0., 'radians'), \
            "the geometry's projection attribute 'rotation' must be 0. in FA."

    headername = _gen_headername()
    CDNOMC = headername

    JPXPAH = FA._FAsoftware_cst['JPXPAH']
    JPXIND = FA._FAsoftware_cst['JPXIND']
    if geometry.rectangular_grid:
        # truncation
        if spectral_geometry is not None:
            truncation = spectral_geometry.truncation
        else:
            # default: linear truncation...
            truncation = truncation_from_gridpoint_dims(geometry.dimensions,
                                                        grid='linear')
        truncation_in_X = truncation['in_X']
        truncation_in_Y = truncation['in_Y']

        # scalars
        if geometry.name == 'regular_lonlat':
            KTYPTR = -11
        else:
            KTYPTR = -1 * truncation_in_X
        if geometry.name in ('lambert', 'polar_stereographic'):
            PSLAPO = (geometry.getcenter()[0].get('radians') -
                      geometry.projection['reference_lon'].get('radians'))
        else:
            PSLAPO = 0.0
        PCLOPO = 0.0
        PSLOPO = 0.0
        if geometry.name == 'academic':
            PCODIL = -1.0
        else:
            PCODIL = 0.0
        if geometry.name == 'regular_lonlat':
            KTRONC = 11
        else:
            KTRONC = truncation_in_Y
        KNLATI = geometry.dimensions['Y']
        KNXLON = geometry.dimensions['X']

        # KNLOPA
        KNLOPA = numpy.zeros(JPXPAH, dtype=numpy.int64)
        KNLOPA[0] = max(0, min(11, truncation_in_X, truncation_in_Y) - 1)
        if geometry.name == 'regular_lonlat':
            corners = geometry.gimme_corners_ij()
            KNLOPA[1] = 0
            KNLOPA[2] = 1 + corners['ll'][0]
            KNLOPA[3] = 1 + corners['ur'][0]
            KNLOPA[4] = 1 + corners['ll'][1]
            KNLOPA[5] = 1 + corners['ur'][1]
            KNLOPA[6] = 8
            KNLOPA[7] = 8
        else:
            corners = geometry.gimme_corners_ij(subzone='CI')
            if geometry.grid['LAMzone'] == 'CIE':
                KNLOPA[1] = 1
            KNLOPA[2] = 1 + corners['ll'][0]
            KNLOPA[3] = 1 + corners['ur'][0]
            KNLOPA[4] = 1 + corners['ll'][1]
            KNLOPA[5] = 1 + corners['ur'][1]
            KNLOPA[6] = geometry.dimensions['X_Iwidth']
            KNLOPA[7] = geometry.dimensions['Y_Iwidth']

        KNOZPA = numpy.zeros(JPXIND, dtype=numpy.int64)

        # PSINLA
        PSINLA = numpy.zeros(max(int((1 + geometry.dimensions['Y']) / 2), 18))
        PSINLA[0] = -1
        if geometry.name == 'regular_lonlat':
            PSINLA[1] = -9
            PSINLA[2] = 0.  # geometry.getcenter()[0].get('radians')  # TOBECHECKED:
            PSINLA[3] = 0.  # geometry.getcenter()[1].get('radians')  # TOBECHECKED:
            PSINLA[4] = geometry.getcenter()[0].get('radians')
            PSINLA[5] = geometry.getcenter()[1].get('radians')
            PSINLA[6] = geometry.grid['X_resolution'].get('radians')
            PSINLA[7] = geometry.grid['Y_resolution'].get('radians')
            PSINLA[8] = geometry.grid['X_resolution'].get('radians') * (geometry.dimensions['X'] - 1)
            PSINLA[9] = geometry.grid['Y_resolution'].get('radians') * (geometry.dimensions['Y'] - 1)
            PSINLA[10] = 0.0
            PSINLA[11] = 0.0
        elif geometry.projected_geometry:
            if geometry.secant_projection:
                raise epygramError("cannot write secant projected" +
                                   " geometries in FA.")
            PSINLA[1] = geometry.projection['reference_lat'].get('cos_sin')[1]
            PSINLA[6] = geometry.grid['X_resolution']
            PSINLA[7] = geometry.grid['Y_resolution']
            PSINLA[8] = geometry.grid['X_resolution'] * (geometry.dimensions['X'] - 1)
            PSINLA[9] = geometry.grid['Y_resolution'] * (geometry.dimensions['Y'] - 1)
            PSINLA[10] = 2. * math.pi / PSINLA[8]
            PSINLA[11] = 2. * math.pi / PSINLA[9]
            PSINLA[2] = geometry.projection['reference_lon'].get('radians')
            PSINLA[3] = geometry.projection['reference_lat'].get('radians')
            PSINLA[4] = geometry.getcenter()[0].get('radians')
            PSINLA[5] = geometry.getcenter()[1].get('radians')
        if geometry.name != 'academic':
            PSINLA[12] = Angle(geometry.ij2ll(*corners['ll'])[0], 'degrees').get('radians')
            PSINLA[13] = Angle(geometry.ij2ll(*corners['ll'])[1], 'degrees').get('radians')
            PSINLA[14] = Angle(geometry.ij2ll(*corners['ur'])[0], 'degrees').get('radians')
            PSINLA[15] = Angle(geometry.ij2ll(*corners['ur'])[1], 'degrees').get('radians')
        else:
            PSINLA[6] = geometry.grid['X_resolution']
            PSINLA[7] = geometry.grid['Y_resolution']
            PSINLA[8] = geometry.grid['X_resolution'] * (geometry.dimensions['X'] - 1)
            PSINLA[9] = geometry.grid['Y_resolution'] * (geometry.dimensions['Y'] - 1)
            PSINLA[10] = 2. * math.pi / PSINLA[8]
            PSINLA[11] = 2. * math.pi / PSINLA[9]

    else:  # global
        if geometry.name == 'reduced_gauss':
            KTYPTR = 1
            PSLAPO = 1.
            PCLOPO = 1.
            PSLOPO = 0.
        elif geometry.name == 'rotated_reduced_gauss':
            KTYPTR = 2
            PSLAPO = geometry.grid['pole_lat'].get('cos_sin')[1]
            PCLOPO = geometry.grid['pole_lon'].get('cos_sin')[0]
            PSLOPO = geometry.grid['pole_lon'].get('cos_sin')[1]
        PCODIL = geometry.grid['dilatation_coef']
        if spectral_geometry is not None:
            KTRONC = spectral_geometry.truncation['max']
        else:
            # default: linear truncation...
            KTRONC = truncation_from_gridpoint_dims(geometry.dimensions,
                                                    grid='linear',
                                                    stretching_coef=geometry.grid['dilatation_coef']
                                                    )['max']
        KNLATI = geometry.dimensions['lat_number']
        KNXLON = geometry.dimensions['max_lon_number']
        KNLOPA = numpy.zeros(JPXPAH, dtype=numpy.int64)
        KNLOPA[0:KNLATI // 2] = geometry.dimensions['lon_number_by_lat'][0:KNLATI // 2]
        KNOZPA = numpy.zeros(JPXIND, dtype=numpy.int64)
        if spectral_geometry is not None:
            KNOZPA[0:KNLATI // 2] = spectral_geometry.truncation['max_zonal_wavenumber_by_lat'][0:KNLATI // 2]
        else:
            spgeom = SpectralGeometry(space='legendre',
                                      truncation={'max':KTRONC})
            KNOZPA[0:KNLATI // 2] = spgeom.trans_inq(geometry.dimensions)[2][0:KNLATI // 2]
        PSINLA = numpy.zeros((1 + geometry.dimensions['lat_number']) // 2,
                             dtype=numpy.float64)
        PSINLA[0:KNLATI // 2] = numpy.array([s.get('cos_sin')[1] for s in
                                            geometry.grid['latitudes'][0:KNLATI // 2]])

    # vertical geometry
    PREFER = FA.reference_pressure
    if geometry.vcoordinate.grid is not None and geometry.vcoordinate.grid != {}:
        Ai = [level[1]['Ai'] for level in geometry.vcoordinate.grid['gridlevels']]
        Bi = [level[1]['Bi'] for level in geometry.vcoordinate.grid['gridlevels']]
        KNIVER = len(Ai) - 1
        PAHYBR = numpy.array(Ai) / FA.reference_pressure
        PBHYBR = numpy.array(Bi)
    else:
        KNIVER = 1
        PAHYBR = numpy.array([0., 0.])
        PBHYBR = numpy.array([0., 1.])
    LDGARD = True
    FA4py.wfacade(CDNOMC,
                  KTYPTR, PSLAPO, PCLOPO, PSLOPO,
                  PCODIL, KTRONC,
                  KNLATI, KNXLON,
                  len(KNLOPA), KNLOPA,
                  len(KNOZPA), KNOZPA,
                  len(PSINLA), PSINLA,
                  KNIVER, PREFER, PAHYBR, PBHYBR,
                  LDGARD)

    return headername


class FA(FileResource):
    """Class implementing all specificities for FA resource format."""

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['FA']),
                default='FA'),
            headername=dict(
                optional=True,
                info="With openmode == 'w', name of an existing header," +
                     " for the new FA to use its geometry."),
            validity=dict(
                type=FieldValidityList,
                optional=True,
                access='rwx',
                info="With openmode == 'w', describes the temporal validity" +
                     " of the resource."),
            default_compression=dict(
                type=FPDict,
                optional=True,
                info="Default compression for writing fields in resource."),
            cdiden=dict(
                optional=True,
                default='unknown',
                access='rwx',
                info="With openmode == 'w', identifies the FA by a keyword," +
                     " usually the model abbreviation."),
            processtype=dict(
                optional=True,
                default='analysis',
                access='rwx',
                info="With openmode == 'w', identifies the processus that" +
                     " produced the resource."),
        )
    )

    # reference pressure coefficient for converting hybrid A coefficients in FA
    reference_pressure = config.FA_default_reference_pressure
    # FA fields dicts
    gribdef = fagribdef
    sfxflddesc = SfxFldDesc_Mod(actual_init=False)

    @classmethod
    def _FAsoft_init(cls):
        """Initialize the FA software maximum dimensions."""
        cls._FAsoftware_cst = dict(zip(('JPXPAH', 'JPXIND', 'JPXGEO', 'JPXNIV'),
                                       FA4py.get_facst()))

    @classmethod
    def field_type(cls, fieldname):
        """
        Get field type, either 'H2D' or the meta-data types registered in
        cls.sfxflddesc
        If field is unknown, supposed to be a H2D.
        """
        ftype = 'H2D'  # default
        if fieldname not in cls.gribdef:  # if field in gribdef: H2D
            if fieldname in cls.sfxflddesc:
                if cls.sfxflddesc.is_metadata(fieldname):
                    ftype = cls.sfxflddesc.get(fieldname).get('type')
            else:
                if cls.sfxflddesc.is_metadata(fieldname):
                    ftype = '?'
        return ftype

    def __init__(self, *args, **kwargs):
        self.isopen = False
        super(FA, self).__init__(*args, **kwargs)
        # Initialization of FA software (if necessary):
        if not hasattr(self, '_FAsoftware_cst'):
            self._FAsoft_init()
        self.fieldscompression = {}
        self._cache_find_re_in_list = {}
        if not self.fmtdelayedopen:
            self.open()

    def open(self,
             geometry=None,
             spectral_geometry=None,
             validity=None,
             openmode=None):
        """
        Opens a FA with ifsaux' FAITOU, and initializes some attributes.

        Actually, as FAITOU needs an existing header with 'w' *openmode*,
        opening file in 'w' *openmode* will require else an existing header to
        which the resource is linked, or to create a header from a geometry via
        *_create_header_from_geometry()* function. This explains the eventual
        need for *geometry/spectral_geometry/validity* in *open()* method. If
        neither headername nor geometry are available, resource is not opened:
        the *open()* will be called again at first writing of a field in
        resource.

        :param geometry: optional, must be a
          :class:`epygram.geometries.Geometry` (or heirs) instance.
        :param spectral_geometry: optional, must be a
          :class:`epygram.geometries.SpectralGeometry` instance.
        :param validity: optional, must be a :class:`epygram.base.FieldValidity` or
          a :class:`epygram.base.FieldValidityList` instance.
        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(FA, self).open(openmode=openmode)

        if self.openmode in ('r', 'a'):
            # FA already exists, including geometry and validity
            if self.openmode in ('r', 'a') and self.headername is None:
                self._attributes['headername'] = _gen_headername()
            if geometry is not None or validity is not None:
                epylog.warning(self.container.abspath + ": FA.open():" +
                               " geometry/validity argument will be ignored" +
                               " with this openmode ('r','a').")
            # open, getting logical unit
            try:
                self._unit = FA4py.wfaitou(self.container.abspath,
                                           'OLD',
                                           self.headername)
            except RuntimeError as e:
                raise IOError(e)
            self.isopen = True
            self.empty = False
            # read info
            self._attributes['cdiden'] = FA4py.wfalsif(self._unit)
            self._read_geometry()
            self._read_validity()
            if self.openmode == 'a':
                if self.default_compression is None:
                    self._attributes['default_compression'] = self._getrunningcompression()
                self._setrunningcompression(**self.default_compression)
        elif self.openmode == 'w':
            if geometry is not None:
                if not isinstance(geometry, Geometry):
                    raise epygramError("geometry must be a Geometry (or heirs) instance.")
                self._attributes['headername'] = _create_header_from_geometry(geometry, spectral_geometry)
            if validity is not None:
                if (not isinstance(validity, FieldValidity)) and (not isinstance(validity, FieldValidityList)):
                    raise epygramError("validity must be a FieldValidity or FieldValidityList instance.")
                if isinstance(validity, FieldValidityList) and len(validity) != 1:
                    raise epygramError("FA can hold only one validity.")
                self._attributes['validity'] = validity
            if self.headername is not None and self.validity is not None:
                # new FA, with an already existing header and validity
                # set geometry from existing header
                self._read_geometry()
                # open
                if os.path.exists(self.container.abspath):
                    if config.protect_unhappy_writes and not self._overwrite:
                        raise IOError('This file already exist: ' + self.container.abspath)
                    else:
                        os.remove(self.container.abspath)
                (self._unit) = FA4py.wfaitou(self.container.abspath,
                                           'NEW',
                                           self.headername)
                self.isopen = True
                self.empty = True
                # set FA date
                if self.validity.get() != self.validity.getbasis():
                    self.processtype = 'forecast'
                self._set_validity()
                # set CDIDEN
                FA4py.wfautif(self._unit, self.cdiden)
                if self.default_compression is None:
                    self._attributes['default_compression'] = config.FA_default_compression
                self._setrunningcompression(**self.default_compression)
            elif self.headername is None or self.validity is None:
                # a header need to be created prior to opening:
                # header definition (then opening) will be done from the
                # geometry taken in the first field to be written in resource
                pass

    def close(self):
        """Closes a FA with ifsaux' FAIRME."""
        if self.isopen:
            try:
                FA4py.wfairme(self._unit, 'KEEP')
            except Exception:
                raise IOError("closing " + self.container.abspath)
            self.isopen = False

################
# ABOUT FIELDS #
################
    def find_fields_in_resource(self, seed=None, fieldtype=[], generic=False):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        :param seed: might be a regular expression, a list of regular expressions
          or *None*. If *None* (default), returns the list of all fields in
          resource.
        :param fieldtype: optional, among ('H2D', 'Misc') or a list of these strings.
          If provided, filters out the fields not of the given types.
        :param generic: if True, returns complete fid's,
          union of {'FORMATname':fieldname} and the according generic fid of
          the fields.
        """
        if isinstance(fieldtype, list):
            fieldtypeslist = fieldtype
        else:
            fieldtypeslist = [fieldtype]
        fieldslist = []
        if seed is None:
            tmplist = self.listfields()
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type_from_file(f) in fieldtypeslist:
                    fieldslist.append(f)
        elif isinstance(seed, str):
            h = (hash(seed), hash(tuple(self.listfields())))
            if h not in self._cache_find_re_in_list:
                self._cache_find_re_in_list[h] = util.find_re_in_list(seed, self.listfields())
            tmplist = self._cache_find_re_in_list[h]
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type_from_file(f) in fieldtypeslist:
                    fieldslist.append(f)
        elif isinstance(seed, list):
            tmplist = []
            for s in seed:
                h = (hash(s), hash(tuple(self.listfields())))
                if h not in self._cache_find_re_in_list:
                    self._cache_find_re_in_list[h] = util.find_re_in_list(s, self.listfields())
                tmplist += self._cache_find_re_in_list[h]
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type_from_file(f) in fieldtypeslist:
                    fieldslist.append(f)
        if fieldslist == []:
            raise epygramError("no field matching: " + str(seed) +
                               " was found in resource " +
                               self.container.abspath)
        if generic:
            fieldslist = [(f, self._hook_generic_fids(f, get_generic_fid(f)))
                          for f in fieldslist]

        return fieldslist

    def listfields(self, **kwargs):
        """
        Returns a list containing the FA identifiers of all the fields of the
        resource.
        """
        return super(FA, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for FA.

        :param complete: - if True method returns a list of {'FA':FA_fid,
                           'generic':generic_fid}
                         - if False method return a list of FA_fid
        """
        records_number = LFI4py.wlfinaf(self._unit)[0]
        LFI4py.wlfipos(self._unit)  # rewind
        fieldslist = []
        for i in range(records_number):
            fieldname = (LFI4py.wlficas(self._unit, True)[0])
            if i >= 7 and fieldname != 'DATX-DES-DONNEES':
                fieldslist.append(fieldname.strip())
            # i >= 7: 7 first fields in LFI are the header ("cadre")
            # 8th field added by P.Marguinaud, DATX-DES-DONNEES, to store dates
            # with 1-second precision, from cy40t1 onwards.
        if complete:
            fieldslist = [{'FA':f, 'generic':get_generic_fid(f)} for f in fieldslist]
            for f in fieldslist:
                f['generic'] = self._hook_generic_fids(f['FA'], f['generic'])
        return fieldslist

    def _hook_generic_fids(self, fid, generic_fid):
        if fid == 'SURFPRESSION':
            # ! hint: SURFPRESSION might be ln(sp) or sp.
            # One way to discriminate is spectralness: spectral SURFPRESSION must be ln(sp)
            if self.fieldencoding(fid)['spectral']:
                generic_fid['parameterNumber'] = 25  # sp = 0 // ln(sp) = 25
        return generic_fid

    def split_UV(self, fieldseed):
        """
        Return two lists of fids corresponding respectively to U and V
        components of wind, given a *fieldseed*.
        """
        fids = self.find_fields_in_resource(fieldseed + '*')
        if fieldseed.startswith('S'):
            Ufid = [f for f in fids if 'WIND.U.PHYS' in f]
            Vfid = [f for f in fids if 'WIND.V.PHYS' in f]
        elif fieldseed[0] in ('P', 'H', 'V'):
            Ufid = [f for f in fids if 'VENT_ZONAL' in f]
            Vfid = [f for f in fids if 'VENT_MERID' in f]
        elif fieldseed.startswith('CLS') and 'VENT' in fids[0]:
            if fieldseed.startswith('CLSVENTNEUTRE'):
                Ufid = [f for f in fids if 'CLSVENTNEUTRE.U' in f]
                Vfid = [f for f in fids if 'CLSVENTNEUTRE.V' in f]
            else:
                Ufid = [f for f in fids if 'VENT.ZONAL' in f]
                Vfid = [f for f in fids if 'VENT.MERIDIEN' in f]
        elif fieldseed.startswith('CLS') and 'RAF' in fids[0]:
            Ufid = [f for f in fids if 'CLSU' in f]
            Vfid = [f for f in fids if 'CLSV' in f]
        else:
            raise NotImplementedError("split_UV: field syntax='" + fieldseed + "'")

        return (sorted(Ufid), sorted(Vfid))

    def sortfields(self):
        """
        Returns a sorted list of fields with regards to their name and nature,
        as a dict of lists.
        """
        re_3D = re.compile(r'(?P<prefix>[A-Z])(?P<level>\d+)(?P<param>[A-Z]+.*)')
        list3D = []
        list2D = []
        # final lists
        list3Dsp = []
        list3Dgp = []
        list2Dsp = []
        list2Dgp = []
        listMisc = []

        params3D = {}
        for f in self.listfields():
            info = self.gribdef.FA2GRIB(f)
            # separate H2D from Misc
            if self._field_type_from_file(f) == 'H2D':
                # separate 3D from 2D
                if info['typeOfFirstFixedSurface'] in (119, 100, 103, 109, 20):
                    re_ok = re_3D.match(f)
                    if re_ok:
                        list3D.append(f)
                        param = re_ok.group('prefix') + re_ok.group('param')
                        if param not in params3D:
                            params3D[param] = []
                        params3D[param].append(f)
                    else:
                        list2D.append(f)
                else:
                    list2D.append(f)
            else:
                listMisc.append(f)
        # separate gp/sp
        for f in list2D:
            if self.fieldencoding(f)['spectral']:
                list2Dsp.append(f)
            else:
                list2Dgp.append(f)
        # sort 2D
        list2Dsp.sort()
        list2Dgp.sort()
        listMisc.sort()
        # sort 3D
        for p in sorted(params3D.keys()):
            if self.fieldencoding(params3D[p][0])['spectral']:
                list3Dsp.extend(sorted(params3D[p]))
            else:
                list3Dgp.extend(sorted(params3D[p]))
        outlists = {'3D spectral fields':list3Dsp,
                    '3D gridpoint fields':list3Dgp,
                    '2D spectral fields':list2Dsp,
                    '2D gridpoint fields':list2Dgp,
                    'Misc-fields':listMisc}

        return outlists

    @FileResource._openbeforedelayed
    def fieldencoding(self, fieldname, update_fieldscompression=False):
        """
        Returns a dict containing info about how the field **fieldname**
        is encoded: spectralness and compression. Interface to ifsaux' FANION.
        If **update_fieldscompression**, store compression info in attribute
        fieldscompression.
        """
        try:
            (LDCOSP,
             KNGRIB,
             KNBITS,
             KSTRON,
             KPUILA) = FA4py.wfanion(self._unit,
                                     fieldname[0:4],
                                     0,
                                     fieldname[4:])[1:6]
        except RuntimeError as e:
            if 'falfilfa4py: Error code -93 was raised' in str(e) or \
               'falfilfa4py: Error code -91 was raised' in str(e):
                raise epygramError(fieldname + ': seems like you try to read a MiscField as a H2DField...')
            else:
                raise e
        encoding = {'spectral':LDCOSP, 'KNGRIB':KNGRIB, 'KNBITS':KNBITS,
                    'KSTRON':KSTRON, 'KPUILA':KPUILA}
        if update_fieldscompression:
            # Save compression in FA
            compression = {'KNGRIB':encoding['KNGRIB'],
                           'KNBPDG':encoding['KNBITS'],
                           'KNBCSP':encoding['KNBITS'],
                           'KSTRON':encoding['KSTRON'],
                           'KPUILA':encoding['KPUILA']}
            self.fieldscompression[fieldname] = compression

        return encoding

    @FileResource._openbeforedelayed
    def readfield(self, fieldname,
                  getdata=True,
                  footprints_proxy_as_builder=config.footprints_proxy_as_builder):
        """
        Reads one field, given its FA name, and returns a Field instance.
        Interface to Fortran routines from 'ifsaux'.

        :param fieldname: FA fieldname
        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        :param footprints_proxy_as_builder: if *True*, uses footprints.proxy
          to build fields.
        """
        if self.openmode == 'w':
            raise epygramError("cannot read fields in resource if with" +
                               " openmode == 'w'.")
        assert fieldname in self.listfields(), ' '.join(["field",
                                                         str(fieldname),
                                                         "not found in resource."])
        if self._field_type_from_file(fieldname) == 'H2D':
            field = self._readH2DField(fieldname,
                                       getdata=getdata,
                                       footprints_proxy_as_builder=footprints_proxy_as_builder)
        else:
            field = self._readMiscField(fieldname, getdata=getdata)
        return field

    @FileResource._openbeforedelayed
    def _readH2DField(self, fieldname,
                      getdata=True,
                      footprints_proxy_as_builder=config.footprints_proxy_as_builder):
        """
        Reads one H2D field, given its FA name, and returns a Field instance.
        Interface to Fortran routine FACILO.

        :param fieldname: FA fieldname
        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        :param footprints_proxy_as_builder: if *True*, uses footprints.proxy
          to build fields.
        """
        if footprints_proxy_as_builder:
            builder = fpx.field
        else:
            builder = H2DField
        encoding = self.fieldencoding(fieldname, update_fieldscompression=True)
        # vertical geometry
        kwargs_vcoord = {'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                         'position_on_grid': self.geometry.vcoordinate.position_on_grid,
                         'grid': self.geometry.vcoordinate.grid,
                         'levels': self.geometry.vcoordinate.levels}
        field_info = self.gribdef.FA2GRIB(fieldname)
        # change from default (FA header) to actual levels
        kwargs_vcoord['typeoffirstfixedsurface'] = field_info.get('typeOfFirstFixedSurface', 0)
        if 'level' in field_info:
            kwargs_vcoord['levels'] = [field_info['level']]
        else:
            kwargs_vcoord['levels'] = [0]
        if 'scaledValueOfFirstFixedSurface' in field_info:
            exp = field_info.get('scaleFactorOfFirstFixedSurface', 0)
            kwargs_vcoord['levels'] = [field_info['scaledValueOfFirstFixedSurface'] * (10 ** -exp)]
        if kwargs_vcoord['typeoffirstfixedsurface'] != 119:  # hybrid-pressure
            kwargs_vcoord.pop('grid', None)
        vcoordinate = VGeometry(**kwargs_vcoord)
        # Prepare field dimensions
        spectral = encoding['spectral']
        if spectral and self.spectral_geometry is not None:
            if 'fourier' in self.spectral_geometry.space:
                # LAM
                gpdims = copy.deepcopy(self.geometry.dimensions)
                gpdims.update({k:v for k, v in self.geometry.grid.items() if 'resolution' in k})
                SPdatasize = self.spectral_geometry.etrans_inq(gpdims)[1]
            elif self.spectral_geometry.space == 'legendre':
                # Global
                # SPdatasize may be stored to avoid calling trans_inq ?
                SPdatasize = self.spectral_geometry.legendre_known_spectraldata_size()
                if SPdatasize is None:
                    # if not, call trans_inq
                    SPdatasize = self.spectral_geometry.trans_inq(self.geometry.dimensions)[1]
                SPdatasize *= 2  # complex coefficients
            datasize = SPdatasize
            spectral_geometry = self.spectral_geometry
        else:
            if self.geometry.rectangular_grid:
                GPdatasize = self.geometry.dimensions['X'] * self.geometry.dimensions['Y']
            else:
                GPdatasize = sum(self.geometry.dimensions['lon_number_by_lat'])
            datasize = GPdatasize
            spectral_geometry = None
        # Make geometry object
        kwargs_geom = dict(name=self.geometry.name,
                           grid=copy.copy(self.geometry.grid),
                           dimensions=self.geometry.dimensions,
                           vcoordinate=vcoordinate,
                           position_on_horizontal_grid=self.geometry.position_on_horizontal_grid,
                           geoid=config.FA_default_geoid)
        if self.geometry.projected_geometry or self.geometry.name == 'academic':
            kwargs_geom['projection'] = self.geometry.projection
        geometry = self.geometry.__class__(**kwargs_geom)
        # Create field
        fid = {self.format:fieldname}
        # Create H2D field
        fid['generic'] = FPDict(get_generic_fid(fieldname))
        cumul = field_info.get('productDefinitionTemplateNumber', None)
        if cumul is None or cumul == 0:
            validity = FieldValidity(basis=self.validity.getbasis(),
                                     term=self.validity.term())
        else:
            validity = self.validity.deepcopy()
            validity.set(statistical_process_on_duration=self.gribdef.FA2GRIB(fieldname).get('typeOfStatisticalProcessing', None))
        # MOCAGE surface fields: different terms can be stored in one file !
        if all([config.FA_allow_MOCAGE_multivalidities,
                fieldname[0:2] in ('SF', 'EM', 'DV'),
                all([c.isdigit() for c in fieldname[2:4]])]
               ):
            term_in_seconds = datetime.timedelta(seconds=3600 * int(fieldname[2:4]))
            validity.set(term=term_in_seconds)
        field = builder(fid=fid,
                        structure=geometry.structure,
                        geometry=geometry,
                        validity=validity,
                        spectral_geometry=spectral_geometry,
                        processtype=self.processtype)
        if 'gauss' in self.geometry.name and config.FA_buffered_gauss_grid:
            # trick: link the gauss lonlat grid so that it can be shared by
            # several geometry objects or fields !
            if not hasattr(self.geometry, '_buffered_gauss_grid'):
                (igrid, jgrid) = self.geometry._allocate_colocation_grid(compressed=False, as_float=True)
                self.geometry._buffered_gauss_grid = {'lons':igrid,
                                                      'lats':jgrid,
                                                      'filled':False}
            field.geometry._buffered_gauss_grid = self.geometry._buffered_gauss_grid
        # Get data if requested
        if getdata:
            if config.spectral_coeff_order == 'model':
                data, masked, masked_value = FA4py.wfacilo(datasize,
                                                           self._unit,
                                                           fieldname[0:4],
                                                           0,
                                                           fieldname[4:],
                                                           spectral)
                data = numpy.array(data)
                if masked:
                    data = numpy.ma.masked_equal(data, masked_value)
            else:
                data = numpy.array(FA4py.wfacile(datasize,
                                                 self._unit,
                                                 fieldname[0:4],
                                                 0,
                                                 fieldname[4:],
                                                 spectral))
            if not field.spectral:
                data = geometry.reshape_data(data)
            field.setdata(data)

        return field

    def _readMiscField(self, fieldname, getdata=True):
        """
        Reads one metadata field, given its FA name, and returns a MiscField instance.
        Interface to Fortran routine FALAIS.

        :param fieldname: FA fieldname
        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        """
        # Create field
        fid = {self.format:fieldname}
        fid['generic'] = FPDict()
        field = MiscField(fid=fid)
        # Get data if requested
        if getdata:
            nature = self.sfxflddesc.nature(fieldname, 'int')
            dim = self.sfxflddesc.dim(fieldname, 1)
            field_length = LFI4py.wlfinfo(self._unit, fieldname)[0]
            data = FA4py.wfalais(self._unit, fieldname, field_length)
            if dim == 0:
                if nature == 'int':
                    dataOut = data.view('int64')[0]
                elif nature == 'str':
                    dataInt = data.view('int64')
                    dataOut = ""
                    for num in dataInt:
                        dataOut += chr(num)
                elif nature == 'bool':
                    dataOut = bool(data.view('int64')[0])
                elif nature == 'float':
                    dataOut = data[0]
                else:
                    raise NotImplementedError("field={}, reading of datatype: {}".format(
                        fieldname, nature))
            else:
                # copy is necessary for garbage collector
                if nature == 'int':
                    dataOut = numpy.copy(data.view('int64')[:])
                elif nature == 'float':
                    dataOut = numpy.copy(data)
                elif nature == 'str':
                    raise NotImplementedError("reading of datatype: " +
                                              nature + " array.")
                    dataOut = numpy.copy(data)
                elif nature == 'bool':
                    dataOut = numpy.copy(data.view('bool')[:])
                else:
                    raise NotImplementedError("field={}, reading of datatype: {} array".format(
                        fieldname, nature))
            field.setdata(dataOut)
        return field

    def readfields(self, requestedfields=None, getdata=True):
        """
        Returns a :class:`epygram.base.FieldSet` containing requested fields
        read in the resource.

        :param requestedfields: might be:\n
          - a regular expression (e.g. 'S\*WIND.[U,V].PHYS')
          - a list of FA fields identifiers with regular expressions (e.g.
            ['SURFTEMPERATURE', 'S0[10-20]WIND.?.PHYS'])
          - if not specified, interpretated as all fields that will be found in
            resource
        :param getdata: optional, if *False*, only metadata are read, the fields
          do not contain data. Default is *True*.
        """
        requestedfields = self.find_fields_in_resource(requestedfields)
        if requestedfields == []:
            raise epygramError("unable to find requested fields in resource.")

        return super(FA, self).readfields(requestedfields, getdata)

    def writefield(self, field, compression=None):
        """
        Write a field in the resource.

        :param field: a :class:`epygram.base.Field` instance or
          :class:`epygram.fields.H2DField`.
        :param compression: optional, a (possibly partial) dict containing
          parameters for field compression (in case of a
          :class:`epygram.fields.H2DField`). Ex: {'KNGRIB': 2, 'KDMOPL': 5,
          'KPUILA': 1, 'KSTRON': 10, 'KNBPDG': 24, 'KNBCSP': 24}
        """
        if self.openmode == 'r':
            raise IOError("cannot write field in a FA with openmode 'r'.")

        if not self.isopen:
            if not isinstance(field, H2DField):
                # FA need a geometry to be open. Maybe this FA has not been
                # given one at opening. For opening it, either write a H2DField
                # in, or call its method open(geometry, validity), geometry
                # being a Geometry (or heirs), validity being a FieldValidity.
                raise epygramError("cannot write a this kind of field on a" +
                                   " non-open FA.")
            if self.validity is None:
                if len(field.validity) != 1:
                    raise epygramError("FA can hold only one validity.")
                self._attributes['validity'] = field.validity
            self._attributes['headername'] = _create_header_from_geometry(field.geometry,
                                                                          field.spectral_geometry)
            self.open()

        if not self.empty and field.fid[self.format] in self.listfields():
            epylog.info("there already is a field with the same name in" +
                        " this FA: overwrite.")

        if isinstance(field, MiscField):
            data = field.getdata()
            if data.shape in ((1,), ()):
                if data.shape == ():
                    data = data.reshape((1,))
                if 'int' in field.datatype.name:
                    dataReal = data.view('float64')
                elif 'str' in field.datatype.name or 'unicode' in field.datatype.name:
                    dataReal = numpy.array([ord(d) for d in data[0]]).view('float64')
                elif 'bool' in field.datatype.name:
                    dataReal = numpy.array([1] if data else [0]).view('float64')
                elif 'float' in field.datatype.name:
                    dataReal = data
                else:
                    raise NotImplementedError("writing of datatype " +
                                              field.datatype.name + ".")
            else:
                try:
                    dataReal = numpy.copy(data.view('float64'))
                except Exception:
                    raise NotImplementedError("writing of datatype " +
                                              field.datatype.__name__ +
                                              " array.")
            FA4py.wfaisan(self._unit,
                          field.fid[self.format],
                          dataReal.size,
                          dataReal)

        elif isinstance(field, H2DField):
            assert (self.geometry.name == field.geometry.name and
                    self.geometry.get_datashape(d4=True, force_dimZ=1, dimT=1) ==
                    field.geometry.get_datashape(d4=True, force_dimZ=1, dimT=1)), \
                "gridpoint geometry incompatibility: a FA can hold only one geometry."
            if field.geometry.vcoordinate.grid and field.geometry.vcoordinate.typeoffirstfixedsurface == 119:
                # tolerant check because of encoding differences in self
                levels = len(self.geometry.vcoordinate.grid['gridlevels'])
                s = self.geometry.vcoordinate.grid['gridlevels']
                f = field.geometry.vcoordinate.grid['gridlevels']
                diffmax = max(numpy.array([s[k][1]['Ai'] - f[k][1]['Ai']
                                           for k in range(levels)]).max(),
                              numpy.array([s[k][1]['Bi'] - f[k][1]['Bi']
                                           for k in range(levels)]).max())
            else:
                diffmax = 0.
            assert diffmax < config.epsilon or not field.geometry.vcoordinate.grid, \
                   "vertical geometry mismatch between field and file."

            if field.spectral_geometry is not None and\
               field.spectral_geometry != self.spectral_geometry:
                # compatibility check
                print(field.spectral_geometry._transforms_lib)
                print(self.spectral_geometry._transforms_lib)
                raise epygramError("spectral geometry incompatibility:" +
                                   " a FA can hold only one geometry.")
            if self.validity.cumulativeduration() is None and field.validity.cumulativeduration() is not None:
                self.validity.set(cumulativeduration=field.validity.cumulativeduration())
                self._set_validity()
            if compression is not None:
                modified_compression = True
            elif field.fid[self.format] in self.fieldscompression:
                compression = self.fieldscompression[field.fid[self.format]]
                modified_compression = True
            else:
                modified_compression = False
                compression = self._getrunningcompression()
            data = field.getdata(d4=True)
            masked = False
            fill_value = 1e20
            if isinstance(data, numpy.ma.core.MaskedArray):
                if compression.get('KNGRIB') > 4:  # GRIB2 compression: deal with bitmap
                    fill_value = data.fill_value
                elif compression.get('KNGRIB') != 0:  # old school compression, need smooth fields
                    fill_value = field.mean()
                data = field.geometry.fill_maskedvalues(data, fill_value=fill_value)
                if compression.get('KNGRIB') > 4:  # GRIB2 compression: if actually masked, activate a bitmap
                    if fill_value in data:
                        masked = True
            data = stretch_array(data.squeeze())
            if compression.get('KNBPDG', config.FA_max_encoding) > config.FA_max_encoding:
                epylog.warning(('FA compression higher than {} ' +
                                '(bits per gridpoint): {} : may be untrustful.').
                               format(config.FA_max_encoding,
                                      compression['KNBPDG']))
            if modified_compression:
                self._setrunningcompression(**compression)
            if config.spectral_coeff_order == 'model':
                FA4py.wfaieno(self._unit,
                              field.fid[self.format][0:4],
                              0,
                              field.fid[self.format][4:],
                              len(data), data,
                              field.spectral,
                              masked,
                              fill_value)
            else:
                # FIXME: next export version CLEANME: when everybody can use faieno (CY41T1_op1 onwards)
                FA4py.wfaienc(self._unit,
                              field.fid[self.format][0:4],
                              0,
                              field.fid[self.format][4:],
                              len(data), data,
                              field.spectral)
            if modified_compression:
                # set back to default
                self._setrunningcompression(**self.default_compression)
            if field.fid[self.format] not in self.fieldscompression:
                self.fieldscompression[field.fid[self.format]] = compression
        if self.empty:
            self.empty = False

    def writefields(self, fieldset, compression=None):
        """
        Write the fields of the *fieldset* in the resource.

        :param fieldset: must be a :class:`epygram.base.FieldSet` instance.
        :param compression: must be a list of compression dicts
          (cf. *writefield()* method), of length equal to the length of the
          *fieldset*, and with the same order.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        if len(fieldset) != len(compression):
            raise epygramError("fieldset and compression must have the same" +
                               " length.")

        # Be sure the first field to be written is an H2DField,
        # for being able to set the header if necessary
        if not self.isopen and not isinstance(fieldset[0], H2DField):
            for f in range(len(fieldset)):
                if isinstance(fieldset[f], H2DField):
                    fieldset.insert(0, fieldset.pop(f))
                    if compression is not None:
                        compression.insert(0, compression.pop(f))
                    break

        if compression is not None:
            # loop separated from the above one,
            # because fieldset is there-above modified
            for f in range(len(fieldset)):
                self.writefield(fieldset[f], compression[f])
        else:
            super(FA, self).writefields(fieldset)

    def rename_field(self, fieldname, new_name):
        """Renames a field "in place"."""
        LFI4py.wlfiren(self._unit, fieldname, new_name)

    def delfield(self, fieldname):
        """Deletes a field from file "in place"."""
        LFI4py.wlfisup(self._unit, fieldname)

    def modify_validity(self, **kwargs):
        """
        Modify the validity of the resource in place.
        All **kwargs** to be passed to self.validity.set(**kwargs)
        """
        self.validity.set(**kwargs)
        self._set_validity()

    @FileResource._openbeforedelayed
    def extractprofile(self, pseudoname, lon=None, lat=None,
                       geometry=None,
                       vertical_coordinate=None,
                       interpolation='nearest',
                       cheap_height=True,
                       external_distance=None):
        """
        Extracts a vertical profile from the FA resource, given its pseudoname
        and the geographic location (*lon*/*lat*) of the profile.

        :param pseudoname: must have syntax: 'K\*PARAMETER',
          K being the kind of surface (S,P,H,V),
          \* being a true star character,
          and PARAMETER being the name of the parameter requested,
          as named in FA.
        :param lon: the longitude of the desired point.
        :param lat: the latitude of the desired point.
          If both None, extract a horizontally-averaged profile.
        :param geometry: can replace *lon*/*lat*, geometry on which to extract
          data. If None, it is built from *lon*/*lat*.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V1DField, as number of GRIB2 norm:
          http://apps.ecmwf.int/codes/grib/format/grib2/ctables/4/5,
          (cf. `epygram.geometries.vertical_coordinates` possible values).
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
        if geometry is None:
            if None in [lon, lat]:
                if not (lon is None and lat is None):
                    raise ValueError('*lon* and *lat* arguments must be ' +
                                     'both None or both not None')
                # mean profile, vertical_coordinate is forgotten
                field3d = self._mk_3dvirtuelfield(pseudoname)
                if field3d.spectral:
                    field3d.sp2gp()
                profG = self.geometry.make_profile_geometry(lon, lat)
                profG.vcoordinate = field3d.geometry.vcoordinate.deepcopy()
                profile = fpx.field(fid=field3d.fid,
                                    structure='V1D',
                                    validity=self.validity.deepcopy(),
                                    geometry=profG,
                                    comment='horizontally-averaged profile')
                data3d = field3d.getdata(d4=True)
                data1d = [data3d[:, i, :, :].mean() for i in range(len(profG.vcoordinate.levels))]
                profile.setdata(data1d)
                return profile
            if self.geometry is None:
                self._read_geometry()
            pointG = self.geometry.make_profile_geometry(lon, lat)
        else:
            if lon is not None or lat is not None:
                raise epygramError("You cannot provide lon or lat when geometry is given")
            if geometry.structure != "V1D":
                raise epygramError("geometry must be a V1D")
            pointG = geometry

        profile = self.extract_subdomain(pseudoname, pointG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         external_distance=external_distance,
                                         cheap_height=cheap_height)

        return profile

    @FileResource._openbeforedelayed
    def extractsection(self, pseudoname, end1=None, end2=None,
                       geometry=None,
                       points_number=None,
                       resolution=None,
                       vertical_coordinate=None,
                       interpolation='linear',
                       cheap_height=True,
                       global_shift_center=None):
        """
        Extracts a vertical section from the FA resource, given its pseudoname
        and the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param pseudoname: must have syntax: 'K\*PARAMETER',
          K being the kind of surface (S,P,H,V),
          \* being a true star character,
          and PARAMETER being the name of the parameter requested, as named in
          FA.
        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param geometry: can replace end1/end2, geometry on which to extract
          data. If None, defaults to
          linearily spaced positions computed from *points_number*.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param vertical_coordinate: defines the requested vertical coordinate of
          the V2DField aka typeOfFirstFixedSurface in GRIB2, (cf.
          `epygram.geometries.vertical_coordinates` possible values).
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
        if geometry is None:
            if None in [end1, end2]:
                raise epygramError("You must give a geometry or end1 *and* end2")
            if self.geometry is None:
                self._read_geometry()
            sectionG = self.geometry.make_section_geometry(end1, end2,
                                                           points_number=points_number,
                                                           resolution=resolution)
        else:
            if end1 is not None or end2 is not None:
                raise epygramError("You cannot provide end1 or end2 when geometry is given")
            if geometry.structure != "V2D":
                raise epygramError("geometry must be a V2D")
            sectionG = geometry

        section = self.extract_subdomain(pseudoname, sectionG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         cheap_height=cheap_height,
                                         global_shift_center=global_shift_center)

        return section

    @FileResource._openbeforedelayed
    def extract_subdomain(self,
                          pseudoname,
                          geometry,
                          vertical_coordinate=None,
                          interpolation='linear',
                          cheap_height=True,
                          external_distance=None,
                          global_shift_center=None):
        """
        Extracts a subdomain from the FA resource, given its fid
        and the geometry to use.

        :param pseudoname: must have syntax: 'K\*PARAMETER',
          K being the kind of surface (S,P,H,V),
          \* being a true star character,
          and PARAMETER being the name of the parameter requested, as named in
          FA.
        :param geometry: is the geometry on which extract data.
                         None to keep the geometry untouched.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. `epygram.geometries.vertical_coordinates`
          possible values).
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
        field3d = self._mk_3dvirtuelfield(pseudoname)
        if field3d.spectral:
            field3d.sp2gp()
        if global_shift_center is not None:
            field3d.global_shift_center(global_shift_center)
        
        if geometry is None or geometry == field3d.geometry:
            subdomain = field3d
            geometry = field3d.geometry
        else:
            subdomain = field3d.extract_subdomain(geometry, interpolation=interpolation)

        # preparation for vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):
            # choose vertical_mean with regards to H/NH
            if 'S001PRESS.DEPART' in self.listfields():
                vertical_mean = 'geometric'
            else:
                vertical_mean = 'arithmetic'
            # surface pressure (hybridP => P,A,H)
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
               vertical_coordinate in (100, 102, 103):
                h_shape = geometry.get_datashape(force_dimZ=2)[1:] #workaround for inconsistency between rectangular and gauss, self.get_datashape(force_dimZ=1)
                Psurf = self.readfield('SURFPRESSION')
                if Psurf.spectral:
                    Psurf.sp2gp()
                ps_transect = numpy.exp(Psurf.getvalue_ll(*geometry.get_lonlat_grid(),
                                                          interpolation=interpolation,
                                                          one=False,
                                                          external_distance=external_distance))
                ps_transect = ps_transect.reshape(h_shape)
                del Psurf
            # P => H necessary profiles
            if vertical_coordinate in (102, 103):
                side_profiles = {'t':'*TEMPERATURE',
                                 'q':'*HUMI.SPECIFI',
                                 'pdep':'*PRESS.DEPART',
                                 'ql':'*CLOUD_WATER',
                                 'qi':'*ICE_CRYSTAL',
                                 'qs':'*SNOW',
                                 'qr':'*RAIN',
                                 'qg':'*GRAUPEL'}
                for p in sorted(side_profiles.keys(), reverse=True):  # reverse to begin by t
                    try:
                        # try to extract profiles for each lon/lat and each parameter
                        if pseudoname == side_profiles[p]:
                            # already extracted as requested profile
                            side_profiles[p] = subdomain
                        else:
                            if cheap_height and p not in ('t', 'q'):
                                raise epygramError()  # to go through "except" instructions below
                            side_profiles[p] = self.extract_subdomain(side_profiles[p], geometry,
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
                          ['q', 'ql', 'qi', 'qr', 'qs', 'qg']])
            if vertical_coordinate == 102:
                try:
                    geopotential = self.readfield('SPECSURFGEOPOTEN')
                except epygramError:
                    geopotential = self.readfield('SURFGEOPOTENTIEL')
                else:
                    geopotential.sp2gp()
                surface_geopotential = geopotential.getvalue_ll(*geometry.get_lonlat_grid(),
                                                                interpolation=interpolation,
                                                                one=False,
                                                                external_distance=external_distance)
                surface_geopotential = surface_geopotential.reshape(h_shape)
                del geopotential
            else:
                surface_geopotential = None

            # effective vertical coords conversion
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
               vertical_coordinate == 100:
                subdomain.geometry.vcoordinate = hybridP2pressure(subdomain.geometry.vcoordinate,
                                                                  ps_transect,
                                                                  vertical_mean,
                                                                  gridposition='mass')
            elif subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
                 vertical_coordinate in (102, 103):
                subdomain.geometry.vcoordinate = hybridP2altitude(subdomain.geometry.vcoordinate,
                                                                  R,
                                                                  side_profiles['t'],
                                                                  ps_transect,
                                                                  vertical_mean,
                                                                  Pdep=side_profiles['pdep'],
                                                                  Phi_surf=surface_geopotential)
            elif subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 100 and \
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

    def _mk_3dvirtuelfield(self, pseudoname):
        """Return a D3VirtualField from pseudoname."""
        fidlist = self.find_fields_in_resource(seed=pseudoname,
                                               fieldtype=['H2D', '3D'],
                                               generic=True)
        if fidlist == []:
            raise epygramError("cannot find profile for " + str(pseudoname) +
                               " in resource.")
        # find the prevailing type of level
        leveltypes = [f[1]['typeOfFirstFixedSurface'] for f in fidlist]
        if len(set(leveltypes)) > 1:
            leveltypes_num = {t:0 for t in set(leveltypes)}
            for t in leveltypes:
                leveltypes_num[t] += 1
            leveltypes = [k for k, v in leveltypes_num.items() if
                          v == max(leveltypes_num.values())]
            if len(leveltypes) > 1:
                raise epygramError("unable to determine type of level" +
                                   " to select.")
        leveltype = leveltypes[0]
        # filter by type of level
        fidlist = [f[0] for f in fidlist if f[1]['typeOfFirstFixedSurface'] == leveltype]

        field3d = fpx.field(fid={'FA':pseudoname},
                            structure='3D',
                            resource=self, resource_fids=fidlist)
        return field3d

###########
# pre-app #
###########

    @FileResource._openbeforedelayed
    def what(self,
             out=sys.stdout,
             details=None,
             sortfields=False,
             **_):
        """
        Writes in file a summary of the contents of the FA.

        :param out: the output open file-like object.
        :param details: 'spectral' if spectralness of fields is requested;
                        'compression' if information about fields compression
                        is requested.
        :param sortfields: **True** if the fields have to be sorted by type.
        """
        for f in self.listfields():
            if self._field_type_from_file(f) == 'H2D':
                first_H2DField = f
                break
        if len(self.listfields()) == 0:
            raise epygramError('empty Resource.')
        firstfield = self.readfield(first_H2DField, getdata=False)
        if not firstfield.spectral and self.spectral_geometry is not None:
            firstfield._attributes['spectral_geometry'] = self.spectral_geometry

        listoffields = self.listfields()
        if sortfields:
            sortedfields = self.sortfields()

        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        out.write("### IDENTIFIER (CDIDEN): " + self.cdiden + "\n")
        out.write("\n")

        FieldValidityList(self.validity).what(out)
        firstfield.what(out,
                        validity=False,
                        vertical_geometry=False,
                        arpifs_var_names=True,
                        fid=False)

        self.geometry.vcoordinate.what(out, levels=False)

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        if sortfields:
            listoffields = []
            for k in sorted(sortedfields.keys()):
                listoffields.append(k)
                listoffields.append('--------------------')
                listoffields.extend(sortedfields[k])
                listoffields.append('--------------------')
            numfields = sum([len(v) for v in sortedfields.values()])
        else:
            numfields = len(listoffields)
        out.write("Number: " + str(numfields) + "\n")
        if details is None:
            write_formatted_fields(out, "Field name")
        elif details == 'spectral':
            write_formatted_fields(out, "Field name", "Spectral")
        elif details == 'compression':
            params = ['KNGRIB', 'KNBITS', 'KSTRON', 'KPUILA']
            width_cp = 8
            compressionline = ""
            for p in params:
                compressionline += '{:^{width}}'.format(p, width=width_cp)
            write_formatted_fields(out, "Field name", "Spectral",
                                   compressionline)
        out.write(separation_line)
        for f in listoffields:
            if details is not None and self._field_type_from_file(f) == 'H2D':
                encoding = self.fieldencoding(f)
                if details == 'spectral':
                    write_formatted_fields(out, f, encoding['spectral'])
                elif details == 'compression':
                    compressionline = ""
                    for p in params:
                        try:
                            compressionline += '{:^{width}}'.format(str(encoding[p]), width=width_cp)
                        except KeyError:
                            compressionline += '{:^{width}}'.format('-', width=width_cp)
                    write_formatted_fields(out, f, encoding['spectral'],
                                           compression=compressionline)
            else:
                write_formatted_fields(out, f)
        out.write(separation_line)

# the FA WAY #
##############
    def _get_header(self, out=sys.stdout, mode='FA'):
        """
        Write the "header" of the resource in **out**.
        :param mode: - if 'LFI', writes the header as the corresponding LFI
                       records.
                     - if **mode**=='FA', writes the header as returned by
                       routines FACIES and FADIEX.
        """
        assert mode in ('FA', 'LFI')
        if mode == 'LFI':
            raise epygramError('that does not work yet ! Soon...')
            self.close()
            from epygram.formats import resource
            lfi = resource(self.container.abspath, 'r', fmt='LFI')
            _level = epylog.getEffectiveLevel()
            epylog.setLevel('WARNING')
            for r in ('CADRE-DIMENSIONS',
                      'CADRE-FRANKSCHMI',
                      'CADRE-REDPOINPOL',
                      'CADRE-FOCOHYBRID',
                      'CADRE-SINLATITUD',
                      'DATE-DES-DONNEES',
                      'DATX-DES-DONNEES'):
                f = lfi.readfield(r)
                out.write(r + '\n')
                out.write(str(f) + '\n')
                out.write('\n')
            epylog.setLevel(_level)
            lfi.close()
            self.open()
        elif mode == 'FA':
            vars_from_facies = ('KTYPTR', 'PSLAPO', 'PCLOPO', 'PSLOPO',
                                'PCODIL', 'KTRONC',
                                'KNLATI', 'KNXLON', 'KNLOPA', 'KNOZPA', 'PSINLA',
                                'KNIVER', 'PREFER', 'PAHYBR', 'PBHYBR')
            zvars = list(zip(vars_from_facies,
                             FA4py.wfacies(self._FAsoftware_cst['JPXPAH'],
                                           self._FAsoftware_cst['JPXIND'],
                                           self._FAsoftware_cst['JPXGEO'],
                                           self._FAsoftware_cst['JPXNIV'],
                                           self.headername)[:-1]))
            for ab in [8, 9, 10]:  # 'KNLOPA', 'KNOZPA', 'PSINLA'
                zvars[ab] = (zvars[ab][0], zvars[ab][1][:zvars[6][1] // 2])  # 6 == KNLATI
            for ab in [13, 14]:  # 'PAHYBR', 'PBHYBR'
                zvars[ab] = (zvars[ab][0], zvars[ab][1][:zvars[11][1] + 1])  # 11 == KNIVER
            for v in zvars:
                out.write(v[0] + '\n')
                out.write(str(v[1]) + '\n')
                out.write('\n')

            KDATEF = FA4py.wfadiex(self._unit)
            out.write('KDATEF\n')
            out.write(str(KDATEF) + '\n')

    def _field_type_from_file(self, fieldname):
        """Return type of the field, based on FANION or FA field dict."""
        try:
            exist = FA4py.wfanion(self._unit,
                                  fieldname[0:4],
                                  0,
                                  fieldname[4:])[0]
        except RuntimeError:
            exist = False
        ftype = self.field_type(fieldname)
        if exist and ftype == '?':
            ftype = 'H2D'  # because fanion fails or answers False for meta-fields
        return ftype

    def _raw_header_get(self):
        vars_from_facies = ('KTYPTR', 'PSLAPO', 'PCLOPO', 'PSLOPO',
                            'PCODIL', 'KTRONC',
                            'KNLATI', 'KNXLON', 'KNLOPA', 'KNOZPA', 'PSINLA',
                            'KNIVER', 'PREFER', 'PAHYBR', 'PBHYBR')
        zvars = list(zip(vars_from_facies,
                         FA4py.wfacies(self._FAsoftware_cst['JPXPAH'],
                                       self._FAsoftware_cst['JPXIND'],
                                       self._FAsoftware_cst['JPXGEO'],
                                       self._FAsoftware_cst['JPXNIV'],
                                       self.headername)[:-1]))
        return zvars

    def _raw_header_set(self, header_vars):
        h = header_vars
        FA4py.wfacade(self.headername,
                      h['KTYPTR'], h['PSLAPO'], h['PCLOPO'], h['PSLOPO'],
                      h['PCODIL'], h['KTRONC'],
                      h['KNLATI'], h['KNXLON'],
                      len(h['KNLOPA']), h['KNLOPA'],
                      len(h['KNOZPA']), h['KNOZPA'],
                      len(h['PSINLA']), h['PSINLA'],
                      h['KNIVER'], h['PREFER'], h['PAHYBR'], h['PBHYBR'],
                      True)

    def _read_geometry(self):
        """
        Reads the geometry in the FA header.
        Interface to Fortran routines from 'ifsaux'.
        """
        (KTYPTR, PSLAPO, PCLOPO, PSLOPO,
         PCODIL, KTRONC,
         KNLATI, KNXLON, KNLOPA, KNOZPA, PSINLA,
         KNIVER, PREFER, PAHYBR, PBHYBR
         ) = FA4py.wfacies(self._FAsoftware_cst['JPXPAH'],
                           self._FAsoftware_cst['JPXIND'],
                           self._FAsoftware_cst['JPXGEO'],
                           self._FAsoftware_cst['JPXNIV'],
                           self.headername)[:-1]
        Ai = [c * PREFER for c in PAHYBR[0:KNIVER + 1]]
        Bi = [c for c in PBHYBR[0:KNIVER + 1]]
        vertical_grid = {'gridlevels': tuple([(i + 1, FPDict({'Ai':Ai[i],
                                                              'Bi':Bi[i]}))
                                              for i in range(len(Ai))]),
                         'ABgrid_position':'flux'}
        kwargs_vcoord = {'typeoffirstfixedsurface':119,
                         'position_on_grid': 'mass',
                         'grid': vertical_grid,
                         'levels': list([i + 1 for i in range(len(Ai) - 1)])
                         }
        vcoordinate_read_in_header = VGeometry(**kwargs_vcoord)
        self.reference_pressure = PREFER

        rectangular_grid = KTYPTR <= 0
        if rectangular_grid:
            LMAP = int(PCODIL) != -1
            # LAM or regular lat/lon
            projected_geometry = ((int(PSINLA[0]) != 0 and
                                   int(PSINLA[1]) != -9) or  # "new" header
                                  (int(PSINLA[0]) == 0 and
                                   int(PSINLA[9]) != -9))  # "old" header
            dimensions = {'X':KNXLON,
                          'Y':KNLATI}
            if projected_geometry:
                geometryclass = ProjectedGeometry
                # LAM (projection)
                dimensions.update({'X_CIzone':KNLOPA[3] - KNLOPA[2] + 1,
                                   'Y_CIzone':KNLOPA[5] - KNLOPA[4] + 1,
                                   'X_Iwidth':KNLOPA[6],
                                   'Y_Iwidth':KNLOPA[7],
                                   'X_Czone':KNLOPA[3] - KNLOPA[2] + 1 - 2 * KNLOPA[6],
                                   'Y_Czone':KNLOPA[5] - KNLOPA[4] + 1 - 2 * KNLOPA[7],
                                   })
                if int(PSINLA[0]) != 0:
                    grid = {'X_resolution':PSINLA[6],
                            'Y_resolution':PSINLA[7]}
                else:
                    grid = {'X_resolution':PSINLA[14],
                            'Y_resolution':PSINLA[15]}
                if KNLOPA[1] == 0 or (dimensions['X'] == dimensions['X_CIzone'] and
                                      dimensions['Y'] == dimensions['Y_CIzone']):  # C+I
                    grid['LAMzone'] = 'CI'
                    dimensions['X'] = dimensions['X_CIzone']
                    dimensions['Y'] = dimensions['Y_CIzone']
                    io, jo = 0, 0
                elif abs(KNLOPA[1]) == 1:  # C+I+E
                    grid['LAMzone'] = 'CIE'
                    io, jo = KNLOPA[2] - 1, KNLOPA[4] - 1
                    dimensions['X_CIoffset'] = io
                    dimensions['Y_CIoffset'] = jo
                if LMAP and int(PSINLA[0]) == 0:
                    # 'old' header : from A. Stanesic (Croatia)
                    projection = {'reference_lon':Angle(PSINLA[7], 'radians'),
                                  'reference_lat':Angle(PSINLA[8], 'radians'),
                                  'rotation':Angle(0., 'radians')}
                    grid.update({'input_lon':Angle(PSINLA[7], 'radians'),
                                 'input_lat':Angle(PSINLA[8], 'radians'),
                                 'input_position':(io + (float(dimensions['X_CIzone']) - 1) / 2.,
                                                   jo + (float(dimensions['Y_CIzone']) - 1) / 2.)})
                    PSINLA[1] = PSINLA[9]
                    PSINLA[0] = 0
                elif LMAP and int(PSINLA[0]) != 0:
                    projection = {'reference_lon':Angle(PSINLA[2], 'radians'),
                                  'reference_lat':Angle(PSINLA[3], 'radians'),
                                  'rotation':Angle(0., 'radians')}
                    grid.update({'input_lon':Angle(PSINLA[4], 'radians'),
                                 'input_lat':Angle(PSINLA[5], 'radians'),
                                 'input_position':(io + (float(dimensions['X_CIzone']) - 1) / 2.,
                                                   jo + (float(dimensions['Y_CIzone']) - 1) / 2.)})

                if abs(PSINLA[1]) <= config.epsilon:
                    geometryname = 'mercator'
                elif 1.0 - abs(PSINLA[1]) <= config.epsilon:
                    geometryname = 'polar_stereographic'
                elif config.epsilon < abs(PSINLA[1]) < 1.0 - config.epsilon:
                    geometryname = 'lambert'
                spectral_space = 'bi-fourier'
                spectral_trunc = {'in_X':KTYPTR * -1,
                                  'in_Y':KTRONC,
                                  'shape':'elliptic'}
                if not LMAP:
                    if dimensions['X'] == 1:
                        spectral_space = 'fourier'
                        spectral_trunc = {'in_Y':KTRONC,
                                          'in_X':KTYPTR * -1}
                    grid.update({'input_lon':1,
                                 'input_lat':1,
                                 'input_position':(0, 0)})
                    projection = {'rotation':Angle(0, 'degrees'),
                                  'reference_dX':grid['X_resolution'],
                                  'reference_dY':grid['X_resolution']}
                    geometryname = 'academic'
                    geometryclass = AcademicGeometry
            elif not projected_geometry:
                # regular lat/lon
                geometryclass = RegLLGeometry
                projection = None
                geometryname = 'regular_lonlat'
                if int(PSINLA[0]) == 0:
                    # 'old' header
                    grid = {'input_lon':Angle(PSINLA[3], 'radians'),
                            'input_lat':Angle(PSINLA[4], 'radians'),
                            'input_position':(0, 0),
                            'X_resolution':Angle(PSINLA[14], 'radians'),
                            'Y_resolution':Angle(PSINLA[15], 'radians')}
                else:
                    grid = {'input_lon':Angle(PSINLA[4], 'radians'),
                            'input_lat':Angle(PSINLA[5], 'radians'),
                            'input_position':((float(dimensions['X']) - 1) / 2.,
                                              (float(dimensions['Y']) - 1) / 2.),
                            'X_resolution':Angle(PSINLA[6], 'radians'),
                            'Y_resolution':Angle(PSINLA[7], 'radians')}
                spectral_space = None
        else:
            # ARPEGE global
            geometryclass = GaussGeometry
            projection = None
            # reconstruction of tables on both hemispheres
            KNLOPA = KNLOPA[:KNLATI // 2]
            KNOZPA = KNOZPA[:KNLATI // 2]
            PSINLA = PSINLA[:KNLATI // 2]
            lon_number_by_lat = [n for n in KNLOPA] + [KNLOPA[-(n + 1)] for n in
                                                       range(0, len(KNLOPA))]
            max_zonal_wavenumber_by_lat = [n for n in KNOZPA] + \
                                          [KNOZPA[-(n + 1)] for n in
                                           range(0, len(KNOZPA))]
            latitudes = [Angle((math.cos(math.asin(sinlat)), sinlat), 'cos_sin')
                         for sinlat in PSINLA] \
                        + [Angle((math.cos(math.asin(PSINLA[-(n + 1)])),
                                  - PSINLA[-(n + 1)]),
                                 'cos_sin')
                           for n in range(0, len(PSINLA))]
            grid = {'dilatation_coef':PCODIL,
                    'latitudes':FPList([l for l in latitudes])
                    }
            if KTYPTR == 1:
                geometryname = 'reduced_gauss'
            elif KTYPTR == 2:
                geometryname = 'rotated_reduced_gauss'
                grid['pole_lat'] = Angle((math.cos(math.asin(PSLAPO)), PSLAPO),
                                         'cos_sin')
                grid['pole_lon'] = Angle((PCLOPO, PSLOPO), 'cos_sin')
            else:
                raise ValueError('wrong value of KTYPTR in FA header: ' + str(KTYPTR))
            dimensions = {'max_lon_number':KNXLON,
                          'lat_number':KNLATI,
                          'lon_number_by_lat':FPList([n for n in
                                                      lon_number_by_lat])
                          }
            spectral_space = 'legendre'
            spectral_trunc = {'max':KTRONC,
                              'shape':'triangular',
                              'max_zonal_wavenumber_by_lat':FPList([k for k in
                                                                    max_zonal_wavenumber_by_lat])
                              }
        kwargs_geom = dict(name=geometryname,
                           grid=grid,
                           dimensions=dimensions,
                           vcoordinate=vcoordinate_read_in_header,
                           position_on_horizontal_grid='center',
                           geoid=config.FA_default_geoid)
        if projection is not None:
            kwargs_geom['projection'] = projection
        self.geometry = geometryclass(**kwargs_geom)
        if spectral_space is not None:
            self.spectral_geometry = SpectralGeometry(space=spectral_space,
                                                      truncation=spectral_trunc)
        else:
            self.spectral_geometry = None

    @FileResource._openbeforedelayed
    def _read_validity(self):
        """
        Reads the validity in the FA header.
        Interface to Fortran routines from 'ifsaux'.
        """
        KDATEF = FA4py.wfadiex(self._unit)
        year = int(KDATEF[0])
        month = int(KDATEF[1])
        day = int(KDATEF[2])
        hour = int(KDATEF[3])
        minute = int(KDATEF[4])
        second = int(KDATEF[13]) - hour * 3600 - minute * 60
        if second >= 60:
            m = second // 60
            second = second % 60
            if m > 60:
                hour += m // 60
            minute = m % 60
        processtype = int(KDATEF[8])
        if processtype == 0:
            self.processtype = 'analysis'
        elif processtype == 1:
            self.processtype = 'initialization'
        elif processtype == 9:
            self.processtype = 'forcings'
        elif processtype == 10:
            self.processtype = 'forecast'
        term_in_seconds = int(KDATEF[14])
        cumulationstart_in_seconds = int(KDATEF[15])
        cumulativeduration_in_seconds = term_in_seconds - cumulationstart_in_seconds
        basis = datetime.datetime(year, month, day, hour, minute, second)
        term = datetime.timedelta(seconds=term_in_seconds)
        cumulativeduration = datetime.timedelta(seconds=cumulativeduration_in_seconds)

        self.validity = FieldValidity(basis=basis, term=term, cumulativeduration=cumulativeduration)

    @FileResource._openbeforedelayed
    def _set_validity(self, termunit='hours'):
        """
        Sets date, hour and the processtype in the resource.
        :param termunit: unit of term, among ('hours, 'days',)
        """
        if not self.isopen:
            raise epygramError("_set_validity must be called after FA is open.")

        basis = self.validity.getbasis()
        KDATEF = numpy.zeros(22, dtype=numpy.int64)
        KDATEF[0] = int(basis.year)
        KDATEF[1] = int(basis.month)
        KDATEF[2] = int(basis.day)
        KDATEF[3] = int(basis.hour)
        KDATEF[4] = int(basis.minute)
        if termunit == 'hours':
            KDATEF[5] = 1
            KDATEF[6] = self.validity.term(fmt='IntHours')
            if self.validity.cumulativeduration() is not None:
                KDATEF[9] = (self.validity.term(fmt='IntHours') -
                             self.validity.cumulativeduration(fmt='IntHours'))
        elif termunit == 'days':
            KDATEF[5] = 2
            KDATEF[6] = self.validity.term('IntHours') // 24
            if self.validity.cumulativeduration() is not None:
                KDATEF[9] = (self.validity.term('IntHours') -
                             self.validity.cumulativeduration('IntHours')) // 24
        else:
            raise NotImplementedError("term unit other than hours/days ?")
        # KDATEF[7] = 0
        if self.processtype == 'analysis':
            processtype = 0
        elif self.processtype == 'initialization':
            processtype = 1
        elif self.processtype == 'forcings':
            processtype = 9
        elif self.processtype == 'forecast':
            processtype = 10
        else:
            processtype = 255  # unknown processtype
        KDATEF[8] = processtype
        # KDATEF[10] = 0
        if not config.FA_fandax:
            FA4py.wfandar(self._unit, KDATEF)
        else:
            # fandax
            KDATEF[11] = 1
            KDATEF[13] = int(basis.second) + \
                         int(basis.minute) * 60 + \
                         int(basis.hour) * 3600
            KDATEF[14] = self.validity.term(fmt='IntSeconds')
            if self.validity.cumulativeduration() is not None:
                KDATEF[15] = (self.validity.term(fmt='IntSeconds') -
                              self.validity.cumulativeduration(fmt='IntSeconds'))
            FA4py.wfandax(self._unit, KDATEF)

    @FileResource._openbeforedelayed
    def _getrunningcompression(self):
        """
        Returns the current compression parameters of the FA (at time of writing).
        Interface to ifsaux' FAVEUR.
        """
        comp = dict()
        (comp['KNGRIB'], comp['KNBPDG'], comp['KNBCSP'], comp['KSTRON'],
         comp['KPUILA'], comp['KDMOPL']) = FA4py.wfaveur(self._unit)

        return comp

    @FileResource._openbeforedelayed
    def _setrunningcompression(self, **kwargs):
        """
        Sets the compression parameters of the FA.
        Interface to FAGOTE (cf. FAGOTE documentation for significance of
        arguments).
        """
        if self.openmode == 'r':
            raise IOError("method _setrunningcompression() can only be" +
                          " called if 'openmode' in('w', 'a').")
        comp = copy.deepcopy(self.default_compression)
        for k in kwargs.keys():
            if k in self.default_compression:
                comp[k] = kwargs[k]
            else:
                raise epygramError("unknown parameter: " + k +
                                   " passed as argument.")
        FA4py.wfagote(self._unit,
                      comp['KNGRIB'],
                      comp['KNBPDG'],
                      comp['KNBCSP'],
                      comp['KSTRON'],
                      comp['KPUILA'],
                      comp['KDMOPL'])
