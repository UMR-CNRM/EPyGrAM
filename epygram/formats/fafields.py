#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
FA fields utilities
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import os
import copy
import re
import six
import io

from bronx.syntax.decorators import nicedeco

import griberies
from epygram import config, epygramError


@nicedeco
def _init_before(mtd):
    """
    Decorator for methods: call method actual_init if not initialized before
    actually calling method.
    """
    def initialized(self, *args, **kwargs):
        if not self.initialized:
            self.actual_init()
        return mtd(self, *args, **kwargs)
    return initialized


class FaGribDef(object):
    """
    Handle FA-related GRIB definition files.
    To add user files:
    a) use env vars GRIBAPI|ECCODES_DEFINITION_PATH or
    b) move it under {config.userlocaldir}/{FaGribDef.dirname}
    """

    re_dirname = re.compile('gribapi\.def\.[\d_]+')
    _official_rootdir = os.path.join(config.installdir, 'data')
    default_grib_edition = 'grib2'
    concepts = ('faFieldName.def', 'faModelName.def', 'faLevelName.def')

    def __init__(self, actual_init=True):
        self.tables = {'grib1':{c:{} for c in self.concepts},
                       'grib2':{c:{} for c in self.concepts}}
        if actual_init:
            self.actual_init()
        else:
            self.initialized = False

    def _find_gribdefs(self, root):
        dirs = []
        for d in sorted(os.listdir(root)):
            if self.re_dirname.match(d):
                d = os.path.join(root, d)
                if os.path.isdir(d):
                    dirs.append(d)
        return dirs

    def actual_init(self):
        # official fagribdef
        defpaths = griberies.get_definition_paths()
        defpaths = self._find_gribdefs(self._official_rootdir) + defpaths
        # complete with local user definitions, if existing
        defpaths.extend(self._find_gribdefs(config.userlocaldir))
        # read gribdef files
        for d in defpaths:
            for grib_edition in ('grib1', 'grib2'):
                for concept in self.concepts:
                    self._read_localConcept(concept, d, grib_edition)
        self.initialized = True

    def _read_localConcept(self, concept, directory, grib_edition=default_grib_edition):
        pathname = os.path.join(directory, grib_edition, 'localConcepts', 'lfpw', concept)
        if os.path.exists(pathname):
            self.tables[grib_edition][concept].update(griberies.read_gribdef(pathname))

    @_init_before
    def FA2GRIB(self, fieldname,
                grib_edition=default_grib_edition,
                include_comments=False,
                fatal=False,
                filter_non_GRIB_keys=True):
        """
        Convert FA field name to GRIB fid, including vertical level identification.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param fatal: if True and fieldname is not retrieved, raise a ValueError;
            else, return a default 255 fid
        :param filter_non_GRIB_keys: filter out the non-GRIb keys that may be
            present in grib def of field
        """
        _re_altitude = re.compile('(?P<ltype>[A-W]+)(?P<level>\d+)(?P<param>.+)')
        if fieldname in self.tables[grib_edition]['faFieldName.def']:
            fid = copy.copy(self.tables[grib_edition]['faFieldName.def'][fieldname])
        elif fieldname.replace(' ', '_') in self.tables[grib_edition]['faFieldName.def']:
            fid = copy.copy(self.tables[grib_edition]['faFieldName.def'][fieldname.replace(' ', '_')])
        else:
            rematch = _re_altitude.match(fieldname)
            if rematch:
                fid = copy.copy(self.tables[grib_edition]['faLevelName.def'].get(
                    rematch.group('ltype'),
                    self.tables[grib_edition]['faLevelName.def']['default']))
                fid.update(
                    self.tables[grib_edition]['faFieldName.def'].get(
                        rematch.group('param'),
                        self.tables[grib_edition]['faFieldName.def']['default']))
                level = int(rematch.group('level'))
                if level == 0:  # formatting issue
                    level = 100000
                if rematch.group('ltype') == 'P':
                    level /= 100
                fid['level'] = level
            else:
                if fatal:
                    raise ValueError("field not found: {}".format(fieldname))
                else:
                    fid = self.tables[grib_edition]['faFieldName.def']['default']
        if not include_comments:
            if '#comment' in fid:
                fid.pop('#comment')
        if filter_non_GRIB_keys:
            for k in ('LSTCUM', 'FMULTM', 'FMULTE'):
                if k in fid:
                    fid.pop(k)
        # productDefinitionTemplateNumber to distinguish between
        # CLSTEMPERATURE and CLSMINI.TEMPERAT / CLSMAXI.TEMPERAT (for instance)
        if 'productDefinitionTemplateNumber' not in fid:
            fid['productDefinitionTemplateNumber'] = 0
        return fid

    @_init_before
    def GRIB2FA(self, gribfid, grib_edition=default_grib_edition):
        """
        Look for a unique matching field in tables.
        ! WARNING ! the unicity might not be ensured depending on the version
        of grib definitions files.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        """
        matching_fields = self.lookup_GRIB(gribfid, grib_edition)
        if len(matching_fields) == 1:
            return matching_fields.keys[0]
        elif len(matching_fields) == 0:
            raise ValueError("No field matching *fid* was found in tables.")
        elif len(matching_fields) > 1:
            raise ValueError("Several fields matching were found. Use method *lookup_GRIB* to refine.")

    @_init_before
    def lookup_FA(self, partial_fieldname,
                  grib_edition=default_grib_edition,
                  include_comments=False):
        """
        Look for all the fields which FA name contain **partial_fieldname**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        """
        fields = {}
        for f, gribfid in self.tables[grib_edition]['faFieldName.def'].items():
            if partial_fieldname in f:
                fields[f] = gribfid
        if not include_comments:
            for f, fid in fields.items():
                fid.pop('#comment')
        return fields

    @_init_before
    def lookup_GRIB(self, partial_fid,
                    grib_edition=default_grib_edition,
                    include_comments=False):
        """
        Look for all the fields which GRIB fid contain **partial_fid**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        """
        if isinstance(partial_fid, six.string_types):
            partial_fid = griberies.parse_GRIBstr_todict(partial_fid)
        partial_fid
        fields = {}
        for f, gribfid in self.tables[grib_edition]['faFieldName.def'].items():
            if set(partial_fid.keys()).issubset(gribfid.keys()):
                ok = True
                for k,v in partial_fid.items():
                    if gribfid[k] != v:
                        ok = False
                        break
                if ok:
                    fields[f] = copy.copy(gribfid)
        if not include_comments:
            for f, fid in fields.items():
                fid.pop('#comment')
        return fields

    def __call__(self, fid,
                 grib_edition=default_grib_edition,
                 include_comments=False):
        """Call methods FA2GRIB or GRIB2FA depending on nature of **fid**."""
        try:
            if isinstance(fid, six.string_types):
                fid = griberies.parse_GRIBstr_todict(fid)
        except SyntaxError:  # fid is a FA fieldname
            convfid = self.FA2GRIB(fid, grib_edition, include_comments)
        else:  # fid is a GRIB fid
            convfid = self.GRIB2FA(fid, grib_edition, include_comments)
        return convfid

    def __contains__(self, fid):
        try:
            if isinstance(fid, six.string_types):
                fid = griberies.parse_GRIBstr_todict(fid)
        except SyntaxError:  # fid is a FA fieldname
            try:
                self.FA2GRIB(fid, fatal=True)
            except ValueError:
                ok = False
            else:
                ok = True
        else:  # fid is a GRIB fid
            ok = len(self.lookup_GRIB(fid)) > 0
        return ok


class SfxFldDesc_Mod(object):
    """Handle fields catalog from sfxflddesc.F90 source file."""
    _pattern = re.compile('"\.(?P<name>[\w\d%]+)\.+(?P<gridtype>\d)\.{2}(?P<type>\w\d)\.{2}(?P<comment>.+)\.{2}(?P<mask>.{20})\.{2}",(&|/)')
    _unit = re.compile('(?P<comment>.*)\((?P<unit>.+)\)')
    _fortran_sourcename = 'sfxflddesc_mod.F90'
    type2nature = {'X':'float', 'L':'bool', 'C':'str', 'N':'int', 'Y':'float',
                   'T':'?T?'}

    def __init__(self, actual_init=True):
        self.table = {}
        if actual_init:
            self.actual_init()
        else:
            self.initialized = False

    def actual_init(self):
        """Read official file and, if existing, user local file."""
        # official one
        filename = os.path.join(config.installdir, 'data',
                                self._fortran_sourcename)
        self.read(filename)
        # local user one
        user_sfxflddesc = os.path.join(config.userlocaldir,
                                       self._fortran_sourcename)
        if os.path.exists(user_sfxflddesc):
            self.read(user_sfxflddesc)
        self.initialized = True

    def read(self, filename):
        """Parse a sfxflddesc_mod.F90 file."""
        with io.open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            m = self._pattern.match(line)
            if m:
                info = {}
                name = m.group('name').replace('.', '')
                info['gridtype'] = int(m.group('gridtype'))
                info['type'] = m.group('type')
                info['comment'] = m.group('comment').replace('.', ' ').strip()
                u = self._unit.match(info['comment'])
                if u:
                    info['unit'] = u.group('unit').replace('.', ' ')
                    info['comment'] = u.group('comment').strip()
                mask = m.group('mask').replace('.', '')
                if mask:
                    info['mask'] = mask
                self.table[name] = info
                self.table['S1D_' + name] = info
                self.table['SFX.' + name] = info

    def update(self, field_dict):
        """Update with a dict of field {fieldname:{}, ...}."""
        self.table.update(field_dict)

    @_init_before
    def nature(self, fieldname):
        """Return type of data in field."""
        return self.type2nature[self.table[fieldname]['type'][0]]

    @_init_before
    def dim(self, fieldname):
        """Return number of dimensions of field."""
        return int(self.table[fieldname]['type'][1])

    @_init_before
    def __getitem__(self, fieldname):
        return self.table[fieldname]

    @_init_before
    def get(self, fieldname, default=None):
        return self.table.get(fieldname, default)

    @_init_before
    def is_metadata(self, fieldname):
        """True if field is not a H2D field but metadata (Misc)."""
        if fieldname in self.table and self.table[fieldname]['type'] not in ('X1', 'X2'):
            return True
        else:
            return False

    @_init_before
    def __contains__(self, fieldname):
        return fieldname in self.table


def find_wind_pair(fieldname):
    """For a wind **fieldname**, find and return the pair."""
    pairs = {'U':'V', 'V':'U',
             'ZONAL':'MERIDIEN', 'MERIDIEN':'ZONAL', 'MERID':'ZONAL'}
    patterns = ['S\d+WIND\.(?P<d>[U,V])\.PHYS',
                '[P,H,V]\d+VENT_(?P<d>(ZONAL)|(MERID)|(MERIDIEN))',
                'CLSVENTNEUTRE.(?P<d>[U,V])',
                'CLSVENT.(?P<d>(ZONAL)|(MERIDIEN))',
                'CLS(?P<d>[U,V]).RAF.MOD.XFU']
    axes = {'U':'x', 'ZONAL':'x',
            'V':'y', 'MERIDIEN':'y', 'MERID':'y'}
    pair = None
    for pattern in patterns:
        re_ok = re.match(pattern, fieldname)
        if re_ok:
            pair = (axes[pairs[re_ok.group(1)]],
                    fieldname.replace(re_ok.group('d'),
                                      pairs[re_ok.group('d')]))
            break
    if pair is None:
        raise epygramError('not a wind field')
    else:
        return pair
