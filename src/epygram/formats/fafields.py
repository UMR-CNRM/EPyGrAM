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
import re
import six
import io

from epygram.extra import griberies
from epygram import config, epygramError


class FaGribDef(griberies.GribDef):
    """
    Handle FA-related GRIB definition files.
    To add user files:
    a) use env vars GRIBAPI|ECCODES_DEFINITION_PATH or
    b) move it under {config.userlocaldir}/{FaGribDef.dirname}
    """

    _re_dirname = re.compile(r'gribapi\.def\.[\d_]+')
    _official_rootdir = os.path.join(config.installdir, 'data')
    _non_GRIB_keys = ('LSTCUM', 'FMULTM', 'FMULTE', 'ZLMULT')

    def __init__(self, actual_init=True,
                 concepts=['faFieldName', 'faModelName', 'faLevelName']):
        super(FaGribDef, self).__init__(actual_init, concepts)

    def _find_gribdefs(self, root):
        dirs = []
        if os.path.exists(root) and os.path.isdir(root):
            for d in sorted(os.listdir(root)):
                if self._re_dirname.match(d):
                    d = os.path.join(root, d)
                    if os.path.isdir(d):
                        dirs.append(d)
        return dirs

    def _actual_init(self):
        """Read definition files."""
        # first, those specified by env var
        defpaths = griberies.get_definition_paths()
        # then those embarked with epygram
        defpaths = self._find_gribdefs(self._official_rootdir) + defpaths
        # complete with local user definitions, if existing
        defpaths = self._find_gribdefs(config.userlocaldir) + defpaths
        # read gribdef files
        for d in defpaths[::-1]:
            for grib_edition in ('grib1', 'grib2'):
                for concept in self._concepts:
                    self._read_localConcept(concept, d, grib_edition)
        self._initialized = True

    def _read_localConcept(self, concept, directory,
                           grib_edition=griberies.GribDef._default_grib_edition):
        pathname = os.path.join(directory, grib_edition, 'localConcepts', 'lfpw', concept + '.def')
        if os.path.exists(pathname):
            self.read(pathname, grib_edition)

    @griberies.init_before
    def FA2GRIB(self, fieldname,
                grib_edition=griberies.GribDef._default_grib_edition,
                include_comments=False,
                fatal=False,
                filter_non_GRIB_keys=True):
        """
        Convert FA field name to GRIB fid, including vertical level identification.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param fatal: if True and fieldname is not retrieved, raise a ValueError;
            else, return a default 255 fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        """
        _re_altitude = re.compile(r'(?P<ltype>[A-W]+)(?P<level>\d+)(?P<param>.+)')
        if fieldname in self.tables[grib_edition]['faFieldName']:
            fid = self._get_def(fieldname, 'faFieldName',
                                grib_edition, include_comments)
        elif fieldname.replace(' ', '_') in self.tables[grib_edition]['faFieldName']:  # FIXME: ?
            fid = self._get_def(fieldname.replace(' ', '_'), 'faFieldName',
                                grib_edition, include_comments)
        else:
            rematch = _re_altitude.match(fieldname)
            if rematch:
                fid = self._get_def(rematch.group('ltype'), 'faLevelName',
                                    grib_edition, include_comments,
                                    fatal=False)
                fid.update(self._get_def(rematch.group('param'), 'faFieldName',
                                         grib_edition, include_comments,
                                         fatal=False))
                level = int(rematch.group('level'))
                if level == 0:  # formatting issue
                    level = 100000
                fid['scaleFactorOfFirstFixedSurface'] = 0
                fid['scaledValueOfFirstFixedSurface'] = level
                if rematch.group('ltype') == 'P':
                    level /= 100
                fid['level'] = level
            else:
                if fatal:
                    raise ValueError("field not found: {}".format(fieldname))
                else:
                    fid = self._get_def('default', 'faFieldName',
                                        grib_edition, include_comments)
                    fid.setdefault('typeOfFirstFixedSurface', 255)
        if filter_non_GRIB_keys:
            fid = self._filter_non_GRIB_keys(fid)
        # productDefinitionTemplateNumber to distinguish between
        # CLSTEMPERATURE and CLSMINI.TEMPERAT / CLSMAXI.TEMPERAT (for instance)
        if 'productDefinitionTemplateNumber' not in fid:
            fid['productDefinitionTemplateNumber'] = 0
        return fid

    @griberies.init_before
    def GRIB2FA(self, gribfid,
                grib_edition=griberies.GribDef._default_grib_edition):
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

    @griberies.init_before
    def lookup_FA(self, partial_fieldname,
                  grib_edition=griberies.GribDef._default_grib_edition,
                  include_comments=False,
                  filter_non_GRIB_keys=True):
        """
        Look for all the fields which FA name contain **partial_fieldname**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        """
        fields = {}
        try:  # first try to get it as an altitude one
            gribfid = self.FA2GRIB(partial_fieldname,
                                   grib_edition=grib_edition,
                                   include_comments=include_comments,
                                   fatal=True,
                                   filter_non_GRIB_keys=filter_non_GRIB_keys)
        except ValueError:
            pass  # field was not found as such; might be partial => finally
        else:
            fields[partial_fieldname] = gribfid
        finally:
            for f, gribfid in self.tables[grib_edition]['faFieldName'].items():
                if partial_fieldname in f:
                    fields[f] = gribfid
            for k, fid in fields.items():
                if not include_comments:
                    fid.pop('#comment', None)
                if filter_non_GRIB_keys:
                    fields[k] = self._filter_non_GRIB_keys(fid)
        return fields

    @griberies.init_before
    def lookup_GRIB(self, partial_fid,
                    grib_edition=griberies.GribDef._default_grib_edition,
                    include_comments=False,
                    filter_non_GRIB_keys=True):
        """
        Look for all the fields which GRIB fid contain **partial_fid**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        """
        return self._lookup_from_kv(partial_fid, 'faFieldName',
                                    grib_edition=grib_edition,
                                    include_comments=include_comments,
                                    filter_non_GRIB_keys=filter_non_GRIB_keys)

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
    _pattern = re.compile(r'"\.(?P<name>[\w\d%]+)\.+(?P<gridtype>\d)\.{2}(?P<type>\w\d)\.{2}(?P<comment>.+)\.{2}(?P<mask>.{20})\.{2}",? ?(&|/)')
    _unit = re.compile(r'(?P<comment>.*)\((?P<unit>.+)\)')
    _fortran_sourcename = 'sfxflddesc_mod.F90'
    type2nature = {'X':'float', 'L':'bool', 'C':'str', 'N':'int', 'Y':'float',
                   'T':'int'}  # FIXME: T: cf mode_write_surf_fa.F90

    def __init__(self, actual_init=True):
        self.table = {}
        if actual_init:
            self._actual_init()
        else:
            self._initialized = False

    def _actual_init(self):
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
        self._initialized = True

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

    @griberies.init_before
    def nature(self, fieldname, default=None):
        """Return type of data in field."""
        try:
            return self.type2nature[self.table[fieldname]['type'][0]]
        except KeyError:
            if default is None:
                raise
            else:
                return default

    @griberies.init_before
    def dim(self, fieldname, default=None):
        """Return number of dimensions of field."""
        try:
            return int(self.table[fieldname]['type'][1])
        except KeyError:
            if default is None:
                raise
            else:
                return default

    @griberies.init_before
    def __getitem__(self, fieldname):
        return self.table[fieldname]

    @griberies.init_before
    def get(self, fieldname, default=None):
        return self.table.get(fieldname, default)

    @griberies.init_before
    def is_metadata(self, fieldname):
        """True if field is not a H2D field but metadata (Misc)."""
        if fieldname in self.table:
            if self.table[fieldname]['type'] in ('X1', 'X2'):
                return False
            else:
                return True
        elif '_FBUF_' in fieldname:
            return True
        else:
            return None

    @griberies.init_before
    def __contains__(self, fieldname):
        return fieldname in self.table


def find_wind_pair(fieldname):
    """For a wind **fieldname**, find and return the pair."""
    pairs = {'U':'V', 'V':'U',
             'ZONAL':'MERIDIEN', 'MERIDIEN':'ZONAL', 'MERID':'ZONAL'}
    patterns = [r'S\d+WIND\.(?P<d>[U,V])\.PHYS',
                r'[P,H,V]\d+VENT_(?P<d>(ZONAL)|(MERID)|(MERIDIEN))',
                r'CLSVENTNEUTRE.(?P<d>[U,V])',
                r'CLSVENT.(?P<d>(ZONAL)|(MERIDIEN))',
                r'CLS(?P<d>[U,V]).RAF.MOD.XFU']
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
