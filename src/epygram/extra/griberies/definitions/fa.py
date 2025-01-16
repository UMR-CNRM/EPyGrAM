#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
FA Grib definitions
"""

import os
import re

from .util import GribDef
from .. import parse_GRIBstr_todict
from ..paths import get_definition_paths
from ...util import init_before
from epygram import config


class FaGribDef(GribDef):
    """
    Handle FA-related GRIB definition files.
    To add user files:
    a) use env var ECCODES_DEFINITION_PATH or
    b) move it under {config.userlocaldir}/{FaGribDef.dirname}
    """

    _re_dirname = re.compile(r'gribapi\.def\.[\d_]+')
    _official_rootdir = os.path.join(config.installdir, 'data')  # TODO: use data
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
        defpaths = get_definition_paths()
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
                           grib_edition=GribDef._default_grib_edition):
        pathname = os.path.join(directory, grib_edition, 'localConcepts', 'lfpw', concept + '.def')
        if os.path.exists(pathname):
            self.read(pathname, grib_edition)

    @init_before
    def FA2GRIB(self, fieldname,
                grib_edition=GribDef._default_grib_edition,
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

    @init_before
    def GRIB2FA(self, gribfid,
                grib_edition=GribDef._default_grib_edition):
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

    @init_before
    def lookup_FA(self, partial_fieldname,
                  grib_edition=GribDef._default_grib_edition,
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

    @init_before
    def lookup_GRIB(self, partial_fid,
                    grib_edition=GribDef._default_grib_edition,
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
            if isinstance(fid, str):
                fid = parse_GRIBstr_todict(fid)
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

fagribdef = FaGribDef(actual_init=False)

