#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities around GRIB definitions files.
"""

import os
import re
import io
import copy

from ...util import init_before


def read_gribdef(filename):
    """Read a grib definition file and return it as a dict."""
    re_name = re.compile(r'("|\')(?P<name>[\w\.\-\_ ]+)("|\')\s*=')
    re_real = r'\+|-?\d*\.\d*e\+|-?\d*'
    re_real_g = r'(?P<real>' + re_real + ')'
    re_int = r'\+|-?\d+'
    re_int_g = r'(?P<int>' + re_int + ')'
    # re_num = '(' + re_real + ')|(' + re_int + ')'
    re_num_g = re_int_g + r'|' + re_real_g
    re_keyvalue = re.compile(r'(?P<key>\w+)\s*=\s*' + re_num_g + r'\s*;')
    # read file
    with io.open(filename, 'r', encoding='utf-8') as f:
        lines_unfold = [l.strip() for l in f.readlines()]
        lines = []
        for line in lines_unfold:
            line = line.replace(';',';\n').replace('{ ','{\n')
            line_split = [l.strip() for l in line.split('\n')]
            lines.extend(line_split)
    # find fields declaration
    dico = {}
    indexes = []
    for i, line in enumerate(lines):
        fmatch = re_name.match(line)
        if fmatch:
            field = fmatch.group('name')
            indexes.append((field, i))
            dico[field] = {}
    # loop on fields
    for (j, (field, istart)) in enumerate(indexes):
        if j + 1 == len(indexes):  # last one
            iend = len(lines)
        else:
            iend = indexes[j + 1][1]
        if istart > 0 and lines[istart - 1].startswith("#"):  # this is a comment
            dico[field]['#comment'] = lines[istart - 1][1:].strip().strip('"')
        for i in range(istart, iend):
            kvmatch = re_keyvalue.match(lines[i])
            if kvmatch:
                if kvmatch.groupdict().get('int'):
                    dico[field][kvmatch.group('key')] = int(kvmatch.group('int'))
                elif kvmatch.groupdict().get('real'):
                    dico[field][kvmatch.group('key')] = float(kvmatch.group('real'))
    return dico


class GribDef(object):

    _default_grib_edition = 'grib2'
    _non_GRIB_keys = []

    def __init__(self, actual_init=True, concepts=[]):
        self._concepts = set(concepts)
        self.tables = {'grib1':{c:{} for c in concepts},
                       'grib2':{c:{} for c in concepts}}
        if actual_init:
            self._actual_init()
        else:
            self._initialized = False

    def _actual_init(self):
        """Read necessary definition files."""
        pass

    def read(self, filename, grib_edition=None):
        """Read a grib def concept file, and update or register it."""
        if grib_edition is None:
            for g in ('grib1', 'grib2'):
                if g in filename:
                    grib_edition = g
        concept = os.path.basename(filename).replace('.def', '')
        if concept in self.tables[grib_edition]:
            self.tables[grib_edition][concept].update(read_gribdef(filename))
        else:
            self.tables[grib_edition][concept] = read_gribdef(filename)

    @classmethod
    def _filter_non_GRIB_keys(cls, fid):
        """Some keys are to be filtered out from gribdef."""
        filtered_out = copy.copy(fid)
        for k in cls._non_GRIB_keys:
            if k in fid:
                filtered_out.pop(k)
        return filtered_out

    @init_before
    def _lookup(self, fid, concept,
                grib_edition=_default_grib_edition,
                include_comments=False,
                filter_non_GRIB_keys=True,
                exact=False):
        """
        Concept equivalence lookup:
          - if **fid** is a **concept** value, get the associated GRIB key/value pairs
          - if **fid is a set of GRIB key/value pairs, get the associated **concept** value(s)

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        :param exact: when **fid** is a concept value,
            - if exact=True, return only the grib def which concept value
            matches exactly **fid**
            - if exact=False, return all the grib def which concept value
            contains **fid**
        """
        try:
            if isinstance(fid, str):
                fid = parse_GRIBstr_todict(fid)
        except SyntaxError:  # fid is a concept value
            retrieved = self._lookup_from_conceptvalue(fid, concept,
                                                       grib_edition=grib_edition,
                                                       include_comments=include_comments,
                                                       filter_non_GRIB_keys=filter_non_GRIB_keys,
                                                       exact=exact)
        else:  # fid is a GRIB fid
            retrieved = self._lookup_from_kv(fid, concept,
                                             grib_edition=grib_edition,
                                             include_comments=include_comments,
                                             filter_non_GRIB_keys=filter_non_GRIB_keys)
        return retrieved

    @init_before
    def _get_def(self, fid, concept,
                 grib_edition=_default_grib_edition,
                 include_comments=False,
                 fatal=True,
                 filter_non_GRIB_keys=True):
        """
        Direct access to a grib definition

        :param include_comments: if a comment is present if grib def, bring it in fid
        :param fatal: if True and fieldname is not retrieved, raise a ValueError;
            else, get 'default' instead of fieldname
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        """
        if fid in self.tables[grib_edition][concept]:
            fid = copy.copy(self.tables[grib_edition][concept][fid])
        else:
            if fatal:
                raise KeyError('{} not found'.format(fid))
            else:
                fid = copy.copy(self.tables[grib_edition][concept]['default'])
        if not include_comments:
            if '#comment' in fid:
                fid.pop('#comment', None)
        if filter_non_GRIB_keys:
            fid = self._filter_non_GRIB_keys(fid)
        return fid

    def _lookup_from_conceptvalue(self, fid, concept,
                                  grib_edition=_default_grib_edition,
                                  include_comments=False,
                                  filter_non_GRIB_keys=True,
                                  exact=False):
        """
        Look for all the fields which **concept** value contain **fid**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        :param exact: when **fid** is a concept value,
            - if exact=True, return only the grib def which concept value
            matches exactly **fid**
            - if exact=False, return all the grib def which concept value
            contains **fid**
        """
        fields = {}
        try:  # first try to get it as exact
            gribfid = self._get_def(fid, concept,
                                    grib_edition=grib_edition,
                                    include_comments=include_comments,
                                    fatal=True,
                                    filter_non_GRIB_keys=filter_non_GRIB_keys)
        except KeyError:
            pass  # field was not found as such; might be partial => finally
        else:
            fields[fid] = gribfid
        finally:
            if not exact:  # complete with all that contains fid
                for f, gribfid in self.tables[grib_edition][concept].items():
                    if fid in f:
                        fields[f] = gribfid
                for gribfid in fields.values():
                    if not include_comments:
                        gribfid.pop('#comment', None)
                    if filter_non_GRIB_keys:
                        gribfid = self._filter_non_GRIB_keys(gribfid)
            else:
                if len(fields) == 1:
                    fields = fields[list(fields.keys())[0]]
        return fields

    @init_before
    def _lookup_from_kv(self, handgrip, concept,
                        grib_edition=_default_grib_edition,
                        include_comments=False,
                        filter_non_GRIB_keys=True):
        """
        Look for all the fields which GRIB def contain **handgrip**,
        returning them identified by their **concept**.

        :param grib_edition: among ('grib1', 'grib2'), the version of GRIB fid
        :param include_comments: if a comment is present if grib def, bring it in fid
        :param filter_non_GRIB_keys: filter out the non-GRIB keys that may be
            present in grib def of field
        """
        if isinstance(handgrip, str):
            handgrip = parse_GRIBstr_todict(handgrip)
        fields = {}
        for f, gribfid in self.tables[grib_edition][concept].items():
            if set(handgrip.keys()).issubset(gribfid.keys()):
                ok = True
                for k,v in handgrip.items():
                    if gribfid[k] != v:
                        ok = False
                        break
                if ok:
                    fields[f] = copy.copy(gribfid)
        for k, fid in fields.items():
            if not include_comments:
                fid.pop('#comment', None)
            if filter_non_GRIB_keys:
                fields[k] = self._filter_non_GRIB_keys(fid)
        return fields

    @init_before
    def known_names_for_concept(self, concept,
                                grib_edition=_default_grib_edition):
        """Get sorted list of names for **concept**."""
        return sorted(self.tables[grib_edition][concept].keys())

    @init_before
    def known_names(self, grib_edition=_default_grib_edition):
        """Get sorted list of names for all concepts."""
        all_names = {}
        concepts = list(self.tables[grib_edition].keys())
        for c in concepts:
            all_names[c] = self.known_names_for_concept(c, grib_edition)
        return all_names

    @init_before
    def known_values_for(self, key,
                         concept=None,
                         grib_edition=_default_grib_edition):
        """Get list of all values present throughout the GribDef for **key**."""
        if concept is None:
            concepts = list(self.tables[grib_edition].keys())
        else:
            concepts = [concept]
        values = set()
        for c in concepts:
            for gribfid in self.tables[grib_edition][c].values():
                if key in gribfid:
                    values.add(gribfid[key])
        return sorted(values)

    def known_values(self,
                     concept=None,
                     grib_edition=_default_grib_edition):
        """Get list of all values present throughout the GribDef for all keys."""
        values = {}
        for key in self._allkeys(grib_edition):
            values[key] = self.known_values_for(key, concept, grib_edition)
        return values

    @init_before
    def _allkeys(self, grib_edition=_default_grib_edition):
        """Get set of all keys present throughout the GribDef."""
        all_keys = set()
        concepts = list(self.tables[grib_edition].keys())
        for c in concepts:
            for gribfid in self.tables[grib_edition][c].values():
                all_keys.update(set(gribfid.keys()))
        if '#comment' in all_keys:
            all_keys.remove('#comment')
        return all_keys

