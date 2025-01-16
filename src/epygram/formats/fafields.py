#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
FA fields utilities
"""

import os
import re
import io

from epygram.extra.util import init_before
from epygram import config, epygramError
from epygram.extra.griberies.definitions.fa import fagribdef


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


def get_generic_fid(fieldname):
    """Return a generic fid from **fieldname** (via FaGribDef)."""
    try:
        fid = fagribdef.FA2GRIB(fieldname,
                                include_comments=False,
                                fatal=True)
    except ValueError:  # not found
        if sfxflddesc.is_metadata(fieldname):
            raise
        else:
            fid = fagribdef.FA2GRIB(fieldname,
                                    include_comments=False,
                                    fatal=False)
    return fid


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

    @init_before
    def nature(self, fieldname, default=None):
        """Return type of data in field."""
        try:
            return self.type2nature[self.table[fieldname]['type'][0]]
        except KeyError:
            if default is None:
                raise
            else:
                return default

    @init_before
    def dim(self, fieldname, default=None):
        """Return number of dimensions of field."""
        try:
            return int(self.table[fieldname]['type'][1])
        except KeyError:
            if default is None:
                raise
            else:
                return default

    @init_before
    def __getitem__(self, fieldname):
        return self.table[fieldname]

    @init_before
    def get(self, fieldname, default=None):
        return self.table.get(fieldname, default)

    @init_before
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

    @init_before
    def __contains__(self, fieldname):
        return fieldname in self.table


sfxflddesc = SfxFldDesc_Mod(actual_init=False)

