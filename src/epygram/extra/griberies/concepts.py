#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
A class for concepts handling.
"""

from .definitions.names import NamesGribDef
from .definitions.fa import FaGribDef
namesgribdef = NamesGribDef()
fagribdef = FaGribDef()
concepts = sorted(list(namesgribdef._concepts) + list(fagribdef._concepts))


class Concept(object):
    def __init__(self, concept):
        self.concept = concept
        if concept in fagribdef._concepts:
            self.gribdef = fagribdef
        else:
            self.gribdef = namesgribdef

    def what(self, fid, grib_edition, exact, format_as_dict=False):
        self.header()
        print()
        fields = self.fields(fid, grib_edition, exact)
        first = True
        for f, fid in fields.items():
            self.niceprint(f, fid, format_as_dict, headline=not first)
            first = False
        print()

    def fields(self, fid, grib_edition, exact):
        if self.concept == 'faFieldName':
            if isinstance(fid, str):
                if exact:
                    try:
                        fields = {fid:self.gribdef.FA2GRIB(fid, fatal=True,
                                                           grib_edition=grib_edition,
                                                           include_comments=True,
                                                           filter_non_GRIB_keys=False)}
                    except ValueError:
                        fields = {}
                else:
                    fields = self.gribdef.lookup_FA(fid,
                                                    include_comments=True,
                                                    grib_edition=grib_edition,
                                                    filter_non_GRIB_keys=False)
            else:
                fields = self.gribdef.lookup_GRIB(fid,
                                                  include_comments=True,
                                                  filter_non_GRIB_keys=False)
        else:
            fields = self.gribdef._lookup(fid, self.concept,
                                          grib_edition=grib_edition,
                                          include_comments=True,
                                          filter_non_GRIB_keys=False,
                                          exact=exact)
        return fields

    def header(self):
        print('{:-^80}'.format('| {} |'.format(self.concept)))

    def niceprint(self, f, fid, format_as_dict=False, headline=True):
        comment = fid.pop('#comment', None)
        filtered_out = self.gribdef._filter_non_GRIB_keys(fid)
        fprint = '> "{}"'.format(f)
        if headline:
            print('-' * len(fprint))
        print(fprint)
        print('=' * len(fprint))
        if comment is not None:
            print("# " + comment)
        if format_as_dict:
            print(fid)
        else:
            for k in sorted(fid.keys()):
                print('{} = {}'.format(k, fid[k]))
        for k in sorted(filtered_out.keys()):
            print('({} = {})'.format(k, filtered_out[k]))
