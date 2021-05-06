#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
from __future__ import print_function, absolute_import, unicode_literals, division
import six

import os
import sys
import argparse

# Automatically set the python path
package_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, package_path)
sys.path.insert(0, os.path.join(package_path, 'site'))
from griberies import parse_GRIBstr_todict
import epygram

namesgribdef = epygram.formats.GRIB.NamesGribDef()
fagribdef = epygram.formats.fafields.FaGribDef()
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
            if isinstance(fid, six.string_types):
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


def main(fid, concept, grib_edition, exact, format_as_dict):
    grib_edition = 'grib{}'.format(grib_edition)
    try:
        if isinstance(fid, six.string_types):
            fid = parse_GRIBstr_todict(fid)
    except SyntaxError:  # fid is a concept value
        if concept is None:
            raise epygram.epygramError("'concept' arg must be given, as fid is not a GRIB fid")
        c = Concept(concept)
        c.what(fid, grib_edition, exact, format_as_dict)
    else:  # fid is a GRIB fid
        if concept is None:
            print_concepts = sorted(concepts)
        else:
            print_concepts = [concept]
        for c in print_concepts:
            c = Concept(c)
            c.what(fid, grib_edition, exact, format_as_dict)


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description='An interface to ask GRIB definitions.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')
    parser.add_argument('fid',
                        help='Complete or partial GRIB fid or fieldname in any concept.')
    parser.add_argument('-c', '--concept',
                        dest='concept',
                        default=None,
                        help='concept, among:({}). Necessary in case fid is a fieldname, optional otherwise.'.format(concepts))
    parser.add_argument('-g', '--grib_edition',
                        dest='grib_edition',
                        default=2,
                        help='1 or 2')
    parser.add_argument('-e', '--exact',
                        action='store_true',
                        default=False,
                        help='Exact: only print fields that match exactly the GRIB fid.',
                        dest='exact')
    parser.add_argument('-d', '--format_as_dict',
                        action='store_true',
                        default=False,
                        help='Format GRIB fids as dict.',
                        dest='format_as_dict')
    args = parser.parse_args()

    # 2. Initializations
    ####################

    # 3. Main
    #########
    main(args.fid,
         args.concept,
         args.grib_edition,
         args.exact,
         args.format_as_dict)
