#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An interface to ask GRIB definitions."""

import argparse

import epygram
from epygram.extra.griberies import parse_GRIBstr_todict
from epygram.extra.griberies.definitions.fa import FaGribDef
from epygram.extra.griberies.definitions.names import NamesGribDef
from epygram.extra.griberies.concepts import Concept
from . import epilog

_description = __doc__

namesgribdef = NamesGribDef()
fagribdef = FaGribDef()
concepts = sorted(list(namesgribdef._concepts) + list(fagribdef._concepts))


def main():
    args = get_args()
    what_the_grib(args.fid,
                  args.concept,
                  args.grib_edition,
                  args.exact,
                  args.format_as_dict)


def what_the_grib(fid, concept, grib_edition, exact, format_as_dict):
    grib_edition = 'grib{}'.format(grib_edition)
    try:
        if isinstance(fid, str):
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


def get_args():

    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
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
    return parser.parse_args()

