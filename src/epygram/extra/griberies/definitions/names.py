#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Handle *name* and *shortName* GRIB definitions.
"""

import os
from .util import GribDef
from ..paths import get_definition_paths


class NamesGribDef(GribDef):
    """Handle *name* and *shortName* GRIB definitions."""

    _non_GRIB_keys = ['is_uerra']

    def __init__(self, actual_init=True,
                 concepts=['name', 'shortName', 'cfName', 'cfVarName']):
        super(NamesGribDef, self).__init__(actual_init, concepts)

    def _actual_init(self):
        """Read definition files."""
        # get definition paths, from env variable
        defpaths = get_definition_paths()
        for d in defpaths[::-1]:
            for grib_edition in ('grib1', 'grib2'):
                for concept in self._concepts:
                    self._readConcept(concept, d, grib_edition)
        self._initialized = True

    def _readConcept(self, concept, directory,
                     grib_edition=GribDef._default_grib_edition):
        pathname = os.path.join(directory, grib_edition, concept + '.def')
        if os.path.exists(pathname):
            self.read(pathname, grib_edition)

    def name(self, fid,
             grib_edition=GribDef._default_grib_edition,
             include_comments=False,
             filter_non_GRIB_keys=True,
             exact=True):
        """
        'name' equivalence lookup:
          - if **fid** is a name, get the associated GRIB key/value pairs
          - if **fid** is a set of GRIB key/value pairs, get the associated name(s)

        Cf. method _lookup() for other optional arguments.
        """
        return self._lookup(fid, 'name',
                            grib_edition=grib_edition,
                            include_comments=include_comments,
                            filter_non_GRIB_keys=filter_non_GRIB_keys,
                            exact=exact)

    def shortName(self, fid,
                  grib_edition=GribDef._default_grib_edition,
                  include_comments=False,
                  filter_non_GRIB_keys=True,
                  exact=True):
        """
        'name' equivalence lookup:
          - if **fid** is a shortName, get the associated GRIB key/value pairs
          - if **fid** is a set of GRIB key/value pairs, get the associated shortName(s)

        Cf. method _lookup() for other optional arguments.
        """
        return self._lookup(fid, 'shortName',
                            grib_edition=grib_edition,
                            include_comments=include_comments,
                            filter_non_GRIB_keys=filter_non_GRIB_keys,
                            exact=exact)

    def cfName(self, fid,
               grib_edition=GribDef._default_grib_edition,
               include_comments=False,
               filter_non_GRIB_keys=True,
               exact=True):
        """
        'name' equivalence lookup:
          - if **fid** is a cfName, get the associated GRIB key/value pairs
          - if **fid** is a set of GRIB key/value pairs, get the associated cfName(s)

        Cf. method _lookup() for other optional arguments.
        """
        return self._lookup(fid, 'cfName',
                            grib_edition=grib_edition,
                            include_comments=include_comments,
                            filter_non_GRIB_keys=filter_non_GRIB_keys,
                            exact=exact)

    def cfVarName(self, fid,
                  grib_edition=GribDef._default_grib_edition,
                  include_comments=False,
                  filter_non_GRIB_keys=True,
                  exact=True):
        """
        'name' equivalence lookup:
          - if **fid** is a cfVarName, get the associated GRIB key/value pairs
          - if **fid** is a set of GRIB key/value pairs, get the associated cfVarName(s)

        Cf. method _lookup() for other optional arguments.
        """
        return self._lookup(fid, 'cfVarName',
                            grib_edition=grib_edition,
                            include_comments=include_comments,
                            filter_non_GRIB_keys=filter_non_GRIB_keys,
                            exact=exact)
