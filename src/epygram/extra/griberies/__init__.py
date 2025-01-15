#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains utilities around GRIB format.
"""

from bronx.syntax.parsing import str2dict


def parse_GRIBstr_todict(strfid):
    """Parse and return a dict GRIB fid from a string."""
    fid = str2dict(strfid, try_convert=int)
    return fid

from . import paths, tables, defaults, definitions
