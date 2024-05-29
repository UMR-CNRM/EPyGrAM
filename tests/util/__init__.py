#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
from __future__ import print_function, absolute_import, division, unicode_literals
import os

file_dir = os.path.dirname(os.path.abspath(__file__))

datadir = os.path.join(file_dir, '..', 'data')

datadir = os.path.normpath(datadir)

delta_assertAlmostEqual = 1e-12
delta_assertAlmostEqual4pyproj = 1e-7  # FIXME: because pyproj not reproducible between python2 and python3
