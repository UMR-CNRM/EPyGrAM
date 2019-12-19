#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains all fields classes.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import glob
import importlib
import os

import footprints

epylog = footprints.loggers.getLogger(__name__)

for module in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'with_*')):
    name, ext = os.path.splitext(module)
    name = os.path.basename(name)
    if ext in ('', '.py'):
        try:
            importlib.import_module('.' + name, __name__)
        except ImportError as e:
            epylog.warning("An error ('" + str(e) + "') occurred when importing the " +
                           name + " plugin; " +
                           "this is certainly due to a missing dependency, some " +
                           "functionalities can be missing.")
