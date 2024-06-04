#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals
import os
# import os
# outside data...
# datadir = os.path.join(os.getenv('TMPDIR',
#                                 os.path.join(os.getenv('HOME'), 'tmp')),
#                       'epygram_test_data')
# or data stored within repository
file_dir = os.path.dirname(__file__)
datadir = os.path.join(file_dir, 'data')
#datadir = './data'
