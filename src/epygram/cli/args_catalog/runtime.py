#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with runtime options
"""

import multiprocessing
from . import _defaults

d = {
    'verbose':[
        '-v', '--verbose',
        dict(action='store_true',
             help="run verbosely. Else, only messages of level Error will be\
                   displayed.",
             default=False)],
    'percentage':[
        '-p', '--percentage',
        dict(action='store_true',
             help="display the percentage done on the run.",
             default=False)],
    'threads_number':[
        '-t', '--threads_number',
        dict(help="number of threads to be run in parallel.",
             type=int,
             default=multiprocessing.cpu_count() // 2)],
                    }

