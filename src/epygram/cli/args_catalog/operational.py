#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Operational arguments
"""

from . import _defaults
from epygram.extra import griberies

d = {
    'suite':[
        '-S', '--suite',
        dict(help="name of the suite, among (oper, dble, research, test).\
                   Defaults to research. Used to build GRIB's \
                   *productionStatusOfProcessedData*.",
             default='research')],
    'typeOfGeneratingProcess':[
        '-g', '--typeOfGeneratingProcess',
        dict(help="GRIB's type of generating process.",
             choices=list(griberies.tables.typeOfGeneratingProcess_dict.keys()),
             default='Forecast')],
    'numod':[
        '-N', '--NUMOD',
        dict(help="model identifier (known as NUMOD at Meteo-France). \
                   A.k.a. 'generatingProcessIdentifier' in GRIB_API. \
                   Default is 255.",
             dest='numod',
             default=None,
             type=int)],
                       }

