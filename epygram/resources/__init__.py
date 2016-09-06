#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Module contains:

- resource for data contained in file (FileResource).
- other resources built on top of FileResource:\n
  - MultiValiditiesResource: join several resources to furnish fields with temporal evolution;
"""

from .FileResource import FileResource
from .MultiValiditiesResource import MultiValiditiesResource, timeresource
from .CombineLevelsResource import CombineLevelsResource
