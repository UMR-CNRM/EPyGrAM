#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Module contains:

- resource for data contained in file (FileResource).
- other resources built on top of FileResource:\n
  - MultiValiditiesResource: join several resources to furnish fields with
    temporal evolution;
  - CombineLevelsResource: from a resource containing 2D fields on adjacent
    levels, emulates a resource that provide 3D fields
- a proxy function to build such special resources 
"""

from .FileResource import FileResource
from .MultiValiditiesResource import MultiValiditiesResource
from .CombineLevelsResource import CombineLevelsResource



def special_resource(filenames, openmode, rtype):
    """
    Factory for special resources, such as MultiValiditiesResource or
    CombineLevelsResource.
    
    *filenames* can be either a filename or a list of
    *openmode*: among 'r', 'w', 'a'
    *rtype': resource type, e.g.:
             'MV' for a MultiValiditiesResource,
             'CL' for a CombineLevelsResource
             'MV+CL' for a composition of both (should be similar to CL+MV)
    """
    from epygram.formats import resource
    from footprints import proxy as fpx

    if '+' in rtype and len(rtype.split('+')) > 2:
        raise NotImplementedError('more than one composition in *rtype*.')

    if rtype.startswith('MV'):
        # step 1
        if isinstance(filenames, str):
            filenames = [filenames]
        resources = [resource(f, openmode, fmtdelayedopen=True) for f in filenames]
        resource = fpx.epyresource(resources=resources,
                                   openmode=openmode,
                                   name='MultiValidities')
        # step 2
        if rtype.endswith('CL'):
            special_resource = fpx.epyresource(resource=resource,
                                               openmode=openmode,
                                               name='CombineLevels')
        else:
            special_resource = resource
    elif rtype.startswith('CL'):
        # step 1
        if isinstance(filenames, list) and len(filenames) > 1:
            assert rtype.endswith('MV'), \
                   '*filenames* should be unique with this *rtype*:' + rtype
        elif isinstance(filenames, str):
            filenames = [filenames]
        resources = [resource(f, openmode, fmtdelayedopen=True) for f in filenames]
        resources = [fpx.epyresource(resource=r,
                                     openmode=openmode,
                                     name='CombineLevels') for r in resources]
        # step 2
        if rtype.endswith('MV'):
            special_resource = fpx.epyresource(resources=resources,
                                               openmode=openmode,
                                               name='MultiValidities')
        else:
            special_resource = resources[0]

    return special_resource
