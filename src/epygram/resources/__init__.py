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
- a proxy function to build such meta resources
"""

from footprints import proxy as fpx

from epygram import epygramError
from epygram.base import Resource


class open_and_close_resource(object):
    """Context manager for automatically open/close resources."""

    def __init__(self, r):
        """
        :param r: resource to open
        """
        self.r = r

    def __enter__(self):
        self.r.open()

    def __exit__(self, t, v, tbk):
        self.r.close()


from .FileResource import FileResource
from .MultiValiditiesResource import MultiValiditiesResource
from .CombineLevelsResource import CombineLevelsResource
from .SubdomainResource import SubdomainResource
from .DiagnosticsResource import DiagnosticsResource


def meta_resource(filenames_or_resources, openmode, rtype):
    """
    Factory for meta resources, such as MultiValiditiesResource or
    CombineLevelsResource.

    :param filenames_or_resources: can be either a filename or a list of,
                                   or a resource or a list of.
    :param openmode: among 'r', 'w', 'a'
    :param rtype: resource type, e.g.:\n
                  - 'MV' for a MultiValiditiesResource,
                  - 'CL' for a CombineLevelsResource
                  - 'MV+CL' for a composition of both (should be similar to CL+MV)
    """
    from epygram.formats import resource

    if '+' in rtype and len(rtype.split('+')) > 2:
        raise NotImplementedError('more than one composition in *rtype*.')

    if not isinstance(filenames_or_resources, list):
        filenames_or_resources = [filenames_or_resources]
    if isinstance(filenames_or_resources[0], str):
        resources = [resource(f, openmode, fmtdelayedopen=True) for f in filenames_or_resources]
    elif isinstance(filenames_or_resources[0], Resource):
        resources = filenames_or_resources
    else:
        raise epygramError('unknown type for *filenames_or_resources*.')

    if rtype.startswith('MV'):
        # step 1
        resource = fpx.epyresource(resources=resources,
                                   openmode=openmode,
                                   name='MultiValidities')
        # step 2
        if rtype.endswith('CL'):
            meta_resource = fpx.epyresource(resource=resource,
                                            openmode=openmode,
                                            name='CombineLevels')
        else:
            meta_resource = resource
    elif rtype.startswith('CL'):
        # step 1
        if len(resources) > 1:
            assert rtype.endswith('MV'), \
                '*filenames* should be unique with this *rtype*:' + rtype
        resources = [fpx.epyresource(resource=r,
                                     openmode=openmode,
                                     name='CombineLevels') for r in resources]
        # step 2
        if rtype.endswith('MV'):
            meta_resource = fpx.epyresource(resources=resources,
                                            openmode=openmode,
                                            name='MultiValidities')
        else:
            meta_resource = resources[0]

    return meta_resource
