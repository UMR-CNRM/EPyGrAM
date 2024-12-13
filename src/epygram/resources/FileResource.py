#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a FileResource.
"""

import os

import footprints
from bronx.syntax.decorators import nicedeco

from epygram.base import Resource
from epygram import epygramError, config

epylog = footprints.loggers.getLogger(__name__)


class FileResource(Resource):
    """Generic abstract class implementing a Resource."""

    _abstract = True
    _collector = ('epyresource', 'dataformat')
    _footprint = dict(
        attr=dict(
            format=dict(
                optional=True,
                info="Format of the resource."),
            filename=dict(
                info="File name (absolute or relative) of the resource."),
            fmtdelayedopen=dict(
                optional=True,
                default=False,
                type=bool,
                info="Opening of the resource delayed (not at time of " +
                     "construction).")
        )
    )

    @nicedeco
    def _openbeforedelayed(mtd):
        """Decorator for Resource: open resource before calling method."""

        def nowopen(self, *args, **kwargs):
            if not self.isopen:
                self.open()
            return mtd(self, *args, **kwargs)

        return nowopen

    _openbeforedelayed = staticmethod(_openbeforedelayed)

    def __init__(self, *args, **kwargs):
        super(Resource, self).__init__(*args, **kwargs)
        self.container = footprints.proxy.container(filename=self.filename)

        if self.openmode in ('r', 'a') and\
           not os.path.exists(self.container.abspath):
            raise IOError(self.container.abspath + " does not exist.")

        if self.openmode in ('r', 'a'):
            assert os.access(self.container.abspath, os.R_OK), \
                'No reading permission for file: ' + self.container.abspath

        # protection against unhappy overwrites...
        if config.protect_unhappy_writes and \
           os.path.exists(self.container.abspath) and self.openmode == 'w':
            self._overwrite = input(' '.join([self.container.abspath,
                                              "will be overwritten:",
                                              "do you want to continue",
                                              "(y/n) ? "])) == 'y'
            if not self._overwrite:
                raise epygramError(self.container.abspath + " already exists.")
        if self.openmode == 'a':
            assert os.access(self.container.abspath, os.W_OK), \
                'No writing permission for file: ' + self.container.abspath
        if self.openmode == 'w':
            if os.path.exists(self.container.abspath):
                assert os.access(self.container.abspath, os.W_OK), \
                    'No overwriting permission for file: ' + self.container.abspath
            else:
                dirpath = os.path.dirname(self.container.abspath)
                assert os.path.exists(dirpath), \
                       'Trying to write into non-existing directory: ' + dirpath
                assert os.access(dirpath, os.W_OK), \
                       'No writing permission in directory: ' + dirpath

    def open(self, openmode=None):
        """
        Opens the resource properly.

        :param openmode: to open with a specific openmode, eventually different from
        the one specified at initialization.
        """
        if openmode is not None:
            self._attributes['openmode'] = openmode


footprints.collectors.get(tag='epyresources').fasttrack = ('format',)
footprints.collectors.get(tag='dataformats').fasttrack = ('format',)
