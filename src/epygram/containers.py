#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the containers for a Resource.
"""

import os

from footprints import FootprintBase

from epygram.util import RecursiveObject


class File(RecursiveObject, FootprintBase):
    """Generic class implementing a File."""

    _collector = ('container', 'file')
    _footprint = dict(
        attr=dict(
            filename=dict(
                alias=('f',),
                info='Relative or absolute pathname.')
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""
        # super to footprints constructor
        super(File, self).__init__(*args, **kwargs)
        # initialise absolute pathname
        self._abspath = os.path.abspath(self.filename)

    @property
    def basename(self):
        """Returns the basename of the file."""
        return os.path.basename(self._abspath)

    @property
    def abspath(self):
        """Returns the absolute path of the file."""
        return self._abspath

    @property
    def absdir(self):
        """Returns the absolute path of the directory of the file."""
        return self._abspath[:-len(self.basename)]
