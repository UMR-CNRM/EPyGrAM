#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to for Miscellaneous fields.
"""

from epygram.base import Field


class MiscField(Field):
    """
    A miscellaneous field class.
    Designed to handle data of various nature (bool, int, float, str...)
    and various dimensions (scalar, 1D/2D).
    """

    @property
    def datatype(self):
        """Returns the data type."""
        return self._data.dtype

    @property
    def shape(self):
        """Returns the data shape."""
        return self._data.shape
