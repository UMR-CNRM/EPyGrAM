#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to handle LFA format.
"""

import os
import numpy
import sys

from falfilfa4py import LFA as LFA4py

import epygram
from epygram import epygramError, config, util
from epygram.resources import FileResource
from epygram.fields import MiscField

__all__ = ['LFA']


class LFA(FileResource):
    """Class implementing all specificities for LFA resource format."""

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['LFA']),
                default='LFA')
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        super(LFA, self).__init__(*args, **kwargs)
        if not LFA4py.wlfatest(self.container.abspath):
            raise IOError("This resource is not a LFA one.")
        if not self.fmtdelayedopen:
            self.open()

    def open(self, openmode=None):
        """
        Opens the LFA in Fortran sense.

        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(LFA, self).open(openmode=openmode)
        if self.openmode in ('r', 'a'):
            # open, getting logical unit
            if self.openmode == 'r':
                self._unit = LFA4py.wlfaouv(self.container.abspath, 'R')
            elif self.openmode == 'a':
                self._unit = LFA4py.wlfaouv(self.container.abspath, 'A')
            self.isopen = True
            self.empty = False
        elif self.openmode == 'w':
            # open
            self._unit = LFA4py.wlfaouv(self.container.abspath, 'W')
            self.isopen = True
            self.empty = True

    def close(self):
        """Closes a LFA properly."""
        if self.isopen:
            try:
                LFA4py.wlfafer(self._unit)
            except Exception:
                raise IOError("closing " + self.container.abspath)
            self.isopen = False
            # Cleanings
            if self.openmode == 'w' and self.empty:
                os.remove(self.container.abspath)

################
# ABOUT FIELDS #
################

    def find_fields_in_resource(self, seed=None):
        """
        Returns a list of the fields from resource whose name match the given
        *seed*.

        :param seed: might be a regular expression, a list of regular expressions
          or *None*. If *None* (default), returns the list of all fields in
          resource.
        """
        if seed is None:
            fieldslist = self.listfields()
        elif isinstance(seed, str):
            fieldslist = util.find_re_in_list(seed, self.listfields())
        elif isinstance(seed, list):
            fieldslist = []
            for s in seed:
                fieldslist += util.find_re_in_list(s, self.listfields())
        if fieldslist == []:
            raise epygramError("no field matching '" + seed +
                               "' was found in resource " +
                               self.container.abspath)

        return fieldslist

    @FileResource._openbeforedelayed
    def readfield(self, fieldname, getdata=True):
        """
        Reads a field in resource.

        :param fieldname: name of the field to be read
        :param getdata: if False, do not read the field data, only metadata.
        """
        field = MiscField(fid={'LFA':fieldname})
        if getdata:
            (fieldtype, fieldlength) = LFA4py.wlfacas(self._unit, fieldname)
            if fieldtype[0] == 'R':
                (data, fieldlength) = LFA4py.wlfalecr(self._unit, fieldname,
                                                      fieldlength)
            elif fieldtype[0] == 'I':
                (data, fieldlength) = LFA4py.wlfaleci(self._unit, fieldname,
                                                      fieldlength)
            elif fieldtype[0] == 'C':
                (data, fieldlength) = LFA4py.wlfalecc(self._unit, fieldname,
                                                      fieldlength,
                                                      config.LFA_maxstrlen)
                data = [data[i].strip().decode() for i in range(fieldlength)]
            field.setdata(numpy.array(data))

        return field

    @FileResource._openbeforedelayed
    def writefield(self, field):
        """Writes a Field in resource."""
        # DEAD-END: we should not need to write with this comdemned format
        raise epygramError("writefield routine has not been tested..." +
                           " If you need to, you might face problems...")

        if not isinstance(field, epygram.base.Field):
            raise epygramError("'field' argument has to be a" +
                               " epygram.base.Field.")
        data = numpy.array(field.getdata())
        if len(data.shape) != 1:
            raise epygramError("LFA can only hold 1D arrays.")

        if data.dtype[0:5] == 'float':
            LFA4py.wlfaecrr(self._unit, field.fid['LFA'], data)
        elif data.dtype[0:3] == 'int':
            LFA4py.wlfaecri(self._unit, field.fid['LFA'], data)
        elif data.dtype[0:3] == 'str':
            LFA4py.wlfaecrc(self._unit, field.fid['LFA'], data)
        else:
            raise epygramError("LFA can only hold float, int or str arrays.")

    def listfields(self):
        """
        Returns a list containing the LFA identifiers of all the fields of
        the resource.
        """
        return super(LFA, self).listfields()

    @FileResource._openbeforedelayed
    def _listfields(self):
        """
        Returns a list containing the names of the fields in LFA.
        """
        (list_length, fieldslist) = LFA4py.wlfalaft(self._unit,
                                                    config.LFA_max_num_fields,
                                                    config.LFA_maxstrlen)
        fieldslist = [fieldslist[i].strip().decode() for i in range(list_length)]
        return fieldslist

    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout,
             sortfields=False,
             **_):
        """
        Writes in file a summary of the contents of the LFA.

        :param out: the output open file-like object
        :param sortfields: **True** if the fields have to be sorted by type.
        """
        firstcolumn_width = 50
        secondcolumn_width = 16
        sepline = '{:-^{width}}'.format('', width=firstcolumn_width +
                                            secondcolumn_width + 1) + '\n'

        listoffields = self.listfields()
        if sortfields:
            listoffields.sort()

        # Write out
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        numfields = len(listoffields)
        out.write("Number: " + str(numfields) + "\n")
        out.write(sepline)
        for f in listoffields:
            out.write(f + "\n")
        out.write(sepline)
