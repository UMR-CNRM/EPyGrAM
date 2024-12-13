#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for GeoPoints format.
"""

import datetime
import numpy
import sys

import footprints
from footprints import FPList, FPDict, proxy as fpx

from epygram import config, epygramError
from epygram.base import FieldSet, FieldValidity
from epygram.resources import FileResource
from epygram.geometries import UnstructuredGeometry, VGeometry
from epygram.fields import H2DField, PointField, D3Field

__all__ = ['GeoPoints']

epylog = footprints.loggers.getLogger(__name__)


class GeoPoints(FileResource):
    """
    Class implementing all specificities for GeoPoints resource format.

    Format:

    <beginning of file>
    #GEO
    #PARAMETER = any name or character string
    #LAT LON DATE TIME VALUE
    #DATA
    32.2671  -12.5501   20140417       21       2.269112E+02
    ...
    <till end of file>
    The set of keys describing each point is defined in its 'columns'
    attribute.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['GeoPoints']),
                default='GeoPoints'),
            parameter=dict(
                optional=True,
                info="The name of the parameter whose data is in the" +
                     " *value* field."),
            columns=dict(
                type=FPList,
                optional=True,
                info="The columns of the geopoints, i.e. the set of" +
                     " keys describing each point."),
            no_header=dict(
                type=bool,
                optional=True,
                default=False,
                info="If True, do not write header (openmode='w')."),
            other_attributes=dict(
                type=FPDict,
                optional=True,
                info="other key:value pairs in header.",
                default=FPDict({})),
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        super(GeoPoints, self).__init__(*args, **kwargs)
        if not self.fmtdelayedopen:
            self.open()

    def open(self,
             parameter=None,
             columns=None,
             openmode=None,
             other_attributes=None):
        """
        Opens a GeoPoint.

        :param parameter: name of the parameter in header (openmode='w').
        :param columns: list of all the items to be written for each point
          (openmode='w').
        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        :param other_attributes: other key:value pairs to be written in header
          (openmode='w').
        """
        super(GeoPoints, self).open(openmode=openmode)

        if self.openmode in ('r', 'a'):
            self._file = open(self.container.abspath, 'r')
            # check format is GeoPoints
            try:
                isgeo = self._file.readline()
            except UnicodeDecodeError:
                isgeo = 'False'
            if isgeo[0:4] != '#GEO':
                raise IOError("this resource is not a GeoPoints one.")
            # read header
            header = {}
            self.headerlength = 1
            headerline = self._file.readline()
            while headerline[0:5] != '#DATA':
                self.headerlength += 1
                if '=' in headerline:
                    metadata = headerline[:-1].replace('#', '').split('=')
                    header[metadata[0].strip()] = metadata[1].strip()
                elif 'FORMAT' in headerline:
                    metadata = headerline[:-1].replace('#', '').split()
                    header[metadata[0].strip()] = metadata[1].strip()
                else:
                    # description of columns
                    self._attributes['columns'] = headerline[1:-1].split()
                headerline = self._file.readline()
            self.headerlength += 1
            for meta in header.keys():
                if meta in ['PARAMETER', ]:
                    self._attributes[meta.lower()] = header[meta]
                else:
                    self.other_attributes[meta] = header[meta]
            # initializations
            self.isopen = True
            self.empty = False
        elif self.openmode == 'w':
            # object initializations
            self.empty = True
            if parameter is not None and self.parameter is None:
                self._attributes['parameter'] = parameter
            if columns is not None and self.columns is None:
                self._attributes['columns'] = columns
            if other_attributes is not None and self.other_attributes is None:
                self._attributes['other_attributes'] = other_attributes
            if len(self.other_attributes) > 0:
                if 'FORMAT' in self.other_attributes.keys():
                    if columns is not None:
                        epylog.warning("Double specification: columns overwritten by FORMAT.")
                    if self.other_attributes['FORMAT'] == 'XYV':
                        self._attributes['columns'] = ['LON', 'LAT', 'VALUE']
                    else:
                        raise NotImplementedError("Not yet. Cf. https://software.ecmwf.int/wiki/display/METV/Geopoints for implementing.")

            if self.parameter is not None and self.columns is not None:
                if not self.no_header:
                    # header initializations
                    self._file = open(self.container.abspath, self.openmode)
                    self._file.write('#GEO\n')
                    for k, v in self.other_attributes.items():
                        self._file.write('#' + k + ' = ' + str(v) + '\n')
                    self._file.write('#PARAMETER = ' + self.parameter + '\n')
                    self._file.write('#')
                    for i in self.columns:
                        self._file.write(i + ' ')
                    self._file.write('\n')
                    self._file.write('#DATA\n')
                    self.isopen = True
                else:
                    self._file = open(self.container.abspath, self.openmode)
                    self.isopen = True
            else:
                epylog.debug("cannot open GeoPoints in 'w' mode without" +
                             " attributes *parameter* and *columns* set." +
                             " Not open yet.")

    def close(self):
        """
        Closes a GeoPoints.
        """
        if hasattr(self, '_file'):
            self._file.close()
        self.isopen = False

    def listfields(self, **kwargs):
        """
        Returns a list containing the LFI identifiers of all the fields of the resource.
        """
        return super(GeoPoints, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for GeoPoints.

        :param complete: - if True method returns a list of {'GeoPoints':GeoPoints_fid, 'generic':generic_fid}
                         - if False method return a list of GeoPoints_fid
        """
        if complete:
            return [{'GeoPoints':self.parameter, 'generic':{}}]
        else:
            return [self.parameter]

################
# ABOUT POINTS #
################

    def countpoints(self):
        """
        Counts and returns the number of points in the resource.
        """
        if self.isopen:
            self._file.close()
            self.isopen = False
            wasopen = True
        n = sum(1 for _ in open(self.container.abspath)) - self.headerlength
        if wasopen:
            if self.openmode == 'w' and not self.empty:
                self._file = open(self.container.abspath, 'a')
                self.isopen = True
            else:
                self._file = open(self.container.abspath, self.openmode)
                self.isopen = True
                if self.openmode == 'r':
                    for _ in range(self.headerlength):
                        self._file.readline()

        return n

    @FileResource._openbeforedelayed
    def readfield(self,
                  parameter='*',
                  footprints_builder=False,
                  as_points=False):
        """
        Reads the GeoPoints.

        :param parameter: is the name of the parameter to read
          (it exists only for compatibility with other formats)
        :param footprints_builder: if *True*, uses footprints.proxy to build
          fields. Defaults to False for performance reasons.
        :param *as_points*: if *True*, returns a collection of points instead
          of a single field.
        """
        if self.openmode in ('w', 'a'):
            raise epygramError("cannot read fields in resource if with" +
                               " openmode == 'w' or 'a'.")
        elif not self.isopen:
            self.open()

        if parameter not in ['*', self.parameter]:
            raise epygramError("This file does not contain the requested parameter.")

        # Rewind
        self._file.close()
        self._file = open(self.container.abspath, self.openmode)
        for _ in range(self.headerlength):
            self._file.readline()

        specialValues_list = {'DATE':[], 'TIME':[], 'LEVEL':[]}  # To contain the date, time and date values
        specialValues_cst = {'DATE':True, 'TIME':True, 'LEVEL':True}  # To know if date, time and level are constant or not
        specialValues_prev = {}  # Previous values of date, time and leve to detect changes in values
        values = []
        latitudes = []
        longitudes = []
        nb = 0
        for line in iter(self._file):
            pointvals = line.split()
            point = dict(zip(self.columns, pointvals))
            # We store date, time and levels as a list of one value if parameter is constant
            # or as a multi-values list if the parameter varies
            for item in ['DATE', 'TIME', 'LEVEL']:
                if item in self.columns:
                    if nb == 0:
                        specialValues_prev[item] = point[item]
                    if specialValues_cst[item] and specialValues_prev[item] != point[item]:
                        specialValues_cst[item] = False
                        specialValues_list[item] = specialValues_list[item] * nb
                    if nb == 0 or not specialValues_cst[item]:
                        specialValues_list[item].append(point[item])
            values.append(float(point['VALUE']))
            latitudes.append(float(point['LAT']))
            longitudes.append(float(point['LON']))
            nb += 1

        field_structure = 'Point' if nb == 1 else 'H2D'

        if footprints_builder:
            field_builder = fpx.field
        else:
            if field_structure == 'Point':
                field_builder = PointField
            else:
                field_builder = H2DField

        if 'LEVEL' in self.columns:
            specialValues_list['LEVEL'] = [float(level)
                                           for level in specialValues_list['LEVEL']]
        vcoordinate = VGeometry(typeoffirstfixedsurface=255,
                                levels=[255] if 'LEVEL' not in self.columns else specialValues_list['LEVEL'])
        geometry = UnstructuredGeometry(name='unstructured',
                                        dimensions={'X':len(values), 'Y':1},
                                        vcoordinate=vcoordinate,
                                        grid={'longitudes':longitudes,
                                              'latitudes':latitudes},
                                        position_on_horizontal_grid='center'
                                        )
        if 'DATE' in self.columns and 'TIME' in self.columns:
            if specialValues_cst['DATE'] and not specialValues_cst['TIME']:
                specialValues_cst['DATE'] = False
                specialValues_list['DATE'] = specialValues_list['DATE'] * nb
            if specialValues_cst['TIME'] and not specialValues_cst['DATE']:
                specialValues_cst['TIME'] = False
                specialValues_list['TIME'] = specialValues_list['TIME'] * nb
        if 'DATE' in self.columns:
            for i in range(len(specialValues_list['DATE'])):
                date = specialValues_list['DATE'][i]
                specialValues_list['DATE'][i] = datetime.datetime(int(date[0:4]),
                                                                  int(date[4:6]),
                                                                  int(date[6:8]))
        if 'TIME' in self.columns:
            for i in range(len(specialValues_list['TIME'])):
                time = specialValues_list['TIME'][i]
                if len(time) in (1, 2):
                    specialValues_list['TIME'][i] = datetime.time(int(time))
                else:
                    specialValues_list['TIME'][i] = datetime.time(int(time[0:2]),
                                                                  int(time[2:4]))
        if 'DATE' in self.columns or 'TIME' in self.columns:
            dates = []
            if 'DATE' in self.columns and 'TIME' in self.columns:
                for i in range(len(specialValues_list['DATE'])):
                    delta = datetime.timedelta(seconds=specialValues_list['TIME'][i].hour * 3600 +
                                                       specialValues_list['TIME'][i].minute * 60 +
                                                       specialValues_list['TIME'][i].second)
                    dates.append(specialValues_list['DATE'][i] + delta)
            elif 'DATE' in self.columns:
                dates = specialValues_list['DATE']
            else:
                dates = specialValues_list['TIME']
            if specialValues_cst['DATE']:
                validity = FieldValidity(date_time=dates[0])
            else:
                validity = FieldValidity(date_time=dates)
        else:
            validity = FieldValidity()

        field = field_builder(structure=field_structure,
                              fid={self.format:self.parameter},
                              geometry=geometry,
                              validity=validity)
        field.setdata(numpy.array(values).reshape(1, len(values)))

        if as_points:
            field = field.as_points()

        return field

    def writefield(self, field,
                   parameter='UNKNOWN',
                   fidkey_for_parameter=None,
                   order='C',
                   llprecision=config.GeoPoints_lonlat_precision,
                   precision=config.GeoPoints_precision,
                   col_width=config.GeoPoints_col_width,
                   subzone=None):
        """
        Write a field in the resource.

        :param field: a :class:`epygram.fields.D3Field` or a
          :class:`epygram.base.FieldSet` of :class:`epygram.fields.PointField`,
          or a dict {'lon':[...], 'lat':[...], 'value':[...], ...}. The keys
          of field that are not in the *columns* are ignored. The keys of the
          columns that are not in the field are filled by zeros.
          At least 'lon' and 'lat' are necessary.
        :param parameter: the name of the parameter, to be given if *field* is a
          dict and if it has not been given at opening.
        :param order: optional, for a rectangular D3Field, whether to flatten
          2D arrays in 'C' (row-major) or 'F' (Fortran, column-major) order.
        :param llprecision: precision (number of digits) in write for lon/lat.
        :param precision: precision (number of digits) in write for other floats.
        :param subzone: optional, among ('C', 'CI').
          Only if *field* is a LAM :class:`epygram.fields.H2DField`,
          extracts the points of the
          resp. C or C+I zone. Default is no subzone, i.e. the whole field.
        """
        if self.openmode == 'r':
            raise IOError("cannot write field in a GeoPoints with openmode 'r'.")
        
        if not self.isopen:
            open_kwargs = {}
            if self.columns:
                columns = self.columns
            else:
                columns = ['LAT', 'LON', 'LEVEL']
            if isinstance(field, D3Field):
                if parameter == 'UNKNOWN':
                    if fidkey_for_parameter is None:
                        open_kwargs['parameter'] = str(field.fid.get(self.format, field.fid))
                    else:
                        open_kwargs['parameter'] = str(field.fid[fidkey_for_parameter])
                if field.validity.get() is not None and not self.columns:
                    columns.extend(['DATE', 'TIME'])
            elif isinstance(field, FieldSet) or isinstance(field, PointField):
                if isinstance(field, PointField):
                    field = FieldSet([field])
                if parameter == 'UNKNOWN':
                    if fidkey_for_parameter is None:
                        open_kwargs['parameter'] = str(field[0].fid.get(self.format, field[0].fid))
                    else:
                        open_kwargs['parameter'] = str(field[0].fid[fidkey_for_parameter])
                if field[0].validity.get() is not None and not self.columns:
                    columns.extend(['DATE', 'TIME'])
            elif isinstance(field, dict) and not self.columns:
                if 'lat' not in field:
                    raise epygramError("'lat' must be in field.keys().")
                if 'lon' not in field:
                    raise epygramError("'lon' must be in field.keys().")
                if 'value' not in field:
                    raise epygramError("'value' must be in field.keys().")
                for i in field.keys():
                    if i.upper() not in ('LAT', 'LON', 'VALUE'):
                        columns.append(i.upper())
            if not self.columns:
                columns.append('VALUE')
            if self.columns is None:
                open_kwargs['columns'] = columns
            self.open(**open_kwargs)

        if isinstance(field, D3Field):
            (lons, lats) = field.geometry.get_lonlat_grid(nb_validities=len(field.validity),
                                                          subzone=subzone,
                                                          d4=True)
            if field.geometry.rectangular_grid:
                lons = lons.flatten(order=order)
                lats = lats.flatten(order=order)
                values = field.getdata(subzone=subzone, d4=True).flatten(order=order)
                if 'LEVEL' in self.columns:
                    levels = field.geometry.get_levels(nb_validities=len(field.validity),
                                                       subzone=subzone,
                                                       d4=True).flatten(order=order)
            else:
                lons = lons.compressed()
                lats = lats.compressed()
                values = field.getdata(subzone=subzone, d4=True).compressed()
                if 'LEVEL' in self.columns:
                    levels = field.geometry.get_levels(nb_validities=len(field.validity),
                                                       subzone=subzone,
                                                       d4=True).flatten(order=order)
            if 'DATE' in self.columns or 'TIME' in self.columns:
                date = field.validity.get(fmt='IntStr')
                if 'TIME' in self.columns:
                    hour = date[8:10]
                if 'DATE' in self.columns:
                    date = date[0:8]
            writebuffer = {'LON':lons, 'LAT':lats, 'VALUE':values}
            if 'LEVEL' in self.columns:
                # level = field.geometry.vcoordinate.get('level', 0)
                writebuffer['LEVEL'] = levels
            for i in self.columns:
                if i == 'DATE':
                    writebuffer[i] = [date] * len(values)
                elif i == 'TIME':
                    writebuffer[i] = [hour] * len(values)
                # elif i == 'LEVEL':
                #     writebuffer[i] = [level] * len(values)
                # others to be implemented here
        elif isinstance(field, FieldSet) or isinstance(field, PointField):
            if isinstance(field, PointField):
                field = FieldSet([field])
            writebuffer = {'LON':[], 'LAT':[], 'VALUE':[], 'DATE':[], 'TIME':[]}
            for pt in field:
                if not isinstance(pt, PointField):
                    raise NotImplementedError("write a FieldSet of fields" +
                                              " that are not PointField.")
                writebuffer['LON'].append(pt.geometry.hlocation['lon'].get('degrees'))
                writebuffer['LAT'].append(pt.geometry.hlocation['lat'].get('degrees'))
                writebuffer['VALUE'].append(float(pt.getdata()))
                if 'DATE' in self.columns or 'TIME' in self.columns:
                    date = pt.validity.get(fmt='IntStr')
                    if 'TIME' in self.columns:
                        writebuffer['TIME'].append(date[8:10])
                    if 'DATE' in self.columns:
                        writebuffer['DATE'].append(date[0:8])
                    # others to be implemented here
        elif isinstance(field, dict):
            writebuffer = {k.upper():v for k, v in field.items()}
        else:
            raise NotImplementedError("*field* of that type.")
        for i in self.columns:
            if i not in writebuffer:
                writebuffer[i] = [0] * len(writebuffer['VALUE'])

        # WRITE
        for pt in range(len(writebuffer['VALUE'])):
            pointstr = ''
            for i in self.columns[:2]:
                pointstr += ' {: .{precision}{type}} '.format(writebuffer[i][pt],
                                                              type='F',
                                                              precision=llprecision)
            for i in self.columns[2:]:
                x = writebuffer[i][pt]
                if type(x) in (float, numpy.float64):
                    x = ' {: .{precision}{type}} '.format(x, type='E',
                                                          precision=precision)
                else:
                    x = '{:^{width}}'.format(x, width=col_width)
                pointstr += x
            self._file.write(pointstr + "\n")

        if self.empty:
            self.empty = False


###########
# pre-app #
###########
    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout, **_):
        """
        Writes in file-like a summary of the contents of the GeoPoints.

        :param out: the output open file-like object.
        """
        points_number = self.countpoints()
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        out.write("#####################" + "\n")
        out.write("# GeoPoints Summary #" + "\n")
        out.write("#####################" + "\n")
        out.write("Parameter (described by VALUE): " + self.parameter + "\n")
        out.write("Number of points              : " + str(points_number) + "\n")
        out.write("Data/meta for each point      : " + str(self.columns) + "\n")
        if 'structure' in self.footprint_attributes:
            out.write("Structure of Points as a Field: " + str(self.structure) + "\n")
