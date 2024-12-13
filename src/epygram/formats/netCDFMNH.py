#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to handle the Meso-NH netCDF format.
"""

import datetime
import os
import copy
import numpy
import math
import re
import sys

import footprints
from footprints import FPDict, proxy as fpx
from bronx.datagrip.misc import read_dict_in_CSV

import netCDF4
import json

from epygram import config, epygramError, util
from epygram.util import Angle, RecursiveObject
from epygram.base import FieldValidity, Field, FieldValidityList
from epygram.resources import FileResource
from epygram.fields import MiscField
from epygram.geometries import VGeometry, AcademicGeometry, ProjectedGeometry
from epygram.geometries.VGeometry import hybridH2altitude, hybridH2pressure

__all__ = ['netCDFMNH']

epylog = footprints.loggers.getLogger(__name__)

gridIndicatorDict = {0:('__unknown__', '__unknown__'),
                     1:('center', 'mass'),
                     2:('center-left', 'mass'),
                     3:('lower-center', 'mass'),
                     4:('center', 'flux'),
                     5:('lower-left', 'mass'),
                     6:('center-left', 'flux'),
                     7:('lower-center', 'flux'),
                     8:('lower-left', 'flux')}

def inquire_field_dict(fieldname):
    """
    Returns the info contained in the netCDFMNH_field_dict for the requested field.
    """
    matching_field = None
    for fd in netCDFMNH._field_dict:
        dictitem = fd['name']
        pattern = re.subn('\.', r'\.', dictitem)[0]  # protect '.'
        pattern = pattern.replace('?', '.')  # change unix '?' to python '.' (any char)
        pattern = pattern.replace('*', '.*')  # change unix '*' to python '.*' (several any char)
        pattern += '(?!.)'
        if re.match(pattern, fieldname):
            matching_field = fd
            break

    if matching_field is None:
        epylog.info("field '" + fieldname + "' is not referenced in Field_Dict_netCDFMNH.")
        matching_field = {'name':fieldname}

    return copy.deepcopy(matching_field)

class empty(RecursiveObject):
    """
    Class to hold some attributes contained in attributes
    normally missing for a Misc field. For example, the
    position on grid can be set in NetCDFMNH Misc field
    whereas this has no meaning (because Misc fields are
    not on a grid).
    """
    pass

class netCDFMNH(FileResource):
    """
    Class implementing all specificities for MesoNH netCDF resource format.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['netCDFMNH']),
                default='netCDFMNH'),
            moveOnMass=dict(
                info="If True, 3d fields are put on mass levels",
                optional=True,
                default=False,
                type=bool)
        )
    )

    # the Field Dictionary gathers info about fields nature
    CSV_field_dictionaries = config.netCDFMNH_field_dictionaries_csv
    # syntax: _field_dict = [{'name':'fieldname1', 'type':'...', ...}, {'name':'fieldname2', 'type':'...', ...}, ...]
    _field_dict = []

    _specialFieldComments = {'CARTESIAN':dict(long_name='CARTESIAN', comment='Logical for cartesian geometry'),
                             'LAT0':dict(long_name='LAT0', comment='Reference latitude for conformal projection', units='degree'),
                             'LON0':dict(long_name='LON0', comment='Reference longitude for conformal projection', units='degree'),
                             'LATORI':dict(long_name='LATORI', units='degree',
                                           comment='Latitude of the point of coordinates x=0, y=0 for conformal projection'),
                             'LATOR':dict(long_name='LATOR', comment='Latitude of 1st mass point', units='degree'),
                             'LONORI':dict(long_name='LONORI', units='degree',
                                           comment='Longitude of the point of coordinates x=0, y=0 for conformal projection'),
                             'LONOR':dict(long_name='LONOR', comment='Longitude of 1st mass point', units='degree'),
                             'RPK':dict(long_name='RPK', comment='Projection parameter for conformal projection'),
                             'BETA':dict(long_name='BETA', comment='Rotation angle for conformal projection', units='degree'),
                             'IMAX':dict(long_name='IMAX', comment='x-dimension of the physical domain'),
                             'JMAX':dict(long_name='JMAX', comment='y-dimension of the physical domain'),
                             'JPHEXT':dict(long_name='JPHEXT', comment='Number of horizontal external points on each side'),
                             'KMAX':dict(long_name='KMAX', comment='z-dimension of the physical domain'),
                             'SLEVE':dict(long_name='SLEVE', comment="Logical for SLEVE coordinate", grid=numpy.int32(4)),
                             'XHAT':dict(standard_name='projection_x_coordinate', long_name='XHAT', units='m',
                                         comment='Position x in the conformal or cartesian plane', _FillValue=9.96920996838687e+36,
                                         valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(2)),
                             'YHAT':dict(standard_name='projection_y_coordinate', long_name='YHAT', units='m',
                                         comment='Position y in the conformal or cartesian plane', _FillValue=9.96920996838687e+36,
                                         valid_min=-1.e+36 , valid_max=1.e+36, grid=numpy.int32(3)),
                             'ZHAT':dict(long_name='ZHAT', units='m', comment='Height level without orography',
                                         _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(4)),
                             'level':dict(long_name='position z in the transformed space', standard_name='', units='m',
                                          axis='Z', positive='up', c_grid_axis_shift=0.,
                                          formula_terms='s: level height: ZTOP orog: ZS',
                                          formula_definition='z(n,k,j,i)=s(k)*(height-orog(j,i))/height+orog(j,i)',
                                          computed_standard_name='altitude'),
                             'level_w':dict(long_name='position z in the transformed space at w location', standard_name='',
                                            units='m', axis='Z', positive='up', c_grid_axis_shift=-0.5,
                                            formula_terms='s: level_w height: ZTOP orog: ZS',
                                            formula_definition='z(n,k,j,i)=s(k)*(height-orog(j,i))/height+orog(j,i)',
                                            computed_standard_name='altitude_at_w_location'),
                             'ZTOP':dict(standard_name='altitude_at_top_of_atmosphere_model', long_name='ZTOP',
                                         units='m', comment='Height of top level', grid=numpy.int32(4)),
                             'ni':dict(long_name='x-dimension of the grid', standard_name='projection_x_coordinate',
                                       units='m', axis='X', c_grid_axis_shift=0.),        
                             'nj':dict(long_name='y-dimension of the grid', standard_name='projection_y_coordinate',
                                       units='m', axis='Y', c_grid_axis_shift=0.),
                             'ni_u':dict(long_name='x-dimension of the grid at u location',
                                         standard_name='projection_x_coordinate_at_u_location',
                                         units='m', axis='X', c_grid_axis_shift=-0.5),
                             'nj_u':dict(long_name='y-dimension of the grid at u location',
                                         standard_name='projection_y_coordinate_at_u_location',
                                         units='m', axis='Y', c_grid_axis_shift=0.),
                             'ni_v':dict(long_name='x-dimension of the grid at v location',
                                         standard_name='projection_x_coordinate_at_v_location',
                                         units='m', axis='X', c_grid_axis_shift=0.),
                             'nj_v':dict(long_name='y-dimension of the grid at v location',
                                         standard_name='projection_y_coordinate_at_v_location',
                                         units='m', axis='Y', c_grid_axis_shift=-0.5),
                             'LAT':dict(long_name='LAT', units='degrees_north', comment='X_Y_latitude',
                                        coordinates='latitude longitude', _FillValue=9.96920996838687e+36,
                                        valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(1)),
                             'LON':dict(long_name='LON', units='degrees_east', comment='X_Y_longitude',
                                        coordinates='latitude longitude', _FillValue=9.96920996838687e+36,
                                        valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(1)),
                             'latitude':dict(standard_name='latitude', long_name='latitude',
                                             units='degrees_north', comment='X_Y_latitude at mass point',
                                             _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(1)),
                             'longitude':dict(standard_name='longitude', long_name='longitude',
                                              units='degrees_east', comment='X_Y_longitude at mass point',
                                              _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(1)),
                             'latitude_u':dict(standard_name='latitude_at_u_location', long_name='latitude at u location',
                                               units='degrees_north', comment='X_Y_latitude at u point',
                                               _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(2)),
                             'longitude_u':dict(standard_name='longitude_at_u_location', long_name='longitude at u location',
                                                units='degrees_east', comment='X_Y_longitude at u point',
                                                _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(2)),
                             'latitude_v':dict(standard_name='latitude_at_v_location', long_name='latitude at v location',
                                               units='degrees_north', comment='X_Y_latitude at v point',
                                               _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(3)),
                             'longitude_v':dict(standard_name='longitude_at_v_location', long_name='longitude at v location',
                                                units='degrees_east', comment='X_Y_longitude at v point',
                                                _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(3)),
                             'latitude_f':dict(standard_name='latitude_at_f_location', long_name='latitude at f location',
                                               units='degrees_north', comment='X_Y_latitude at f point',
                                               _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(5)),
                             'longitude_f':dict(standard_name='longitude_at_f_location', long_name='longitude at f location',
                                                units='degrees_east', comment='X_Y_longitude at f point',
                                                _FillValue=9.96920996838687e+36, valid_min=-1.e+36, valid_max=1.e+36, grid=numpy.int32(5)),
                             #'DTMOD':dict(long_name='DTMOD', comment='Time and date of model beginning', calendar='standard'),
                             'DTCUR':dict(standard_name='time', long_name='DTCUR', comment='Current time and date', calendar='standard'),
                             'DTEXP':dict(long_name='DTEXP', comment='Time and date of experiment beginning', calendar='standard'),
                             #'DTSEG':dict(long_name='DTSEG', comment='Time and date of segment beginning', calendar='standard'),
                             'time':dict(long_name='time axis', standard_name='time', axis='T', calendar='standard'),
                            }

    @classmethod
    def _read_field_dict(cls, fd_abspath):
        """Reads the CSV fields dictionary of the format."""
        field_dict, file_priority = read_dict_in_CSV(fd_abspath)
        if file_priority == 'main':
            cls._field_dict = field_dict
        elif file_priority == 'underwrite':
            for fd in field_dict:
                found = False
                for cfd in cls._field_dict:
                    if fd['name'] == cfd['name']:
                        found = True
                        break
                if not found:
                    cls._field_dict.append(fd)
        elif file_priority == 'overwrite':
            for cfd in cls._field_dict:
                found = False
                for fd in field_dict:
                    if fd['name'] == cfd['name']:
                        found = True
                        break
                if not found:
                    field_dict.append(cfd)
            cls._field_dict = field_dict

    def __init__(self, *args, **kwargs):
        self.isopen = False
        self.geometry = None
        self.validity = None

        # At creation of the first netCDFMNH, initialize netCDFMNH._field_dict
        if self._field_dict == []:
            self._read_field_dict(self.CSV_field_dictionaries['default'])
            if os.path.exists(self.CSV_field_dictionaries['user']):
                self._read_field_dict(self.CSV_field_dictionaries['user'])

        super(netCDFMNH, self).__init__(*args, **kwargs)
        if self.openmode in ('r', 'a'):
            try:
                guess = netCDF4.Dataset(self.container.abspath, self.openmode)
                if not all([v in guess.variables for v in ['IMAX', 'JMAX', 'CARTESIAN']]):
                    raise IOError("this resource is not a netCDFMNH one.")
            except (RuntimeError, UnicodeEncodeError):
                raise IOError("this resource is not a netCDF one.")
            else:
                guess.close()

        if not self.fmtdelayedopen:
            self.open()

    def open(self, openmode=None):
        """
        Opens the file.

        :param openmode: optional, to open with a specific openmode, eventually
                         different from the one specified at initialization.
        """

        super(netCDFMNH, self).open(openmode=openmode)
        self._nc = netCDF4.Dataset(self.container.abspath, self.openmode)
        self.isopen = True
        self.empty = self.isopen == 'w'

    def close(self):
        """Closes a Meso-NH netCDF file."""
        if self.isopen:
            self._nc.close()
        self.isopen = False

################
# ABOUT FIELDS #
################
    def _get_type(self, fieldname):
        """Returns the type of field from the dimensions"""
        dimensions = self._nc[fieldname].dimensions
        if len(dimensions) > 0 and dimensions[0] == 'time':
            has_time = True
            dimensions = dimensions[1:]
        else:
            has_time = False
        if len(dimensions) > 0 and dimensions[0] in ['level', 'level_w']:
            has_level = True
            dimensions = dimensions[1:]
        else:
            has_level = False
        if dimensions in [('nj', 'ni'), ('nj_u', 'ni_u'),
                          ('nj_v', 'ni_v'), ('nj_v', 'ni_u')]:
            return '3D' if has_level else 'H2D'
        elif dimensions in [('nj',), ('ni',), ('nj_u',), ('ni_u',), ('nj_v',), ('ni_v'),
                            ('size1', 'ni_u'), ('size1', 'ni_v'), ('size1', 'ni')]:
            return 'V2D' if has_level else 'H1D'
        elif len(dimensions) == 0:
            return 'V1D' if has_level else 'Point'
        else:
            return 'Misc'

    def _get_generic_fid(self, fieldname):
        """Return a generic fid from fieldname (via Field Dict)."""
        fid = inquire_field_dict(fieldname)
        if self._get_type(fieldname) in ['3D', 'V2D', 'V1D']:
            if self.geometry is None:
                self._read_geometry()
            fid['typeOfFirstFixedSurface'] = self.geometry.vcoordinate.typeoffirstfixedsurface
        fid.pop('name')

        return fid
    
    def find_fields_in_resource(self, seed=None, fieldtype=[], generic=False):
        """
        Returns a list of the fields from resource whose identifier match the given seed.

        :param seed: might be:\n
          - a regular expressions,
          - a list of regular expressions
          - *None*. If *None* (default), returns the list of all fields in resource.
        :param fieldtype: optional, among ('H2D', 'Misc') or a list of these strings.
          If provided, filters out the fields not of the given types.
        :param generic: if True, returns a list of tuples (fieldname, generic fid) of
          the fields.
        """
        if isinstance(fieldtype, list):
            fieldtypeslist = list(fieldtype)
        else:
            fieldtypeslist = [fieldtype]
        fieldslist = []
        def fill_fieldslist(tmplist):
            for f in tmplist:
                if fieldtypeslist == [] or self._get_type(f) in fieldtypeslist:
                    fieldslist.append(f)
        if seed is None:
            tmplist = self.listfields()
            fill_fieldslist(tmplist)
        elif isinstance(seed, str):
            tmplist = util.find_re_in_list(seed, self.listfields())
            fill_fieldslist(tmplist)
        elif isinstance(seed, list):
            tmplist = []
            for s in seed:
                tmplist += self.find_fields_in_resource(seed=s)
            fill_fieldslist(tmplist)
        else:
            raise epygramError("seed must be a list, None or a string")
        if fieldslist == []:
            raise epygramError("no field matching '" + str(seed) + "' was found in resource " + self.container.abspath)

        if generic:
            fieldslist = [(f, self._get_generic_fid(f)) for f in fieldslist]

        return fieldslist

    def listfields(self, **kwargs):
        """
        Returns a list containing the netCDFMNH identifiers of all the fields of the
        resource.
        """
        return super(netCDFMNH, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method.

        :param complete: - if True method returns a list of {'netCDFMNH':fid,
                           'generic':generic_fid}
                         - if False method return a list of fid
        """
        fieldslist = [f for f in self._nc.variables.keys() if f not in self._specialFieldComments]

        if complete:
            return [{'netCDFMNH':f, 'generic':self._get_generic_fid(f)} for f in fieldslist]
        else:
            return fieldslist

    def sortfields(self):
        """
        Returns a sorted list of fields with regards to their name and nature,
        as a dict of lists.
        """
        listMisc = []
        list3D = []
        list2D = []

        for field in self.listfields():
            info = self._get_type(field)
            if info in ['V2D', 'H2D']:
                list2D.append(field)
            elif info == '3D':
                list3D.append(field)
            else:
                listMisc.append(field)

        # sort
        list2D.sort()
        list3D.sort()
        listMisc.sort()

        outlists = {'3D fields':list(set(list3D)),
                    '2D fields':list(set(list2D)),
                    'Misc-fields':list(set(listMisc))}

        return outlists

    @FileResource._openbeforedelayed
    def readfield(self, fieldidentifier, getdata=True):
        """
        Reads one field, given its name and returns a Field instance.

        :param fieldidentifier: field name.
        :param getdata: optional, if *False*, only metadata are read, the field do not contain data.
                        Default is *True*.
        """
        if not isinstance(fieldidentifier, str):
            raise epygramError("fieldidentifier of a netCDFMNH file is a string.")
        if fieldidentifier in self._specialFieldComments:
            raise epygramError("This is a special field that is not allowed do be read as field")

        # Get field info
        var = self._nc.variables[fieldidentifier]
        field_info = inquire_field_dict(fieldidentifier)
        field_info['type'] = self._get_type(fieldidentifier)
        if field_info['type'] in ['H1D', 'V1D', 'V2D', 'H2D', '3D']:
            if self.geometry is None:
                self._read_geometry()
            if self.validity is None:
                self._read_validity()

            # Make geometry object
            kwargs_geom = dict(name=self.geometry.name,
                               grid=copy.deepcopy(self.geometry.grid),
                               dimensions=copy.deepcopy(self.geometry.dimensions),
                               geoid=config.netCDFMNH_default_geoid,
                               projection=copy.deepcopy(self.geometry.projection)  # Also used for academic geometries
                               )

            if self.geometry.vcoordinate is not None:
                # vertical geometry
                kwargs_vcoord = {'typeoffirstfixedsurface': self.geometry.vcoordinate.typeoffirstfixedsurface,
                                 'position_on_grid': self.geometry.vcoordinate.position_on_grid,
                                 'grid': copy.copy(self.geometry.vcoordinate.grid),
                                 'levels': copy.copy(self.geometry.vcoordinate.levels)}
                if field_info['type'] == 'H2D' and 'level' not in field_info:
                    field_info['level'] = 0
                for k in field_info:
                    if k == 'typeOfFirstFixedSurface':
                        kwargs_vcoord['typeoffirstfixedsurface'] = field_info[k]
                    elif k == 'level':
                        kwargs_vcoord['levels'] = [field_info[k]]
                if field_info['type'] not in ('V1D', 'V2D', '3D'):
                    kwargs_vcoord.pop('grid', None)

        # Get field metadata
        if 'grid' in var.ncattrs():
            (h, v) = gridIndicatorDict[var.grid]
        else:
            h = v = None
        gridIndicator = {'vertical':v, 'horizontal':h}
        comment = {}
        for a in var.ncattrs():
            if a not in ['grid', 'c_grid_dynamic_range']:
                if isinstance(var.getncattr(a), numpy.float32):  # pb with json and float32
                    comment.update({a:numpy.float64(var.getncattr(a))})
                elif isinstance(var.getncattr(a), numpy.int32):  # pb with json and int32
                    comment.update({a:numpy.int64(var.getncattr(a))})
                elif isinstance(var.getncattr(a), numpy.ndarray):  # pb with json and numpy arrays
                    comment.update({a:numpy.float64(var.getncattr(a)).tolist()})
                else:
                    comment.update({a:var.getncattr(a)})
        comment = json.dumps(comment)
        if comment == '{}':
            comment = None

        # Create field
        if field_info['type'] in ['H1D', 'V1D', 'V2D', 'H2D', '3D']:
            # Create H2D field
            fid = {self.format: fieldidentifier,
                   'generic':FPDict(self._get_generic_fid(fieldidentifier))
                   }
            kwargs_geom['position_on_horizontal_grid'] = gridIndicator['horizontal']
            if field_info['type'] in ['V1D', 'V2D', '3D'] and gridIndicator['vertical'] == 'flux' and self.moveOnMass:
                kwargs_vcoord['position_on_grid'] = 'mass'
            else:
                kwargs_vcoord['position_on_grid'] = gridIndicator['vertical']
            kwargs_geom['vcoordinate'] = VGeometry(**kwargs_vcoord)
            geometry = self.geometry.__class__(**kwargs_geom)
            if 'time' in var.dimensions:
                validity = self.validity.deepcopy()
            else:
                validity = FieldValidity()
            field = fpx.field(fid=fid,
                              structure=geometry.structure,
                              geometry=geometry, validity=validity,
                              processtype='forecast', comment=comment)
        elif field_info['type'] == 'Misc':
            # Create Misc field
            fid = {self.format: fieldidentifier,
                   'generic': FPDict()
                   }
            geometry = empty()
            geometry.position_on_horizontal_grid = gridIndicator['horizontal']
            geometry.vcoordinate = empty()
            geometry.vcoordinate.position_on_grid = gridIndicator['vertical']
            field = MiscField(fid=fid, comment=comment)
            if 'time' in var.dimensions:
                if self.validity is None:
                    self._read_validity()
                field.validity = FieldValidityList(self.validity.deepcopy())
            field.geometry = geometry
        if getdata:
            data = var[...]
            if field_info['type'] in('H1D', 'H2D'):
                # Only one horizontal level
                data = geometry.reshape_data(data.flatten())
            elif field_info['type'] in ('V1D', 'V2D', '3D'):
                # 3D data
                data = data.reshape((len(self.geometry.vcoordinate.grid['gridlevels']) + 1,
                                    self.geometry.dimensions['X'] * self.geometry.dimensions['Y']))
                data = geometry.reshape_data(data, 'Z')
                if gridIndicator['vertical'] == 'flux' and self.moveOnMass:
                    data[:-1, :] = 0.5 * (data[:-1, :] + data[1:, :]) #last level kept untouched
            field.setdata(data)

        return field

    def readfields(self, requestedfields=None, getdata=True):
        """
        Returns a :class:`epygram.base.FieldSet` containing requested fields read in the resource.

        :param requestedfields: might be \n
          - a string or a regular expression (e.g. 'R?T' or 'RCT')
          - a list of strings or regular expressions (e.g. ['COVER???', 'RVT'])
          - if not specified, interpretated as all fields that will be found in resource
        :param getdata: optional, if *False*, only metadata are read, the fields do not contain data.
                        Default is *True*.
        """

        requestedfields = self.find_fields_in_resource(requestedfields)
        if requestedfields == []:
            raise epygramError("unable to find requested fields in resource.")

        return super(netCDFMNH, self).readfields(requestedfields, getdata)

    @FileResource._openbeforedelayed
    def writefield(self, field, miscdims=None):
        """
        Write a field in the resource.

        :param field: the field to write to the resource.
        :param miscdims: to force dimensions (useful only for misc fields)
        """
        if not isinstance(field, Field):
            raise epygramError("*field* must be a Field instance.")
        if not self.format in field.fid:
            raise epygramError("field fid must contain the " + self.format + " key")
        varname = field.fid[self.format]
        if not isinstance(varname, str):
            raise epygramError("field fid of a netCDFMNH file is a string.")
        if varname in self._specialFieldComments:
            raise epygramError("This fid (" + varname + \
                               ") is reserved for special fields holding geometry or validity")
        if varname in self.listfields():
            raise epygramError("there already is a field with the same name in this file")
        if hasattr(field, 'spectral') and field.spectral:
            raise epygramError("Spectral fields cannot be written on this file")
        
        dims = []
        time_dim = self._write_validity_from_field(field)
        if time_dim is not None:
            dims.append(time_dim)
        geom_dims = self._write_geometry_from_field(field)
        if geom_dims is not None:
            dims.extend(geom_dims)
        if isinstance(field, MiscField):
            if miscdims is not None:
                dims = miscdims
            elif len(dims) != len(field.getdata().shape):
                #Guessed dimensions do not fit
                kind = 'char' if numpy.array(field.getdata()).dtype.kind in {'U', 'S'} else 'size'
                dims = [kind + str(d) for d in field.getdata().shape]
            for d in dims:
                if d not in self._nc.dimensions:
                    try:
                        self._nc.createDimension(d, int(d[4:]))
                    except ValueError:
                        epylog.error("This dimension is unknown (if it is a valid standard dimension, " +
                                     "please write first a known field using this dimension)")
                        raise
        comment = {}
        if field.comment not in (None, ''):
            try:
                d = json.loads(field.comment)
                if not isinstance(d, dict):
                    raise ValueError("We only deal with dict")
                comment = d
            except ValueError:
                comment['comment'] = field.comment
        data = field.getdata().squeeze()
        self._nc.createVariable(varname, data.dtype, dimensions=dims,
                                fill_value=comment.get('_FillValue'))
        for k, v in comment.items():
            if not isinstance(k, str):
                raise ValueError("Key must be a string")
            if k != '_FillValue':
                self._nc.variables[varname].setncattr(k, v)
        self._nc.variables[varname][...] = data
        if field.units not in (None, ''):
            self._nc.variables[varname].units = field.units
        (h, v) = field.geometry.position_on_horizontal_grid, field.geometry.vcoordinate.position_on_grid
        if h is not None and v is not None:
            gridIndicator = [key for key, value in gridIndicatorDict.items() if value == (h, v)][0]
            self._nc.variables[varname].grid = numpy.int32(gridIndicator)

#    def rename_field(self, fid, new_fid):
#        """Renames a field "in place"."""
#        wlfi.wlfiren(self._unit, fid if self.true3d else fid[0], new_fid if self.true3d else new_fid[0])
#
#    def delfield(self, fid):
#        """Deletes a field from file "in place"."""
#        wlfi.wlfisup(self._unit, fid if self.true3d else fid[0])
#
    @FileResource._openbeforedelayed
    def extractprofile(self, fid, lon=None, lat=None,
                       geometry=None,
                       vertical_coordinate=None,
                       interpolation='nearest',
                       external_distance=None,
                       cheap_height=None):
        """
        Extracts a vertical profile from the netCDFMNH resource, given its fid
        and the geographic location (*lon*/*lat*) of the profile.

        :param fid: must have syntax: 'PARAMETER' PARAMETER being the name of the
                    parameter requested, as named in the file.
        :param lon: the longitude of the desired point.
        :param lat: the latitude of the desired point.
        :param geometry: the geometry on which extract data. If None, it is built from
          lon/lat.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V1DField (cf. :module:`epygram.geometries` coordinate possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile at requested lon/lat from the fields grid:\n
          - if 'nearest' (default), extracts profile at the horizontal nearest neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear spline interpolation;
          - if 'cubic', computes profile with horizontal cubic spline interpolation.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = |target_value - external_field.data|
        :param cheap_height: has no effect (compatibity with FA format)
        """
        if geometry is None:
            if None in [lon, lat]:
                raise epygramError("You must give a geometry or lon *and* lat")
            if self.geometry is None:
                self._read_geometry()
            pointG = self.geometry.make_profile_geometry(lon, lat)
        else:
            if lon is not None or lat is not None:
                raise epygramError("You cannot provide lon or lat when geometry is given")
            if geometry.structure != "V1D":
                raise epygramError("geometry must be a V1D")
            pointG = geometry

        profile = self.extract_subdomain(fid, pointG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate)

        return profile

    @FileResource._openbeforedelayed
    def extractsection(self, fid, end1=None, end2=None,
                       geometry=None,
                       points_number=None,
                       resolution=None,
                       vertical_coordinate=None,
                       interpolation='linear',
                       cheap_height=None,
                       global_shift_center=None):
        """
        Extracts a vertical section from the netCDFMNH resource, given its fid
        and the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param fid: must have syntax: 'PARAMETER' PARAMETER being the name of the
                    parameter requested, as named in the file.
        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param geometry: the geometry on which extract data. If None, defaults to
          linearily spaced positions computed from  *points_number*.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. :module:`epygram.geometries` coordinate possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: has no effect (compatibity with FA format)
        :param global_shift_center: for global lon/lat grids, shift the center by the
            requested angle (in degrees). Enables a [0,360] grid
            to be shifted to a [-180,180] grid, for instance (with -180 argument).
        """
        if geometry is None:
            if None in [end1, end2]:
                raise epygramError("You must give a geometry or end1 *and* end2")
            if self.geometry is None:
                self._read_geometry()
            sectionG = self.geometry.make_section_geometry(end1, end2,
                                                           points_number=points_number,
                                                           resolution=resolution)
        else:
            if end1 is not None or end2 is not None:
                raise epygramError("You cannot provide end1 or end2 when geometry is given")
            if geometry.structure != "V2D":
                raise epygramError("geometry must be a V2D")
            sectionG = geometry

        section = self.extract_subdomain(fid, sectionG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         global_shift_center=global_shift_center)

        return section

    @FileResource._openbeforedelayed
    def extract_subdomain(self, fid, geometry,
                          vertical_coordinate=None,
                          interpolation='linear',
                          exclude_extralevels=True,
                          cheap_height=None,
                          global_shift_center=None):
        """
        Extracts a subdomain from the netCDFMNH resource, given its fid
        and the geometry to use.

        :param fid: must have syntax: 'PARAMETER' PARAMETER being the name of the
                    parameter requested, as named in the file.
        :param geometry: the geometry on which extract data.
                         None to keep the geometry untouched.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. :module:`epygram.geometries` coordinate possible values).
        :param interpolation defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: has no effect (compatibity with FA format)
        :param global_shift_center: for global lon/lat grids, shift the center by the
            requested angle (in degrees). Enables a [0,360] grid
            to be shifted to a [-180,180] grid, for instance (with -180 argument).
        """
        field3d = self.readfield(fid)
        if global_shift_center is not None:
            field3d.global_shift_center(global_shift_center)

        if geometry is None or geometry == field3d.geometry:
            subdomain = field3d
            if exclude_extralevels:
                subdomain = subdomain.extract_physicallevels()
            geometry = subdomain.geometry
        elif geometry == field3d.geometry.make_physicallevels_geometry():
            subdomain = field3d.extract_physicallevels()
            geometry = subdomain.geometry
        else:
            subdomain = field3d.extract_subdomain(geometry, interpolation=interpolation)
            if exclude_extralevels:
                subdomain = subdomain.extract_physicallevels()

        # vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 118 and \
               vertical_coordinate in (102, 103):
                zsfield = self.readfield('ZS')
                zs_values = zsfield.getvalue_ll(*geometry.get_lonlat_grid(),
                                                interpolation=interpolation, one=False)
                zs_values = zs_values.reshape(geometry.get_datashape(force_dimZ=1))
                subdomain.geometry.vcoordinate = hybridH2altitude(subdomain.geometry.vcoordinate,
                                                                  zs_values,
                                                                  gridposition=subdomain.geometry.vcoordinate.position_on_grid,
                                                                  conv2height=(vertical_coordinate == 103))
            elif subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 118 and \
                   vertical_coordinate == 100:
                pid = ['PABSM', 'PABST']
                try:
                    P = self.extract_subdomain(pid[0], geometry, interpolation=interpolation)
                except:
                    P = self.extract_subdomain(pid[1], geometry, interpolation=interpolation)
#                    
#                try:
#                    P3d = self.readfield('PABSM')
#                except:
#                    P3d = self.readfield('PABST')
#                P = P3d.extract_subdomain(geometry, interpolation=interpolation)
                subdomain.geometry.vcoordinate = hybridH2pressure(subdomain.geometry.vcoordinate,
                                                                  P.getdata(),
                                                                  P.geometry.vcoordinate.position_on_grid)
            else:
                raise NotImplementedError("this vertical coordinate conversion.")

        return subdomain

###########
# pre-app #
###########
    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout,
             details=False,
             sortfields=False,
             **_):
        """
        Writes in file a summary of the contents of the netCDFMNH file.

        :param out: the output open file-like object
        :param sortfields: **True** if the fields have to be sorted by type.
        """
        firstcolumn_width = 50
        secondcolumn_width = 16
        sepline = '{:-^{width}}'.format('', width=firstcolumn_width + secondcolumn_width + 1) + '\n'

        first_MTOField = [f for f in self.listfields() if self._get_type(f) in ['V1D', 'H2D', '3D']][0]
        firstfield = self.readfield(first_MTOField, getdata=False)

        listfields = self.listfields()
        listoffields = listfields
        if sortfields:
            sortedfields = self.sortfields()

        def write_formatted(dest, label, value):
            dest.write('{:<{width}}'.format(label, width=firstcolumn_width) +
                       ':' +
                       '{:>{width}}'.format(str(value), width=secondcolumn_width) +
                       '\n')
        def write_formatted_col(dest, label, value):
            dest.write('{:>{width}}'.format(label, width=firstcolumn_width) +
                       ':' +
                       '{:>{width}}'.format(str(value), width=secondcolumn_width) +
                       '\n')
        def write_formatted_fields(dest, label, gridIndicator=None, comment=None):
            if gridIndicator is None and comment is None:
                dest.write('{:<{width}}'.format(label, width=20) +
                           '\n')
            else:
                dest.write('{:<{width}}'.format(label, width=20) +
                           ':' +
                           '{:^{width}}'.format(str(gridIndicator), width=10) +
                           ':' +
                           comment +
                           '\n')
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")

        firstfield.what(out, vertical_geometry=False, fid=False)

        if self.geometry.vcoordinate is not None:
            self.geometry.vcoordinate.what(out, levels=False)

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        if sortfields:
            listoffields = []
            for k in sorted(sortedfields.keys()):
                listoffields.append(k)
                listoffields.append('--------------------')
                listoffields.extend(sortedfields[k])
                listoffields.append('--------------------')
            numfields = sum([len(v) for v in sortedfields.values()])
        else:
            numfields = len(listoffields)
        out.write("Number: " + str(numfields) + "\n")
        if details is not None:
            write_formatted_fields(out, "Field name", "Grid ind.", "Comment")
        else:
            write_formatted_fields(out, "Field name")
        out.write(sepline)
        done = []
        for f in listoffields:
            if f not in done:
                if details:
                    field = self.readfield(f)
                    if hasattr(field, 'geometry'):
                        gridIndicator = {('__unknown__', '__unknown__'):0,
                                         ('center', 'mass'):1,
                                         ('center-left', 'mass'):2,
                                         ('lower-center', 'mass'):3,
                                         ('center', 'flux'):4,
                                         ('lower-left', 'mass'):5,
                                         ('center-left', 'flux'):6,
                                         ('lower-center', 'flux'):7,
                                         ('lower-left', 'flux'):8}[field.geometry.position_on_horizontal_grid,
                                                                   field.geometry.vcoordinate.position_on_grid]
                    else:
                        gridIndicator = '-'
                    write_formatted_fields(out, f, gridIndicator, field.comment)
                else:
                    write_formatted_fields(out, f)
                done.append(f)
        out.write(sepline)

# the netCDFMNH WAY #
###############
    def get_dimensions(self, fid):
        """
        Returns the dimensions associated with this field identifier
        as a tuple of strings.
        :param fid: field identifier
        """
        return self._nc[fid].dimensions

    @staticmethod
    def _get_latin1_latin2_lambert(lat0, rpk):
        def k(latin2):
            latin1 = lat0
            m1 = math.cos(math.radians(latin1))
            m2 = math.cos(math.radians(latin2))
            t1 = math.tan(math.pi / 4. - math.radians(latin1) / 2.)
            t2 = math.tan(math.pi / 4. - math.radians(latin2) / 2.)
            return (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2)) - rpk
        try:
            import scipy.optimize as op
            latin2 = Angle(op.fsolve(k, math.degrees(2 * math.asin(rpk)) - lat0)[0],
                           'degrees')
            latin1 = Angle(float(lat0), 'degrees')
        except Exception:
            def solve(function, x0):
                """A solver adapted to this problem. Do not try to use it elsewhere!"""
                x1 = x0 + 1.
                x2 = x0
                y1 = function(x1)
                y2 = function(x2)
                while math.fabs(y2) > 10.E-20 and math.fabs(y2) != math.fabs(y1):
                    x1, x2 = x2, x2 - (x2 - x1) / (y2 - y1) * y2
                    y1, y2 = y2, function(x2)
                return x2
            latin2 = Angle(solve(k, math.degrees(2 * math.asin(rpk)) - lat0),
                           'degrees')
            latin1 = Angle(float(lat0), 'degrees')
        return (latin1, latin2)

    @FileResource._openbeforedelayed
    def _read_geometry(self):
        """
        Reads the geometry in the netCDFMNH articles.
        """
        def v(name): return self._nc.variables[name][...]
        def vdefault(name, default): return self._nc.variables[name][...] if name in self._nc.variables else default
        listnames = self._nc.variables.keys()
        
        if 'CARTESIAN' in listnames:
            cartesian = v('CARTESIAN')
        else:
            cartesian = False
        imax = int(v('IMAX'))
        jmax = int(v('JMAX'))
        jphext = int(vdefault('JPHEXT', 1))
        xhat = v('XHAT')
        yhat = v('YHAT')
        if 'KMAX' in listnames:
            kmax = int(v('KMAX'))
            if kmax > 1:
                zhat = v('ZHAT')
                kmax += 2
        else:
            kmax = 0
        dimensions = {'X':1 if (cartesian and imax == 1) else imax + 2 * jphext,
                      'Y':1 if (cartesian and jmax == 1) else jmax + 2 * jphext,
                      'X_CIzone':imax,
                      'Y_CIzone':jmax,
                      'X_Iwidth':0,
                      'Y_Iwidth':0,
                      'X_Czone':imax,
                      'Y_Czone':jmax,
                      'X_CIoffset':jphext,
                      'Y_CIoffset':jphext
                      }

        if cartesian:
            geometryclass = AcademicGeometry
            lat0 = v('LAT0')
            lon0 = v('LON0')
            
            grid = {'X_resolution':xhat[1] - xhat[0],
                    'Y_resolution':yhat[1] - yhat[0],
                    'LAMzone':'CIE',
                    'latitude':Angle(float(lat0), 'degrees'),
                    'longitude':Angle(float(lon0), 'degrees'),
                    'input_lon':1,
                    'input_lat':1,
                    'input_position':(0, 0)
                    }
            projection = {'rotation':Angle(0., 'degrees'),
                          'reference_dX':grid['X_resolution'],
                          'reference_dY':grid['X_resolution']}
            geometryname = 'academic'
            kwargs_geom = dict(structure='3D',
                               name=geometryname,
                               grid=grid,
                               dimensions=dimensions,
                               projection=projection,
                               geoid=config.netCDFMNH_default_geoid,
                               )
        else:
            geometryclass = ProjectedGeometry
            lat0 = v('LAT0')
            lon0 = v('LON0')
            xres, yres = xhat[1] - xhat[0], yhat[1] - yhat[0]
            if 'latitude' in listnames and 'longitude' in listnames:
                input_position = (0, 0)
                lat1 = v('latitude')[input_position[0], input_position[1]]
                lon1 = v('longitude')[input_position[0], input_position[1]]
            elif 'LATORI' in listnames and 'LONORI' in listnames:
                lat1 = v('LATORI')
                lon1 = v('LONORI')
                xoffset, yoffset = .5 * (xhat[0] + xhat[1]) / xres, .5 * (yhat[0] + yhat[1]) / yres
                input_position = (jphext - 1 - xoffset, jphext - 1 - yoffset) #deduced from run output
            else:
                lat1 = v('LATOR')
                lon1 = v('LONOR')
                input_position = (0, 0)
            rpk = v('RPK')
            beta = v('BETA')

            projection = {'reference_lon':Angle(float(lon0), 'degrees'),
                          'rotation': Angle(float(beta), 'degrees')
                          }
            if abs(rpk - math.sin(math.radians(lat0))) <= config.epsilon:
                # non secant
                projection['reference_lat'] = Angle(float(lat0), 'degrees')
            else:
                if abs(rpk) in [0., 1.]:
                    # mercator or polar stereographic: one secant latitude
                    projection['reference_lat'] = Angle(float(numpy.copysign(90, lat0)), 'degrees')
                    projection['secant_lat'] = Angle(float(lat0), 'degrees')
                else:
                    # lambert: two secant latitudes
                    latin1, latin2 = self._get_latin1_latin2_lambert(lat0, rpk)
                    projection['secant_lat1'] = latin1
                    projection['secant_lat2'] = latin2
            grid = {'X_resolution':xres,
                    'Y_resolution':yres,
                    'LAMzone':'CIE',
                    'input_lon':Angle(float(lon1), 'degrees'),
                    'input_lat':Angle(float(lat1), 'degrees'),
                    'input_position':input_position,
                    }
            if abs(rpk) <= config.epsilon:
                geometryname = 'mercator'
            elif abs(1 - abs(rpk)) <= config.epsilon:
                geometryname = 'polar_stereographic'
            else:
                geometryname = 'lambert'

            kwargs_geom = dict(name=geometryname,
                               grid=grid,
                               dimensions=dimensions,
                               geoid=config.netCDFMNH_default_geoid,
                               projection=projection
                               )

        if kmax > 1:
            H = zhat[-1]
            if 'SLEVE' in listnames:
                sleve = v('SLEVE')
            else:
                sleve = False
            Ai = [c for c in zhat[0:kmax + 2]][1:]
            Bi = [1 - c / H for c in zhat[0:kmax + 2]][1:]
            grid = {'gridlevels': tuple([(i + 1, FPDict({'Ai':Ai[i], 'Bi':Bi[i]})) for
                                         i in range(len(Ai))]),
                    'ABgrid_position':'flux'}
            kwargs_vcoord = {'typeoffirstfixedsurface':118 if not sleve else 255,
                             'position_on_grid': 'mass',
                             'grid': grid,
                             'levels': list([i for i in range(len(Ai) + 1)])
                             }
        else:
            kwargs_vcoord = {'typeoffirstfixedsurface': 255,
                             'position_on_grid': '__unknown__',
                             'levels':[255]}
        kwargs_geom['position_on_horizontal_grid'] = 'center'
        kwargs_geom['vcoordinate'] = VGeometry(**kwargs_vcoord)
        self.geometry = geometryclass(**kwargs_geom)

    @FileResource._openbeforedelayed
    def _read_validity(self):
        """Reads the validity in the netCDFMNH articles."""
        listnames = self._nc.variables.keys()
        def todate(fieldname):
            var = self._nc.variables[fieldname]
            if var.units.startswith("seconds since ") and \
               var.units.endswith(" +0:00"):
                value = datetime.datetime.strptime(var.units, "seconds since %Y-%m-%d %H:%M:%S +0:00")
                value += datetime.timedelta(seconds=float(var[...]))
            else:
                raise epygramError("Unknown unit for time: " + var.units)
            return value
        kwargs = {}
        if 'DTEXP' in listnames:
            kwargs['basis'] = todate('DTEXP')
            if 'DTCUR' in listnames:
                kwargs['term'] = todate('DTCUR') - kwargs['basis']
        elif 'DTCUR' in listnames:
            kwargs['date_time'] = todate('DTCUR')
        kwargs['cumulativeduration'] = datetime.timedelta(seconds=0)
        self.validity = FieldValidity(**kwargs)

    def _write_special_records(self, records, field):
        """
        Write special records that control geometry and validity
        """
        jphext = records.get('JPHEXT', None) #not available when writing special fields representing the validity
        for name, value in records.items():
            meta = self._specialFieldComments[name]
            if name in ['DTEXP', 'DTCUR', 'time']:
                assert isinstance(value, tuple)
                assert len(value) == 2
                assert isinstance(value[0], datetime.datetime) and isinstance(value[1], datetime.timedelta)
                units = value[0].strftime("seconds since %Y-%m-%d %H:%M:%S +0:00")
                if name in self._nc.variables:
                    if self._nc.variables[name].units != units:
                        raise epygramError(name + " already exists in file with different value")
                meta['units'] = units
                value = value[1].total_seconds()
            elif name in ['level']:
                meta['c_grid_dynamic_range'] = '2:' + str(len(value) - 1)
            elif name in ['level_w']:
                meta['c_grid_dynamic_range'] = '2:' + str(len(value))
            elif name in ['ni', 'ni_v']:
                meta['c_grid_dynamic_range'] = str(jphext + 1) + ':' + str(len(value) - jphext)
            elif name in ['ni_u']:
                meta['c_grid_dynamic_range'] = str(jphext + 1) + ':' + str(len(value))
            elif name in ['nj', 'nj_u']:
                meta['c_grid_dynamic_range'] = str(jphext + 1) + ':' + str(len(value) - jphext)
            elif name in ['nj_v']:
                meta['c_grid_dynamic_range'] = str(jphext + 1) + ':' + str(len(value))
            if numpy.array(value).dtype == bool:
                value = numpy.where(numpy.array(value), numpy.array(1, dtype=numpy.int8),
                                                        numpy.array(0, dtype=numpy.int8))
            if name in self._nc.variables:
                #Variable already stored in file, we check equality
                recordedValue = self._nc.variables[name][...]
                if name in ['CARTESIAN', 'SLEVE'] + ['IMAX', 'JMAX', 'KMAX', 'JPHEST']:
                    #boolean or integer
                    check = self._nc.variables[name][...] == value
                elif name in ['LAT0', 'LON0', 'LATOR', 'LATORI', 'LONOR', 'LONORI',
                              'RPK', 'BETA', 'ZHAT', 'DTMOD', 'DTCUR']:
                    # Float comparisons
                    special_rpk = (name == 'RPK' and
                                   field.geometry.secant_projection and
                                   field.geometry.name not in ['mercator', 'polar_stereographic'])
                    if special_rpk:
                        # In geometries, when seant, we store latin1 and latin2 which are computed from RPK
                        # Computation is not exact, computing back RPK from latin1 and latin2 does not give exactly the same result
                        latin1_field = field.geometry.projection['secant_lat1'].get('degrees')
                        latin2_field = field.geometry.projection['secant_lat2'].get('degrees')
                        latin1, latin2 = self._get_latin1_latin2_lambert(self._nc.variables['LAT0'][...], self._nc.variables['RPK'][...])
                        check = numpy.all(util.nearlyEqualArray([latin1_field, latin2_field],
                                                                [latin1.get('degrees'), latin2.get('degrees')])) or \
                                util.nearlyEqual(value, recordedValue)
                    else:
                        check = numpy.all(util.nearlyEqualArray(value, self._nc.variables[name][...]))
                elif name in ['XHAT', 'YHAT']:
                    # We check deltaX and deltaY because XHAT and YHAT can be computed in two different ways (prep_ideal or prep_pgd)
                    check = util.nearlyEqualArray(value[1] - value[0],
                                                  recordedValue[1] - recordedValue[0])
                else:
                    check = numpy.all(value == recordedValue)
                if not check:
                    raise epygramError(name + " already exists in file with different value")
            else:
                #Variable not in file, we write it
                dimensions = dict(time=('time', ),
                                  XHAT=('ni_u', ), YHAT=('nj_v', ),
                                  ZHAT=('level_w', ), level=('level',), level_w=('level_w',),
                                  LAT=('nj', 'ni'), LON=('nj', 'ni'),
                                  ni=('ni',), nj=('nj',), ni_u=('ni_u',), nj_u=('nj_u',),
                                  ni_v=('ni_v',), nj_v=('nj_v',),
                                  latitude=('nj', 'ni'), longitude=('nj', 'ni'),
                                  latitude_u=('nj_u', 'ni_u'), longitude_u=('nj_u', 'ni_u'),
                                  latitude_v=('nj_v', 'ni_v'), longitude_v=('nj_v', 'ni_v'),
                                  latitude_f=('nj_v', 'ni_u'), longitude_f=('nj_v', 'ni_u')).get(name, tuple())
                self._nc.createVariable(name, numpy.array(value).dtype,
                                        dimensions=dimensions, fill_value=meta.get('_FillValue'))
                self._nc[name][...] = value
                for k, v in meta.items():
                    if k != '_FillValue':
                        self._nc[name].setncattr(k, v)

    def _write_validity_from_field(self, field):
        """
        Write special records needed to represent the validity
        returns the dimension to use to represent the field
        """
        specialFieldValues = dict()
        dim = None
        if hasattr(field, 'validity') and field.validity is not None:
            if len(field.validity) != 1:
                raise epygramError("netCDFMNH can hold only one validity.")
            basis = field.validity.getbasis()
            if basis is not None:
                dim = 'time'
                reference = datetime.datetime(basis.year, basis.month, basis.day)
                specialFieldValues['DTEXP'] = (reference, basis - reference)
                validity = field.validity.get()
                if validity is not None:
                    specialFieldValues['DTCUR'] = (reference, validity - reference)
                    specialFieldValues['time'] = (reference, validity - reference)
                else:
                    specialFieldValues['time'] = (reference, basis - reference)
                if not dim in self._nc.dimensions:
                    self._nc.createDimension(dim, 1)
        self._write_special_records(specialFieldValues, field)
        return dim
    
    def _write_geometry_from_field(self, field):
        """
        Write special records needed to represent the geometry
        returns the dimension to use for writing the field
        """
        def f2m(f):
            m = numpy.empty_like(f)
            m[:-1] = 0.5 * (f[:-1] + f[1:])
            m[-1] = f[-1] + 0.5 * (f[-1] - f[-2])
            return m
        specialFieldValues = dict()
        g = field.geometry
        dim = []
        
        #Vertical grid
        if g.vcoordinate is not None and \
            hasattr(g.vcoordinate, 'grid') and \
           'gridlevels' in g.vcoordinate.grid and \
           g.structure in ['3D', 'H2D', 'V1D', 'H1D', 'V2D']:
            specialFieldValues['SLEVE'] = g.vcoordinate.typeoffirstfixedsurface != 118
            kmax = len(g.vcoordinate.grid['gridlevels']) + 1
            if kmax > 1:
                kmax -= 2
            specialFieldValues['KMAX'] = numpy.int32(kmax)
            Ai = [level[1]['Ai'] for level in g.vcoordinate.grid['gridlevels']]
            Ai = [-Ai[1]] + Ai
            specialFieldValues['ZHAT'] = numpy.array(Ai)
            if not 'level' in self._nc.dimensions:
                self._nc.createDimension('level', kmax + 2)
            if not 'level_w' in self._nc.dimensions:
                self._nc.createDimension('level_w', kmax + 2)
            specialFieldValues['level_w'] = specialFieldValues['ZHAT']
            specialFieldValues['level'] = f2m(specialFieldValues['level_w'])
            specialFieldValues['ZTOP'] = specialFieldValues['ZHAT'][-1]
            if g.vcoordinate.grid['ABgrid_position'] != 'flux':
                raise epygramError("Don't know how to deal with ABgrid_position!='flux'")
            if field.geometry.vcoordinate.position_on_grid == 'flux':
                dim.append('level_w')
            else:
                dim.append('level')
        
        #Horizontal grid
        if hasattr(g, 'dimensions') and hasattr(g, 'grid'):
            jphext = g.dimensions['X_CIoffset']
            specialFieldValues['JPHEXT'] = jphext
            assert g.dimensions['Y_CIoffset'] == jphext, 'JPHEXT must have the same value on both horizontal dimensions'
            #Meso-NH format is not consistent. dimensions['X']/dimensions['Y'] include
            #the jphext points except if it has a length of one
            #dimX/dimY always contain the jphext points
            dimX = (1 + 2 * jphext) if g.dimensions['X'] == 1 else g.dimensions['X']
            dimY = (1 + 2 * jphext) if g.dimensions['Y'] == 1 else g.dimensions['Y']
            specialFieldValues['IMAX'] = numpy.int32(dimX - 2 * jphext)
            specialFieldValues['JMAX'] = numpy.int32(dimY - 2 * jphext)
            specialFieldValues['XHAT'] = numpy.arange(1 / 2. - jphext,
                                                      dimX - jphext,
                                                      1.) * g.grid['X_resolution']
            for d in ['ni', 'ni_u', 'ni_v']:
                if d not in self._nc.dimensions:
                    self._nc.createDimension(d, dimX)    
            specialFieldValues['ni_u'] = specialFieldValues['XHAT']
            specialFieldValues['ni'] = f2m(specialFieldValues['ni_u'])
            specialFieldValues['ni_v'] = specialFieldValues['ni']
            specialFieldValues['YHAT'] = numpy.arange(1 / 2. - jphext,
                                                      dimY - jphext,
                                                      1.) * g.grid['Y_resolution']
            for d in ['nj', 'nj_u', 'nj_v']:
                if d not in self._nc.dimensions:
                    self._nc.createDimension(d, dimY)    
            specialFieldValues['nj_v'] = specialFieldValues['YHAT']
            specialFieldValues['nj'] = f2m(specialFieldValues['nj_v'])
            specialFieldValues['nj_u'] = specialFieldValues['nj']
            if g.structure in ('Point', 'V1D'):
                #Meso-NH file does not contain horizontal dimensions for V1D fields
                pass
            elif g.structure in ('H1D', 'V2D'):
                #Meso-NH file does not contain y-dimension for V2D fields
                dim.extend({'center':('ni',), #1, 4
                            'center-left':('ni_u',), #2, 6
                            'lower-center':('ni_v',), #3, 7
                            'lower-left':('ni_u'), #5, 8
                           }[field.geometry.position_on_horizontal_grid])
                if 'size1' not in self._nc.dimensions:
                    self._nc.createDimension('size1', 1)
            else:
                dim.extend({'center':('nj', 'ni'), #1, 4
                            'center-left':('nj_u', 'ni_u'), #2, 6
                            'lower-center':('nj_v', 'ni_v'), #3, 7
                            'lower-left':('nj_v', 'ni_u'), #5, 8
                           }[field.geometry.position_on_horizontal_grid])
            specialFieldValues['LON'],  specialFieldValues['LAT'] = g.get_lonlat_grid(position='center')
            specialFieldValues['longitude'] = specialFieldValues['LON']
            specialFieldValues['latitude'] = specialFieldValues['LAT']
            specialFieldValues['longitude_u'], specialFieldValues['latitude_u'] = g.get_lonlat_grid(position='center-left')
            specialFieldValues['longitude_v'], specialFieldValues['latitude_v'] = g.get_lonlat_grid(position='lower-center')
            specialFieldValues['longitude_f'], specialFieldValues['latitude_f'] = g.get_lonlat_grid(position='lower-left')

        #Projection
        if hasattr(g, 'name') and g.name == 'academic' and hasattr(g, 'grid'):
            specialFieldValues['BETA'] = 0.
            specialFieldValues['CARTESIAN'] = True
            if 'latitude' in g.grid:
                specialFieldValues['LAT0'] = g.grid['latitude'].get('degrees')
            else:
                specialFieldValues['LAT0'] = 0.
            if 'longitude' in g.grid:
                specialFieldValues['LON0'] = g.grid['longitude'].get('degrees')
            else:
                specialFieldValues['LON0'] = 0.
        elif hasattr(g, 'name') and g.name != 'academic' and hasattr(g, 'projection') and hasattr(g, 'grid'):
            specialFieldValues['BETA'] = g.projection['rotation'].get('degrees')
            specialFieldValues['CARTESIAN'] = False
            specialFieldValues['LON0'] = g.projection['reference_lon'].get('degrees')
            if not g.secant_projection:
                specialFieldValues['LAT0'] = g.projection['reference_lat'].get('degrees')
                specialFieldValues['RPK'] = math.sin(g.projection['reference_lat'].get('radians'))
            else:
                if g.name in ['mercator', 'polar_stereographic']:
                    specialFieldValues['LAT0'] = g.projection['secant_lat'].get('degrees')
                    if g.name == 'mercator':
                        specialFieldValues['RPK'] = 0.
                    else:
                        specialFieldValues['RPK'] = numpy.copysign(1, specialFieldValues['LAT0'])
                else:
                    latin1 = g.projection['secant_lat1'].get('degrees')
                    latin2 = g.projection['secant_lat2'].get('degrees')
                    m1 = math.cos(math.radians(latin1))
                    m2 = math.cos(math.radians(latin2))
                    t1 = math.tan(math.pi / 4. - math.radians(latin1) / 2.)
                    t2 = math.tan(math.pi / 4. - math.radians(latin2) / 2.)
                    specialFieldValues['LAT0'] = latin1
                    specialFieldValues['RPK'] = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
            if g.grid['input_position'] == (jphext - 1, jphext - 1):
                specialFieldValues['LONORI'] = g.grid['input_lon'].get('degrees')
                specialFieldValues['LATORI'] = g.grid['input_lat'].get('degrees')
            if g.grid['input_position'] == (0, 0):
                specialFieldValues['LONOR'] = g.grid['input_lon'].get('degrees')
                specialFieldValues['LATOR'] = g.grid['input_lat'].get('degrees')
            if any([k not in specialFieldValues for k in ['LONORI', 'LATORI', 'LONOR', 'LATOR']]):
                lon, lat = g.get_lonlat_grid(position='center')
                if 'LONOR' not in specialFieldValues:
                    specialFieldValues['LONOR'] = lon[0, 0]
                    specialFieldValues['LATOR'] = lat[0, 0]
                if 'LONORI' not in specialFieldValues:
                    specialFieldValues['LONORI'] = lon[jphext - 1, jphext - 1]
                    specialFieldValues['LATORI'] = lat[jphext - 1, jphext - 1]

        self._write_special_records(specialFieldValues, field)
        return dim



