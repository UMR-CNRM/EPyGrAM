#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains classes for netCDF4 resource.
"""

__all__ = ['netCDF']

import datetime
from dateutil import parser as dtparser
import copy
import json
import sys

import footprints
from footprints import proxy as fpx, FPDict

from epygram import config, epygramError, util
from epygram.base import FieldValidity
from epygram.resources import FileResource
from epygram.fields import H2DField
from epygram.util import stretch_array

import netCDF4

epylog = footprints.loggers.getLogger(__name__)



class netCDF(FileResource):
    """Class implementing all specificities for netCDF (4) resource format."""

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['netCDF']),
                default='netCDF'),
            behaviour=dict(
                info="Describes how fields are defined in resource.",
                type=FPDict,
                optional=True,
                default=config.netCDF_default_behaviour)
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        self.isopen = False
        super(netCDF, self).__init__(*args, **kwargs)

        if self.openmode in ('r', 'a'):
            try:
                guess = netCDF4.Dataset(self.container.abspath, self.openmode)
            except RuntimeError:
                raise IOError("this resource is not a netCDF one.")
            else:
                guess.close()
        behaviour = copy.copy(config.netCDF_default_behaviour)
        behaviour.update(self.behaviour)
        self._attributes['behaviour'] = behaviour
        if not self.fmtdelayedopen:
            self.open()

    def open(self, openmode=None):
        """
        Opens a netCDF and initializes some attributes.

        - *openmode*: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """

        super(netCDF, self).open(openmode=openmode)
        self._nc = netCDF4.Dataset(self.container.abspath, self.openmode)
        self.isopen = True

    def close(self):
        """
        Closes a netCDF.
        """

        if hasattr(self, '_nc') and self._nc._isopen:
            self._nc.close()
        self.isopen = False

    def variables_number(self):
        """Return the number of variables in resource."""
        return len(self._variables)

    def find_fields_in_resource(self, seed=None, generic=False, **kwargs):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        Args: \n
        - *seed*: might be a regular expression, a list of regular expressions
          or *None*. If *None* (default), returns the list of all fields in
          resource.
        - *fieldtype*: optional, among ('H2D', 'Misc') or a list of these strings.
          If provided, filters out the fields not of the given types.
        - *generic*: if True, returns complete fid's,
          union of {'FORMATname':fieldname} and the according generic fid of
          the fields.
        """

        if seed == None:
            fieldslist = self.listfields()
        elif isinstance(seed, str):
            fieldslist = util.find_re_in_list(seed, self.listfields())
        elif isinstance(seed, list):
            fieldslist = []
            for s in seed:
                fieldslist += util.find_re_in_list(s, self.listfields())
        if fieldslist == []:
            raise epygramError("no field matching '" + seed + \
                               "' was found in resource " + \
                               self.container.abspath)
        if generic:
            raise NotImplementedError("not yet.")

        return fieldslist

    def _listfields(self):
        """Returns the fid list of the fields inside the resource."""
        return self._variables.keys()

    @FileResource._openbeforedelayed
    def readfield(self, fid,
                  getdata=True):
        """
        Reads one field, given its netCDF name, and returns a Field instance.
        
        Args: \n
        - *fid*: netCDF field identifier
        - *getdata*: if *False*, only metadata are read, the field do not
          contain data.
        """

        if self.openmode == 'w':
            raise epygramError("cannot read fields in resource if with" + \
                               " openmode == 'w'.")
        assert fid in self.listfields(), ' '.join(["field",
                                                   fid,
                                                   "not found in resource."])
        field_kwargs = {'fid':{'netCDF':fid}}

        # geometry
        dimensions = {d:len(self._dimensions[d]) for d in self._variables[fid].dimensions}
        geometryname = 'unstructured'
        if set(self._variables[fid].dimensions) == set(self.behaviour['H2D_dimensions_names']):
            read_as_miscfield = False
            # this is a H2D field
            structure = 'H2D'
            lons = self._variables[self.behaviour['variable_name_for_longitudes']][:, :]
            lats = self._variables[self.behaviour['variable_name_for_latitudes']][:, :]
            lons_dim = self._variables[self.behaviour['variable_name_for_longitudes']].dimensions[:]
            var_dim = [self._variables[fid].dimensions[i:i + len(lons_dim)]
                       for i in range(0, len(self._variables[fid].dimensions) - len(lons_dim) + 1)]
            assert lons_dim in var_dim, \
                   "lons/lats and variable " + fid + " dimensions mismatch"
            grid = {'longitudes':lons,
                    'latitudes':lats}
        else:
            epylog.warning("unable to assume geometry of field. Read as MiscField.")
            read_as_miscfield = True

        # validity
        if self.behaviour.get('variable_name_for_validity') in self._variables[fid].dimensions:
            # temporal dimension
            _validity = self._variables[self.behaviour['variable_name_for_validity']]
            raise NotImplementedError('temporal dimension of field: not yet !')
        elif 'validity' in self._variables[fid].ncattrs():
            # validity stored as an attribute of variable (as in writefield !)
            try:
                _validity = json.loads(self._variables[fid].validity)
                _validity['basis'] = dtparser.parse(_validity['basis'])
                _validity['date_time'] = dtparser.parse(_validity['date_time'])
                if _validity.get('cumulativeduration') is not None:
                    _validity['cumulativeduration'] = datetime.timedelta(seconds=float(_validity['cumulativeduration']))
            except (KeyError, ValueError):
                epylog.warning("unable to decode validity attribute.")
                validity = FieldValidity()
            else:
                validity = FieldValidity(**_validity)
        else:
            validity = FieldValidity()

        # build field
        if not read_as_miscfield:
            field_kwargs['validity'] = validity
            kwargs_geom = {'structure':structure,
                           'name':geometryname,
                           'grid':grid,
                           'dimensions':dimensions,
                           'vcoordinate':{'structure':'V',
                                          'typeoffirstfixedsurface':255,
                                          'levels':[0]},
                           'position_on_horizontal_grid':'center'}
            geometry = fpx.geometry(**kwargs_geom)
            field_kwargs['geometry'] = geometry
            field_kwargs['structure'] = structure
        comment = {}
        for a in self._variables[fid].ncattrs():
            if read_as_miscfield or (not read_as_miscfield and a != 'validity'):
                comment.update({a:self._variables[fid].getncattr(a)})
        comment = json.dumps(comment)
        if comment != '{}':
            field_kwargs['comment'] = comment
        field = fpx.field(**field_kwargs)
        if getdata:
            field.setdata(self._variables[fid][...])

        return field

    def writefield(self, field, compression=4, metadata={}):
        """
        Write a field in resource.
        Args:\n
        - *compression* ranges from 1 (low compression, fast writing)
          to 9 (high compression, slow writing). 0 is no compression.
        - *metadata* can be filled by any meta-data, that will be stored
          as attribute of the netCDF variable.
        """

        def check_or_add_dim(k, size=None):
            if size is None:
                size = field.geometry.dimensions[k]
            if k not in self._dimensions:
                self._nc.createDimension(k, size=size)
            else:
                assert len(self._dimensions[k]) == size, \
                       "dimensions mismatch: " + k + ": " + \
                       str(self._dimensions[k]) + " != " + str(size)

        vartype = 'f8'
        if isinstance(field, H2DField):
            # dimensions
            if field.geometry.rectangular_grid:
                # rectangular grid case
                dims = (k for k in field.geometry.dimensions.keys() if len(k) == 1)
                dims = sorted(dims, reverse=(not self.behaviour['transpose_data_ordering']))
                for k in dims:
                    check_or_add_dim(k)
            elif not self.behaviour['flatten_non_rectangular_grids']:
                fill_value = -999999.9
                if 'gauss' in field.geometry.name:
                    dims = ('lat_number', 'max_lon_number')
                    for k in dims:
                        check_or_add_dim(k)
                else:
                    raise NotImplementedError("grid not rectangular nor a gauss one.")
            else:
                # gauss case
                if 'gauss' in field.geometry.name:
                    check_or_add_dim('lat_number')
                    check_or_add_dim('max_lon_number')
                    gn = sum(field.geometry.dimensions['lon_number_by_lat'])
                    check_or_add_dim('gridpoints_number', size=gn)
                    self._nc.lon_number_by_lat = field.geometry.dimensions['lon_number_by_lat']
                    dims = ('gridpoints_number',)
                else:
                    raise NotImplementedError("grid not rectangular nor a gauss one.")

            # geometry (lons/lats)
            (lons, lats) = field.geometry.get_lonlat_grid()
            if self.behaviour['transpose_data_ordering']:
                lons = lons.transpose()
                lats = lats.transpose()
            if self.behaviour['variable_name_for_longitudes'] in self._variables:
                if field.geometry.rectangular_grid or not self.behaviour['flatten_non_rectangular_grids']:
                    lons_ok = lons.shape == self._variables[self.behaviour['variable_name_for_longitudes']].shape
                else:
                    # gauss case, flattened
                    gn = sum(field.geometry.dimensions['lon_number_by_lat'])
                    lons_ok = (gn,) == self._variables[self.behaviour['variable_name_for_longitudes']].shape
                assert lons_ok, "dimensions mismatch: lons grid."
            else:
                lons_var = self._nc.createVariable(self.behaviour['variable_name_for_longitudes'],
                                                   vartype, dims)
                if field.geometry.rectangular_grid:
                    lons_var[...] = lons
                elif not self.behaviour['flatten_non_rectangular_grids']:
                    lons_var.missing_value = fill_value
                    lons_var[...] = lons.filled(fill_value)
                else:
                    lons_var[...] = stretch_array(lons)
            if self.behaviour['variable_name_for_latitudes'] in self._variables:
                if field.geometry.rectangular_grid or not self.behaviour['flatten_non_rectangular_grids']:
                    lats_ok = lats.shape == self._variables[self.behaviour['variable_name_for_latitudes']].shape
                else:
                    # gauss case, flattened
                    gn = sum(field.geometry.dimensions['lon_number_by_lat'])
                    lats_ok = (gn,) == self._variables[self.behaviour['variable_name_for_latitudes']].shape
                assert lats_ok, "dimensions mismatch: lats grid."
            else:
                lats_var = self._nc.createVariable(self.behaviour['variable_name_for_latitudes'],
                                                   vartype, dims)
                if field.geometry.rectangular_grid:
                    lats_var[...] = lats
                elif not self.behaviour['flatten_non_rectangular_grids']:
                    lats_var.missing_value = fill_value
                    lats_var[...] = lats.filled(fill_value)
                else:
                    lats_var[...] = stretch_array(lats)

            # validity
            if len(field.validity) == 1:
                validity = {'basis':None, 'date_time':None}
                if field.validity[0].getbasis() is not None:
                    validity['basis'] = field.validity[0].getbasis().isoformat()
                if field.validity[0].get() is not None:
                    validity['date_time'] = field.validity[0].get().isoformat()
                if field.validity[0].cumulativeduration() is not None:
                    validity['cumulativeduration'] = str(field.validity[0].cumulativeduration().total_seconds())
            elif len(field.validity) > 1:
                raise NotImplementedError("not yet !")
                #TODO: create a 'time' dimension
                #if 'time' not in self._dimensions:
                #    self._nc.createDimension('time', size=None)
                #    for v in field.validity:
                #        self._nc.dimensions['time'].append(v.get('IntStr'))
                #else:
                #    # check that time dimension is compatible ?
                #    raise NotImplementedError("not yet !")

            # create variable
            var = self._nc.createVariable(util.linearize2str(field.fid.get('netCDF', field.fid)),
                                          vartype, dims,
                                          zlib=bool(compression), complevel=compression)
            # set metadata
            if len(field.validity) == 1:
                var.validity = json.dumps(validity)
            if field.comment is not None:
                metadata.update(comment=field.comment)
            for k, v in metadata.items():
                setattr(var, k, v)
            # set data
            if self.behaviour['transpose_data_ordering']:
                data = field.data.transpose()
            else:
                data = field.data
            if field.geometry.rectangular_grid:
                var[...] = data
            elif not self.behaviour['flatten_non_rectangular_grids']:
                var.missing_value = fill_value
                var[...] = data.filled(fill_value)
            else:
                var[...] = stretch_array(data)
        else:
            raise NotImplementedError("not yet !")

    def behave(self, **kwargs):
        """
        Set-up the given arguments in self.behaviour, for the purpose of
        building fields from netCDF.
        """

        self.behaviour.update(kwargs)

    def what(self, out=sys.stdout):
        """Writes in file a summary of the contents of the GRIB."""

        # adapted from http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
        def ncdump(nc, out):
            '''
            ncdump outputs dimensions, variables and their attribute information.
            The information is similar to that of NCAR's ncdump utility.
            ncdump requires a valid instance of Dataset.
        
            Parameters
            ----------
            nc : netCDF4.Dataset
                A netCDF4 dateset object
            verb : Boolean
                whether or not nc_attrs, nc_dims, and nc_vars are printed
        
            Returns
            -------
            nc_attrs : list
                A Python list of the NetCDF file global attributes
            nc_dims : list
                A Python list of the NetCDF file dimensions
            nc_vars : list
                A Python list of the NetCDF file variables
            '''

            def outwrite(*items):
                items = list(items)
                stritem = items.pop(0)
                for i in items:
                    stritem += ' ' + str(i)
                out.write(stritem + '\n')

            def print_ncattr(key):
                """
                Prints the NetCDF file attributes for a given key
        
                Parameters
                ----------
                key : unicode
                    a valid netCDF4.Dataset.variables key
                """
                try:
                    outwrite("\t\ttype:", repr(nc.variables[key].dtype))
                    for ncattr in nc.variables[key].ncattrs():
                        outwrite('\t\t%s:' % ncattr, \
                                 repr(nc.variables[key].getncattr(ncattr)))
                except KeyError:
                    outwrite("\t\tWARNING: %s does not contain variable attributes" % key)

            # NetCDF global attributes
            nc_attrs = nc.ncattrs()
            outwrite("NetCDF Global Attributes:")
            for nc_attr in nc_attrs:
                outwrite('\t%s:' % nc_attr, repr(nc.getncattr(nc_attr)))
            nc_dims = [dim for dim in nc.dimensions]  # list of nc dimensions
            # Dimension shape information.
            outwrite("NetCDF dimension information:")
            for dim in nc_dims:
                outwrite("\tName:", dim)
                outwrite("\t\tsize:", len(nc.dimensions[dim]))
                #print_ncattr(dim)
            # Variable information.
            nc_vars = [var for var in nc.variables]  # list of nc variables
            outwrite("NetCDF variable information:")
            for var in nc_vars:
                if var not in nc_dims:
                    outwrite('\tName:', var)
                    outwrite("\t\tdimensions:", nc.variables[var].dimensions)
                    outwrite("\t\tsize:", nc.variables[var].size)
                    print_ncattr(var)
            return nc_attrs, nc_dims, nc_vars

        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        ncdump(self._nc, out)

    @property
    @FileResource._openbeforedelayed
    def _dimensions(self):
        return self._nc.dimensions

    @property
    @FileResource._openbeforedelayed
    def _variables(self):
        return self._nc.variables

