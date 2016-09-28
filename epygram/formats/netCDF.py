#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains classes for netCDF4 resource.
"""

__all__ = ['netCDF']

import copy
import json
import sys
import numpy
from collections import OrderedDict

import netCDF4

import footprints
from footprints import proxy as fpx, FPDict, FPList

from epygram import config, epygramError, util
from epygram.base import FieldValidity, FieldValidityList
from epygram.resources import FileResource
from epygram.util import stretch_array, Angle, nearlyEqual

epylog = footprints.loggers.getLogger(__name__)


_typeoffirstfixedsurface_dict = {'altitude':102,
                                 'height':103,
                                 'hybrid-pressure':119,
                                 'hybrid-height':118,
                                 'pressure':100}
_typeoffirstfixedsurface_dict_inv = {v:k for k, v in _typeoffirstfixedsurface_dict.items()}


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
    def ncinfo_field(self, fid):
        """
        Get info about the field (dimensions and meta-data of the netCDF variable).
        
        Args: \n
        - *fid*: netCDF field identifier
        """

        assert fid in self.listfields(), 'field: ' + fid + ' not in resource.'
        dimensions = OrderedDict()
        for d in self._variables[fid].dimensions:
            dimensions[d] = len(self._dimensions[d])
        metadata = {a:getattr(self._variables[fid], a) for a in self._variables[fid].ncattrs()}

        return {'dimensions':dimensions,
                'metadata':metadata}


    @FileResource._openbeforedelayed
    def readfield(self, fid,
                  getdata=True,
                  only=None,
                  adhoc_behaviour=None):
        """
        Reads one field, given its netCDF name, and returns a Field instance.
        
        Args: \n
        - *fid*: netCDF field identifier
        - *getdata*: if *False*, only metadata are read, the field do not
          contain data.
        - *only*: to specify indexes [0 ... n-1] of specific dimensions,
          e.g. {'time':5,} to select only the 6th term of time dimension.
        - *adhoc_behaviour*: to specify "on the fly" a behaviour (usual
          dimensions or grids).
        """

        # 0. initialization
        assert self.openmode != 'w', \
               "cannot read fields in resource if with openmode == 'w'."
        assert fid in self.listfields(), \
               ' '.join(["field", fid, "not found in resource."])
        only = util.ifNone_emptydict(only)
        adhoc_behaviour = util.ifNone_emptydict(adhoc_behaviour)
        field_kwargs = {'fid':{'netCDF':fid}}
        variable = self._variables[fid]
        behaviour = self.behaviour.copy()
        behaviour.update(adhoc_behaviour)
        return_Yaxis = False

        # 1.1 identify usual dimensions
        variable_dimensions = {d:len(self._dimensions[d]) for d in variable.dimensions}
        for d in variable_dimensions.keys():
            for sd in config.netCDF_standard_dimensions:
                # if behaviour is not explicitly given,
                # try to find out who is "d" among the standard dimensions
                if not sd in behaviour.keys() and d in config.netCDF_usualnames_for_standard_dimensions[sd]:
                    behaviour[sd] = d
        dims_dict_n2e = {}
        for d in variable_dimensions.keys():
            for k in config.netCDF_standard_dimensions:
                if d == behaviour.get(k):
                    dims_dict_n2e[d] = k
        dims_dict_e2n = {v:k for k, v in dims_dict_n2e.items()}

        # 1.2 try to identify grids
        for f in self.listfields():
            for sd in config.netCDF_standard_dimensions:
                sg = sd.replace('dimension', 'grid')
                # if behaviour is not explicitly given,
                # try to find out who is "f" among the standard grids
                if not sg in behaviour.keys() and f in config.netCDF_usualnames_for_standard_dimensions[sd] \
                or f == behaviour.get(sd):
                    behaviour[sg] = f

        # 2. time
        validity = FieldValidity()
        def get_validity(T_varname):
            if not T_varname in self._variables.keys():
                raise epygramError('unable to find T_grid in variables.')
            T = self._variables[T_varname][:]
            time_unit = self._variables[T_varname].units
            T = netCDF4.num2date(T, time_unit)
            validity = FieldValidityList()
            validity.pop()
            basis = netCDF4.num2date(0, time_unit)
            for v in T:
                validity.append(FieldValidity(date_time=v, basis=basis))
            return validity
        if 'T_dimension' in dims_dict_e2n.keys():
            # field has a time dimension
            var_corresponding_to_T_grid = behaviour.get('T_grid', False)
            validity = get_validity(var_corresponding_to_T_grid)
            if dims_dict_e2n['T_dimension'] in only.keys():
                validity = validity[only[dims_dict_e2n['T_dimension']]]
        elif any([t in self._variables.keys() for t in config.netCDF_usualnames_for_standard_dimensions['T_dimension']]):
            # look for a time variable
            T_varnames = [t for t in config.netCDF_usualnames_for_standard_dimensions['T_dimension'] if t in self._variables.keys()]
            if len(T_varnames) == 1:
                _v = get_validity(T_varnames[0])
                if len(_v) == 1:
                    validity = _v
        field_kwargs['validity'] = validity

        # 3. GEOMETRY
        # ===========
        kwargs_geom = {}
        # 3.1 identify the structure
        keys = copy.copy(variable_dimensions.keys())
        for k in only.keys():
            if k in keys:
                keys.remove(k)
            else:
                raise ValueError("dimension: " + k + " from 'only' not in field variable.")
        if 'T_dimension' in dims_dict_e2n.keys() and dims_dict_e2n['T_dimension'] not in only.keys():
            keys.remove(dims_dict_e2n['T_dimension'])
        squeezed_variables = [dims_dict_n2e.get(k)
                              for k in keys
                              if variable_dimensions[k] != 1]
        H2D = set(squeezed_variables) == set(['X_dimension',
                                              'Y_dimension'])
        D3 = set(squeezed_variables) == set(['X_dimension',
                                             'Y_dimension',
                                             'Z_dimension'])
        V2D = set(squeezed_variables) == set(['N_dimension',
                                              'Z_dimension'])
        H1D = set(squeezed_variables) == set(['N_dimension'])
        V1D = set(squeezed_variables) == set(['Z_dimension'])
        points = set(squeezed_variables) == set([])
        if not any([D3, H2D, V2D, H1D, V1D, points]):
            raise epygramError("unable to guess structure of the field: " \
                               + str(variable_dimensions.keys())\
                               + "refine behaviour dimensions or filter dimensions with 'only'.")
        else:
            if D3:
                structure = '3D'
            elif H2D:
                structure = 'H2D'
            elif V2D:
                structure = 'V2D'
            elif H1D:
                structure = 'H1D'
                raise NotImplementedError('H1D: not yet.')  #TODO:
            elif V1D:
                structure = 'V1D'
            elif points:
                structure = 'Point'
            kwargs_geom['structure'] = structure

        # 3.2 vertical geometry (default)
        default_kwargs_vcoord = {'structure':'V',
                                 'typeoffirstfixedsurface':255,
                                 'position_on_grid':'mass',
                                 'grid':{'gridlevels': []},
                                 'levels':[0]}
        # TODO: complete with field dict when we have one
        # + fid generic ?
        kwargs_vcoord = default_kwargs_vcoord

        # 3.3 Specific parts
        # 3.3.1 dimensions
        dimensions = {}
        kwargs_geom['name'] = 'unstructured'
        if D3 or H2D:
            dimensions['X'] = variable_dimensions[dims_dict_e2n['X_dimension']]
            dimensions['Y'] = variable_dimensions[dims_dict_e2n['Y_dimension']]
        if V2D or H1D:
            dimensions['X'] = variable_dimensions[dims_dict_e2n['N_dimension']]
            dimensions['Y'] = 1
        if V1D or points:
            dimensions['X'] = 1
            dimensions['Y'] = 1
        # 3.3.2 vertical part
        if D3 or V1D or V2D:
            var_corresponding_to_Z_grid = behaviour.get('Z_grid', False)
            #assert var_corresponding_to_Z_grid in self._variables.keys(), \
            #       'unable to find Z_grid in variables.'
            levels = None
            if var_corresponding_to_Z_grid in self._variables.keys():
                if hasattr(self._variables[var_corresponding_to_Z_grid], 'standard_name') \
                   and self._variables[var_corresponding_to_Z_grid].standard_name in ('atmosphere_hybrid_sigma_pressure_coordinate',
                                                                                      'atmosphere_hybrid_height_coordinate'):
                    formula_terms = self._variables[var_corresponding_to_Z_grid].formula_terms.split(' ')
                    if 'a:' in formula_terms and 'p0:' in formula_terms:
                        a_name = formula_terms[formula_terms.index('a:') + 1]
                        p0_name = formula_terms[formula_terms.index('p0:') + 1]
                        a = self._variables[a_name][:] * self._variables[p0_name][:]
                    elif 'ap:' in formula_terms:
                        a_name = formula_terms[formula_terms.index('ap:') + 1]
                        a = self._variables[a_name][:]
                    b_name = formula_terms[formula_terms.index('b:') + 1]
                    b = self._variables[b_name][:]
                    gridlevels = [(i + 1, {'Ai':a[i], 'Bi':b[i]}) for i in range(len(a))]
                    levels = [i + 1 for i in range(len(a))]
                else:
                    gridlevels = self._variables[var_corresponding_to_Z_grid][:]
                if hasattr(self._variables[behaviour['Z_grid']], 'standard_name'):
                    kwargs_vcoord['typeoffirstfixedsurface'] = _typeoffirstfixedsurface_dict.get(self._variables[behaviour['Z_grid']].standard_name, 255)
                # TODO: complete the reading of variable units to convert
                if hasattr(self._variables[behaviour['Z_grid']], 'units'):
                    if self._variables[behaviour['Z_grid']].units == 'km':
                        gridlevels = gridlevels * 1000.  # get back to metres
            else:
                gridlevels = range(1, variable_dimensions[dims_dict_e2n['Z_dimension']] + 1)
                kwargs_vcoord['typeoffirstfixedsurface'] = 255
            kwargs_vcoord['grid']['gridlevels'] = [p for p in gridlevels]  # footprints purpose
            if levels is None:
                kwargs_vcoord['levels'] = kwargs_vcoord['grid']['gridlevels']  # could it be else ?
            else:
                kwargs_vcoord['levels'] = levels

        # 3.3.3 horizontal part
        # find grid in variables
        if H2D or D3:
            var_corresponding_to_X_grid = behaviour.get('X_grid', False)
            if not var_corresponding_to_X_grid in self._variables.keys():
                raise epygramError('unable to find X_grid in variables.')
            var_corresponding_to_Y_grid = behaviour.get('Y_grid', False)
            if not var_corresponding_to_Y_grid in self._variables.keys():
                raise epygramError('unable to find Y_grid in variables.')
            else:
                if hasattr(self._variables[var_corresponding_to_Y_grid], 'standard_name') \
                and self._variable[var_corresponding_to_Y_grid].standard_name == 'projection_y_coordinate':
                    behaviour['grid_is_lonlat'] = False
                elif 'lat' in var_corresponding_to_Y_grid.lower() \
                and 'grid_is_lonlat' not in behaviour.keys():
                    behaviour['grid_is_lonlat'] = True
            if len(self._variables[var_corresponding_to_X_grid].dimensions) == 1 \
            and len(self._variables[var_corresponding_to_Y_grid].dimensions) == 1:
                X = self._variables[var_corresponding_to_X_grid][:]
                Y = self._variables[var_corresponding_to_Y_grid][:]
                Xgrid = numpy.ones((Y.size, X.size)) * X
                Ygrid = (numpy.ones((Y.size, X.size)).transpose() * Y).transpose()
            elif len(self._variables[var_corresponding_to_X_grid].dimensions) == 2 \
            and len(self._variables[var_corresponding_to_Y_grid].dimensions) == 2:
                Xgrid = self._variables[var_corresponding_to_X_grid][:, :]
                Ygrid = self._variables[var_corresponding_to_Y_grid][:, :]
            if Ygrid[0, 0] > Ygrid[-1, 0]:
                return_Yaxis = True
                Ygrid = Ygrid[::-1, :]

            # projection or grid
            if hasattr(variable, 'grid_mapping'):
                # geometry described as "grid_mapping" meta-data
                gm = variable.grid_mapping
                grid_mapping = self._variables[gm]
                if grid_mapping.grid_mapping_name in ('lambert_conformal_conic',):
                    if (hasattr(self._variables[variable.grid_mapping], 'resolution') \
                    or not behaviour.get('grid_is_lonlat', False)):
                        # if resolution is either in grid_mapping attributes or in the grid itself
                        if grid_mapping.grid_mapping_name == 'lambert_conformal_conic':
                            kwargs_geom['name'] = 'lambert'
                            if hasattr(grid_mapping, 'resolution'):
                                Xresolution = Yresolution = grid_mapping.resolution
                            else:
                                Xresolution = abs(Xgrid[0, 0] - Xgrid[0, 1])
                                Yresolution = abs(Ygrid[0, 0] - Ygrid[1, 0])
                            grid = {'input_lon':Angle(grid_mapping.longitude_of_central_meridian, 'degrees'),
                                    'input_lat':Angle(grid_mapping.latitude_of_projection_origin, 'degrees'),
                                    'input_position':((float(dimensions['X']) - 1) / 2.,
                                                      (float(dimensions['Y']) - 1) / 2.),
                                    'X_resolution':Xresolution,
                                    'Y_resolution':Yresolution,
                                    'LAMzone':None}
                            kwargs_geom['projection'] = {'reference_lon':Angle(grid_mapping.longitude_of_central_meridian, 'degrees'),
                                                         'reference_lat':Angle(grid_mapping.standard_parallel, 'degrees'),
                                                         'rotation':Angle(0., 'radians')}
                            if hasattr(grid_mapping, 'standard_meridian'):
                                kwargs_geom['projection']['reference_lon'] = Angle(grid_mapping.standard_meridian, 'degrees')
                        else:
                            raise NotImplementedError('grid_mapping.grid_mapping_name == ' + grid_mapping.grid_mapping_name)
                    else:
                        # no resolution available: grid mapping is useless
                        gm = None
                elif 'gauss' in grid_mapping.grid_mapping_name.lower():
                    # NOTE: this is a (good) approximation actually, the true latitudes are the roots of Legendre polynoms
                    if hasattr(grid_mapping, 'latitudes'):
                        latitudes = self._variables[grid_mapping.latitudes.split(' ')[1]][:]
                    else:
                        raise NotImplementedError('(re-)computation of Gauss latitudes (not in file metadata.')
                    grid = {'latitudes':FPList([Angle(l, 'degrees') for l in latitudes]),
                            'dilatation_coef':1.}
                    if hasattr(grid_mapping, 'lon_number_by_lat'):
                        if not isinstance(grid_mapping.lon_number_by_lat, int):
                            kwargs_geom['name'] = 'regular_gauss'
                            lon_number_by_lat = self._variables[grid_mapping.lon_number_by_lat.split(' ')[1]][:]
                        else:
                            lon_number_by_lat = [dimensions['X'] for _ in range(dimensions['Y'])]
                        if hasattr(grid_mapping, 'pole_lon'):
                            kwargs_geom['name'] = 'rotated_reduced_gauss'
                            grid['pole_lon'] = Angle(grid_mapping.pole_lon, 'degrees')
                            grid['pole_lat'] = Angle(grid_mapping.pole_lat, 'degrees')
                            if hasattr(grid_mapping, 'dilatation_coef'):
                                grid['dilatation_coef'] = grid_mapping.dilatation_coef
                        else:
                            kwargs_geom['name'] = 'reduced_gauss'
                    dimensions = {'max_lon_number':int(max(lon_number_by_lat)),
                                  'lat_number':len(latitudes),
                                  'lon_number_by_lat':FPList([int(n) for n in
                                                      lon_number_by_lat])}
                    return_Yaxis = False
                else:
                    raise NotImplementedError('grid_mapping.grid_mapping_name == ' + grid_mapping.grid_mapping_name)
            else:
                # grid only in variables
                if behaviour.get('grid_is_lonlat', False):
                    grid = {'longitudes':Xgrid,
                            'latitudes':Ygrid}
                else:
                    # grid is not lon/lat and no other metadata available : Academic
                    kwargs_geom['name'] = 'academic'
                    grid = {'LAMzone':None,
                            'X_resolution':abs(Xgrid[0, 1] - Xgrid[0, 0]),
                            'Y_resolution':abs(Ygrid[1, 0] - Ygrid[0, 0])}
                    #grid = {'X':Xgrid,
                    #        'Y':Ygrid}
        elif V1D or V2D or points:
            var_corresponding_to_X_grid = behaviour.get('X_grid', False)
            if not var_corresponding_to_X_grid in self._variables.keys():
                if points or V1D:
                    lon = ['_']
                else:
                    raise epygramError('unable to find X_grid in variables.')
            else:
                lon = self._variables[var_corresponding_to_X_grid][:]
            var_corresponding_to_Y_grid = behaviour.get('Y_grid', False)
            if not var_corresponding_to_Y_grid in self._variables.keys():
                if points or V1D:
                    lat = ['_']
                else:
                    raise epygramError('unable to find Y_grid in variables.')
            else:
                lat = self._variables[var_corresponding_to_Y_grid][:]
            grid = {'longitudes':lon,
                    'latitudes':lat,
                    'LAMzone':None}

        # 3.4 build geometry
        vcoordinate = fpx.geometry(**kwargs_vcoord)
        kwargs_geom['grid'] = grid
        kwargs_geom['dimensions'] = dimensions
        kwargs_geom['vcoordinate'] = vcoordinate
        kwargs_geom['position_on_horizontal_grid'] = 'center'
        geometry = fpx.geometry(**kwargs_geom)

        # 4. build field
        field_kwargs['geometry'] = geometry
        field_kwargs['structure'] = kwargs_geom['structure']
        comment = {}
        for a in variable.ncattrs():
            if a != 'validity':
                if isinstance(variable.getncattr(a), numpy.float32):  # pb with json and float32
                    comment.update({a:numpy.float64(variable.getncattr(a))})
                else:
                    comment.update({a:variable.getncattr(a)})
        comment = json.dumps(comment)
        if comment != '{}':
            field_kwargs['comment'] = comment
        field = fpx.field(**field_kwargs)
        if getdata:
            if only:
                n = len(variable.dimensions)
                buffdata = variable
                for k, i in only.items():
                    d = variable.dimensions.index(k)
                    buffdata = util.restrain_to_index_i_of_dim_d(buffdata, i, d, n=n)
            else:
                buffdata = variable[...]
            # check there is no leftover unknown dimension
            field_dim_num = 1 if len(field.validity) > 1 else 0
            if field.structure != 'Point':
                field_dim_num += [int(c) for c in field.structure if c.isdigit()][0]
            assert field_dim_num == len(buffdata.squeeze().shape), \
                   ' '.join(['shape of field and identified usual dimensions',
                             'do not match: use *only* to filter or',
                             '*adhoc_behaviour* to identify dimensions'])
            if len(buffdata.shape) > 1:
                # re-shuffle to have data indexes in order (t,z,y,x)
                positions = []
                shp4D = [1, 1, 1, 1]
                if 'T_dimension' in dims_dict_e2n.keys():
                    idx = variable.dimensions.index(dims_dict_e2n['T_dimension'])
                    positions.append(idx)
                    shp4D[0] = buffdata.shape[idx]
                if  'Z_dimension' in dims_dict_e2n.keys():
                    idx = variable.dimensions.index(dims_dict_e2n['Z_dimension'])
                    positions.append(idx)
                    shp4D[1] = buffdata.shape[idx]
                if  'Y_dimension' in dims_dict_e2n.keys():
                    idx = variable.dimensions.index(dims_dict_e2n['Y_dimension'])
                    positions.append(idx)
                    shp4D[2] = buffdata.shape[idx]
                if  'X_dimension' in dims_dict_e2n.keys():
                    idx = variable.dimensions.index(dims_dict_e2n['X_dimension'])
                    positions.append(idx)
                    shp4D[3] = buffdata.shape[idx]
                elif  'N_dimension' in dims_dict_e2n.keys():
                    idx = variable.dimensions.index(dims_dict_e2n['N_dimension'])
                    positions.append(idx)
                    shp4D[3] = buffdata.shape[idx]
                for d in variable.dimensions:
                    # whatever the order of these, they must have been filtered and dimension 1 (only)
                    if d not in dims_dict_e2n.values():
                        positions.append(variable.dimensions.index(d))
                shp4D = tuple(shp4D)
                buffdata = buffdata.transpose(*positions).squeeze()
                if isinstance(buffdata, numpy.ma.masked_array):
                    data = numpy.ma.zeros(shp4D)
                else:
                    data = numpy.empty(shp4D)
                data[...] = buffdata.reshape(data.shape)
                if return_Yaxis:
                    data[...] = data[:, :, ::-1, :]
            else:
                data = buffdata
            field.setdata(data)

        return field

    def writefield(self, field, compression=4, metadata=None):
        """
        Write a field in resource.
        Args:\n
        - *compression* ranges from 1 (low compression, fast writing)
          to 9 (high compression, slow writing). 0 is no compression.
        - *metadata*: dict, can be filled by any meta-data, that will be stored
          as attribute of the netCDF variable.
        """

        metadata = util.ifNone_emptydict(metadata)
        vartype = 'f8'
        fill_value = -999999.9
        def check_or_add_dim(d, d_in_field=None, size=None):
            if size is None:
                if d_in_field is None:
                    d_in_field = d
                size = field.geometry.dimensions[d_in_field]
            if d not in self._dimensions:
                self._nc.createDimension(d, size=size)
            else:
                assert len(self._dimensions[d]) == size, \
                       "dimensions mismatch: " + d + ": " + \
                       str(self._dimensions[d]) + " != " + str(size)
        def check_or_add_variable(varname, vartype,
                                  dimensions=(),
                                  **kwargs):
            if unicode(varname) not in self._variables.keys():
                var = self._nc.createVariable(varname, vartype,
                                              dimensions=dimensions,
                                              **kwargs)
                status = 'created'
            else:
                assert self._variables[varname].dtype == vartype, \
                       ' '.join(['variable', varname,
                                 'already exist with other type:',
                                 self._variables[varname].dtype])
                if isinstance(dimensions, str):
                    dimensions = (dimensions,)
                assert self._variables[varname].dimensions == tuple(dimensions), \
                       ' '.join(['variable', varname,
                                 'already exist with other dimensions:',
                                 str(self._variables[varname].dimensions)])
                var = self._variables[varname]
                status = 'match'
            return var, status

        assert field.fid.has_key('netCDF')
        assert not field.spectral

        # 1. dimensions
        T = Y = X = G = N = None
        dims = []
        # time
        if len(field.validity) > 1:
            # default
            T = config.netCDF_usualnames_for_standard_dimensions['T_dimension'][0]
            # or any existing identified time dimension
            T = {'found':v for v in self._dimensions
                 if (v in config.netCDF_usualnames_for_standard_dimensions['T_dimension']
                     and len(self._dimensions[v]) == len(field.validity))}.get('found', T)
            # or specified behaviour
            T = self.behaviour.get('T_dimension', T)
            check_or_add_dim(T, size=len(field.validity))
        # vertical part
        # default
        Z = config.netCDF_usualnames_for_standard_dimensions['Z_dimension'][0]
        # or any existing identified time dimension
        Z = {'found':v for v in self._dimensions
                 if (v in config.netCDF_usualnames_for_standard_dimensions['Z_dimension']
                     and len(self._dimensions[v]) == len(field.geometry.vcoordinate.levels))}.get('found', Z)
        # or specified behaviour
        Z = self.behaviour.get('Z_dimension', Z)
        if 'gridlevels' in field.geometry.vcoordinate.grid.keys():
            Z_gridsize = max(len(field.geometry.vcoordinate.grid['gridlevels']), 1)
            if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
                Z_gridsize -= 1
        else:
            Z_gridsize = 1
        if Z_gridsize > 1:
            check_or_add_dim(Z, size=Z_gridsize)
        # horizontal
        if field.geometry.rectangular_grid:
            if field.geometry.dimensions['Y'] > 1 and field.geometry.dimensions['X'] > 1:
                Y = self.behaviour.get('Y_dimension',
                                       config.netCDF_usualnames_for_standard_dimensions['Y_dimension'][0])
                check_or_add_dim(Y, d_in_field='Y')
                X = self.behaviour.get('X_dimension',
                                       config.netCDF_usualnames_for_standard_dimensions['X_dimension'][0])
                check_or_add_dim(X, d_in_field='X')
            elif field.geometry.dimensions['X'] > 1:
                N = self.behaviour.get('N_dimension',
                                       config.netCDF_usualnames_for_standard_dimensions['N_dimension'][0])
                check_or_add_dim(N, d_in_field='X')
        elif 'gauss' in field.geometry.name:
            Y = self.behaviour.get('Y_dimension', 'latitude')
            check_or_add_dim(Y, d_in_field='lat_number')
            X = self.behaviour.get('X_dimension', 'longitude')
            check_or_add_dim(X, d_in_field='max_lon_number')
            if self.behaviour['flatten_non_rectangular_grids']:
                _gpn = sum(field.geometry.dimensions['lon_number_by_lat'])
                G = 'gridpoints_number'
                check_or_add_dim(G, size=_gpn)
        else:
            raise NotImplementedError("grid not rectangular nor a gauss one.")

        # 2. validity
        #TODO: deal with unlimited time dimension ?
        if field.validity[0] != FieldValidity():
            tgrid = config.netCDF_usualnames_for_standard_dimensions['T_dimension'][0]
            tgrid = {'found':v for v in self._variables
                     if v in config.netCDF_usualnames_for_standard_dimensions['T_dimension']}.get('found', tgrid)
            tgrid = self.behaviour.get('T_grid', tgrid)
            if len(field.validity) == 1:
                _, _status = check_or_add_variable(tgrid, int)
            else:
                _, _status = check_or_add_variable(tgrid, int, T)
                dims.append(tgrid)
            datetime0 = field.validity[0].getbasis().isoformat(sep=' ')
            datetimes = [int((dt.get() - field.validity[0].getbasis()).total_seconds()) for dt in field.validity]
            if _status == 'created':
                self._variables[tgrid][:] = datetimes
                self._variables[tgrid].units = ' '.join(['seconds', 'since', datetime0])
            else:
                assert (self._variables[tgrid][:] == datetimes).all(), \
                       ' '.join(['variable', tgrid, 'mismatch.'])

        # 3. geometry
        # 3.1 vertical part
        if len(field.geometry.vcoordinate.levels) > 1:
            dims.append(Z)
        if Z_gridsize > 1:
            zgridname = config.netCDF_usualnames_for_standard_dimensions['Z_dimension'][0]
            zgridname = {'found':v for v in self._variables
                         if v in config.netCDF_usualnames_for_standard_dimensions['Z_dimension']}.get('found', zgridname)
            zgridname = self.behaviour.get('Z_grid', zgridname)
            if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
                ZP1 = Z + '+1'
                check_or_add_dim(ZP1, size=Z_gridsize + 1)
                zgrid, _status = check_or_add_variable(zgridname, int)
                if _status == 'created':
                    if field.geometry.vcoordinate.typeoffirstfixedsurface == 119:
                        zgrid.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
                        zgrid.positive = "down"
                        zgrid.formula_terms = "ap: hybrid_coef_A b: hybrid_coef_B ps: surface_air_pressure"
                        check_or_add_variable('hybrid_coef_A', vartype, ZP1)
                        self._variables['hybrid_coef_A'][:] = [iab[1]['Ai'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                        check_or_add_variable('hybrid_coef_B', vartype, ZP1)
                        self._variables['hybrid_coef_B'][:] = [iab[1]['Bi'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                    elif field.geometry.vcoordinate.typeoffirstfixedsurface == 118:
                        # TOBECHECKED:
                        zgrid.standard_name = "atmosphere_hybrid_height_coordinate"
                        zgrid.positive = "up"
                        zgrid.formula_terms = "a: hybrid_coef_A b: hybrid_coef_B orog: orography"
                        check_or_add_variable('hybrid_coef_A', vartype, ZP1)
                        self._variables['hybrid_coef_A'][:] = [iab[1]['Ai'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                        check_or_add_variable('hybrid_coef_B', vartype, ZP1)
                        self._variables['hybrid_coef_B'][:] = [iab[1]['Bi'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                else:
                    epylog.info('assume 118/119 type vertical grid matches.')
            else:
                if len(numpy.shape(field.geometry.vcoordinate.grid['gridlevels'])) > 1:
                    dims_Z = [d for d in [Z, Y, X, G, N] if d is not None]
                else:
                    dims_Z = Z
                zgrid, _status = check_or_add_variable(zgridname, vartype, dims_Z)
                u = {102:'m', 103:'m', 100:'hPa'}.get(field.geometry.vcoordinate.typeoffirstfixedsurface, None)
                if u is not None:
                    zgrid.units = u
                if _status == 'created':
                    zgrid[:] = field.geometry.vcoordinate.grid['gridlevels']
                else:
                    assert zgrid[:].all() == numpy.array(field.geometry.vcoordinate.grid['gridlevels']).all(), \
                           ' '.join(['variable', zgrid, 'mismatch.'])
            if _typeoffirstfixedsurface_dict_inv.get(field.geometry.vcoordinate.typeoffirstfixedsurface, False):
                zgrid.short_name = _typeoffirstfixedsurface_dict_inv[field.geometry.vcoordinate.typeoffirstfixedsurface]
        # 3.2 grid (lonlat)
        dims_lonlat = []
        (lons, lats) = field.geometry.get_lonlat_grid()
        if self.behaviour['flatten_non_rectangular_grids'] \
        and not field.geometry.rectangular_grid:
            dims_lonlat.append(G)
            dims.append(G)
            lons = stretch_array(lons)
            lats = stretch_array(lats)
        elif 'gauss' in field.geometry.name or field.geometry.dimensions.get('Y', 0) > 1:  # both Y and X dimensions
            dims_lonlat.extend([Y, X])
            dims.extend([Y, X])
        elif field.geometry.dimensions['X'] > 1:  # only X == N
            dims_lonlat.append(N)
            dims.append(N)
        # else: pass (single point or profile)
        if isinstance(lons, numpy.ma.masked_array):
            lons = lons.filled(fill_value)
            lats = lats.filled(fill_value)
        else:
            fill_value = None
        try:
            _ = float(stretch_array(lons)[0])
        except ValueError:
            no_lonlat = True
        else:
            no_lonlat = False
        if not no_lonlat:
            lons_var, _status = check_or_add_variable('longitude', vartype, dims_lonlat, fill_value=fill_value)
            lats_var, _status = check_or_add_variable('latitude', vartype, dims_lonlat, fill_value=fill_value)
            if _status == 'match':
                epylog.info('assume lons/lats match.')
            else:
                lons_var[...] = lons[...]
                lats_var[...] = lats[...]
        # 3.3 meta-data
        if field.geometry.dimensions.get('Y', field.geometry.dimensions.get('lat_number', 0)) > 1:
            if all([k in field.geometry.name for k in ('reduced', 'gauss')]):
                # reduced Gauss case
                meta = 'Gauss_grid'
                _, _status = check_or_add_variable(meta, int)
                if _status == 'created':
                    self._variables[meta].grid_mapping_name = "gauss_grid"
                    self._variables[meta].lon_number_by_lat = 'var: lon_number_by_lat'
                    check_or_add_variable('lon_number_by_lat', int, Y)
                    self._variables['lon_number_by_lat'][:] = field.geometry.dimensions['lon_number_by_lat']
                    self._variables[meta].latitudes = 'var: gauss_latitudes'
                    check_or_add_variable('gauss_latitudes', float, Y)
                    self._variables['gauss_latitudes'][:] = [l.get('degrees') for l in field.geometry.grid['latitudes']]
                    if 'pole_lon' in field.geometry.grid.keys():
                        self._variables[meta].pole_lon = field.geometry.grid['pole_lon'].get('degrees')
                        self._variables[meta].pole_lat = field.geometry.grid['pole_lat'].get('degrees')
                    if 'dilatation_coef' in field.geometry.grid.keys():
                        self._variables[meta].dilatation_coef = field.geometry.grid['dilatation_coef']
                else:
                    epylog.info('assume Gauss grid parameters match.')
            elif field.geometry.projected_geometry:
                # projections
                if field.geometry.name == 'lambert':
                    meta = 'Lambert_Conformal'
                    _, _status = check_or_add_variable(meta, int)
                    if _status == 'created':
                        self._variables[meta].grid_mapping_name = 'lambert_conformal_conic'
                        if field.geometry.grid['X_resolution'] != field.geometry.grid['Y_resolution']:
                            raise NotImplementedError('anisotropic resolution')
                        else:
                            self._variables[meta].resolution = field.geometry.grid['X_resolution']
                        if field.geometry.grid.get('LAMzone'):
                            _lon_cen, _lat_cen = field.geometry.ij2ll(float(field.geometry.dimensions['X'] - 1.) / 2.,
                                                                      float(field.geometry.dimensions['Y'] - 1.) / 2.)
                        else:
                            _lon_cen = field.geometry._center_lon.get('degrees')
                            _lat_cen = field.geometry._center_lat.get('degrees')
                        self._variables[meta].longitude_of_central_meridian = _lon_cen
                        self._variables[meta].latitude_of_projection_origin = _lat_cen
                        if not nearlyEqual(field.geometry._center_lon.get('degrees'),
                                           field.geometry.projection['reference_lon'].get('degrees')):
                            epylog.warning('center_lon != reference_lon (tilting) is not "on the cards" in CF convention 1.6')
                            self._variables[meta].standard_meridian = field.geometry.projection['reference_lon'].get('degrees')
                        if field.geometry.secant_projection:
                            self._variables[meta].standard_parallel = [field.geometry.projection['secant_lat1'].get('degrees'),
                                                                                      field.geometry.projection['secant_lat2'].get('degrees')]
                        else:
                            self._variables[meta].standard_parallel = field.geometry.projection['reference_lat'].get('degrees')
                    else:
                        epylog.info('assume projection parameters match.')
                else:
                    raise NotImplementedError('field.geometry.name == ' + field.geometry.name)
            else:
                meta = False
        else:
            meta = False

        # 4. Variable
        varname = field.fid['netCDF'].replace('.', config.netCDF_replace_dot_in_variable_names)
        _, _status = check_or_add_variable(varname, vartype, dims,
                                           zlib=bool(compression),
                                           complevel=compression,
                                           fill_value=fill_value)
        if meta:
            self._variables[varname].grid_mapping = meta
        if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
            self._variables[varname].vertical_grid = zgridname
        if self.behaviour['flatten_non_rectangular_grids'] \
        and not field.geometry.rectangular_grid:
            data = stretch_array(field.getdata())
        else:
            data = field.getdata()
        if isinstance(data, numpy.ma.masked_array):
            data = data.filled(fill_value)
        if _status == 'match':
            epylog.info('overwrite data in variable ' + varname)
        self._variables[varname][...] = data

        # 5. metadata
        for k, v in metadata.items():
            self._nc.setncattr(k, v)

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
                # print_ncattr(dim)
            # Variable information.
            nc_vars = [var for var in nc.variables]  # list of nc variables
            outwrite("NetCDF variable information:")
            for var in nc_vars:
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

