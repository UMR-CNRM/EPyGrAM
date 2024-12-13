#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains classes for netCDF4 resource.
"""

import copy
import json
import sys
import numpy
from collections import OrderedDict
import datetime
import re

import netCDF4

import footprints
from footprints import proxy as fpx, FPDict, FPList
from bronx.syntax.arrays import stretch_array

from epygram import config, epygramError, util
from epygram.base import FieldValidity, FieldValidityList
from epygram.resources import FileResource
from epygram.util import Angle
from epygram.geometries import (VGeometry, ProjectedGeometry, GaussGeometry, 
                                RegLLGeometry, AcademicGeometry, UnstructuredGeometry)

__all__ = ['netCDF']

epylog = footprints.loggers.getLogger(__name__)


_typeoffirstfixedsurface_dict = {'altitude':102,
                                 'height':103,
                                 'atmosphere_hybrid_sigma_pressure_coordinate':119,
                                 'atmosphere_hybrid_height_coordinate':118,
                                 'pressure':100}
_typeoffirstfixedsurface_short_dict = _typeoffirstfixedsurface_dict.copy()
_typeoffirstfixedsurface_short_dict.update({'hybrid-pressure':119,
                                            'hybrid-height':118})
_typeoffirstfixedsurface_dict_inv = {v:k for k, v in _typeoffirstfixedsurface_dict.items()}
_typeoffirstfixedsurface_short_dict_inv = {v:k for k, v in _typeoffirstfixedsurface_short_dict.items()}

_proj_dict = {'lambert':'lambert_conformal_conic',
              'mercator':'mercator',
              'polar_stereographic':'polar_stereographic',
              'space_view':'vertical_perspective'}
_proj_dict_inv = {v:k for k, v in _proj_dict.items()}


def _default_numpy2json(obj):
    """Helper to encode numpy arrays to json."""
    if type(obj).__module__ == numpy.__name__:
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif any([isinstance(obj, t) for t in [numpy.int8, numpy.int16, numpy.int32, numpy.int64]]):
            return int(obj)
        elif any([isinstance(obj, t) for t in [numpy.float32, numpy.float64]]):
            return float(obj)
        else:
            return obj.item()
    raise TypeError('Unknown type:', type(obj))


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
        self.isopen = False
        super(netCDF, self).__init__(*args, **kwargs)
        if self.openmode in ('r', 'a'):
            try:
                guess = netCDF4.Dataset(self.container.abspath, self.openmode)
            except (RuntimeError, UnicodeEncodeError):
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

        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(netCDF, self).open(openmode=openmode)
        self._nc = netCDF4.Dataset(self.container.abspath, self.openmode)
        self.isopen = True
        if self.openmode == 'w':
            self.set_default_global_attributes()

    def close(self):
        """Closes a netCDF."""
        if hasattr(self, '_nc') and self._nc._isopen:
            self._nc.close()
        self.isopen = False

    def variables_number(self):
        """Return the number of variables in resource."""
        return len(self._variables)

    def find_fields_in_resource(self, seed=None, generic=False, **_):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        :param seed: might be a regular expression, a list of regular expressions
          or *None*. If *None* (default), returns the list of all fields in
          resource.
        :param fieldtype: optional, among ('H2D', 'Misc') or a list of these strings.
          If provided, filters out the fields not of the given types.
        :param generic: if True, returns complete fid's,
          union of {'FORMATname':fieldname} and the according generic fid of
          the fields.
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
        if generic:
            raise NotImplementedError("not yet.")

        return fieldslist

    @FileResource._openbeforedelayed
    def ncinfo_field(self, fid):
        """
        Get info about the field (dimensions and meta-data of the netCDF variable).

        :param fid: netCDF field identifier
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

        :param fid: netCDF field identifier
        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        :param only: to specify indexes [0 ... n-1] of specific dimensions,
          e.g. {'time':5,} to select only the 6th term of time dimension.
        :param adhoc_behaviour: to specify "on the fly" a behaviour (usual
          dimensions or grids, ...).
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

        # 1.1 identify usual dimensions
        all_dimensions_e2n = {}
        for d in self._dimensions.keys():
            for sd in config.netCDF_standard_dimensions:
                if d in config.netCDF_usualnames_for_standard_dimensions[sd]:
                    all_dimensions_e2n[sd] = d
        variable_dimensions = {d:len(self._dimensions[d]) for d in variable.dimensions}
        for d in variable_dimensions.keys():
            for sd in config.netCDF_standard_dimensions:
                # if behaviour is not explicitly given,
                # try to find out who is "d" among the standard dimensions
                if sd not in behaviour and d in config.netCDF_usualnames_for_standard_dimensions[sd]:
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
                if sg not in behaviour and f in config.netCDF_usualnames_for_standard_dimensions[sd] or \
                   f == behaviour.get(sd):
                    behaviour[sg] = f

        # 2. time
        def get_validity(T_varname):
            try:
                validity = FieldValidityList()
                validity.pop()
                if T_varname not in self._listfields():
                    raise epygramError('unable to find T_grid in variables.')
                T = self._variables[T_varname][:]
                if len(self._variables[T_varname].dimensions) == 0:
                    T = [T]
                time_unit = getattr(self._variables[T_varname], 'units', '')
                if re.match(r'(hours|seconds|days|minutes)\s+since.+$', time_unit):
                    T = netCDF4.num2date(T, time_unit).squeeze().reshape((len(T),))
                    T = [datetime.datetime(*t.timetuple()[:6]) for t in T]  # FIXME: not sure of that for dates older than julian/gregorian calendar
                    basis = netCDF4.num2date(0, time_unit)
                    if basis.year <= 1582:
                        epylog.warning('suspicion of inconsistency of julian/gregorian dates')
                    basis = datetime.datetime(*basis.timetuple()[:6])  # FIXME: not sure of that for dates older than julian/gregorian calendar
                    for v in T:
                        validity.append(FieldValidity(date_time=v, basis=basis))
                else:
                    epylog.warning('temporal unit is not CF1.6-compliant, cannot decode.')
                    for v in T:
                        validity.append(FieldValidity())
            except Exception:
                epylog.warning("Failed to read validity, but ignore !")
                validity = FieldValidityList()
            return validity
        if 'T_dimension' in dims_dict_e2n:
            # field has a time dimension
            var_corresponding_to_T_grid = behaviour.get('T_grid', False)
            validity = get_validity(var_corresponding_to_T_grid)
            if dims_dict_e2n['T_dimension'] in only:
                validity = validity[only[dims_dict_e2n['T_dimension']]]
        elif any([t in self._listfields()
                  for t in config.netCDF_usualnames_for_standard_dimensions['T_dimension']]):
            # look for a time variable
            T_varnames = [t for t in config.netCDF_usualnames_for_standard_dimensions['T_dimension'] if t in self._listfields()]
            _v = get_validity(T_varnames[0])
            if len(_v) == 1:
                validity = _v
            else:
                validity = FieldValidity()
        else:
            validity = FieldValidity()
        field_kwargs['validity'] = validity

        # 3. GEOMETRY
        # ===========
        kwargs_geom = {}
        kwargs_geom['position_on_horizontal_grid'] = 'center'
        # 3.1 identify the structure
        keys = copy.copy(list(variable_dimensions.keys()))
        for k in only.keys():
            if k in keys:
                keys.remove(k)
            else:
                raise ValueError("dimension: " + k + " from 'only' not in field variable.")
        if 'T_dimension' in dims_dict_e2n and dims_dict_e2n['T_dimension'] not in only:
            keys.remove(dims_dict_e2n['T_dimension'])
        squeezed_variables = [dims_dict_n2e.get(k)
                              for k in keys
                              if variable_dimensions[k] != 1]
        H2D = set(squeezed_variables) == set(['X_dimension',
                                              'Y_dimension']) \
              or (set(squeezed_variables) == set(['N_dimension']) and
                  all([d in all_dimensions_e2n for d in ['X_dimension',  # flattened grids
                                                         'Y_dimension']])) \
              or (set(squeezed_variables) == set(['N_dimension']) \
                  and behaviour.get('H1D_is_H2D_unstructured', False))  # or 2D unstructured grids
        D3 = set(squeezed_variables) == set(['X_dimension',
                                             'Y_dimension',
                                             'Z_dimension']) \
             or (set(squeezed_variables) == set(['N_dimension',
                                                 'Z_dimension']) and
                 all([d in all_dimensions_e2n for d in ['X_dimension',  # flattened grids
                                                        'Y_dimension']])) \
             or (set(squeezed_variables) == set(['N_dimension',
                                                 'Z_dimension']) \
                 and behaviour.get('H1D_is_H2D_unstructured', False))  # or 2D unstructured grids
        V2D = set(squeezed_variables) == set(['N_dimension',
                                              'Z_dimension']) and not D3
        H1D = set(squeezed_variables) == set(['N_dimension']) and not H2D
        V1D = set(squeezed_variables) == set(['Z_dimension'])
        points = set(squeezed_variables) == set([])
        if not any([D3, H2D, V2D, H1D, V1D, points]):
            for d in set(variable_dimensions.keys()) - set(dims_dict_n2e.keys()):
                epylog.error(" ".join(["dimension:",
                                       d,
                                       "has not been identified as a usual",
                                       "(T,Z,Y,X) dimension. Please specify",
                                       "with readfield() argument",
                                       "adhoc_behaviour={'T_dimension':'" + d + "'}",
                                       "for instance or",
                                       "my_resource.behave(T_dimension=" + d + ")",
                                       "or complete",
                                       "config.netCDF_usualnames_for_standard_dimensions",
                                       "in $HOME/.epygram/userconfig.py"]))
            raise epygramError(" ".join(["unable to guess structure of field:",
                                         str(list(variable_dimensions.keys())),
                                         "=> refine behaviour dimensions or",
                                         "filter dimensions with 'only'."]))
        elif H1D:
            raise NotImplementedError('H1D: not yet.')  # TODO:

        # 3.2 vertical geometry (default)
        default_kwargs_vcoord = {'typeoffirstfixedsurface':255,
                                 'position_on_grid':'mass',
                                 'grid':{'gridlevels': []},
                                 'levels':[0]}
        kwargs_vcoord = default_kwargs_vcoord

        # 3.3 Specific parts
        # 3.3.1 dimensions
        dimensions = {}
        kwargs_geom['name'] = 'unstructured'
        if V2D or H1D:
            dimensions['X'] = variable_dimensions[dims_dict_e2n['N_dimension']]
            dimensions['Y'] = 1
        if V1D or points:
            dimensions['X'] = 1
            dimensions['Y'] = 1
        if D3 or H2D:
            if set(['X_dimension', 'Y_dimension']).issubset(set(dims_dict_e2n.keys())):
                flattened = False
            elif 'N_dimension' in dims_dict_e2n:
                flattened = True
            else:
                raise epygramError('unable to find grid dimensions.')
            if not flattened:
                dimensions['X'] = variable_dimensions[dims_dict_e2n['X_dimension']]
                dimensions['Y'] = variable_dimensions[dims_dict_e2n['Y_dimension']]
            else:  # flattened
                if 'X_dimension' in all_dimensions_e2n \
                   and 'Y_dimension' in all_dimensions_e2n:
                    dimensions['X'] = len(self._dimensions[all_dimensions_e2n['X_dimension']])
                    dimensions['Y'] = len(self._dimensions[all_dimensions_e2n['Y_dimension']])
                elif behaviour.get('H1D_is_H2D_unstructured', False):
                    dimensions['X'] = variable_dimensions[dims_dict_e2n['N_dimension']]
                    dimensions['Y'] = 1
                else:
                    assert 'X_dimension' in all_dimensions_e2n, \
                           ' '.join(["unable to find X_dimension of field:",
                                     "please specify with readfield() argument",
                                     "adhoc_behaviour={'X_dimension':'" + d + "'}",
                                     "for instance or",
                                     "my_resource.behave(X_dimension=" + d + ")",
                                     "or complete",
                                     "config.netCDF_usualnames_for_standard_dimensions",
                                     "in $HOME/.epygram/userconfig.py"])
                    assert 'Y_dimension' in all_dimensions_e2n, \
                           ' '.join(["unable to find Y_dimension of field:",
                                     "please specify with readfield() argument",
                                     "adhoc_behaviour={'Y_dimension':'" + d + "'}",
                                     "for instance or",
                                     "my_resource.behave(Y_dimension=" + d + ")",
                                     "or complete",
                                     "config.netCDF_usualnames_for_standard_dimensions",
                                     "in $HOME/.epygram/userconfig.py"])
        # 3.3.2 vertical part
        if D3 or V1D or V2D:
            var_corresponding_to_Z_grid = behaviour.get('Z_grid', dims_dict_e2n['Z_dimension'])
            # assert var_corresponding_to_Z_grid in self._variables.keys(), \
            #        'unable to find Z_grid in variables.'
            levels = None
            if var_corresponding_to_Z_grid in self._listfields():
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
                    gridlevels = tuple([(i + 1, FPDict({'Ai':a[i], 'Bi':b[i]})) for i in range(len(a))])
                    levels = [i + 1 for i in range(variable_dimensions[dims_dict_e2n['Z_dimension']])]
                    if len(gridlevels) == len(levels):
                        kwargs_vcoord['grid']['ABgrid_position'] = 'mass'
                    else:
                        kwargs_vcoord['grid']['ABgrid_position'] = 'flux'
                else:
                    gridlevels = self._variables[var_corresponding_to_Z_grid][:]
                if hasattr(self._variables[behaviour['Z_grid']], 'standard_name'):
                    kwargs_vcoord['typeoffirstfixedsurface'] = _typeoffirstfixedsurface_dict.get(self._variables[behaviour['Z_grid']].standard_name, 255)
                elif hasattr(self._variables[behaviour['Z_grid']], 'short_name'):
                    kwargs_vcoord['typeoffirstfixedsurface'] = _typeoffirstfixedsurface_short_dict.get(self._variables[behaviour['Z_grid']].short_name, 255)
                # TODO: complete the reading of variable units to convert
                if hasattr(self._variables[behaviour['Z_grid']], 'units'):
                    if self._variables[behaviour['Z_grid']].units == 'km':
                        gridlevels = gridlevels * 1000.  # get back to metres
            else:
                gridlevels = list(range(1, variable_dimensions[dims_dict_e2n['Z_dimension']] + 1))
                kwargs_vcoord['typeoffirstfixedsurface'] = 255
            kwargs_vcoord['grid']['gridlevels'] = [p for p in gridlevels]  # footprints purpose
            if levels is None:
                kwargs_vcoord['levels'] = kwargs_vcoord['grid']['gridlevels']  # could it be else ?
            else:
                kwargs_vcoord['levels'] = levels

        # 3.3.3 horizontal part
        # find grid in variables
        geometryclass = UnstructuredGeometry #Default geometry
        if H2D or D3:
            def find_grid_in_variables():
                var_corresponding_to_X_grid = behaviour.get('X_grid', False)
                if var_corresponding_to_X_grid not in self._listfields():
                    epylog.error(" ".join(["unable to find X_grid in variables.",
                                           "Please specify with readfield()",
                                           "argument",
                                           "adhoc_behaviour={'X_grid':'name_of_the_variable'}",
                                           "for instance or",
                                           "my_resource.behave(X_grid='name_of_the_variable')"]))
                    raise epygramError('unable to find X_grid in variables.')
                var_corresponding_to_Y_grid = behaviour.get('Y_grid', False)
                if var_corresponding_to_Y_grid not in self._listfields():
                    epylog.error(" ".join(["unable to find Y_grid in variables.",
                                           "Please specify with readfield()",
                                           "argument",
                                           "adhoc_behaviour={'Y_grid':'name_of_the_variable'}",
                                           "for instance or",
                                           "my_resource.behave(Y_grid='name_of_the_variable')"]))
                    raise epygramError('unable to find Y_grid in variables.')
                else:
                    if hasattr(self._variables[var_corresponding_to_Y_grid], 'standard_name') and \
                       self._variables[var_corresponding_to_Y_grid].standard_name == 'projection_y_coordinate' and \
                       self._variables[var_corresponding_to_X_grid].standard_name == 'projection_x_coordinate':
                        behaviour['grid_is_lonlat'] = False
                    elif 'lat' in var_corresponding_to_Y_grid.lower() and \
                         'lon' in var_corresponding_to_X_grid.lower() and \
                         'grid_is_lonlat' not in behaviour:
                        behaviour['grid_is_lonlat'] = True
                if len(self._variables[var_corresponding_to_X_grid].dimensions) == 1 and \
                   len(self._variables[var_corresponding_to_Y_grid].dimensions) == 1:
                    # case of a flat grid
                    X = self._variables[var_corresponding_to_X_grid][:]
                    Y = self._variables[var_corresponding_to_Y_grid][:]
                    if flattened:
                        if len(X) == dimensions.get('X') * dimensions.get('Y'):
                            # case of a H2D field with flattened grid
                            Xgrid = X.reshape((dimensions['Y'], dimensions['X']))
                            Ygrid = Y.reshape((dimensions['Y'], dimensions['X']))
                        elif behaviour.get('H1D_is_H2D_unstructured', False):
                            # case of a H2D unstructured field
                            X = self._variables[var_corresponding_to_X_grid][:]
                            Y = self._variables[var_corresponding_to_Y_grid][:]
                            Xgrid = X.reshape((1, len(X)))
                            Ygrid = Y.reshape((1, len(Y)))
                        else:
                            raise epygramError('unable to reconstruct 2D grid.')
                    else:
                        # case of a regular grid where X is constant on a column
                        # and Y constant on a row: reconstruct 2D
                        Xgrid = numpy.ones((Y.size, X.size)) * X
                        Ygrid = (numpy.ones((Y.size, X.size)).transpose() * Y).transpose()
                elif len(self._variables[var_corresponding_to_X_grid].dimensions) == 2 and \
                     len(self._variables[var_corresponding_to_Y_grid].dimensions) == 2:
                    Xgrid = self._variables[var_corresponding_to_X_grid][:, :]
                    Ygrid = self._variables[var_corresponding_to_Y_grid][:, :]
                elif len(self._variables[var_corresponding_to_X_grid].dimensions) == 3 and \
                     len(self._variables[var_corresponding_to_Y_grid].dimensions) == 3:
                    # In this case, we check that X and Y are constant on Z axis
                    Xgrid = self._variables[var_corresponding_to_X_grid][:, :, :]
                    Ygrid = self._variables[var_corresponding_to_Y_grid][:, :, :]
                    if not all([numpy.all(Xgrid[0] == Xgrid[i] for i in range(len(Xgrid)))]):
                        raise epygramError('X coordinate must be constant on the vertical')
                    if not all([numpy.all(Ygrid[0] == Ygrid[i] for i in range(len(Ygrid)))]):
                        raise epygramError('Y coordinate must be constant on the vertical')
                    Xgrid = Xgrid[0]
                    Ygrid = Ygrid[0]
                else:
                    raise epygramError('Unknown case for X and Y dimensions')
                if Ygrid[0, 0] > Ygrid[-1, 0] and not behaviour.get('reverse_Yaxis'):
                    epylog.warning("Ygrid seems to be reversed; shouldn't behaviour['reverse_Yaxis'] be True ?")
                elif behaviour.get('reverse_Yaxis'):
                    Ygrid = Ygrid[::-1, :]
                    Xgrid = Xgrid[::-1, :]

                return Xgrid, Ygrid

            # projection or grid
            def grid_mapping_ok():
                ok = False
                if hasattr(variable, 'grid_mapping'):
                    if variable.grid_mapping in self._variables:
                        if hasattr(self._variables[variable.grid_mapping], 'grid_mapping_name'):
                            gmm = self._variables[variable.grid_mapping].grid_mapping_name
                            if gmm in ('lambert_conformal_conic',
                                       'mercator',
                                       'polar_stereographic',
                                       'latitude_longitude') or 'gauss' in gmm.lower():
                                ok = True
                    elif variable.grid_mapping == 'GeosCoordinateSystem':
                        if 'geos' in self._variables:
                            if hasattr(self._variables['geos'], 'grid_mapping_name'):
                                gmm = self._variables['geos'].grid_mapping_name
                                if gmm == 'vertical_perspective':
                                    ok = True
                return ok

            if grid_mapping_ok():
                # geometry described as "grid_mapping" meta-data
                if (variable.grid_mapping == 'GeosCoordinateSystem' and
                    'geos' in self._variables
                    and hasattr(self._variables['geos'], 'grid_mapping_name')
                    and self._variables['geos'].grid_mapping_name == 'vertical_perspective'):
                    gm = 'geos'
                else:
                    gm = variable.grid_mapping
                grid_mapping = self._variables[gm]
                if hasattr(grid_mapping, 'ellipsoid'):
                    kwargs_geom['geoid'] = {'ellps':grid_mapping.ellipsoid}
                elif hasattr(grid_mapping, 'earth_radius'):
                    kwargs_geom['geoid'] = {'a':grid_mapping.earth_radius,
                                            'b':grid_mapping.earth_radius}
                elif hasattr(grid_mapping, 'semi_major_axis') and hasattr(grid_mapping, 'semi_minor_axis'):
                    kwargs_geom['geoid'] = {'a':grid_mapping.semi_major_axis,
                                            'b':grid_mapping.semi_minor_axis}
                elif hasattr(grid_mapping, 'semi_major_axis') and hasattr(grid_mapping, 'inverse_flattening'):
                    kwargs_geom['geoid'] = {'a':grid_mapping.semi_major_axis,
                                            'rf':grid_mapping.inverse_flattening}
                else:
                    kwargs_geom['geoid'] = config.default_geoid
                if hasattr(grid_mapping, 'position_on_horizontal_grid'):
                    kwargs_geom['position_on_horizontal_grid'] = grid_mapping.position_on_horizontal_grid
                if grid_mapping.grid_mapping_name in ('lambert_conformal_conic',
                                                      'mercator',
                                                      'polar_stereographic',
                                                      'vertical_perspective'):
                    if (hasattr(grid_mapping, 'x_resolution') or
                        not behaviour.get('grid_is_lonlat', False)):
                        # if resolution is either in grid_mapping attributes or in the grid itself
                        geometryclass = ProjectedGeometry
                        kwargs_geom['name'] = _proj_dict_inv[grid_mapping.grid_mapping_name]
                        if hasattr(grid_mapping, 'x_resolution'):
                            Xresolution = grid_mapping.x_resolution
                            Yresolution = grid_mapping.y_resolution
                        else:
                            Xgrid, Ygrid = find_grid_in_variables()
                            if 1 in Xgrid.shape and behaviour.get('H1D_is_H2D_unstructured', False):
                                raise epygramError('unable to retrieve both X_resolution and Y_resolution from a 1D list of points.')
                            else:
                                Xresolution = abs(Xgrid[0, 0] - Xgrid[0, 1])
                                Yresolution = abs(Ygrid[0, 0] - Ygrid[1, 0])
                        grid = {'X_resolution':Xresolution,
                                'Y_resolution':Yresolution,
                                'LAMzone':None}
                        import pyproj
                        if kwargs_geom['name'] == 'lambert':
                            kwargs_geom['projection'] = {'reference_lon':Angle(grid_mapping.longitude_of_central_meridian, 'degrees'),
                                                         'rotation':Angle(0., 'degrees')}
                            if hasattr(grid_mapping, 'rotation'):
                                kwargs_geom['projection']['rotation'] = Angle(grid_mapping.rotation, 'degrees')
                            if isinstance(grid_mapping.standard_parallel, numpy.ndarray):
                                s1, s2 = grid_mapping.standard_parallel
                                kwargs_geom['projection']['secant_lat1'] = Angle(s1, 'degrees')
                                kwargs_geom['projection']['secant_lat2'] = Angle(s2, 'degrees')
                            else:
                                r = grid_mapping.standard_parallel
                                kwargs_geom['projection']['reference_lat'] = Angle(r, 'degrees')
                                s1 = s2 = r
                            fe = grid_mapping.false_easting
                            fn = grid_mapping.false_northing
                            reference_lat = grid_mapping.latitude_of_projection_origin

                            # compute x_0, y_0...
                            p = pyproj.Proj(proj='lcc',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_1=s1, lat_2=s2,
                                            **kwargs_geom['geoid'])
                            dx, dy = p(kwargs_geom['projection']['reference_lon'].get('degrees'),
                                       reference_lat)
                            # ... for getting center coords from false_easting, false_northing
                            p = pyproj.Proj(proj='lcc',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_1=s1, lat_2=s2,
                                            x_0=-dx, y_0=-dy,
                                            **kwargs_geom['geoid'])
                            ll00 = p(-fe, -fn, inverse=True)
                            del p
                            grid['input_lon'] = Angle(ll00[0], 'degrees')
                            grid['input_lat'] = Angle(ll00[1], 'degrees')
                            grid['input_position'] = (0, 0)
                        elif kwargs_geom['name'] == 'mercator':
                            kwargs_geom['projection'] = {'reference_lon':Angle(grid_mapping.longitude_of_central_meridian, 'degrees'),
                                                         'rotation':Angle(0., 'degrees')}
                            if hasattr(grid_mapping, 'rotation'):
                                kwargs_geom['projection']['rotation'] = Angle(grid_mapping.rotation, 'degrees')
                            kwargs_geom['projection']['reference_lat'] = Angle(0., 'degrees')
                            if grid_mapping.standard_parallel != 0.:
                                lat_ts = grid_mapping.standard_parallel
                                kwargs_geom['projection']['secant_lat'] = Angle(lat_ts, 'degrees')
                            else:
                                lat_ts = 0.
                            fe = grid_mapping.false_easting
                            fn = grid_mapping.false_northing
                            # compute x_0, y_0...
                            p = pyproj.Proj(proj='merc',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_ts=lat_ts,
                                            **kwargs_geom['geoid'])
                            dx, dy = p(kwargs_geom['projection']['reference_lon'].get('degrees'),
                                       0.)
                            # ... for getting center coords from false_easting, false_northing
                            p = pyproj.Proj(proj='merc',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_ts=lat_ts,
                                            x_0=-dx, y_0=-dy,
                                            **kwargs_geom['geoid'])
                            ll00 = p(-fe, -fn, inverse=True)
                            del p
                            grid['input_lon'] = Angle(ll00[0], 'degrees')
                            grid['input_lat'] = Angle(ll00[1], 'degrees')
                            grid['input_position'] = (0, 0)
                        elif kwargs_geom['name'] == 'polar_stereographic':
                            kwargs_geom['projection'] = {'reference_lon':Angle(grid_mapping.straight_vertical_longitude_from_pole, 'degrees'),
                                                         'rotation':Angle(0., 'degrees')}
                            if hasattr(grid_mapping, 'rotation'):
                                kwargs_geom['projection']['rotation'] = Angle(grid_mapping.rotation, 'degrees')
                            kwargs_geom['projection']['reference_lat'] = Angle(grid_mapping.latitude_of_projection_origin, 'degrees')
                            # CF-1.4
                            if hasattr(grid_mapping, 'scale_factor_at_projection_origin'):
                                if grid_mapping.scale_factor_at_projection_origin == 1.:
                                    lat_ts = grid_mapping.latitude_of_projection_origin
                                else:
                                    raise NotImplementedError("CF-1.4 encoding of secant parallel using 'scale_factor_at_projection_origin'.")
                            # CF-1.6
                            elif hasattr(grid_mapping, 'standard_parallel'):
                                lat_ts = grid_mapping.standard_parallel
                                if grid_mapping.standard_parallel != grid_mapping.latitude_of_projection_origin:
                                    kwargs_geom['projection']['secant_lat'] = Angle(lat_ts, 'degrees')
                            fe = grid_mapping.false_easting
                            fn = grid_mapping.false_northing
                            # compute x_0, y_0...
                            p = pyproj.Proj(proj='stere',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_0=kwargs_geom['projection']['reference_lat'].get('degrees'),
                                            lat_ts=lat_ts,
                                            **kwargs_geom['geoid'])
                            dx, dy = p(kwargs_geom['projection']['reference_lon'].get('degrees'),
                                       kwargs_geom['projection']['reference_lat'].get('degrees'),)
                            # ... for getting center coords from false_easting, false_northing
                            p = pyproj.Proj(proj='stere',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            lat_0=kwargs_geom['projection']['reference_lat'].get('degrees'),
                                            lat_ts=lat_ts,
                                            x_0=-dx, y_0=-dy,
                                            **kwargs_geom['geoid'])
                            ll00 = p(-fe, -fn, inverse=True)
                            del p
                            grid['input_lon'] = Angle(ll00[0], 'degrees')
                            grid['input_lat'] = Angle(ll00[1], 'degrees')
                            grid['input_position'] = (0, 0)
                        elif kwargs_geom['name'] == 'space_view':
                            kwargs_geom['projection'] = {'reference_lon':Angle(grid_mapping.longitude_of_central_meridian, 'degrees'),
                                                         'rotation':Angle(0., 'degrees')}
                            if hasattr(grid_mapping, 'rotation'):
                                kwargs_geom['projection']['rotation'] = Angle(grid_mapping.rotation, 'degrees')
                            r = grid_mapping.latitude_of_projection_origin
                            kwargs_geom['projection']['satellite_lat'] = Angle(r, 'degrees')
                            kwargs_geom['projection']['satellite_lon'] = Angle(grid_mapping.longitude_of_central_meridian, 'degrees')
                            kwargs_geom['projection']['satellite_height'] = grid_mapping.height_above_earth
                            fe = grid_mapping.false_easting
                            fn = grid_mapping.false_northing
                            reference_lat = grid_mapping.standard_parallel

                            # compute x_0, y_0...
                            p = pyproj.Proj(proj='geos',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            h=kwargs_geom['projection']['satellite_height'],
                                            **kwargs_geom['geoid'])
                            dx, dy = p(kwargs_geom['projection']['reference_lon'].get('degrees'),
                                       reference_lat)
                            # ... for getting center coords from false_easting, false_northing
                            p = pyproj.Proj(proj='geos',
                                            lon_0=kwargs_geom['projection']['reference_lon'].get('degrees'),
                                            h=kwargs_geom['projection']['satellite_height'],
                                            x_0=-dx, y_0=-dy,
                                            **kwargs_geom['geoid'])
                            ll00 = p(-fe, -fn, inverse=True)
                            del p
                            grid['input_lon'] = Angle(ll00[0], 'degrees')
                            grid['input_lat'] = Angle(ll00[1], 'degrees')
                            grid['input_position'] = (0, 0)
                    else:
                        # no resolution available: grid mapping is useless
                        gm = None
                elif 'gauss' in grid_mapping.grid_mapping_name.lower():
                    if hasattr(grid_mapping, 'latitudes'):
                        latitudes = self._variables[grid_mapping.latitudes.split(' ')[1]][:]
                    else:
                        # NOTE: this is a (good) approximation actually, the true latitudes are the roots of Legendre polynoms
                        raise NotImplementedError('(re-)computation of Gauss latitudes (not in file metadata).')
                    grid = {'latitudes':FPList([Angle(l, 'degrees') for l in latitudes]),
                            'dilatation_coef':1.}
                    if hasattr(grid_mapping, 'lon_number_by_lat'):
                        geometryclass = GaussGeometry
                        if isinstance(grid_mapping.lon_number_by_lat, str):
                            lon_number_by_lat = self._variables[grid_mapping.lon_number_by_lat.split(' ')[1]][:]
                        else:
                            kwargs_geom['name'] = 'regular_gauss'
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
                elif grid_mapping.grid_mapping_name == 'latitude_longitude':
                    # try to find out longitude, latitude arrays
                    for f in config.netCDF_usualnames_for_lonlat_grids['X']:
                        if f in self.listfields():
                            behaviour['X_grid'] = f
                            break
                    for f in config.netCDF_usualnames_for_lonlat_grids['Y']:
                        if f in self.listfields():
                            behaviour['Y_grid'] = f
                            break
                    Xgrid, Ygrid = find_grid_in_variables()
                    grid = {'longitudes':Xgrid,
                            'latitudes':Ygrid}
                    if ((hasattr(self._variables[variable.grid_mapping],
                                 'x_resolution') and
                         hasattr(self._variables[variable.grid_mapping],
                                 'y_resolution') or
                         hasattr(self._variables[variable.grid_mapping],
                                 'resolution')) and
                        len(Xgrid.shape) == 2):
                            # then this is a regular lon lat
                            geometryclass = RegLLGeometry
                            kwargs_geom['name'] = 'regular_lonlat'
                            if hasattr(self._variables[variable.grid_mapping],
                                       'x_resolution'):
                                x_res = grid_mapping.x_resolution
                            else:
                                x_res = grid_mapping.resolution
                            if hasattr(self._variables[variable.grid_mapping],
                                       'y_resolution'):
                                y_res = grid_mapping.y_resolution
                            else:
                                y_res = grid_mapping.resolution
                            grid = {'input_lon':Angle(Xgrid[0, 0], 'degrees'),
                                    'input_lat':Angle(Ygrid[0, 0], 'degrees'),
                                    'input_position':(0, 0),
                                    'X_resolution':Angle(x_res, 'degrees'),
                                    'Y_resolution':Angle(y_res, 'degrees')}
                else:
                    raise NotImplementedError('grid_mapping.grid_mapping_name == ' +
                                              grid_mapping.grid_mapping_name)
            else:
                if hasattr(variable, 'grid_mapping'):
                    epylog.info('grid_mapping ignored: unknown case')
                # grid only in variables
                Xgrid, Ygrid = find_grid_in_variables()
                if behaviour.get('grid_is_lonlat', False):
                    regll = False
                    if len(Xgrid.shape) == 2 and len(Ygrid.shape) == 2:
                        loose_epsilon = 1e-12  # config.epsilon was potentially too narrow
                        x_res = Xgrid[:,1:] - Xgrid[:,:-1]
                        y_res = Ygrid[1:,:] - Ygrid[:-1,:]
                        if Xgrid.shape[0] == 1:  # dim 1 on Y
                            y_res = x_res
                        elif Ygrid.shape[1] == 1:  # dim 1 on X
                            x_res = y_res
                        x_regular = numpy.all(x_res - x_res[0, 0] <= loose_epsilon)
                        y_regular = numpy.all(y_res - y_res[0, 0] <= loose_epsilon)
                        regll = x_regular and y_regular
                    if regll:
                        geometryclass = RegLLGeometry
                        kwargs_geom['name'] = 'regular_lonlat'
                        grid = {'input_lon': Angle(Xgrid[0, 0], 'degrees'),
                                'input_lat': Angle(Ygrid[0, 0], 'degrees'),
                                'input_position': (0, 0),
                                'X_resolution': Angle(x_res[0, 0], 'degrees'),
                                'Y_resolution': Angle(y_res[0, 0], 'degrees')}
                    else:
                        grid = {'longitudes':Xgrid,
                                'latitudes':Ygrid}
                else:
                    # grid is not lon/lat and no other metadata available : Academic
                    geometryclass = AcademicGeometry
                    kwargs_geom['name'] = 'academic'
                    grid = {'LAMzone':None,
                            'X_resolution':abs(Xgrid[0, 1] - Xgrid[0, 0]),
                            'Y_resolution':abs(Ygrid[1, 0] - Ygrid[0, 0]),
                            'input_lat':1.,
                            'input_lon':1.,
                            'input_position':(0, 0)}
                    kwargs_geom['projection'] = {'reference_dX':grid['X_resolution'],
                                                 'reference_dY':grid['Y_resolution'],
                                                 'rotation': Angle(0, 'degrees')}
        elif V1D or V2D or points:
            var_corresponding_to_X_grid = behaviour.get('X_grid', False)
            if var_corresponding_to_X_grid not in self._listfields():
                if points or V1D:
                    lon = ['_']
                else:
                    raise epygramError('unable to find X_grid in variables.')
            else:
                lon = self._variables[var_corresponding_to_X_grid][:]
            var_corresponding_to_Y_grid = behaviour.get('Y_grid', False)
            if var_corresponding_to_Y_grid not in self._listfields():
                if points or V1D:
                    lat = ['_']
                else:
                    raise epygramError('unable to find Y_grid in variables.')
            else:
                lat = self._variables[var_corresponding_to_Y_grid][:]
            grid = {'longitudes':lon,
                    'latitudes':lat}

        # 3.4 build geometry
        vcoordinate = VGeometry(**kwargs_vcoord)
        kwargs_geom['grid'] = grid
        kwargs_geom['dimensions'] = dimensions
        kwargs_geom['vcoordinate'] = vcoordinate
        geometry = geometryclass(**kwargs_geom)

        # 4. build field
        field_kwargs['geometry'] = geometry
        field_kwargs['structure'] = geometry.structure
        comment = {}
        for a in variable.ncattrs():
            if a != 'validity':
                # CLEANME: all these conversions might be useless since the use of _default_numpy2json ?
                if isinstance(variable.getncattr(a), numpy.float32):  # pb with json and float32
                    comment.update({a:numpy.float64(variable.getncattr(a))})
                elif isinstance(variable.getncattr(a), numpy.int32):  # pb with json and int32
                    comment.update({a:numpy.int64(variable.getncattr(a))})
                elif isinstance(variable.getncattr(a), numpy.ndarray):  # pb with json and numpy arrays
                    comment.update({a:numpy.float64(variable.getncattr(a)).tolist()})
                elif (isinstance(variable.getncattr(a), numpy.int16) or
                      isinstance(variable.getncattr(a), numpy.uint16)):  # pb with json and int16
                    comment.update({a:numpy.int64(variable.getncattr(a))})
                elif isinstance(variable.getncattr(a), numpy.int8):  # pb with json and int8
                    comment.update({a:numpy.int64(variable.getncattr(a))})
                else:
                    comment.update({a:variable.getncattr(a)})
        comment = json.dumps(comment, default=_default_numpy2json)
        if comment != '{}':
            field_kwargs['comment'] = comment
        field = fpx.field(**field_kwargs)
        if getdata:
            if only:
                n = len(variable.dimensions)
                buffdata = variable
                for k, i in only.items():
                    d = variable.dimensions.index(k)
                    if isinstance(i, slice):
                        ibuff = [util.restrain_to_index_i_of_dim_d(buffdata, j, d, n=n) for j in
                                 range(i.start, i.stop, i.step)]
                        buffdata = numpy.concatenate(ibuff, axis=d)
                    else:
                        buffdata = util.restrain_to_index_i_of_dim_d(buffdata, i, d, n=n)
            else:
                buffdata = variable[...]
            # check there is no leftover unknown dimension
            field_dim_num = 1 if len(field.validity) > 1 else 0
            if field.structure != 'Point':
                field_dim_num += [int(c) for c in field.structure if c.isdigit()][0]
                if (H2D or D3) and flattened:
                    field_dim_num -= 1
            assert field_dim_num == len(buffdata.squeeze().shape), \
                   ' '.join(['shape of field and identified usual dimensions',
                             'do not match: use *only* to filter or',
                             '*adhoc_behaviour* to identify dimensions'])
            # re-shuffle to have data indexes in order (t,z,y,x)
            positions = []
            shp4D = [1, 1, 1, 1]
            if 'T_dimension' in dims_dict_e2n:
                idx = variable.dimensions.index(dims_dict_e2n['T_dimension'])
                positions.append(idx)
                shp4D[0] = buffdata.shape[idx]
            if 'Z_dimension' in dims_dict_e2n:
                idx = variable.dimensions.index(dims_dict_e2n['Z_dimension'])
                positions.append(idx)
                shp4D[1] = buffdata.shape[idx]
            if 'Y_dimension' in dims_dict_e2n:
                idx = variable.dimensions.index(dims_dict_e2n['Y_dimension'])
                positions.append(idx)
                shp4D[2] = buffdata.shape[idx]
            if 'X_dimension' in dims_dict_e2n:
                idx = variable.dimensions.index(dims_dict_e2n['X_dimension'])
                positions.append(idx)
                shp4D[3] = buffdata.shape[idx]
            elif 'N_dimension' in dims_dict_e2n:
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
            if (H2D or D3) and flattened:
                if len(buffdata.shape) == 2:
                    if D3:
                        first_dimension = 'Z'
                    else:
                        first_dimension = 'T'
                else:
                    first_dimension = None
                data = geometry.reshape_data(buffdata, first_dimension=first_dimension, d4=True)
            else:
                data[...] = buffdata.reshape(data.shape)
            if behaviour.get('reverse_Yaxis'):
                data[...] = data[:, :, ::-1, :]
            field.setdata(data)

        return field

    def writefield(self, field,
                   compression=4,
                   metadata=None,
                   adhoc_behaviour=None,
                   fill_value=config.netCDF_default_variables_fill_value):
        """
        Write a field in resource.

        :param field: the :class:`~epygram.base.Field` object to write
        :param compression ranges from 1 (low compression, fast writing)
          to 9 (high compression, slow writing). 0 is no compression.
        :param metadata: dict, can be filled by any meta-data, that will be stored
          as attribute of the netCDF variable.
        :param adhoc_behaviour: to specify "on the fly" a behaviour (usual
          dimensions or grids, ...).
        """
        metadata = util.ifNone_emptydict(metadata)
        adhoc_behaviour = util.ifNone_emptydict(adhoc_behaviour)
        behaviour = self.behaviour.copy()
        behaviour.update(adhoc_behaviour)

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
            if str(varname) not in self._listfields():
                var = self._nc.createVariable(varname, vartype,
                                              dimensions=dimensions,
                                              **kwargs)
                status = 'created'
            else:
                assert self._variables[varname].dtype == vartype, \
                    ' '.join(['variable', varname,
                              'already exist with other type:',
                              str(self._variables[varname].dtype)])
                if isinstance(dimensions, str):
                    dimensions = (dimensions,)
                assert self._variables[varname].dimensions == tuple(dimensions), \
                       ' '.join(['variable', varname,
                                 'already exist with other dimensions:',
                                 str(self._variables[varname].dimensions)])
                var = self._variables[varname]
                status = 'match'
            return var, status

        assert 'netCDF' in field.fid
        assert not field.spectral

        # 1. dimensions
        T = Y = X = G = N = None
        dims = []
        # time
        if len(field.validity) > 1 or behaviour.get('force_a_T_dimension', False):
            # default
            T = config.netCDF_usualnames_for_standard_dimensions['T_dimension'][0]
            # or any existing identified time dimension
            T = {'found':v for v in self._dimensions
                 if (v in config.netCDF_usualnames_for_standard_dimensions['T_dimension'] and
                     len(self._dimensions[v]) == len(field.validity))}.get('found', T)
            # or specified behaviour
            T = behaviour.get('T_dimension', T)
            check_or_add_dim(T, size=len(field.validity))
        # vertical part
        # default
        Z = config.netCDF_usualnames_for_standard_dimensions['Z_dimension'][0]
        # or any existing identified time dimension
        Z = {'found':v for v in self._dimensions
                 if (v in config.netCDF_usualnames_for_standard_dimensions['Z_dimension'] and
                     len(self._dimensions[v]) == len(field.geometry.vcoordinate.levels))}.get('found', Z)
        # or specified behaviour
        Z = behaviour.get('Z_dimension', Z)
        if 'gridlevels' in field.geometry.vcoordinate.grid:
            Z_gridsize = max(len(field.geometry.vcoordinate.grid['gridlevels']), 1)
            if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119) and \
               field.geometry.vcoordinate.grid.get('ABgrid_position') != \
               field.geometry.vcoordinate.position_on_grid:
                Z_gridsize -= 1
        else:
            Z_gridsize = len(field.geometry.vcoordinate.levels)  # 1
        if Z_gridsize > 1:
            check_or_add_dim(Z, size=Z_gridsize)
        # horizontal
        if behaviour.get('flatten_horizontal_grids', False):
            _gpn = field.geometry.gridpoints_number
            G = 'gridpoints_number'
            check_or_add_dim(G, size=_gpn)
        if field.geometry.rectangular_grid:
            if field.geometry.dimensions['Y'] > 1 and field.geometry.dimensions['X'] > 1:
                Y = behaviour.get('Y_dimension',
                                  config.netCDF_usualnames_for_standard_dimensions['Y_dimension'][0])
                check_or_add_dim(Y, d_in_field='Y')
                X = behaviour.get('X_dimension',
                                  config.netCDF_usualnames_for_standard_dimensions['X_dimension'][0])
                check_or_add_dim(X, d_in_field='X')
            elif field.geometry.dimensions['X'] > 1:
                N = behaviour.get('N_dimension',
                                  config.netCDF_usualnames_for_standard_dimensions['N_dimension'][0])
                check_or_add_dim(N, d_in_field='X')
        elif 'gauss' in field.geometry.name:
            Y = behaviour.get('Y_dimension', 'latitude')
            check_or_add_dim(Y, d_in_field='lat_number')
            X = behaviour.get('X_dimension', 'longitude')
            check_or_add_dim(X, d_in_field='max_lon_number')
        else:
            raise NotImplementedError("grid not rectangular nor a gauss one.")

        # 2. validity
        # TODO: deal with unlimited time dimension ?
        if field.validity[0] != FieldValidity():
            tgrid = config.netCDF_usualnames_for_standard_dimensions['T_dimension'][0]
            tgrid = {'found':v for v in self._variables
                     if v in config.netCDF_usualnames_for_standard_dimensions['T_dimension']}.get('found', tgrid)
            tgrid = behaviour.get('T_grid', tgrid)
            if len(field.validity) > 1 or behaviour.get('force_a_T_dimension', False):
                _, _status = check_or_add_variable(tgrid, float, T)
                dims.append(T)  # FIXME: T instead of tgrid ?
            else:
                _, _status = check_or_add_variable(tgrid, float)
            datetime0 = field.validity[0].getbasis().isoformat(sep=str(' '))
            datetimes = [int((dt.get() - field.validity[0].getbasis()).total_seconds()) for dt in field.validity]
            if _status == 'created':
                self._variables[tgrid][:] = datetimes
                self._variables[tgrid].units = ' '.join(['seconds', 'since', datetime0])
            else:
                units = self._variables[tgrid].units.split(sep=' ')
                units_base = datetime.datetime.strptime(' '.join(units[2:4]), '%Y-%m-%d %H:%M:%S')
                if field.validity[0].getbasis() != units_base:
                    basediff = field.validity[0].getbasis() - units_base
                    datetimes = [dt + int(basediff.total_seconds()) for dt in datetimes]
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
            zgridname = behaviour.get('Z_grid', zgridname)
            if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
                if field.geometry.vcoordinate.grid.get('ABgrid_position') == 'mass' and \
                   field.geometry.vcoordinate.position_on_grid == 'mass':
                    ZP1 = Z
                    check_or_add_dim(ZP1, size=Z_gridsize)
                else:
                    ZP1 = Z + '+1'
                    check_or_add_dim(ZP1, size=Z_gridsize + 1)
                zgrid, _status = check_or_add_variable(zgridname, int)
                if _status == 'created':
                    zgrid.standard_name = _typeoffirstfixedsurface_dict_inv[field.geometry.vcoordinate.typeoffirstfixedsurface]
                    if field.geometry.vcoordinate.typeoffirstfixedsurface == 119:
                        zgrid.positive = "down"
                        zgrid.formula_terms = "ap: hybrid_coef_A b: hybrid_coef_B ps: surface_air_pressure"
                        check_or_add_variable('hybrid_coef_A',
                                              behaviour.get('meta_var_type',
                                                            config.netCDF_default_metavariables_dtype),
                                              ZP1)
                        self._variables['hybrid_coef_A'][:] = [iab[1]['Ai'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                        check_or_add_variable('hybrid_coef_B',
                                              behaviour.get('meta_var_type',
                                                            config.netCDF_default_metavariables_dtype),
                                              ZP1)
                        self._variables['hybrid_coef_B'][:] = [iab[1]['Bi'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                    elif field.geometry.vcoordinate.typeoffirstfixedsurface == 118:
                        # TOBECHECKED:
                        zgrid.positive = "up"
                        zgrid.formula_terms = "a: hybrid_coef_A b: hybrid_coef_B orog: orography"
                        check_or_add_variable('hybrid_coef_A',
                                              behaviour.get('meta_var_type',
                                                            config.netCDF_default_metavariables_dtype),
                                              ZP1)
                        self._variables['hybrid_coef_A'][:] = [iab[1]['Ai'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                        check_or_add_variable('hybrid_coef_B',
                                              behaviour.get('meta_var_type',
                                                            config.netCDF_default_metavariables_dtype),
                                              ZP1)
                        self._variables['hybrid_coef_B'][:] = [iab[1]['Bi'] for iab in field.geometry.vcoordinate.grid['gridlevels']]
                else:
                    epylog.info('assume 118/119 type vertical grid matches.')
            else:
                if len(numpy.shape(field.geometry.vcoordinate.levels)) > 1:
                    dims_Z = [d for d in [Z, Y, X, G, N] if d is not None]
                else:
                    dims_Z = Z
                zgrid, _status = check_or_add_variable(zgridname,
                                                       behaviour.get('vartype',
                                                                     config.netCDF_default_variables_dtype),
                                                       dims_Z)
                u = {102:'m', 103:'m', 100:'hPa'}.get(field.geometry.vcoordinate.typeoffirstfixedsurface, None)
                if u is not None:
                    zgrid.units = u
                if _status == 'created':
                    zgrid[:] = field.geometry.vcoordinate.levels
                else:
                    assert zgrid[:].all() == numpy.array(field.geometry.vcoordinate.levels).all(), \
                        ' '.join(['variable', zgrid, 'mismatch.'])
            if _typeoffirstfixedsurface_short_dict_inv.get(field.geometry.vcoordinate.typeoffirstfixedsurface, False):
                zgrid.short_name = _typeoffirstfixedsurface_short_dict_inv[field.geometry.vcoordinate.typeoffirstfixedsurface]
        # 3.2 grid (lonlat)
        dims_lonlat = []
        (lons, lats) = field.geometry.get_lonlat_grid()
        if behaviour.get('flatten_horizontal_grids'):
            dims_lonlat.append(G)
            dims.append(G)
            lons = stretch_array(lons)
            lats = stretch_array(lats)
        elif field.geometry.dimensions.get('Y', field.geometry.dimensions.get('lat_number', 0)) > 1:  # both Y and X dimensions
            dims_lonlat.extend([Y, X])
            dims.extend([Y, X])
        elif field.geometry.dimensions['X'] > 1:  # only X ==> N
            dims_lonlat.append(N)
            dims.append(N)
        # else: pass (single point or profile)
        if isinstance(lons, numpy.ma.masked_array):
            lonlat_fill_value = fill_value
            lons = lons.filled(lonlat_fill_value)
            lats = lats.filled(lonlat_fill_value)
        else:
            lonlat_fill_value = None
        try:
            _ = float(stretch_array(lons)[0])
        except ValueError:
            write_lonlat_grid = False
        else:
            write_lonlat_grid = behaviour.get('write_lonlat_grid', True)
        if write_lonlat_grid:
            lons_var, _status = check_or_add_variable(behaviour.get('X_grid', 'longitude'),
                                                      behaviour.get('vartype',
                                                                     config.netCDF_default_variables_dtype),
                                                      dims_lonlat,
                                                      fill_value=lonlat_fill_value)
            lats_var, _status = check_or_add_variable(behaviour.get('Y_grid', 'latitude'),
                                                      behaviour.get('vartype',
                                                                     config.netCDF_default_variables_dtype),
                                                      dims_lonlat,
                                                      fill_value=lonlat_fill_value)
            if _status == 'match':
                epylog.info('assume lons/lats match.')
            else:
                lons_var[...] = lons[...]
                lons_var.units = 'degrees'
                lats_var[...] = lats[...]
                lats_var.units = 'degrees'

        # 3.3 meta-data
        def set_ellipsoid(meta):
            if 'ellps' in field.geometry.geoid:
                self._variables[meta].ellipsoid = field.geometry.geoid['ellps']
            elif field.geometry.geoid.get('a', False) == field.geometry.geoid.get('b', True):
                self._variables[meta].earth_radius = field.geometry.geoid['a']
            elif field.geometry.geoid.get('a', False) and field.geometry.geoid.get('b', False):
                self._variables[meta].semi_major_axis = field.geometry.geoid['a']
                self._variables[meta].semi_minor_axis = field.geometry.geoid['b']
            elif field.geometry.geoid.get('a', False) and field.geometry.geoid.get('rf', False):
                self._variables[meta].semi_major_axis = field.geometry.geoid['a']
                self._variables[meta].inverse_flattening = field.geometry.geoid['rf']
            else:
                raise NotImplementedError('this kind of geoid:' + str(field.geometry.geoid))
        if field.geometry.dimensions.get('Y', field.geometry.dimensions.get('lat_number', 0)) > 1:
            if 'gauss' in field.geometry.name:
                # reduced Gauss case
                meta = 'Gauss_grid'
                _, _status = check_or_add_variable(meta, int)
                if _status == 'created':
                    self._variables[meta].grid_mapping_name = field.geometry.name + "_grid"
                    set_ellipsoid(meta)
                    if 'reduced' in field.geometry.name:
                        self._variables[meta].lon_number_by_lat = 'var: lon_number_by_lat'
                        check_or_add_variable('lon_number_by_lat', int, Y)
                        self._variables['lon_number_by_lat'][:] = field.geometry.dimensions['lon_number_by_lat']
                    self._variables[meta].latitudes = 'var: gauss_latitudes'
                    check_or_add_variable('gauss_latitudes', float, Y)
                    self._variables['gauss_latitudes'][:] = [l.get('degrees') for l in field.geometry.grid['latitudes']]
                    if 'pole_lon' in field.geometry.grid:
                        self._variables[meta].pole_lon = field.geometry.grid['pole_lon'].get('degrees')
                        self._variables[meta].pole_lat = field.geometry.grid['pole_lat'].get('degrees')
                    if 'dilatation_coef' in field.geometry.grid:
                        self._variables[meta].dilatation_coef = field.geometry.grid['dilatation_coef']
                else:
                    epylog.info('assume Gauss grid parameters match.')
            elif field.geometry.projected_geometry:
                # projections
                if field.geometry.name in ('lambert', 'mercator', 'polar_stereographic'):
                    meta = 'Projection_parameters'
                    _, _status = check_or_add_variable(meta, int)
                    if _status == 'created':
                        self._variables[meta].grid_mapping_name = _proj_dict[field.geometry.name]
                        set_ellipsoid(meta)
                        if field.geometry.position_on_horizontal_grid != 'center':
                            self._variables[meta].position_on_horizontal_grid = field.geometry.position_on_horizontal_grid
                        self._variables[meta].x_resolution = field.geometry.grid['X_resolution']
                        self._variables[meta].y_resolution = field.geometry.grid['Y_resolution']
                        if field.geometry.grid.get('LAMzone'):
                            _lon_cen, _lat_cen = field.geometry.ij2ll(float(field.geometry.dimensions['X'] - 1.) / 2.,
                                                                      float(field.geometry.dimensions['Y'] - 1.) / 2.)
                        else:
                            _lon_cen = field.geometry._center_lon.get('degrees')
                            _lat_cen = field.geometry._center_lat.get('degrees')
                        if field.geometry.name == 'lambert':
                            if field.geometry.secant_projection:
                                std_parallel = [field.geometry.projection['secant_lat1'].get('degrees'),
                                                field.geometry.projection['secant_lat2'].get('degrees')]
                                latitude_of_projection_origin = (std_parallel[0] + std_parallel[1]) / 2.
                            else:
                                std_parallel = field.geometry.projection['reference_lat'].get('degrees')
                                latitude_of_projection_origin = std_parallel
                            x00, y00 = field.geometry.ij2xy(0, 0)
                            x0, y0 = field.geometry.ll2xy(field.geometry.projection['reference_lon'].get('degrees'),
                                                          latitude_of_projection_origin)
                            (dx, dy) = (x00 - x0, y00 - y0)
                            self._variables[meta].longitude_of_central_meridian = field.geometry.projection['reference_lon'].get('degrees')
                            self._variables[meta].latitude_of_projection_origin = latitude_of_projection_origin
                            self._variables[meta].standard_parallel = std_parallel
                            self._variables[meta].false_easting = -dx
                            self._variables[meta].false_northing = -dy
                        elif field.geometry.name == 'mercator':
                            if field.geometry.secant_projection:
                                std_parallel = field.geometry.projection['secant_lat'].get('degrees')
                            else:
                                std_parallel = field.geometry.projection['reference_lat'].get('degrees')
                            x00, y00 = field.geometry.ij2xy(0, 0)
                            x0, y0 = field.geometry.ll2xy(field.geometry.projection['reference_lon'].get('degrees'), 0.)
                            (dx, dy) = (x00 - x0, y00 - y0)
                            self._variables[meta].longitude_of_central_meridian = field.geometry.projection['reference_lon'].get('degrees')
                            self._variables[meta].standard_parallel = std_parallel
                            self._variables[meta].false_easting = -dx
                            self._variables[meta].false_northing = -dy
                        elif field.geometry.name == 'polar_stereographic':
                            if field.geometry.secant_projection:
                                std_parallel = field.geometry.projection['secant_lat'].get('degrees')
                            else:
                                std_parallel = field.geometry.projection['reference_lat'].get('degrees')
                            x00, y00 = field.geometry.ij2xy(0, 0)
                            x0, y0 = field.geometry.ll2xy(field.geometry.projection['reference_lon'].get('degrees'),
                                                          field.geometry.projection['reference_lat'].get('degrees'))
                            (dx, dy) = (x00 - x0, y00 - y0)
                            self._variables[meta].straight_vertical_longitude_from_pole = field.geometry.projection['reference_lon'].get('degrees')
                            self._variables[meta].latitude_of_projection_origin = field.geometry.projection['reference_lat'].get('degrees')
                            self._variables[meta].standard_parallel = std_parallel
                            self._variables[meta].false_easting = -dx
                            self._variables[meta].false_northing = -dy
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
        _, _status = check_or_add_variable(varname,
                                           behaviour.get('vartype',
                                                         config.netCDF_default_variables_dtype),
                                           dims,
                                           zlib=bool(compression),
                                           complevel=compression,
                                           fill_value=fill_value)
        if meta:
            self._variables[varname].grid_mapping = meta
        if field.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
            self._variables[varname].vertical_grid = zgridname
        data = field.getdata(d4=True)
        if isinstance(data, numpy.ma.masked_array):
            if 'gauss' in field.geometry.name:
                data = field.geometry.fill_maskedvalues(data)
            else:
                data = data.filled(fill_value)
        if behaviour.get('flatten_horizontal_grids'):
            data = field.geometry.horizontally_flattened(data)
        data = data.squeeze()
        if _status == 'match':
            epylog.info('overwrite data in variable ' + varname)
        self._variables[varname][...] = data
        if field.units not in (None, ''):
            self._variables[varname].units = field.units

        # 5. metadata
        for k, v in metadata.items():
            self._variables[varname].setncattr(k, v)

    def set_default_global_attributes(self):
        """
        Set default global attributes (those from
        config.netCDF_default_global_attributes).
        """
        self._nc.setncatts(config.netCDF_default_global_attributes)

    def set_global_attributes(self, **attributes):
        """Set the given global attributes."""
        self._nc.setncatts(attributes)

    def behave(self, **kwargs):
        """
        Set-up the given arguments in self.behaviour, for the purpose of
        building fields from netCDF.
        """
        self.behaviour.update(kwargs)

    # adapted from http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    def ncdump(self, out):
        """
        Outputs dimensions, variables and their attribute information.
        The information is similar to that of NCAR's ncdump utility.

        :param out: IO where nc_attrs, nc_dims, and nc_vars are printed

        Returns
        -------
        nc_attrs : list
            A Python list of the NetCDF file global attributes
        nc_dims : list
            A Python list of the NetCDF file dimensions
        nc_vars : list
            A Python list of the NetCDF file variables
        """
        nc = self._nc
        
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
                    outwrite('\t\t%s:' % ncattr,
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

    def what(self, out=sys.stdout, **_):
        """Writes in file a summary of the contents of the GRIB."""
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        self.ncdump(out)

    def _listfields(self):
        """Returns the fid list of the fields inside the resource."""
        return list(self._variables.keys())

    @property
    @FileResource._openbeforedelayed
    def _dimensions(self):
        return self._nc.dimensions

    @property
    @FileResource._openbeforedelayed
    def _variables(self):
        return self._nc.variables
