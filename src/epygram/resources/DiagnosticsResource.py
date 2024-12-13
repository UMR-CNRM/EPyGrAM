#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
Contains the class that handle a DiagnosticsResource.
This resource exposes fields diagnosed from the low level resource fields
in addition to the fields already present in it.

Some diagnostics need approximation or options. Here are the lists:
List of approximations:
    - approx_ff: see diag_ff_from_U_V_W
    - approx_gIsConstant: see diag_equiv_hgeopot_hgeom
    - approx_R: see diag_R_from_q
    - approx_qt: see diag_qt_from_q

List of options:
    - center: see readfields
    - pos_P_from_sigma: see diag_P_from_hybridP
    - vertical_mean_P_from_sigma: see diag_P_from_hybridP and diag_Hgeopot_from_T_q
    - cheap_height: see diag_Hgeopot_from_T_q

The possibility to compute fields depend on the available fields in the low-
level resource. Some circumstances have not been tested and are likely to
produce strange behavior. Among them:
    - if 3D fields exist in the low level resource at the same time as equivalent 2D fields
    - if fields are defined with several validities (especially for subdomain extraction)
    - if fields are present on several geometries in the same low-level resource

Fields are identified using their generic fid which is the GRIB code. Some codes
are not well defined. Here are some explanations on how they are used in this module:
    - graupel and hail are not difficult to use because GRIB code does not make difference
                       between mixing-ratio and specific contents for this categories
                       We choose specific content of graupel to be parameter 201
                                 mixing ratio of graupel to be 253
                                 specific content of hail to be 31
                                 mixing ratio of hail to be 252
    - R is defined here as being 'discipline':0, 'parameterCategory':1, 'parameterNumber':254
    - qt is defined here as being 'discipline':0, 'parameterCategory':1, 'parameterNumber':251

Diagnostics that could be easily added:
    - conversions (both ways) between specific contents and mixing-ratios
    - conversions to and from relative humidity
    - computation of wet density
More tricky diagnostics:
    - LWP computation; normally feasible on hybrid levels because in this case we have a description
      of the whole atmosphere. Can we put typeoffirstfixedsurface in asked fid and return an fid
      without it? This could be the way to specify the kind of levels to use in the computation
      and abort if this is different from 118 and 119.
"""

import numpy

from footprints import FPDict, proxy as fpx
from bronx.meteo import constants
from bronx.meteo.conversion import q2R

from epygram.base import Resource, FieldSet
from epygram import epygramError
from epygram import profiles
from epygram.geometries.VGeometry import hybridP2pressure, hybridH2altitude

##########
# Temporary extension, waiting for bronx update
constants.P0 = 100000.


class DiagnosticsResource(Resource):
    """Class implementing a DiagnosticsResource."""

    _collector = ('resource_modificator', 'epyresource')
    _footprint = dict(
        attr=dict(
            resource=dict(
                type=Resource,
                info="Low level resource"),
            name=dict(
                values=set(['Diagnostics'])),
            approx=dict(
                type=FPDict,
                optional=True,
                default=FPDict(),
                info="dictionary used to define the degree of approximation to use."),
            options=dict(
                type=FPDict,
                optional=True,
                default=FPDict(),
                info="dictionary used to define some options.")
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(Resource, self).__init__(*args, **kwargs)

        if self.resource.openmode != self.openmode:
            raise epygramError("The low level resource must be opened using the same mode as this high-level resource.")

        self.format = "Diagnostics"
        self._diagMethod = None
        self._cache = {}  # cache for various method, key is the method name

    def open(self):
        """Opens the low level resource"""

        if not self.resource.isopen:
            self.resource.open()

    def close(self):
        """Closes the low level resource."""
        try:
            self.resource.close()
        except Exception:
            pass

    def find_fields_in_resource(self, seed=None, generic=False, **kwargs):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        Args: \n
        - *seed*: might be a 'handgrip', i.e. a dict where you can store all
          requested keys,
          e.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl', 'level':850},
          a list of handgrips or *None*. If *None* (default), returns the list
          of all fields in resource.
        - *generic*: if True, returns a list of tuples (fid, fid) of
          the fields (for mimetism with other formats).
        """
        if seed is None or isinstance(seed, dict):
            fieldslist = self.listfields(select=seed)
        elif isinstance(seed, list):
            fieldslist = []
            for s in seed:
                fieldslist.extend(self.listfields(select=s))
        else:
            raise epygramError("unknown type for seed: " + str(type(seed)))
        if fieldslist == []:
            raise epygramError("no field matching '" + str(seed) +
                               "' was found in resource ")
        if generic:
            fieldslist = [(fieldslist[i], fieldslist[i]) for i in range(len(fieldslist))]

        return fieldslist

    @staticmethod
    def _hash_dict(d):
        """Returns a tuple of key/value"""
        return tuple([(k, d[k]) for k in sorted(d.keys())])

    def _create_list(self):
        """Creates the list of available fields associated with the original fids."""

        if self._diagMethod is None:
            # self._diagMethod is a dictionary that associate to each available field the way to produce it
            self._diagMethod_done = False  # Boolean to inform the different diag that all available fields are not yet known
            self._diagMethod = {}
            # First: fields directly available in the low level resource
            for fid in self.resource.listfields(complete=True):
                self._diagMethod[self._hash_dict(fid['generic'])] = {'generic':dict(fid['generic']),
                                                                     'method':None,
                                                                     'original_fid':fid[self.resource.format]}
            # Second: list fields that can be produced by the different diag
            #         must be done several times to deal with dependencies
            #         list is complete when do not evolve during iteration
            previousLength = -1
            diags = [method for method in dir(self) if method.startswith('diag_') and callable(getattr(self, method))]
            while len(self._diagMethod) > previousLength:
                previousLength = len(self._diagMethod)
                for method in diags:  # loop over the different diags
                    for fid in getattr(self, method)(dolist=True):  # Loop over the fields the diag can produce
                        hash_fid = self._hash_dict(fid['generic'])
                        if hash_fid not in self._diagMethod:
                            # We keep the first solution found to select the solution with the least conversions
                            self._diagMethod[hash_fid] = {'generic':dict(fid['generic']), 'method':method}
            self._diagMethod_done = True  # All fields are known

        return self._diagMethod.values()

    def listfields(self, onlykey=None, select=None, complete=False):
        """Lists the available fields."""

        fidlist = [v['generic'] for v in self._create_list()]
        if select is not None:
            fidlist = [f for f in fidlist if all([(k in f and f[k] == select[k]) for k in select.keys()])]
        if onlykey is not None:
            if isinstance(onlykey, str):
                fidlist = [f[onlykey] for f in fidlist]
            elif isinstance(onlykey, tuple):
                fidlist = [{k:f[k] for k in onlykey} for f in fidlist]
        if complete:
            fidlist = [{'generic': fid, self.format:fid} for fid in fidlist]
        return fidlist

    def sortfields(self, sortingkey, onlykey=None):
        """
        Returns a sorted list of fields with regards to the given *sortingkey*
        of their fid, as a dict of lists.

        Argument *onlykey* can be specified as a string or a tuple of strings,
        so that only specified keys of the fid will returned.
        """
        # Taken from GRIB.py

        sortedfields = {}
        listoffields = self.listfields()
        onlykeylistoffields = self.listfields(onlykey=onlykey)
        for f in range(len(listoffields)):
            try:
                category = listoffields[f][sortingkey]
            except KeyError:
                category = 'None'
            field = onlykeylistoffields[f]
            if category in sortedfields:  # FIXME: ? sortedfields is empty...
                sortedfields[category].append(field)
            else:
                sortedfields[category] = [field]

        return sortedfields

    def readfield(self, handgrip, getdata=True):
        """Read the field in the low level resource and join the levels."""

        result = self.readfields(handgrip=handgrip, getdata=getdata)
        if len(result) != 1:
            raise epygramError(str(len(result)) + " field(s) have been found, one and only one expected (handgrib=" + str(handgrip) + ").")
        return result[0]

    def readfields(self, handgrip, getdata=True):
        """
        Read the field in the low level resource and join the levels.
        Options:
        - center: if True, center the field on the grid (shuman transformation)
        """
        fieldset = FieldSet()
        fieldsGroupedByMethod = {}
        original_fids = []
        for fid in self.listfields(select=handgrip):
            found = None
            for v in self._create_list():
                if fid == v['generic']:
                    if found is not None:
                        raise epygramError("Internal error...")
                    else:
                        found = v['generic']
                        fieldsGroupedByMethod[v['method']] = fieldsGroupedByMethod.get(v['method'], []) + [found]
                        if v['method'] is None:
                            original_fids.append(v['original_fid'])
            if found is None:
                raise epygramError("Internal error....")
        for method, fids in fieldsGroupedByMethod.items():
            if method is None:
                for original_fid in original_fids:
                    for field in self.resource.readfields(original_fid, getdata=getdata):
                        field.fid[self.format] = field.fid['generic']
                        if self.options.get("center", False):
                            field.center()
                        fieldset.append(field)
            else:
                fieldset.extend(getattr(self, method)(fids, getdata=getdata))

        return fieldset

    def writefield(self, *args, **kwargs):
        """Write field."""
        self.resource.writefield(*args, **kwargs)

    def writefields(self, *args, **kwargs):
        """Write fields."""
        self.resource.writefields(*args, **kwargs)

    def _get_field3d_from_handrgrip(self, handgrip):
        """
        Returns the requested 3D field described by the handgrip

        :param handgrip:  MUST define the parameter and the type of levels
        """
        if len(self.find_fields_in_resource(handgrip)) > 1:
            field3d = fpx.field(fid={self.format:handgrip},
                                structure='3D',
                                resource=self, resource_fids=[handgrip])
        else:
            field3d = self.readfield(handgrip)

        return field3d

    def extractprofile(self, handgrip, lon=None, lat=None,
                       geometry=None,
                       vertical_coordinate=None,
                       interpolation='nearest'):
        """
        Extracts a vertical profile from the resource, given its handgrip
        and the geographic location (*lon*/*lat*) of the profile.

        :param handgrip: MUST define the parameter and the type of levels
        :param lon: is the longitude of the desired point.
        :param lat: is the latitude of the desired point.
        :param geometry: is the geometry on which extract data.
                         If None, it is built from lon/lat.
        :param vertical_coordinate: defines the requested vertical coordinate of the
                                    V1DField (as number of GRIB2 norm:
                                    http://apps.ecmwf.int/codes/grib/format/grib2/ctables/4/5).
        :param interpolation: defines the interpolation function used to compute
                              the profile at requested lon/lat from the fields grid:
                                  - if 'nearest' (default), extracts profile at the
                                    horizontal nearest neighboring gridpoint;
                                  - if 'linear', computes profile with horizontal
                                    linear spline interpolation;
                                  - if 'cubic', computes profile with horizontal
                                    cubic spline interpolation.
        """

        field3d = self._get_field3d_from_handrgrip(handgrip)

        if geometry is None:
            if None in [lon, lat]:
                raise epygramError("You must give a geometry or lon *and* lat")
            pointG = field3d.geometry.make_profile_geometry(lon, lat)
        else:
            if lon is not None or lat is not None:
                raise epygramError("You cannot provide lon or lat when geometry is given")
            if geometry.structure != "V1D":
                raise epygramError("geometry must be a V1D")
            pointG = geometry

        profile = self.extract_subdomain(handgrip, pointG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         field3d=field3d)

        return profile

    def extractsection(self, handgrip, end1=None, end2=None,
                       geometry=None, points_number=None,
                       resolution=None, vertical_coordinate=None,
                       interpolation='linear',
                       global_shift_center=None):
        """
        Extracts a vertical section from the resource, given its handgrip
        and the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param handgrip: MUST define the parameter and the type of levels
        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param geometry: is the geometry on which extract data. If None, defaults to
                         linearily spaced positions computed from  *points_number*.
        :param points_number: defines the total number of horizontal points of the
                              section (including ends). If None, defaults to a number
                              computed from the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
                           field. If None, defaults to the horizontal resolution
                           of the field.
        :param vertical_coordinate: defines the requested vertical coordinate of the
                                    V2DField (cf. :module:`epygram.geometries`
                                    coordinate possible values).
        :param interpolation: defines the interpolation function used to compute
                              the profile points locations from the fields grid: \n
                                - if 'nearest', each horizontal point of the section is
                                  taken as the horizontal nearest neighboring gridpoint;
                                - if 'linear' (default), each horizontal point of the section is
                                  computed with linear spline interpolation;
                                - if 'cubic', each horizontal point of the section is
                                  computed with linear spline interpolation.
        :param global_shift_center: for global lon/lat grids, shift the center by the
            requested angle (in degrees). Enables a [0,360] grid
            to be shifted to a [-180,180] grid, for instance (with -180 argument).
        """
        field3d = self._get_field3d_from_handrgrip(handgrip)
        if global_shift_center is not None:
            field3d.global_shift_center(global_shift_center)
        if geometry is None:
            if None in [end1, end2]:
                raise epygramError("You must give a geometry or end1 *and* end2")
            sectionG = field3d.geometry.make_section_geometry(end1, end2,
                                                              points_number=points_number,
                                                              resolution=resolution)
        else:
            if end1 is not None or end2 is not None:
                raise epygramError("You cannot provide end1 or end2 when geometry is given")
            if geometry.structure != "V2D":
                raise epygramError("geometry must be a V2D")
            sectionG = geometry

        section = self.extract_subdomain(handgrip, sectionG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         field3d=field3d)

        return section

    def extract_subdomain(self, handgrip, geometry, vertical_coordinate=None,
                          interpolation='linear', exclude_extralevels=True, field3d=None):
        """
        Extracts a subdomain from the resource, given its handgrip
        and the geometry to use.

        :param handgrip: MUST define the parameter and the type of levels
        :param geometry: is the geometry on which extract data.
        :param vertical_coordinate: defines the requested vertical coordinate of the
                                    V2DField (cf. :module:`epygram.geometries`
                                    coordinate possible values).
        :param interpolation: defines the interpolation function used to compute
                              the profile points locations from the fields grid: \n
                                - if 'nearest', each horizontal point of the section is
                                  taken as the horizontal nearest neighboring gridpoint;
                                - if 'linear' (default), each horizontal point of the section is
                                  computed with linear spline interpolation;
                                - if 'cubic', each horizontal point of the section is
                                  computed with linear spline interpolation.
        :param exclude_extralevels: if True, not physical levels are removed
        """

        if field3d is None:
            field3d = self._get_field3d_from_handrgrip(handgrip)

        if field3d.spectral:
            field3d.sp2gp()
        subdomain = field3d.extract_subdomain(geometry, interpolation=interpolation)
        if exclude_extralevels:
            subdomain = subdomain.extract_physicallevels()

        # preparation for vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):

            vertical_fid = {100:{'discipline':0, 'parameterCategory':3, 'parameterNumber':0},
                            102:{'discipline':0, 'parameterCategory':3, 'parameterNumber':6},
                            103:{'discipline':0, 'parameterCategory':3, 'parameterNumber':6}
                            }
            if vertical_coordinate not in vertical_fid:
                raise NotImplementedError("this vertical coordinate conversion.")
            vertical_fid = vertical_fid[vertical_coordinate]
            if 'typeOfFirstFixedSurface' in handgrip:
                vertical_fid['typeOfFirstFixedSurface'] = handgrip['typeOfFirstFixedSurface']

            vertical_field = self._get_field3d_from_handrgrip(vertical_fid)
            if vertical_field.spectral:
                vertical_field.sp2gp()
            levels = vertical_field.extract_subdomain(subdomain.geometry, interpolation=interpolation)
            if exclude_extralevels:
                levels = levels.extract_physicallevels()

            if vertical_coordinate == 100:
                pass
            elif vertical_coordinate == 103:
                surface_height = self.readfield({'discipline':2, 'parameterCategory':0, 'parameterNumber':7})
                surface_geom = geometry.deepcopy()
                surface_geom.vcoordinate = surface_height.geometry.vcoordinate
                surface_height = surface_height.extract_subdomain(surface_geom, interpolation=interpolation)
                if exclude_extralevels:
                    surface_height = surface_height.extract_physicallevels()
                if surface_height.spectral:
                    surface_height.sp2gp()
                levels.setdata(levels.getdata() - surface_height.getdata())

            subdomain.use_field_as_vcoord(levels, force_kind=vertical_coordinate)

        return subdomain

    @staticmethod
    def _common_levels(parameters):
        """Look for common levels between several sets of generic fids"""
        # we must exclude unknown levels
        setList = []
        for p in parameters:
            pSet = set()
            for fid in p:
                f = fid.copy()
                f.pop('parameterCategory', None)
                f.pop('discipline', None)
                f.pop('parameterNumber', None)
                f.pop('productDefinitionTemplateNumber', None)
                if (('level' not in f) or f['level'] != 255) and \
                   (('typeOfFirstFixedSurface' not in f) or f['typeOfFirstFixedSurface'] != 255):
                    pSet.add(DiagnosticsResource._hash_dict(f))
            setList.append(pSet)
        result = []
        for r in set.intersection(*setList):
            result.append({k:v for (k, v) in r})
        return result

    def _check_compatibilty(self, fid1=None, fid2=None, field1=None, field2=None):
        """Checks if two fields are compatibles (validity, horizontal geometry, point position)"""
        if field1 is None and fid1 is None:
            raise epygramError("field1 and fid1 cannot be None at the same time")
        if field2 is None and fid2 is None:
            raise epygramError("field2 and fid2 cannot be None at the same time")

        if field1 is not None and fid1 is not None:
            if field1.fid[self.format] != fid1:
                raise epygramError("fids are not equal")
        if field2 is not None and fid2 is not None:
            if field2.fid[self.format] != fid2:
                raise epygramError("fids are not equal")

        if field1 is not None and fid1 is None:
            fid1 = field1.fid[self.format]
        if field2 is not None and fid2 is None:
            fid2 = field2.fid[self.format]

        if 'check_compatibilty' not in self._cache:
            self._cache['check_compatibilty'] = {}

        hfid1 = self._hash_dict(fid1)
        hfid2 = self._hash_dict(fid2)
        if not (hfid1, hfid2) in self._cache['check_compatibilty']:
            if field1 is None:
                field1 = self._get_field(fid1, fid1, getdata=False)
            if field2 is None:
                field2 = self._get_field(fid2, fid2, getdata=False)
            ok = True
            if field1.validity != field2.validity:
                ok = False
            if not field1.geometry.eq_Hgeom(field2.geometry):
                ok = False
            if 'gridlevels' in field1.geometry.vcoordinate.grid and \
               'gridlevels' in field2.geometry.vcoordinate.grid:
                if field1.geometry.vcoordinate.position_on_grid != field2.geometry.vcoordinate.position_on_grid:
                    ok = False
            self._cache['check_compatibilty'][(hfid1, hfid2)] = ok
        return self._cache['check_compatibilty'][(hfid1, hfid2)]

    def _diag_checks(self, fids, dolist):
        """Common checks for all diags"""
        if dolist and len(fids) != 0:
            raise epygramError("if dolist, no fid must be given")
        if len(fids) == 0 and not dolist:
            raise epygramError("Nothing to do")
        if not dolist:
            if not all([fid in self.listfields() for fid in fids]):
                raise epygramError("Some fids are not available")

    def _simple_HDiag(self, inputs, output):
        """Common list algorithm for horizontal diag"""

        result = []
        try:
            parameters = [self.find_fields_in_resource(seed) for seed in inputs]
        except:
            return result
        for level in self._common_levels(parameters):
            InputFids = []
            for seed in inputs:
                fid = level.copy()
                fid.update(seed)
                InputFids.append(fid)
            if all([self._check_compatibilty(InputFids[0], InputFids[i + 1]) for i in range(len(InputFids) - 1)]):
                fid = level.copy()
                fid.update(output)
                result.append({'generic':fid})
        return result

    def _get_field(self, param, level, getdata, deepcopy=False):
        """Wrapper to simplify field reading"""
        myfid = level.copy()
        myfid.update(param)
        if 'level' in level:
            field = self.readfield(myfid, getdata=getdata)
        else:
            # We must select the field that do not have level in its fid
            fields = self.readfields(myfid, getdata=getdata)
            for fid in [field.fid for field in fields]:
                if 'level' in fid[self.format]:
                    fields.remove(field.fid)
            if len(fields) != 1:
                raise epygramError("Zero or more than 1 field correspond to fid")
            field = fields[0]
        if deepcopy:
            field = field.deepcopy()
        if field.spectral and getdata:
            field.sp2gp()
        return field

    def diag_equiv_geopot_hgeopot(self, *args, **kwargs):
        equiv = ({'discipline':0, 'parameterCategory':3, 'parameterNumber':4},  # Geopotential
                 {'discipline':0, 'parameterCategory':3, 'parameterNumber':5},  # Geopotential height
                 True, constants.g0, None)  # Equivalent with approximation, factor is g
        return self._diag_equivalence(equiv, *args, **kwargs)

    def diag_equiv_hgeopot_hgeom(self, *args, **kwargs):
        equiv = ({'discipline':0, 'parameterCategory':3, 'parameterNumber':5},  # Geopotential height
                 {'discipline':0, 'parameterCategory':3, 'parameterNumber':6},  # Geometric height
                 self.approx.get('approx_gIsConstant', 0) == 1, 1, None)  # Equivalent with approximation, no factor
        return self._diag_equivalence(equiv, *args, **kwargs)

    def diag_equiv_model_terrain_height(self, *args, **kwargs):
        equiv = ({'discipline':2, 'parameterCategory':0, 'parameterNumber':7},  # Model terrain height
                 {'discipline':0, 'parameterCategory':3, 'parameterNumber':6, 'typeOfFirstFixedSurface':1},  # Geometric height
                 True, 1, None)  # Always equivalent, no factor
        return self._diag_equivalence(equiv, *args, **kwargs)

    def diag_equiv_P_logP(self, *args, **kwargs):
        equiv = ({'discipline':0, 'parameterCategory':3, 'parameterNumber':25},  # log(P)
                 {'discipline':0, 'parameterCategory':3, 'parameterNumber':0},  # P
                 True, 1, numpy.log)  # Always equivalent, no factor
        return self._diag_equivalence(equiv, *args, **kwargs)

    def _diag_equivalence(self, equiv, fids=[], dolist=False, getdata=True):
        """
        Some fields can receive multiple grib code
        equiv is a tuple, whose components are fid1, fid2, ok, factor, operation
            ok to enable or not this conversion
            factor and operation are such that read(fid1)=operation(factor*read(fid2))
        """
        if dolist:
            result = []
            if equiv[2]:
                for e in [(equiv[0], equiv[1]), (equiv[1], equiv[0])]:
                    try:
                        fidlist = self.find_fields_in_resource(e[0])
                    except:
                        fidlist = []
                    for fid in fidlist:
                        fid = fid.copy()
                        fid.update(e[1])
                        result.append({'generic':fid})
            return result
        else:
            fieldset = FieldSet()
            for fid in fids:
                if equiv[2]:
                    for e in [(equiv[0], equiv[1], equiv[3], equiv[4]), (equiv[1], equiv[0], 1. / equiv[3], {numpy.log:numpy.exp,
                                                                                                             numpy.exp:numpy.log,
                                                                                                             None:None}[equiv[4]])]:
                        if all([k in fid and fid[k] == e[0][k] for k in e[0]]):
                            myfid = fid.copy()
                            myfid.update(e[1])
                            field = self.readfield(myfid, getdata=getdata)
                            if e[2] != 1 or e[3] is not None:
                                for fmt in [fmt for fmt in field.fid.keys() if fmt not in ['generic', self.format]]:
                                    field.fid.pop(fmt)
                            fid_generic = field.fid.get('generic', {})
                            fid_generic.update(e[0])
                            fid_format = field.fid.get(self.format, {})
                            fid_format.update(e[0])
                            field.fid['generic'] = fid_generic
                            field.fid[self.format] = fid_format
                            if getdata:
                                if field.spectral:
                                    field.sp2gp()
                                field.setdata(field.getdata() * e[2])
                                if e[3] is not None:
                                    field.setdata(e[3](field.getdata()))
                            fieldset.append(field)
            return fieldset

    def diag_P_from_Plevels(self, fids=[], dolist=False, getdata=True):
        """Returns pressure on levels with constant pressure."""
        self._diag_checks(fids, dolist)
        fid_P = {'discipline': 0, 'parameterCategory': 3, 'parameterNumber': 0, 'productDefinitionTemplateNumber':0}
        fid_field = {'typeOfFirstFixedSurface': 100}

        if dolist:
            # We need a cache to know in which order fids have been found
            # to not try to use a fid depending on the current diag
            if 'P_from_Plevels' not in self._cache:
                self._cache['P_from_Plevels'] = {'levels':set(), 'fids':[]}
            try:
                fidlist = self.find_fields_in_resource(fid_field)
            except:
                fidlist = []
            for fid in fidlist:
                if fid['level'] not in self._cache['P_from_Plevels']['levels']:
                    self._cache['P_from_Plevels']['levels'].add(fid['level'])
                    if fid not in self._cache['P_from_Plevels']['fids']:
                        self._cache['P_from_Plevels']['fids'].append(fid)
            return [{'generic':dict(level=level, typeOfFirstFixedSurface=100, **fid_P)} for level in self._cache['P_from_Plevels']['levels']]
        else:
            fieldset = FieldSet()
            for fid in fids:
                fidlist = []
                for one_fid in self._cache['P_from_Plevels']['fids']:
                    fidlist.extend(self.find_fields_in_resource(one_fid))
                for f in fidlist:
                    if not all([f[k] == v for k, v in fid_P.iteritems()]):
                        readfid = f
                        break
                field = self.readfield(readfid, getdata=getdata)
                for fmt in list(field.fid.keys()):
                    field.fid.pop(fmt)
                field.fid[self.format] = dict(level=fid['level'], typeOfFirstFixedSurface=100, **fid_P)
                field.fid['generic'] = dict(level=fid['level'], typeOfFirstFixedSurface=100, **fid_P)
                if len(field.geometry.vcoordinate.levels) != 1:
                    raise NotImplementedError("Pb could be easily resolved here")
                    # To solve pb, one must set a new numpy array for data instead of using field.getdata()
                del field.geometry.vcoordinate.levels[:]
                field.geometry.vcoordinate.levels.append(fid['level'])
                if getdata:
                    if field.spectral:
                        field.sp2gp()
                    data = field.getdata()
                    data[...] = fid['level'] * 100  # convert hPa into Pa
                    field.setdata(data)
                fieldset.append(field)
            return fieldset

    def diag_theta_from_T_P(self, fids=[], dolist=False, getdata=True):
        """Computes theta from T and P (no approximation)"""
        self._diag_checks(fids, dolist)
        fid_T = {'discipline':0, 'parameterCategory':0, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_P = {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_Theta = {'discipline':0, 'parameterCategory':0, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}

        if dolist:
            return self._simple_HDiag([fid_T, fid_P], fid_Theta)
        else:
            fieldset = FieldSet()
            for fid in fids:
                T = self._get_field(fid_T, fid, getdata, True)
                if T.spectral and getdata:
                    T.sp2gp()
                P = self._get_field(fid_P, fid, getdata)
                if not self._check_compatibilty(field1=T, field2=P):
                    raise epygramError("Internal error: T and P fields are not compatible")
                if getdata:
                    pass
                    T.setdata(T.getdata() * (constants.P0 / P.getdata()) ** (constants.Rd / constants.Cpd))
                myfid_Theta = fid.copy()
                myfid_Theta.update(fid_Theta)
                for fmt in list(T.fid.keys()):
                    T.fid.pop(fmt)
                T.fid[self.format] = myfid_Theta
                T.fid['generic'] = myfid_Theta
                fieldset.append(T)
            return fieldset

    def diag_thetaV_from_theta_qv_qt(self, fids=[], dolist=False, getdata=True):
        """Computes theta from T and P (no approximation)"""
        self._diag_checks(fids, dolist)
        fid_theta = {'discipline':0, 'parameterCategory':0, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}
        fid_qv = {'discipline':0, 'parameterCategory':1, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_qt = {'discipline':0, 'parameterCategory':1, 'parameterNumber':251, 'productDefinitionTemplateNumber':0}
        fid_ThetaV = {'discipline':0, 'parameterCategory':0, 'parameterNumber':15, 'productDefinitionTemplateNumber':0}

        if dolist:
            return self._simple_HDiag([fid_theta, fid_qv, fid_qt], fid_ThetaV)
        else:
            fieldset = FieldSet()
            for fid in fids:
                theta = self._get_field(fid_theta, fid, getdata, True)
                qv = self._get_field(fid_qv, fid, getdata)
                qt = self._get_field(fid_qt, fid, getdata)
                if not self._check_compatibilty(field1=theta, field2=qv):
                    raise epygramError("Internal error: theta and qv fields are not compatible")
                if not self._check_compatibilty(field1=qt, field2=qv):
                    raise epygramError("Internal error: qt and qv fields are not compatible")
                if getdata:
                    theta.setdata(theta.getdata() * (1 + qv.getdata() * constants.Rv / constants.Rd - qt.getdata()))
                myfid_ThetaV = fid.copy()
                myfid_ThetaV.update(fid_ThetaV)
                for fmt in list(theta.fid.keys()):
                    theta.fid.pop(fmt)
                theta.fid[self.format] = myfid_ThetaV
                theta.fid['generic'] = myfid_ThetaV
                fieldset.append(theta)
            return fieldset

    def diag_T_from_Theta_P(self, fids=[], dolist=False, getdata=True):
        """Computes T from theta and P (no approximation)"""
        self._diag_checks(fids, dolist)
        fid_T = {'discipline':0, 'parameterCategory':0, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_P = {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_Theta = {'discipline':0, 'parameterCategory':0, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}

        if dolist:
            return self._simple_HDiag([fid_P, fid_Theta], fid_T)
        else:
            fieldset = FieldSet()
            for fid in fids:
                Theta = self._get_field(fid_Theta, fid, getdata, True)
                P = self._get_field(fid_P, fid, getdata)
                if not self._check_compatibilty(field1=Theta, field2=P):
                    raise epygramError("Internal error: Theta and P fields are not compatible")
                if getdata:
                    Theta.setdata(Theta.getdata() * (P.getdata() / constants.P0) ** (constants.Rd / constants.Cpd))
                myfid_Theta = fid.copy()
                myfid_Theta.update(fid_T)
                for fmt in list(Theta.fid.keys()):
                    Theta.fid.pop(fmt)
                Theta.fid[self.format] = myfid_Theta
                Theta.fid['generic'] = myfid_Theta
                fieldset.append(Theta)
            return fieldset

    def diag_qt_from_q(self, fids=[], dolist=False, getdata=True):
        """
        Computes qt from the hydrometeor contents
        Approximation:
        - approx_qt: 0 to not approximate (we need qv, qc, qr, qi, qs, qg and qh)
                     1 we only need qv
                     2 we only need qv, qc and qi
                     3 we only need qv, qc, qi, qr, qs and qg
        """
        self._diag_checks(fids, dolist)
        fid_qv = {'discipline':0, 'parameterCategory':1, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_qc = {'discipline':0, 'parameterCategory':1, 'parameterNumber':83, 'productDefinitionTemplateNumber':0}
        fid_qr = {'discipline':0, 'parameterCategory':1, 'parameterNumber':85, 'productDefinitionTemplateNumber':0}
        fid_qi = {'discipline':0, 'parameterCategory':1, 'parameterNumber':84, 'productDefinitionTemplateNumber':0}
        fid_qs = {'discipline':0, 'parameterCategory':1, 'parameterNumber':86, 'productDefinitionTemplateNumber':0}
        fid_qg = {'discipline':0, 'parameterCategory':1, 'parameterNumber':201, 'productDefinitionTemplateNumber':0}
        fid_qh = {'discipline':0, 'parameterCategory':1, 'parameterNumber':31, 'productDefinitionTemplateNumber':0}
        fid_qt = {'discipline':0, 'parameterCategory':1, 'parameterNumber':251, 'productDefinitionTemplateNumber':0}

        approx_qt = self.approx.get('approx_qt', 0)

        if dolist:
            if approx_qt == 0:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qr, fid_qi,
                                           fid_qs, fid_qg, fid_qh], fid_qt)
            elif approx_qt == 1:
                return self._simple_HDiag([fid_qv], fid_qt)
            elif approx_qt == 2:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qi], fid_qt)
            elif approx_qt == 3:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qr, fid_qi,
                                           fid_qs, fid_qg], fid_qt)
            else:
                raise ValueError("Unknown value for approx_R: " + str(approx_qt))
        else:
            fieldset = FieldSet()
            for fid in fids:
                qv = self._get_field(fid_qv, fid, getdata, True)
                if getdata:
                    if qv.spectral:
                        qv.gp2sp()
                    hydrometeors = {}
                    if approx_qt in [0, 2, 3]:
                        hydrometeors['ql'] = self._get_field(fid_qc, fid, getdata).getdata()
                        hydrometeors['qi'] = self._get_field(fid_qi, fid, getdata).getdata()
                    if approx_qt in [0, 3]:
                        hydrometeors['qr'] = self._get_field(fid_qr, fid, getdata).getdata()
                        hydrometeors['qs'] = self._get_field(fid_qs, fid, getdata).getdata()
                        hydrometeors['qg'] = self._get_field(fid_qg, fid, getdata).getdata()
                    if approx_qt == 0:
                        hydrometeors['qh'] = self._get_field(fid_qh, fid, getdata).getdata()
                    qt = qv.getdata()
                    for h in hydrometeors.values():
                        qt += h
                    qv.setdata(qt)
                myfid_qt = fid.copy()
                myfid_qt.update(fid_qt)
                for fmt in list(qv.fid.keys()):
                    qv.fid.pop(fmt)
                qv.fid[self.format] = myfid_qt
                qv.fid['generic'] = myfid_qt
                fieldset.append(qv)
            return fieldset

    def diag_qt_from_r(self, fids=[], dolist=False, getdata=True):
        """
        Computes qt from the hydrometeor mixing ratio
        Approximation:
        - approx_qt: 0 to not approximate (we need qv, qc, qr, qi, qs, qg and qh)
                     1 we only need qv
                     2 we only need qv, qc and qi
                     3 we only need qv, qc, qi, qr, qs and qg
        """
        self._diag_checks(fids, dolist)
        fid_rv = {'discipline':0, 'parameterCategory':1, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}
        fid_rc = {'discipline':0, 'parameterCategory':1, 'parameterNumber':22, 'productDefinitionTemplateNumber':0}
        fid_rr = {'discipline':0, 'parameterCategory':1, 'parameterNumber':24, 'productDefinitionTemplateNumber':0}
        fid_ri = {'discipline':0, 'parameterCategory':1, 'parameterNumber':23, 'productDefinitionTemplateNumber':0}
        fid_rs = {'discipline':0, 'parameterCategory':1, 'parameterNumber':25, 'productDefinitionTemplateNumber':0}
        fid_rg = {'discipline':0, 'parameterCategory':1, 'parameterNumber':253, 'productDefinitionTemplateNumber':0}
        fid_rh = {'discipline':0, 'parameterCategory':1, 'parameterNumber':252, 'productDefinitionTemplateNumber':0}
        fid_qt = {'discipline':0, 'parameterCategory':1, 'parameterNumber':251, 'productDefinitionTemplateNumber':0}

        approx_qt = self.approx.get('approx_qt', 0)

        if dolist:
            if approx_qt == 0:
                return self._simple_HDiag([fid_rv, fid_rc, fid_rr, fid_ri,
                                           fid_rs, fid_rg, fid_rh], fid_qt)
            elif approx_qt == 1:
                return self._simple_HDiag([fid_rv], fid_qt)
            elif approx_qt == 2:
                return self._simple_HDiag([fid_rv, fid_rc, fid_ri], fid_qt)
            elif approx_qt == 3:
                return self._simple_HDiag([fid_rv, fid_rc, fid_rr, fid_ri,
                                           fid_rs, fid_rg], fid_qt)
            else:
                raise ValueError("Unknown value for approx_R: " + str(approx_qt))
        else:
            fieldset = FieldSet()
            for fid in fids:
                rv = self._get_field(fid_rv, fid, getdata, True)
                if getdata:
                    if rv.spectral:
                        rv.gp2sp()
                    hydrometeors = {}
                    if approx_qt in [0, 2, 3]:
                        hydrometeors['rl'] = self._get_field(fid_rc, fid, getdata).getdata()
                        hydrometeors['ri'] = self._get_field(fid_ri, fid, getdata).getdata()
                    if approx_qt in [0, 3]:
                        hydrometeors['rr'] = self._get_field(fid_rr, fid, getdata).getdata()
                        hydrometeors['rs'] = self._get_field(fid_rs, fid, getdata).getdata()
                        hydrometeors['rg'] = self._get_field(fid_rg, fid, getdata).getdata()
                    if approx_qt == 0:
                        hydrometeors['rh'] = self._get_field(fid_rh, fid, getdata).getdata()
                    qt = rv.getdata()  # mixing ratio
                    for h in hydrometeors.values():
                        qt += h  # mixing ratio
                    qt = qt / (1 + qt)  # convert mixing ratio into specific content
                    rv.setdata(qt)
                myfid_qt = fid.copy()
                myfid_qt.update(fid_qt)
                for fmt in list(rv.fid.keys()):
                    rv.fid.pop(fmt)
                rv.fid[self.format] = myfid_qt
                rv.fid['generic'] = myfid_qt
                fieldset.append(rv)
            return fieldset

    def _diag_q_from_r_qt(self, fid_r, fid_q, fids, dolist, getdata):
        fid_qt = {'discipline':0, 'parameterCategory':1, 'parameterNumber':251, 'productDefinitionTemplateNumber':0}
        if dolist:
            return self._simple_HDiag([fid_r, fid_qt], fid_q)
        else:
            fieldset = FieldSet()
            for fid in fids:
                qt = self._get_field(fid_qt, fid, getdata, True)
                r = self._get_field(fid_r, fid, getdata)
                if not self._check_compatibilty(field1=qt, field2=r):
                    raise epygramError("Internal error: qt and r fields are not compatible")
                if getdata:
                    qt.setdata(r.getdata() * (1. - qt.getdata()))
                myfid_Theta = fid.copy()
                myfid_Theta.update(fid_q)
                for fmt in list(qt.fid.keys()):
                    qt.fid.pop(fmt)
                qt.fid[self.format] = myfid_Theta
                qt.fid['generic'] = myfid_Theta
                fieldset.append(qt)
            return fieldset

    def diag_qv_from_rv(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qc_from_rc(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':22, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':83, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qr_from_rr(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':24, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':85, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qi_from_ri(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':23, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':84, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qs_from_rs(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':25, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':86, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qg_from_rg(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':253, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':201, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_qh_from_rh(self, fids=[], dolist=False, getdata=True):
        fid_r = {'discipline':0, 'parameterCategory':1, 'parameterNumber':252, 'productDefinitionTemplateNumber':0}
        fid_q = {'discipline':0, 'parameterCategory':1, 'parameterNumber':31, 'productDefinitionTemplateNumber':0}
        return self._diag_q_from_r_qt(fid_r, fid_q, fids, dolist, getdata)

    def diag_R_from_q(self, fids=[], dolist=False, getdata=True):
        """
        Computes R from the hydrometeor contents
        Approximation:
        - approx_R: 0 to not approximate (we need qv, qc, qr, qi, qs, qg and qh)
                    1 we only need qv
                    2 we only need qv, qc and qi
                    3 we only need qv, qc, qi, qr, qs and qg
        """
        self._diag_checks(fids, dolist)
        fid_qv = {'discipline':0, 'parameterCategory':1, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        fid_qc = {'discipline':0, 'parameterCategory':1, 'parameterNumber':83, 'productDefinitionTemplateNumber':0}
        fid_qr = {'discipline':0, 'parameterCategory':1, 'parameterNumber':85, 'productDefinitionTemplateNumber':0}
        fid_qi = {'discipline':0, 'parameterCategory':1, 'parameterNumber':84, 'productDefinitionTemplateNumber':0}
        fid_qs = {'discipline':0, 'parameterCategory':1, 'parameterNumber':86, 'productDefinitionTemplateNumber':0}
        fid_qg = {'discipline':0, 'parameterCategory':1, 'parameterNumber':201, 'productDefinitionTemplateNumber':0}
        fid_qh = {'discipline':0, 'parameterCategory':1, 'parameterNumber':31, 'productDefinitionTemplateNumber':0}
        fid_R = {'discipline':0, 'parameterCategory':1, 'parameterNumber':254, 'productDefinitionTemplateNumber':0}

        approx_R = self.approx.get('approx_R', 0)

        if dolist:
            if approx_R == 0:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qr, fid_qi,
                                           fid_qs, fid_qg, fid_qh], fid_R)
            elif approx_R == 1:
                return self._simple_HDiag([fid_qv], fid_R)
            elif approx_R == 2:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qi], fid_R)
            elif approx_R == 3:
                return self._simple_HDiag([fid_qv, fid_qc, fid_qr, fid_qi,
                                           fid_qs, fid_qg], fid_R)
            else:
                raise ValueError("Unknown value for approx_R: " + str(approx_R))
        else:
            fieldset = FieldSet()
            for fid in fids:
                qv = self._get_field(fid_qv, fid, getdata, True)
                if getdata:
                    if qv.spectral:
                        qv.gp2sp()
                    hydrometeors = {}
                    if approx_R in [0, 2, 3]:
                        hydrometeors['ql'] = self._get_field(fid_qc, fid, getdata).getdata()
                        hydrometeors['qi'] = self._get_field(fid_qi, fid, getdata).getdata()
                    if approx_R in [0, 3]:
                        hydrometeors['qr'] = self._get_field(fid_qr, fid, getdata).getdata()
                        hydrometeors['qs'] = self._get_field(fid_qs, fid, getdata).getdata()
                        hydrometeors['qg'] = self._get_field(fid_qg, fid, getdata).getdata()
                    if approx_R == 0:
                        hydrometeors['qh'] = self._get_field(fid_qh, fid, getdata).getdata()
                    R = q2R(qv.getdata(), **hydrometeors)
                    qv.setdata(R)
                myfid_R = fid.copy()
                myfid_R.update(fid_R)
                for fmt in list(qv.fid.keys()):
                    qv.fid.pop(fmt)
                qv.fid[self.format] = myfid_R
                qv.fid['generic'] = myfid_R
                fieldset.append(qv)
            return fieldset

    def diag_ff_from_U_V_W(self, fids=[], dolist=False, getdata=True):
        """
        Computes the wind speed from its components.
        *approx_ff*: 0: we need the 3 components
                     1: we only need the 2 horizontal components
        """
        self._diag_checks(fids, dolist)
        fid_U = {'discipline':0, 'parameterCategory':2, 'parameterNumber':2, 'productDefinitionTemplateNumber':0}
        fid_V = {'discipline':0, 'parameterCategory':2, 'parameterNumber':3, 'productDefinitionTemplateNumber':0}
        fid_W = {'discipline':0, 'parameterCategory':2, 'parameterNumber':9, 'productDefinitionTemplateNumber':0}
        fid_ff = {'discipline':0, 'parameterCategory':2, 'parameterNumber':1, 'productDefinitionTemplateNumber':0}
        approx_ff = self.approx.get('approx_ff', 0)

        if dolist:
            if approx_ff == 0:
                return self._simple_HDiag([fid_U, fid_V, fid_W], fid_ff)
            else:
                return self._simple_HDiag([fid_U, fid_V], fid_ff)
        else:
            fieldset = FieldSet()
            for fid in fids:
                U = self._get_field(fid_U, fid, getdata, True)
                V = self._get_field(fid_V, fid, getdata)
                if not self._check_compatibilty(field1=U, field2=V):
                    raise epygramError("Internal error: U and V fields are not compatible")
                if getdata:
                    U.setdata(U.getdata() ** 2 + V.getdata() ** 2)
                if approx_ff == 0:
                    W = self._get_field(fid_W, fid, getdata)
                    if not self._check_compatibilty(field1=U, field2=W):
                        raise epygramError("Internal error or approx needed: U and W fields are not compatible")
                    if getdata:
                        U.setdata(U.getdata() + W.getdata() ** 2)
                if getdata:
                    U.setdata(numpy.sqrt(U.getdata()))
                myfid_ff = fid.copy()
                myfid_ff.update(fid_ff)
                for fmt in list(U.fid.keys()):
                    U.fid.pop(fmt)
                U.fid[self.format] = myfid_ff
                U.fid['generic'] = myfid_ff
                fieldset.append(U)
            return fieldset

    def _hybridP(self):
        """
        Helper method that returns information to deal with hybrid-P coordinates
        """
        fid_Ps = [{'discipline':0, 'parameterCategory':3, 'parameterNumber':25, 'typeOfFirstFixedSurface':1, 'productDefinitionTemplateNumber':0},
                  {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'typeOfFirstFixedSurface':1, 'productDefinitionTemplateNumber':0}]

        fid_for_vcoord = None
        if 'hybridP' not in self._cache:
            fid_for_vcoord = [fid[self.resource.format] for fid in self.resource.listfields(complete=True)
                              if (fid['generic'].get('typeOfFirstFixedSurface', 255) == 119 and
                                  fid['generic'].get('parameterNumber', 255) != 255)] #to be able to read it from a CombineLevels resource
            if len(fid_for_vcoord) != 0:
                hybridP_geometry = self.resource.readfield(fid_for_vcoord[0], getdata=False).geometry.vcoordinate
                A = [level[1]['Ai'] for level in hybridP_geometry.grid['gridlevels']][1:]
                B = [level[1]['Bi'] for level in hybridP_geometry.grid['gridlevels']][1:]
                if not len(A) == len(B):
                    raise ValueError("A, B must have the same size.")
                for fid in fid_Ps:
                    try:
                        fid_Ps_complete = self.find_fields_in_resource(fid_Ps)
                        break
                    except:
                        fid_Ps_complete = None
                if fid_Ps_complete is None:
                    result = None
                else:
                    field_3d = len(hybridP_geometry.levels) > 1  # True if the first field on hybrid-P level is 3D
                    result = dict(A=A, B=B, ABgrid_position=hybridP_geometry.grid['ABgrid_position'],
                                  fid_for_vcoord=fid_for_vcoord[0], field_3d=field_3d,
                                  vgeometry=hybridP_geometry)
            else:
                result = None

            if self._diagMethod_done or result is not None:
                self._cache['hybridP'] = result
        else:
            result = self._cache['hybridP']

        return result

    def diag_P_from_hybridP(self, fids=[], dolist=False, getdata=True):
        """
        Computes P on hybrid-P levels
        Options:
        - pos_P_from_sigma: 'mass' to compute P on mass levels, 'flux' to return P on flux levels (defaults to 'mass')
        - vertical_mean_P_from_sigma: defines the kind of averaging done on the vertical
               to compute half-levels from full-levels, or inverse: 'geometric' or 'arithmetic'
        """
        self._diag_checks(fids, dolist)
        fid_Ps = {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'typeOfFirstFixedSurface':1, 'productDefinitionTemplateNumber':0}
        fid_P = {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'productDefinitionTemplateNumber':0}
        option_pos = self.options.get('pos_P_from_sigma', 'mass')
        option_vertical_mean = self.options.get('vertical_mean_P_from_sigma', None)
        if option_pos not in ['mass', 'flux']:
            raise epygramError("pos_P_from_sigma != 'mass' or 'flux'.")
        if getdata and option_vertical_mean is None and not dolist:
            raise epygramError("You must define the vertical_mean_P_from_sigma option to use this diagnostic")

        hybridP = self._hybridP()

        if dolist:
            if hybridP is not None:
                A = hybridP['A']
                if hybridP['field_3d']:
                    # First field on hybrid-P level is 3D, so we put the 3D P field in the list
                    result = [{'generic':dict(typeOfFirstFixedSurface=119, **fid_P)}]
                else:
                    # First field on hybrid-P level is 2D, so we put the 2D P fields in the list
                    result = [{'generic':dict(typeOfFirstFixedSurface=119, level=l, **fid_P)} for l in range(1, len(A) + 1)]
            else:
                result = []
            return result
        else:
            Ps = self.readfield(fid_Ps, getdata=getdata)
            if Ps.spectral and getdata:
                Ps.sp2gp()
            if getdata and Ps.getdata().flatten()[0] < 100:
                # Ps contains lnsp and not sp
                # This is a workaround for a bug in FA
                Ps.setdata(numpy.exp(Ps.getdata()))
            if getdata:
                pressure = hybridP2pressure(hybridP['vgeometry'], Ps.getdata(), option_vertical_mean,
                                            gridposition=option_pos).levels
                pressure = numpy.array(pressure) * 100
            fieldset = FieldSet()
            for fid in fids:
                P = self.resource.readfield(hybridP['fid_for_vcoord'], getdata=False)
                geometry = P.geometry.deepcopy()
                if 'level' in fid:
                    geometry.vcoordinate.levels = [fid['level']]
                else:
                    A = hybridP['A']
                    geometry.vcoordinate.levels = list(range(1, len(A) + 1))
                geometry.vcoordinate.position_on_grid = option_pos
                myfid_P = fid.copy()
                myfid_P.update(fid_P)
                myfid_P = {self.format:myfid_P, 'generic':myfid_P}
                P = fpx.field(structure=geometry.structure,
                              geometry=geometry,
                              validity=Ps.validity,
                              processtype=Ps.processtype,
                              fid=myfid_P)
                if getdata:
                    if 'level' in fid:
                        P.setdata(pressure[fid['level'] - 1])
                    else:
                        P.setdata(pressure)
                fieldset.append(P)
            return fieldset

    def _hybridH(self):
        """
        Helper method that returns information to deal with hybrid-H coordinates
        """
        fid_Zs = {'discipline':2, 'parameterCategory':0, 'parameterNumber':7, 'productDefinitionTemplateNumber':0}

        fid_for_vcoord = None
        if 'hybridH' not in self._cache:
            fid_for_vcoord = [fid[self.resource.format] for fid in self.resource.listfields(complete=True) if fid['generic'].get('typeOfFirstFixedSurface', 255) == 118]
            if len(fid_for_vcoord) != 0:
                hybridH_geometry = self.resource.readfield(fid_for_vcoord[0], getdata=False).geometry.vcoordinate
                A = [level[1]['Ai'] for level in hybridH_geometry.grid['gridlevels']][1:]
                B = [level[1]['Bi'] for level in hybridH_geometry.grid['gridlevels']][1:]
                if not len(A) == len(B):
                    raise ValueError("A, B must have the same size.")
                try:
                    fid_Zs = self.find_fields_in_resource(fid_Zs)
                except:
                    fid_Zs = None
                if fid_Zs is None:
                    result = None
                else:
                    field_3d = len(hybridH_geometry.levels) > 1  # True if the first field on hybrid-H level is 3D
                    result = dict(A=A, B=B, ABgrid_position=hybridH_geometry.grid['ABgrid_position'],
                                  fid_for_vcoord=fid_for_vcoord[0], field_3d=field_3d,
                                  vgeometry=hybridH_geometry)
            else:
                result = None

            if self._diagMethod_done or result is not None:
                self._cache['hybridH'] = result
        else:
            result = self._cache['hybridH']

        return result

    def diag_H_from_hybridH(self, fids=[], dolist=False, getdata=True):
        """
        Computes H on hybrid-H levels
        Options:
        - pos_H_from_sigma: 'mass' to compute H on mass levels, 'flux' to return H on flux levels (defaults to 'mass')
        """
        self._diag_checks(fids, dolist)
        fid_Zs = {'discipline':2, 'parameterCategory':0, 'parameterNumber':7, 'productDefinitionTemplateNumber':0}
        fid_Z = {'discipline':0, 'parameterCategory':3, 'parameterNumber':5, 'productDefinitionTemplateNumber':0}
        option_pos = self.options.get('pos_H_from_sigma', 'mass')
        if option_pos not in ['mass', 'flux']:
            raise epygramError("pos_H_from_sigma != 'mass' or 'flux'.")

        hybridH = self._hybridH()

        if dolist:
            if hybridH is not None:
                A = hybridH['A']
                if hybridH['field_3d']:
                    # First field on hybrid-H level is 3D, so we put the 3D H field in the list
                    result = [{'generic':dict(typeOfFirstFixedSurface=118, **fid_Z)}]
                else:
                    # First field on hybrid-H level is 2D, so we put the 2D H fields in the list
                    result = [{'generic':dict(typeOfFirstFixedSurface=118, level=l, **fid_Z)} for l in range(0, len(A) + 2)]
            else:
                result = []
            return result
        else:
            Zs = self.readfield(fid_Zs, getdata=getdata)
            if getdata:
                A = hybridH['A']
                del hybridH['vgeometry'].levels[:]
                hybridH['vgeometry'].levels.extend(range(0, len(A) + 2))
                if Zs.spectral:
                    Zs.sp2gp()
                altitude = hybridH2altitude(hybridH['vgeometry'], Zs.getdata(),
                                            gridposition=option_pos,
                                            conv2height=False).levels
            fieldset = FieldSet()
            for fid in fids:
                Z = self.resource.readfield(hybridH['fid_for_vcoord'], getdata=False)
                geometry = Z.geometry.deepcopy()
                geometry.position_on_horizontal_grid = Zs.geometry.position_on_horizontal_grid
                if 'level' in fid:
                    geometry.vcoordinate.levels = [fid['level']]
                else:
                    geometry.vcoordinate.levels = list(range(0, len(A) + 2))
                geometry.vcoordinate.position_on_grid = option_pos
                myfid_Z = fid.copy()
                myfid_Z.update(fid_Z)
                myfid_Z = {self.format:myfid_Z, 'generic':myfid_Z}
                Z = fpx.field(structure=geometry.structure,
                              geometry=geometry,
                              validity=Zs.validity,
                              processtype=Zs.processtype,
                              fid=myfid_Z)
                if getdata:
                    if 'level' in fid:
                        Z.setdata(altitude[fid['level']])
                    else:
                        Z.setdata(altitude)
                fieldset.append(Z)
            return fieldset


class DiagnosticsAROMEResource(DiagnosticsResource):
    """Class implementing a DiagnosticsResource for AROME."""

    _collector = ('resource_modificator', 'epyresource')
    _footprint = dict(
        attr=dict(
            name=dict(
                values=set(['DiagnosticsAROME']))
        )
    )

    def diag_equiv_model_terrain_height_AROME(self, *args, **kwargs):
        equiv = ({'discipline':0, 'parameterCategory':193, 'parameterNumber':5, 'typeOfFirstFixedSurface':1},
                 {'discipline':0, 'parameterCategory':3, 'parameterNumber':4, 'typeOfFirstFixedSurface':1},
                 True, 1, None)  # Always equivalent, no factor
        return self._diag_equivalence(equiv, *args, **kwargs)
    
    def diag_Hgeopot_from_T_q(self, fids=[], dolist=False, getdata=True):
        """
        Computes Z on hybrid-P levels
        Options:
        - cheap_height: 0 to compute R using all species and to use pressure departure
                        1 only use qv to compute R (even if other species are present) and do not use pressure departure
        - vertical_mean_P_from_sigma: defines the kind of averaging done on the vertical
               to compute half-levels from full-levels, or inverse: 'geometric' or 'arithmetic'
        Approximation:
        - approx_R is used, see. diag_R_from_q
        """
        self._diag_checks(fids, dolist)
        fid_Ps = {'discipline':0, 'parameterCategory':3, 'parameterNumber':0, 'typeOfFirstFixedSurface':1, 'productDefinitionTemplateNumber':0}
        fid_T = {'discipline':0, 'parameterCategory':0, 'parameterNumber':0, 'typeOfFirstFixedSurface':119, 'productDefinitionTemplateNumber':0}
        fid_Pdep = {'discipline':0, 'parameterCategory':193, 'parameterNumber':4, 'typeOfFirstFixedSurface':119, 'productDefinitionTemplateNumber':0}
        fid_R = {'discipline':0, 'parameterCategory':1, 'parameterNumber':254, 'typeOfFirstFixedSurface': 119, 'productDefinitionTemplateNumber':0}
        fid_PhiSurf = {'discipline':0, 'parameterCategory':3, 'parameterNumber':4, 'typeOfFirstFixedSurface': 1, 'productDefinitionTemplateNumber':0}
        fid_HGeopot = {'discipline':0, 'parameterCategory':3, 'parameterNumber':5, 'typeOfFirstFixedSurface':119, 'productDefinitionTemplateNumber':0}

        option_vertical_mean = self.options.get('vertical_mean_P_from_sigma', None)
        option_cheap_height = self.options.get('cheap_height', 0)

        hybridP = self._hybridP()
        A = hybridP['A']

        if dolist:
            if hybridP is not None:
                params = [fid_T, fid_R]
                if option_cheap_height != 1:
                    params.append(fid_Pdep)
                result = self._simple_HDiag(params, fid_HGeopot)  # list of fids for HGeopot that we can actually produce
                field_3d = not any(['level' in r['generic'] for r in result])
                if field_3d:
                    result = [r for r in result if 'level' not in r]  # we limit to 3D results
                else:
                    if set([r['generic']['level'] for r in result if 'level' in r['generic']]) != set(range(1, len(A) + 1)):
                        result = []
                    else:
                        result = [r for r in result if 'level' in r['generic']]  # we limit to 2D results
            else:
                result = []
            return result
        else:
            Ps = self.readfield(fid_Ps, getdata=getdata)

            if getdata:
                if 'HGeopot' not in self._cache:
                    if Ps.spectral:
                        Ps.sp2gp()
                    if Ps.getdata().flatten()[0] < 100:
                        # Ps contains lnsp and not sp
                        # This is a workaround for a bug in FA
                        Ps.setdata(numpy.exp(Ps.getdata()))

                    # Pressure on mass levels
                    pi = hybridP2pressure(hybridP['vgeometry'], Ps.getdata(), option_vertical_mean,
                                          gridposition='mass').levels
                    pi = numpy.array(pi) * 100

                    field_3d = not any(['level' in fid for fid in fids])

                    # R, T
                    if field_3d:
                        R = self.readfield(fid_R).getdata()
                        T = self.readfield(fid_T)
                        if T.spectral:
                            T.sp2gp()
                        T = T.getdata()
                    else:
                        R = numpy.empty_like(pi)
                        T = numpy.empty_like(pi)
                        for k in range(1, len(A) + 1):
                            R[k - 1] = self.readfield(dict(fid_R, **dict(level=k))).getdata()
                            Tl = self.readfield(dict(fid_T, **dict(level=k)))
                            if Tl.spectral:
                                Tl.sp2gp()
                            T[k - 1] = Tl.getdata()
                            del Tl

                    # Pressure departure
                    if option_cheap_height != 1:
                        if field_3d:
                            Pdep = self.readfield(fid_Pdep)
                            if Pdep.spectral:
                                Pdep.sp2gp()
                            Pdep = Pdep.getdata()
                        else:
                            Pdep = numpy.empty_like(pi)
                            for k in range(1, len(A) + 1):
                                Pdepl = self.readfield(dict(fid_Pdep, **dict(level=k)))
                                if Pdepl.spectral:
                                    Pdepl.sp2gp()
                                Pdep[k - 1] = Pdepl.getdata()
                                del Pdepl
                    else:
                        Pdep = None

                    Phi_surf = self.readfield(fid_PhiSurf, getdata=getdata)
                    if Phi_surf.spectral:
                        Phi_surf.sp2gp()
                    self._cache['HGeopot'] = profiles.pressure2altitude(R, T,
                                                                        vertical_mean=option_vertical_mean,
                                                                        pi=pi,
                                                                        Phi_surf=Phi_surf.getdata(),
                                                                        Pdep=Pdep)
                HGeopotVal = self._cache['HGeopot']

            fieldset = FieldSet()
            for fid in fids:
                HGeopot = self.resource.readfield(hybridP['fid_for_vcoord'], getdata=False)
                geometry = HGeopot.geometry.deepcopy()
                if 'level' in fid:
                    geometry.vcoordinate.levels = [fid['level']]
                else:
                    geometry.vcoordinate.levels = list(range(1, len(A) + 1))
                myfid_HGeopot = fid.copy()
                myfid_HGeopot.update(fid_HGeopot)
                myfid_HGeopot = {self.format:myfid_HGeopot, 'generic':myfid_HGeopot}
                HGeopot = fpx.field(structure=geometry.structure,
                                    geometry=geometry,
                                    validity=Ps.validity,
                                    processtype=Ps.processtype,
                                    fid=myfid_HGeopot)
                if getdata:
                    if 'level' in fid:
                        HGeopot.setdata(HGeopotVal[fid['level'] - 1])
                    else:
                        HGeopot.setdata(HGeopotVal)
                fieldset.append(HGeopot)

            return fieldset
