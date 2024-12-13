#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for FA format.
"""

import numpy
import re
import sys

import footprints
from footprints import proxy as fpx
from bronx.meteo.conversion import q2R

from epygram import epygramError, util
from epygram.util import separation_line, write_formatted_fields
from epygram.base import FieldSet, FieldValidityList
from epygram.resources import FileResource
from epygram.geometries.VGeometry import (hybridP2pressure, hybridP2altitude,
                                          pressure2altitude)
from epygram.formats.FA import FA, get_generic_fid

__all__ = ['FA3d']

epylog = footprints.loggers.getLogger(__name__)


class FA3d(FileResource):
    """Class implementing a direct 3D access to a FA resource format."""

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['FA']),
                default='FA'),
            true3d=dict(
                values=set([True]),
                info="Open the FA file in a 3d mode")
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False

        super(FA3d, self).__init__(*args, **kwargs)

        self.resource = FA(*args, **kwargs)
        if not self.fmtdelayedopen:
            self.open()

        self._cache_find_re_in_list = {}
        self._3dnames = [(re.compile(r'P[0-9][0-9][0-9][0-9][0-9]'), 'P-----'),
                         (re.compile(r'S[0-9][0-9][0-9]'), 'S---'),
                         (re.compile(r'H[0-9][0-9][0-9][0-9][0-9]'), 'H-----'),
                         (re.compile(r'V[0-9][0-9][0-9]'), 'V---'),
                         (re.compile(r'T[0-9][0-9][0-9]'), 'T---'),
                         (re.compile(r'KT[0-9][0-9][0-9]'), 'KT---')]
        self._cont = {}

    def open(self, openmode=None):
        """
        Opens a FA file using the normal FA interface.

        :param openmode: optional, to open with a specific openmode, eventually
          different from the one specified at initialization.
        """
        super(FA3d, self).open(openmode=openmode)

        if self.openmode not in ('r', 'a'):
            raise epygramError("FA files with true3d option cannot be opened in 'w' mode")
        self.resource.open()
        self.isopen = True
        self.geometry = self.resource.geometry
        self.validity = self.resource.validity

    def close(self):
        """Closes the underlying FA file."""
        if self.isopen:
            self.resource.close()

################
# ABOUT FIELDS #
################
    def find_fields_in_resource(self, seed=None, fieldtype=[], generic=False):
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
        if isinstance(fieldtype, list):
            fieldtypeslist = fieldtype
        else:
            fieldtypeslist = [fieldtype]
        fieldslist = []
        if seed is None:
            tmplist = self.listfields()
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type(f) in fieldtypeslist:
                    fieldslist.append(f)
        elif isinstance(seed, str):
            h = (hash(seed), hash(tuple(self.listfields())))
            if h not in self._cache_find_re_in_list:
                self._cache_find_re_in_list[h] = util.find_re_in_list(seed, self.listfields())
            tmplist = self._cache_find_re_in_list[h]
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type(f) in fieldtypeslist:
                    fieldslist.append(f)
        elif isinstance(seed, list):
            tmplist = []
            for s in seed:
                h = (hash(s), hash(tuple(self.listfields())))
                if h not in self._cache_find_re_in_list:
                    self._cache_find_re_in_list[h] = util.find_re_in_list(s, self.listfields())
                tmplist += self._cache_find_re_in_list[h]
            for f in tmplist:
                if fieldtypeslist == [] or\
                   self._field_type(f) in fieldtypeslist:
                    fieldslist.append(f)
        if fieldslist == []:
            raise epygramError("no field matching: " + str(seed) +
                               " was found in resource " +
                               self.container.abspath)
        if generic:
            fieldslist = [(f, get_generic_fid(f)) for f in fieldslist]

        return fieldslist

    def listfields(self, **kwargs):
        """
        Returns a list containing the identifiers of all the fields of the
        resource.
        """
        return super(FA3d, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for FA3d.

        :param complete: - if True method returns a list of {'FA':FA_fid,
                           'generic':generic_fid}
                         - if False method return a list of FA_fid
        """
        fieldslist = []
        for fid in self.resource.listfields(complete=complete):
            if complete:
                fidFA = fid['FA']
            else:
                fidFA = fid
            fidFA3d = fidFA
            for c, f in self._3dnames:
                fidFA3d = c.sub(f, fidFA3d)
            if fidFA3d not in fieldslist:
                fieldslist.append(fidFA3d)
                self._cont[fidFA3d] = []
            self._cont[fidFA3d].append(fidFA)

        if complete:
            completefidList = []
            for fidFA3d, fidFA in self._cont.items():
                generic = get_generic_fid(fidFA[0])
                if len(fidFA) != 1:
                    del generic['level']
                completefidList.append({'FA':fidFA3d, 'generic':generic})
            return completefidList
        else:
            return fieldslist

    def split_UV(self, fieldseed):
        """
        Return two lists of fids corresponding respectively to U and V
        components of wind, given a *fieldseed*.
        """
        fids = self.find_fields_in_resource(fieldseed + '*')
        if fieldseed.startswith('S'):
            Ufid = [f for f in fids if 'WIND.U.PHYS' in f]
            Vfid = [f for f in fids if 'WIND.V.PHYS' in f]
        elif fieldseed[0] in ('P', 'H', 'V'):
            Ufid = [f for f in fids if 'VENT_ZONAL' in f]
            Vfid = [f for f in fids if 'VENT_MERID' in f]
        elif fieldseed.startswith('CLS') and 'VENT' in fids[0]:
            if fieldseed.startswith('CLSVENTNEUTRE'):
                Ufid = [f for f in fids if 'CLSVENTNEUTRE.U' in f]
                Vfid = [f for f in fids if 'CLSVENTNEUTRE.V' in f]
            else:
                Ufid = [f for f in fids if 'VENT.ZONAL' in f]
                Vfid = [f for f in fids if 'VENT.MERIDIEN' in f]
        elif fieldseed.startswith('CLS') and 'RAF' in fids[0]:
            Ufid = [f for f in fids if 'CLSU' in f]
            Vfid = [f for f in fids if 'CLSV' in f]
        else:
            raise NotImplementedError("split_UV: field syntax='" + fieldseed + "'")

        return (sorted(Ufid), sorted(Vfid))

    def sortfields(self):
        """
        Returns a sorted list of fields with regards to their name and nature,
        as a dict of lists.
        """
        re_3D = re.compile(r'(?P<prefix>[A-Z])(?P<level>\d+)(?P<param>[A-Z]+.*)')
        list3D = []
        list2D = []
        # final lists
        list3Dsp = []
        list3Dgp = []
        list2Dsp = []
        list2Dgp = []
        listMisc = []

        params3D = {}
        for f in self.listfields():
            info = FA.gribdef.FA2GRIB(self._cont[f][0])
            # separate H2D from Misc
            if self._field_type(f) == 'H2D':
                # separate 3D from 2D
                if info['typeOfFirstFixedSurface'] in (119, 100, 103, 109, 20):
                    re_ok = re_3D.match(f)
                    if re_ok:
                        list3D.append(f)
                        param = re_ok.group('prefix') + re_ok.group('param')
                        if param not in params3D:
                            params3D[param] = []
                        params3D[param].append(f)
                    else:
                        list2D.append(f)
                else:
                    list2D.append(f)
            elif self._field_type(f) == '3D':
                list3D.append(f)
                params3D[f] = [f]
            else:
                listMisc.append(f)
        # separate gp/sp
        for f in list2D:
            if self.resource.fieldencoding(self._cont[f][0])['spectral']:
                list2Dsp.append(f)
            else:
                list2Dgp.append(f)
        # sort 2D
        list2Dsp.sort()
        list2Dgp.sort()
        listMisc.sort()
        # sort 3D
        for p in sorted(params3D.keys()):
            if self.resource.fieldencoding(self._cont[params3D[p][0]][0])['spectral']:
                list3Dsp.extend(sorted(params3D[p]))
            else:
                list3Dgp.extend(sorted(params3D[p]))
        outlists = {'3D spectral fields':list3Dsp,
                    '3D gridpoint fields':list3Dgp,
                    '2D spectral fields':list2Dsp,
                    '2D gridpoint fields':list2Dgp,
                    'Misc-fields':listMisc}

        return outlists

    def _field_type(self, fieldname):
        """Return type of the field."""
        if fieldname in self.resource.listfields():
            return self.resource.field_type(fieldname)
        else:
            return '3D'

    @FileResource._openbeforedelayed
    def readfield(self, fieldname,
                  getdata=True):
        """
        Reads one field, given its FA name, and returns a Field instance.
        Interface to Fortran routines from 'ifsaux'.

        :param fieldname: FA3d fieldname
        :param getdata: if *False*, only metadata are read, the field do not
          contain data.
        """
        if self.openmode == 'w':
            raise epygramError("cannot read fields in resource if with" +
                               " openmode == 'w'.")
        assert fieldname in self.listfields(), ' '.join(["field",
                                                         str(fieldname),
                                                         "not found in resource."])
        # Get field info
        ftype = self._field_type(fieldname)
        if ftype in ['H2D', 'Misc']:
            field = self.resource.readfield(fieldname, getdata=getdata)
        else:
            fieldset = FieldSet()
            for fid in self._cont[fieldname]:
                fieldset.append(self.resource.readfield(fid))
            assert all([len(f.geometry.vcoordinate.levels) == 1 for f in fieldset]), "Internal error"
            levels = [f.geometry.vcoordinate.levels[0] for f in fieldset]

            vcoordinate = fieldset[0].geometry.vcoordinate
            del vcoordinate.levels[:]
            vcoordinate.levels.extend(sorted(levels))

            kwargs_field = fieldset[0]._attributes
            kwargs_field['structure'] = '3D'
            kwargs_field['geometry'] = fieldset[0].geometry.deepcopy()

            field = fpx.field(**kwargs_field)

            field.fid['FA'] = fieldname
            if 'generic' in field.fid:
                del field.fid['generic']['level']

            if getdata:
                indexes = sorted(range(len(levels)), key=lambda k: levels[k])
                field.setdata(numpy.array([fieldset[i].getdata() for i in indexes]))

        return field

    def readfields(self, requestedfields=None, getdata=True):
        """
        Returns a :class:`epygram.base.FieldSet` containing requested fields
        read in the resource.

        :param requestedfields: might be:\n
          - a regular expression (e.g. 'S\*WIND.[U,V].PHYS')
          - a list of FA fields identifiers with regular expressions (e.g.
            ['SURFTEMPERATURE', 'S0[10-20]WIND.?.PHYS'])
          - if not specified, interpretated as all fields that will be found in
            resource
        :param getdata: optional, if *False*, only metadata are read, the fields
          do not contain data. Default is *True*.
        """
        requestedfields = self.find_fields_in_resource(requestedfields)
        if requestedfields == []:
            raise epygramError("unable to find requested fields in resource.")

        return super(FA3d, self).readfields(requestedfields, getdata)

    def writefield(self, field, compression=None):
        raise epygramError("FA files with true3d option cannot be opened in 'w' mode")

    def writefields(self, fieldset, compression=None):
        raise epygramError("FA files with true3d option cannot be opened in 'w' mode")

    @FileResource._openbeforedelayed
    def extractprofile(self, fieldname, lon=None, lat=None,
                       geometry=None,
                       vertical_coordinate=None,
                       interpolation='nearest',
                       cheap_height=True,
                       external_distance=None):
        """
        Extracts a vertical profile from the FA3d resource, given its fieldname
        and the geographic location (*lon*/*lat*) of the profile.

        :param fieldname: name of the field.
        :param lon: the longitude of the desired point.
        :param lat: the latitude of the desired point.
          If both None, extract a horizontally-averaged profile.
        :param geometry: can replace *lon*/*lat*, geometry on which to extract
          data. If None, it is built from *lon*/*lat*.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V1DField, as number of GRIB2 norm:
          http://apps.ecmwf.int/codes/grib/format/grib2/ctables/4/5,
          (cf. `epygram.geometries.vertical_coordinates` possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile at requested lon/lat from the fields grid:\n
          - if 'nearest' (default), extracts profile at the horizontal nearest neighboring gridpoint;
          - if 'linear', computes profile with horizontal linear spline interpolation;
          - if 'cubic', computes profile with horizontal cubic spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
        :param external_distance: can be a dict containing the target point value
          and an external field on the same grid as self, to which the distance
          is computed within the 4 horizontally nearest points; e.g.
          {'target_value':4810, 'external_field':an_H2DField_with_same_geometry}.
          If so, the nearest point is selected with
          distance = abs(target_value - external_field.data)
        """
        if geometry is None:
            if None in [lon, lat]:
                if not (lon is None and lat is None):
                    raise ValueError('*lon* and *lat* arguments must be ' +
                                     'both None or both not None')
                # mean profile, vertical_coordinate is forgotten
                field3d = self.readfield(fieldname)
                if field3d.spectral:
                    field3d.sp2gp()
                profG = self.geometry.make_profile_geometry(lon, lat)
                profG.vcoordinate = field3d.geometry.vcoordinate.deepcopy()
                profile = fpx.field(fid=field3d.fid,
                                    structure='V1D',
                                    validity=self.validity.deepcopy(),
                                    geometry=profG,
                                    comment='horizontally-averaged profile')
                data3d = field3d.getdata(d4=True)
                data1d = [data3d[:, i, :, :].mean() for i in range(len(profG.vcoordinate.levels))]
                profile.setdata(data1d)
                return profile
            if self.geometry is None:
                self._read_geometry()
            pointG = self.geometry.make_profile_geometry(lon, lat)
        else:
            if lon is not None or lat is not None:
                raise epygramError("You cannot provide lon or lat when geometry is given")
            if geometry.structure != "V1D":
                raise epygramError("geometry must be a V1D")
            pointG = geometry

        profile = self.extract_subdomain(fieldname, pointG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         external_distance=external_distance,
                                         cheap_height=cheap_height)

        return profile

    @FileResource._openbeforedelayed
    def extractsection(self, fieldname, end1=None, end2=None,
                       geometry=None,
                       points_number=None,
                       resolution=None,
                       vertical_coordinate=None,
                       interpolation='linear',
                       cheap_height=True,
                       global_shift_center=None):
        """
        Extracts a vertical section from the FA resource, given its fieldoname
        and the geographic (lon/lat) coordinates of its ends.
        The section is returned as a V2DField.

        :param fieldname: name of the field.
        :param end1: must be a tuple (lon, lat).
        :param end2: must be a tuple (lon, lat).
        :param geometry: can replace end1/end2, geometry on which to extract
          data. If None, defaults to
          linearily spaced positions computed from *points_number*.
        :param points_number: defines the total number of horizontal points of the
          section (including ends). If None, defaults to a number computed from
          the *ends* and the *resolution*.
        :param resolution: defines the horizontal resolution to be given to the
          field. If None, defaults to the horizontal resolution of the field.
        :param vertical_coordinate: defines the requested vertical coordinate of
          the V2DField aka typeOfFirstFixedSurface in GRIB2, (cf.
          `epygram.geometries.vertical_coordinates` possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
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

        section = self.extract_subdomain(fieldname, sectionG,
                                         interpolation=interpolation,
                                         vertical_coordinate=vertical_coordinate,
                                         cheap_height=cheap_height,
                                         global_shift_center=global_shift_center)

        return section

    @FileResource._openbeforedelayed
    def extract_subdomain(self, fieldname, geometry, vertical_coordinate=None,
                          interpolation='linear', cheap_height=True,
                          external_distance=None,
                          global_shift_center=None):
        """
        Extracts a subdomain from the FA resource, given its fieldname
        and the geometry to use.

        :param fieldname: name of the field.
        :param geometry: is the geometry on which extract data.
                         None to keep the geometry untouched.
        :param vertical_coordinate: defines the requested vertical coordinate of the
          V2DField (cf. `epygram.geometries.vertical_coordinates`
          possible values).
        :param interpolation: defines the interpolation function used to compute
          the profile points locations from the fields grid: \n
          - if 'nearest', each horizontal point of the section is
            taken as the horizontal nearest neighboring gridpoint;
          - if 'linear' (default), each horizontal point of the section is
            computed with linear spline interpolation;
          - if 'cubic', each horizontal point of the section is
            computed with linear spline interpolation.
        :param cheap_height: if True and *vertical_coordinate* among
          ('altitude', 'height'), the computation of heights is done without
          taking hydrometeors into account (in R computation) nor NH Pressure
          departure (Non-Hydrostatic data). Computation therefore faster.
        :param global_shift_center: for global lon/lat grids, shift the center by the
            requested angle (in degrees). Enables a [0,360] grid
            to be shifted to a [-180,180] grid, for instance (with -180 argument).
        """
        field3d = self.readfield(fieldname)
        if field3d.spectral:
            field3d.sp2gp()
        if global_shift_center is not None:
            field3d.global_shift_center(global_shift_center)
            
        if geometry is None or geometry == field3d.geometry:
            subdomain = field3d
            geometry = field3d.geometry
        else:
            subdomain = field3d.extract_subdomain(geometry,
                                                  interpolation=interpolation)

        # preparation for vertical coords conversion
        if vertical_coordinate not in (None, subdomain.geometry.vcoordinate.typeoffirstfixedsurface):
            # choose vertical_mean with regards to H/NH
            if 'S---PRESS.DEPART' in self.listfields():
                vertical_mean = 'geometric'
            else:
                vertical_mean = 'arithmetic'
            # surface pressure (hybridP => P,A,H)
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
               vertical_coordinate in (100, 102, 103):
                h_shape = geometry.get_datashape(force_dimZ=1)
                Psurf = self.readfield('SURFPRESSION')
                if Psurf.spectral:
                    Psurf.sp2gp()
                ps_transect = numpy.exp(Psurf.getvalue_ll(*geometry.get_lonlat_grid(),
                                                          interpolation=interpolation,
                                                          one=False,
                                                          external_distance=external_distance))
                ps_transect = ps_transect.reshape(h_shape)
                del Psurf
            # P => H necessary profiles
            if vertical_coordinate in (102, 103):
                prefix = fieldname[:fieldname.rfind('-') + 1]
                side_profiles = {'t':prefix + 'TEMPERATURE',
                                 'q':prefix + 'HUMI.SPECIFI',
                                 'pdep':prefix + 'PRESS.DEPART',
                                 'ql':prefix + 'CLOUD_WATER',
                                 'qi':prefix + 'ICE_CRYSTAL',
                                 'qs':prefix + 'SNOW',
                                 'qr':prefix + 'RAIN',
                                 'qg':prefix + 'GRAUPEL'}
                for p in sorted(side_profiles.keys(), reverse=True):  # reverse to begin by t
                    try:
                        # try to extract profiles for each lon/lat and each parameter
                        if fieldname == side_profiles[p]:
                            # already extracted as requested profile
                            side_profiles[p] = subdomain
                        else:
                            if cheap_height and p not in ('t', 'q'):
                                raise epygramError()  # to go through "except" instructions below
                            side_profiles[p] = self.extract_subdomain(side_profiles[p], geometry,
                                                                      interpolation=interpolation,
                                                                      external_distance=external_distance)
                        side_profiles[p] = side_profiles[p].getdata()
                    except epygramError:
                        # fields not present in file
                        if p in ('t', 'q'):
                            raise epygramError("Temperature and Specific" +
                                               " Humidity must be in" +
                                               " resource.")
                        else:
                            side_profiles[p] = numpy.zeros(side_profiles['t'].shape)
                R = q2R(*[side_profiles[p] for p in
                          ['q', 'ql', 'qi', 'qr', 'qs', 'qg']])
            if vertical_coordinate == 102:
                try:
                    geopotential = self.readfield('SPECSURFGEOPOTEN')
                except epygramError:
                    geopotential = self.readfield('SURFGEOPOTENTIEL')
                else:
                    geopotential.sp2gp()
                surface_geopotential = geopotential.getvalue_ll(*geometry.get_lonlat_grid(),
                                                                interpolation=interpolation,
                                                                one=False,
                                                                external_distance=external_distance)
                surface_geopotential = surface_geopotential.reshape(h_shape)
                del geopotential
            else:
                surface_geopotential = None

            # effective vertical coords conversion
            if subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
               vertical_coordinate == 100:
                subdomain.geometry.vcoordinate = hybridP2pressure(subdomain.geometry.vcoordinate,
                                                                  ps_transect,
                                                                  vertical_mean,
                                                                  gridposition='mass')
            elif subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 119 and \
                 vertical_coordinate in (102, 103):
                subdomain.geometry.vcoordinate = hybridP2altitude(subdomain.geometry.vcoordinate,
                                                                  R,
                                                                  side_profiles['t'],
                                                                  ps_transect,
                                                                  vertical_mean,
                                                                  Pdep=side_profiles['pdep'],
                                                                  Phi_surf=surface_geopotential)
            elif subdomain.geometry.vcoordinate.typeoffirstfixedsurface == 100 and \
                 vertical_coordinate in (102, 103):
                subdomain.geometry.vcoordinate = pressure2altitude(subdomain.geometry.vcoordinate,
                                                                   R,
                                                                   side_profiles['t'],
                                                                   vertical_mean,
                                                                   Pdep=side_profiles['pdep'],
                                                                   Phi_surf=surface_geopotential)
            else:
                raise NotImplementedError("this vertical coordinate" +
                                          " conversion.")

        return subdomain

###########
# pre-app #
###########

    @FileResource._openbeforedelayed
    def what(self,
             out=sys.stdout,
             details=None,
             sortfields=False,
             **_):
        """
        Writes in file a summary of the contents of the FA.

        :param out: the output open file-like object.
        :param details: 'spectral' if spectralness of fields is requested;
                        'compression' if information about fields compression
                        is requested.
        :param sortfields: **True** if the fields have to be sorted by type.
        """
        for f in self.listfields():
            if self._field_type(f) == 'H2D':
                first_H2DField = f
                break
        if len(self.listfields()) == 0:
            raise epygramError('empty Resource.')
        firstfield = self.readfield(first_H2DField, getdata=False)
        if not firstfield.spectral and self.spectral_geometry is not None:
            firstfield._attributes['spectral_geometry'] = self.spectral_geometry

        listoffields = self.listfields()
        if sortfields:
            sortedfields = self.sortfields()

        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        out.write("### IDENTIFIER (CDIDEN): " + self.resource.cdiden + "\n")
        out.write("\n")

        FieldValidityList(self.validity).what(out)
        firstfield.what(out,
                        validity=False,
                        vertical_geometry=False,
                        arpifs_var_names=True,
                        fid=False)

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
        if details is None:
            write_formatted_fields(out, "Field name")
        elif details == 'spectral':
            write_formatted_fields(out, "Field name", "Spectral")
        elif details == 'compression':
            params = ['KNGRIB', 'KNBITS', 'KSTRON', 'KPUILA']
            width_cp = 8
            compressionline = ""
            for p in params:
                compressionline += '{:^{width}}'.format(p, width=width_cp)
            write_formatted_fields(out, "Field name", "Spectral",
                                   compressionline)
        out.write(separation_line)
        for f in listoffields:
            if details is not None and self._field_type(f) in ['H2D', '3D']:
                encoding = self.resource.fieldencoding(self._cont[f][0])
                if details == 'spectral':
                    write_formatted_fields(out, f, encoding['spectral'])
                elif details == 'compression':
                    compressionline = ""
                    for p in params:
                        try:
                            compressionline += '{:^{width}}'.format(str(encoding[p]), width=width_cp)
                        except KeyError:
                            compressionline += '{:^{width}}'.format('-', width=width_cp)
                    write_formatted_fields(out, f, encoding['spectral'],
                                           compression=compressionline)
            else:
                write_formatted_fields(out, f)
        out.write(separation_line)
