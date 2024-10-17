#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for DDH-in-LFA format.
"""

import math
import numpy
import datetime
import sys

import footprints
from footprints import FPDict, FPList
from bronx.meteo.constants import g0

from  falfilfa4py import LFA as LFA4py

from epygram import config, epygramError
from epygram.util import Angle, write_formatted, separation_line
from epygram.base import FieldValidity, FieldValidityList, FieldSet
from .LFA import LFA
from epygram.geometries import UnstructuredGeometry, VGeometry
from epygram.fields import V1DField, PointField, MiscField
from epygram.resources import FileResource
from epygram import profiles

__all__ = ['DDHLFA']

epylog = footprints.loggers.getLogger(__name__)


class DDHLFA(LFA):
    """
    Class implementing all specificities for DDHLFA resource format.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['DDHLFA']),
                default='DDHLFA'),
            validity=dict(
                type=FieldValidityList,
                access='rwx',
                optional=True,
                info="Describes the temporal validity of the resource."),
            domains=dict(
                type=FPDict,
                optional=True,
                info="Describes the domains covered by the resource."),
            levels=dict(
                type=FPDict,
                optional=True,
                info="Number of levels for variables/tendencies ('VT') and\
                      fluxes ('F')."),
            xpid=dict(
                optional=True,
                info="Experiment identifier.")
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        super(LFA, self).__init__(*args, **kwargs)
        if not LFA4py.wlfatest(self.container.abspath):
            raise IOError("This resource is not a LFA one.")
        if self.openmode in ('r', 'a'):
            guess = LFA(filename=self.container.abspath, openmode=self.openmode)
            if not all([v in guess.listfields() for v in ['ECHEANCE', 'DATE', 'DOCFICHIER',
                                                          'VPP0', 'VPP1']]):
                raise IOError("this resource is not a DDHLFA one.")
            guess.close()
        if self.openmode != 'r':
            raise NotImplementedError("openmode != 'r' for DDHLFA.")
        elif self.openmode == 'r':
            if self.domains is not None or self.validity is not None:
                epylog.warning(self.container.abspath +
                               ": DDHLFA.__init__(): domains/validity" +
                               " argument will be ignored with" +
                               " this openmode ('r').")
        if not self.fmtdelayedopen:
            self.open()

    def open(self, openmode=None):
        """
        Opens the DDHLFA in Fortran sense, and initializes domains, validity
        and vertical geometry.

        :param openmode: optional, to open with a specific openmode, eventually
                         different from the one specified at initialization.
        """
        super(DDHLFA, self).open(openmode=openmode)

        if self.openmode == 'r':
            # read info
            self._attributes['xpid'] = super(DDHLFA, self).readfield('INDICE EXPERIENCE').getdata()[0]
            self._read_geometry()
            self._read_validity()
        else:
            raise NotImplementedError("openmode != 'r' for DDHLFA.")

################
# ABOUT FIELDS #
################

    def find_fields_in_resource(self,
                                seed=None,
                                fieldtype=None,
                                generic=False,
                                **_):
        """
        Returns a list of the fields from resource whose name match the given
        *seed*.

        :param seed: might be a regular expression, a list of regular expressions
          or *None*. If *None* (default), returns the list of all fields in
          resource.
        :param fieldtype: optional, among ('V1D', 'Point', 'Misc'). If provided,
          filters out the fields not of the given type.
        """
        fieldslist = super(DDHLFA, self).find_fields_in_resource(seed=seed)
        if generic:
            raise NotImplementedError("not yet.")
        if fieldtype is not None:
            fs = FieldSet()
            for f in fieldslist:
                fld = self.readfield(f, getdata=False)
                if isinstance(fld, FieldSet):
                    fs.append(fld[0])
                else:
                    fs.append(fld)
            if fieldtype == 'V1D':
                fieldslist = [f.fid[self.format] for f in fs
                              if isinstance(f, V1DField)]
            elif fieldtype == 'Point':
                fieldslist = [f.fid[self.format] for f in fs
                              if isinstance(f, PointField)]
            elif fieldtype == 'Misc':
                fieldslist = [f.fid[self.format] for f in fs
                              if isinstance(f, MiscField)]
            else:
                raise epygramError("unknown fieldtype for DDHLFA: " + fieldtype)

        return fieldslist

    @FileResource._openbeforedelayed
    def readfield(self, fieldname,
                  getdata=True,
                  footprints_proxy_as_builder=config.footprints_proxy_as_builder):
        """
        Reads a field in resource.

        - Documentation fields ('INDICE EXPERIENCE', 'DATE', 'DOCFICHIER',
          'ECHEANCE', 'DOCDnnn') are returned as :class:`epygram.formats.LFA`
          returns.
        - Profile/surface fields are returned as a
          :class:`epygram.base.FieldSet` of 1D/Point fields, one for each
          domain.

        :param getdata: if False, do not read data but only metadata
        :param footprints_proxy_as_builder: if *True*, uses footprints.proxy
               to build fields and geometry.
        """
        if footprints_proxy_as_builder:
            field_builder = footprints.proxy.field
        else:
            if fieldname[0] in ('S', 'G'):
                field_builder = PointField
            else:
                field_builder = V1DField
        field_from_LFA = super(DDHLFA, self).readfield(fieldname,
                                                       getdata=getdata)
        if fieldname in ('INDICE EXPERIENCE', 'DATE', 'DOCFICHIER', 'ECHEANCE')\
           or fieldname[0:4] == 'DOCD':
            # documentation
            field_from_LFA.setfid({self.format:field_from_LFA.fid['LFA']})
            toreturn = field_from_LFA
        else:
            # true field
            validity = FieldValidity(basis=self.validity.getbasis())
            # distinction between initial, final and integrated fields
            if fieldname[0:2] == 'SV':
                validity.set(date_time=self.validity.get())
            elif fieldname[0:2] == 'SF':
                validity.set(date_time=self.validity.get(),
                             cumulativeduration=self.validity.cumulativeduration())
            elif fieldname[0] in ('V', 'S'):
                if fieldname[-1] == '0':
                    validity.set(date_time=self.validity.getbasis())
                elif fieldname[-1] == '1':
                    validity.set(date_time=self.validity.get())
            elif fieldname[0] in ('T', 'F', 'G') or fieldname == 'PPP':
                validity.set(date_time=self.validity.get(),
                             cumulativeduration=self.validity.cumulativeduration())
            toreturn = FieldSet()
            if fieldname[0] in ('S', 'G'):
                # surface fields
                for d in range(self.domains['number']):
                    if getdata:
                        value = field_from_LFA.getdata()[d]
                    domain = self.domains['geometry'][d]
                    pgeometry = UnstructuredGeometry(name='DDH:' + domain['type'],
                                                     grid={'DDH_domain':domain},
                                                     vcoordinate=VGeometry(
                                                         typeoffirstfixedsurface=1,
                                                         levels=FPList([1])),
                                                     dimensions={'X':1, 'Y':1})
                    field = field_builder(structure='Point',
                                          fid={self.format:fieldname},
                                          geometry=pgeometry,
                                          validity=validity,
                                          processtype=self.processtype)
                    if getdata:
                        field.setdata(value)
                    toreturn.append(field)
            else:
                # profile fields
                if not getdata:
                    fieldlevels = (LFA4py.wlfacas(self._unit, fieldname)[1] //
                                   self.domains['number'])
                else:
                    fieldlevels = (len(field_from_LFA.getdata()) //
                                   self.domains['number'])
                if fieldlevels == self.levels['VT']:
                    position_on_grid = 'mass'
                    # gridposition = 'mass'
                    # gridsize = self.levels['VT']
                elif fieldlevels == self.levels['F']:
                    position_on_grid = 'flux'
                    # gridposition = 'flux'
                    # gridsize = self.levels['F']
                for d in range(self.domains['number']):
                    domain = self.domains['geometry'][d]
                    if fieldlevels == self.levels['VT']:
                        if validity.term(fmt='IntHours') == 0:
                            pressure_vertical_grid = self.domains['vertical_grid'][d]['masslevels_pressure_init']
                        else:
                            pressure_vertical_grid = self.domains['vertical_grid'][d]['masslevels_pressure_term']
                    elif fieldlevels == self.levels['F']:
                        if validity.term(fmt='IntHours') == 0:
                            pressure_vertical_grid = self.domains['vertical_grid'][d]['fluxlevels_pressure_init']
                        else:
                            pressure_vertical_grid = self.domains['vertical_grid'][d]['fluxlevels_pressure_term']
                    vcoordinate = VGeometry(typeoffirstfixedsurface=100,
                                            levels=FPList(pressure_vertical_grid / 100.),
                                            # grid={'gridposition':gridposition,  # TODO: ?
                                            #       'gridlevels':pressure_vertical_grid},  # TODO: ?
                                            position_on_grid=position_on_grid)
                    vgeometry = UnstructuredGeometry(name='DDH:' + domain['type'],
                                                     grid={'DDH_domain':domain},
                                                     vcoordinate=vcoordinate,
                                                     dimensions={'X':1, 'Y':1})
                    if getdata:
                        profile = field_from_LFA.getdata()[d * fieldlevels:(d + 1) * fieldlevels]
                    field = field_builder(structure='V1D',
                                          fid={self.format:fieldname, 'generic':FPDict()},
                                          geometry=vgeometry,
                                          validity=validity,
                                          processtype=self.processtype)
                    if getdata:
                        field.setdata(profile)
                    toreturn.append(field)
        return toreturn

    def readfields(self, requestedfields, **kwargs):
        """
        Inactivation of readfields because readfield already returns a FieldSet.
        """
        raise epygramError("readfields() method disabled for DDHLFA.")

#################
# DDH interface #
#################

    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout, sortfields=False, **_):
        """
        Writes in file a summary of the contents of the DDHLFA.

        :param out: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        :param sortfields: **True** if the fields have to be sorted by type.
        """
        domains = self.domains
        mass_vertical_levels = len(self.readfield('VPP0')[0].geometry.vcoordinate.levels)
        listoffields = self.listfields()
        validity = self.validity

        # Write out
        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")
        out.write("### IDENTIFIER (INDICE EXPERIENCE): " + self.xpid + "\n")
        out.write("\n")

        out.write("################\n")
        out.write("### VALIDITY ###\n")
        out.write("################\n")
        write_formatted(out, "Validity", validity.get())
        write_formatted(out, "Basis", validity.getbasis())
        write_formatted(out, "Term", validity.term())
        write_formatted(out, "Duration for cumulative quantities",
                        validity.cumulativeduration())
        write_formatted(out, "Kind of producting process", self.processtype)
        out.write(separation_line)
        out.write("\n")

        out.write("###############\n")
        out.write("### DOMAINS ###\n")
        out.write("###############\n")
        out.write("Domain(s): " + str(domains['number']) + " " +
                  domains['geometry'][0]['type'] + "\n")
        for d in range(domains['number']):
            geom = domains['geometry'][d]
            if geom['type'] != 'globe':
                out.write("# " + str(d + 1) + ":\n")
            if geom['type'] == 'ij_point':
                write_formatted(out, "Longitude", geom['lon'].get('degrees'))
                write_formatted(out, "Latitude", geom['lat'].get('degrees'))
                write_formatted(out, "Point index (X)", geom['jlon'])
                write_formatted(out, "Point index (Y)", geom['jgl'])
            elif geom['type'] == 'point':
                write_formatted(out, "Longitude", geom['lon'].get('degrees'))
                write_formatted(out, "Latitude", geom['lat'].get('degrees'))
                write_formatted(out, "Point index (X)", geom['jlon'])
                write_formatted(out, "Point index (Y)", geom['jgl'])
                write_formatted(out, "User Longitude",
                                geom['user_lon'].get('degrees'))
                write_formatted(out, "User Latitude",
                                geom['user_lat'].get('degrees'))
            elif geom['type'] in ('quadrilateral', 'rectangle'):
                write_formatted(out, "Longitude 1", geom['lon1'].get('degrees'))
                write_formatted(out, "Latitude 1", geom['lat1'].get('degrees'))
                write_formatted(out, "Longitude 2", geom['lon2'].get('degrees'))
                write_formatted(out, "Latitude 2", geom['lat2'].get('degrees'))
                write_formatted(out, "Longitude 3", geom['lon3'].get('degrees'))
                write_formatted(out, "Latitude 3", geom['lat3'].get('degrees'))
                write_formatted(out, "Longitude 4", geom['lon4'].get('degrees'))
                write_formatted(out, "Latitude 4", geom['lat4'].get('degrees'))
            elif geom['type'] == 'zonal_bands':
                write_formatted(out, "Center Latitude",
                                geom['lat'].get('degrees'))
        out.write(separation_line)
        out.write("\n")

        out.write("#########################\n")
        out.write("### VERTICAL GEOMETRY ###\n")
        out.write("#########################\n")
        write_formatted(out, "Number of vertical levels for variables",
                        mass_vertical_levels)
        write_formatted(out, "Number of vertical levels for fluxes",
                        mass_vertical_levels + 1)
        out.write(separation_line)
        out.write("\n")

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        numfields = len(listoffields)
        out.write("Number: " + str(numfields) + "\n")
        out.write(separation_line)
        if sortfields:
            out.write("Documentation:\n")
            out.write("-------------\n")
            docfields = [f for f in listoffields if
                         f in ('INDICE EXPERIENCE', 'DATE', 'DOCFICHIER',
                               'ECHEANCE') or f[0:4] == 'DOCD']
            for f in docfields:
                out.write(f + "\n")
            out.write("-----------------\n")
            out.write("Vertical profiles:\n")
            out.write("-----------------\n")
            profilefields = [f for f in listoffields
                             if f[0] not in ('S', 'G') and f not in docfields]
            for f in profilefields:
                out.write(f + "\n")
            out.write("--------------\n")
            out.write("Surface fields:\n")
            out.write("--------------\n")
            surfacefields = [f for f in listoffields
                             if f not in docfields and f not in profilefields]
            for f in surfacefields:
                out.write(f + "\n")
        else:
            for f in listoffields:
                out.write(f + "\n")
        out.write(separation_line)

    @FileResource._openbeforedelayed
    def _read_geometry(self):
        """
        Reads DOCFICHIER and DOCDnnn fields, to fill-in domains and levels
        attributes. Cf. J.M. Piriou's documentation about DDH.
        """
        domains = {}
        DOCFICHIER = self.readfield('DOCFICHIER')
        vpp_init = super(DDHLFA, self).readfield('VPP0')
        vpp_term = super(DDHLFA, self).readfield('VPP1')

        NFLEV = DOCFICHIER.getdata()[5]
        self._attributes['levels'] = {'VT':NFLEV, 'F':NFLEV + 1}

        if DOCFICHIER.getdata()[0] == 1:
            domains['type'] = 'limited_area_domains'
        elif DOCFICHIER.getdata()[0] == 5:
            domains['type'] = 'global_domain'
        elif DOCFICHIER.getdata()[0] == 6:
            domains['type'] = 'zonal_bands'
        else:
            raise NotImplementedError("DOCFICHIER[0] not among (1,5,6).")
        domains['number'] = DOCFICHIER.getdata()[14]
        domains['geometry'] = []
        domains['vertical_grid'] = []
        if self.xpid == 'ARPE':
            vertical_mean = 'arithmetic'
        else:
            vertical_mean = 'geometric'
        for d in range(1, domains['number'] + 1):
            DOCD = self.readfield('DOCD' + '{:0>{width}}'.format(str(d),
                                                                 width=3))
            DOCD_data = DOCD.getdata()
            geom = {}
            if DOCD_data[10] == 1.:
                geom['type'] = 'ij_point'
                geom['lon'] = Angle(DOCD_data[2], 'radians')
                geom['lat'] = Angle(math.asin(DOCD_data[3]), 'radians')
                geom['jlon'] = int(DOCD_data[4])
                geom['jlgl'] = int(DOCD_data[5])
            elif DOCD_data[10] == 4.:
                geom['type'] = 'point'
                geom['lon'] = Angle(DOCD_data[2], 'radians')
                geom['lat'] = Angle(math.asin(DOCD_data[3]), 'radians')
                geom['jlon'] = int(DOCD_data[4])
                geom['jgl'] = int(DOCD_data[5])
                geom['user_lon'] = Angle(DOCD_data[6], 'radians')
                geom['user_lat'] = Angle(math.asin(DOCD_data[7]), 'radians')
            elif DOCD_data[10] == 2.:
                geom['type'] = 'quadrilateral'
                geom['lon1'] = Angle(DOCD_data[2], 'radians')
                geom['lat1'] = Angle(math.asin(DOCD_data[3]), 'radians')
                geom['lon2'] = Angle(DOCD_data[4], 'radians')
                geom['lat2'] = Angle(math.asin(DOCD_data[5]), 'radians')
                geom['lon3'] = Angle(DOCD_data[6], 'radians')
                geom['lat3'] = Angle(math.asin(DOCD_data[7]), 'radians')
                geom['lon4'] = Angle(DOCD_data[8], 'radians')
                geom['lat4'] = Angle(math.asin(DOCD_data[9]), 'radians')
            elif DOCD_data[10] == 3.:
                geom['type'] = 'rectangle'
                geom['lon1'] = Angle(DOCD_data[2], 'radians')
                geom['lat1'] = Angle(math.asin(DOCD_data[3]), 'radians')
                geom['lon2'] = Angle(DOCD_data[4], 'radians')
                geom['lat2'] = Angle(math.asin(DOCD_data[5]), 'radians')
                geom['lon3'] = Angle(DOCD_data[6], 'radians')
                geom['lat3'] = Angle(math.asin(DOCD_data[7]), 'radians')
                geom['lon4'] = Angle(DOCD_data[8], 'radians')
                geom['lat4'] = Angle(math.asin(DOCD_data[9]), 'radians')
            elif DOCD_data[10] == 5.:
                geom['type'] = 'globe'
            elif DOCD_data[10] == 6.:
                geom['type'] = 'zonal_bands'
                geom['lat'] = Angle(math.asin(DOCD_data[3]), 'radians')
                geom['band'] = DOCD_data[1]
                geom['bands_number'] = DOCD_data[2]
            else:
                raise ValueError('unknown value: ' + str(DOCD_data[10]) +
                                 ' for DOCD' + '{:0>{width}}'.format(str(d),
                                                                     width=3))
            # vertical grid
            dm1 = d - 1
            vgrid = {}
            fluxlevels_pressure_init = vpp_init.getdata()[dm1 * self.levels['VT']:(dm1 + 1) * self.levels['VT']
                                                          ] * g0
            fluxlevels_pressure_init = [sum(fluxlevels_pressure_init[0:i])
                                        for i in range(1, len(fluxlevels_pressure_init) + 1)]
            vgrid['fluxlevels_pressure_init'] = numpy.array([0.] + fluxlevels_pressure_init)
            vgrid['masslevels_pressure_init'] = profiles.flux2masspressures(fluxlevels_pressure_init,
                                                                            vertical_mean)
            fluxlevels_pressure_term = vpp_term.getdata()[dm1 * self.levels['VT']:(dm1 + 1) * self.levels['VT']
                                                          ] * g0
            fluxlevels_pressure_term = [sum(fluxlevels_pressure_term[0:i])
                                        for i in range(1, len(fluxlevels_pressure_term) + 1)]
            vgrid['fluxlevels_pressure_term'] = numpy.array([0.] + fluxlevels_pressure_term)
            vgrid['masslevels_pressure_term'] = profiles.flux2masspressures(fluxlevels_pressure_term,
                                                                            vertical_mean)
            domains['geometry'].append(geom)
            domains['vertical_grid'].append(vgrid)

        self._attributes['domains'] = domains

    @FileResource._openbeforedelayed
    def _read_validity(self):
        """
        Reads DATE and ECHEANCE fields, to fill-in a FieldValidity object.
        Cf. J.M. Piriou's documentation about DDH.
        """
        DATE = self.readfield('DATE')
        DATE_data = DATE.getdata()
        year = DATE_data[0]
        month = DATE_data[1]
        day = DATE_data[2]
        hour = DATE_data[3]
        minute = DATE_data[4]
        processtype = DATE_data[8]
        if processtype == 0:
            self.processtype = 'analysis'
        elif processtype == 1:
            self.processtype = 'initialization'
        elif processtype == 10:
            self.processtype = 'forecast'
        # cumulated_timesteps = DATE_data[9]

        basis = datetime.datetime(year, month, day, hour, minute)

        ECHEANCE = self.readfield('ECHEANCE')
        term = datetime.timedelta(seconds=int(ECHEANCE.getdata()[0]))
        cumulativeduration = term

        self.validity = FieldValidity(basis=basis,
                                      term=term,
                                      cumulativeduration=cumulativeduration)
