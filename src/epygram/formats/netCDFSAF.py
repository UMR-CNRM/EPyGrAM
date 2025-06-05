#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to handle HDF5 format used by SAF products.
"""

import datetime
import sys
import numpy

import footprints
from footprints import FPDict

import netCDF4
from epygram import epygramError, util
from epygram.util import Angle
from epygram.base import FieldSet, FieldValidity, Field
from epygram.geometries import VGeometry, ProjectedGeometry
from epygram.resources import FileResource
from epygram.fields import H2DField

__all__ = ['netCDFSAF']

epylog = footprints.loggers.getLogger(__name__)


class netCDFSAF(FileResource):
    """
    Class implementing all specificities for netCDFSAF resource format.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['netCDFSAF']),
                default='netCDFSAF'),
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        self.geometry = None
        self.validity = None

        super(netCDFSAF, self).__init__(*args, **kwargs)

        if not self.fmtdelayedopen:
            self.open()

    def open(self):
        """Opens a netCDFSAF, and initializes some attributes."""
        # Opening
        if self.openmode in ('r'):
            self.nc = netCDF4.Dataset(self.container.abspath, 'r')
            if not all(hasattr(self.nc, k) for k in ['nominal_product_time',
                                                     'cgms_projection',
                                                     'region_name',
                                                     'sub-satellite_longitude',
                                                     'satellite_identifier']):
                self.nc.close()
                raise IOError('netCDF file is not a netCDFSAF file')
            self.isopen = True
        else:
            raise NotImplementedError("netCDFSAF is only implemented for reading")

        # Reading of metadata
        if self.geometry is None:
            self._read_dategeom()

    def close(self):
        """Closes a netCDFSAF file."""
        if self.isopen:
            self.isopen = False
            # Cleanings
            self.nc.close()
            del self.nc

################
# ABOUT FIELDS #
################
    def find_fields_in_resource(self, seed=None, fieldtype=[], generic=False):
        """
        Returns a list of the fields from resource whose identifier match the given seed.

        Args: \n
        - *seed*: might be a tuple of regular expressions, a list of regular expressions tuples or *None*.
                  If *None* (default), returns the list of all fields in resource.
        - *fieldtype*: optional, among ('H2D', 'Misc') or a list of these strings.
          If provided, filters out the fields not of the given types.
        - *generic*: if True, returns a list of tuples (fieldname, generic fid) of
          the fields.
        """

        if type(fieldtype) == type(list()):
            fieldtypeslist = list(fieldtype)
        else:
            fieldtypeslist = [fieldtype]
        fieldslist = []
        if seed is None:
            fieldslist = self.listfields()
        elif isinstance(seed, str):
            fieldslist = util.find_re_in_list(seed, self.listfields())
        elif isinstance(seed, list):
            fieldslist = []
            for s in seed:
                fieldslist += util.find_re_in_list(s, self.listfields())
        if fieldtypeslist != [] and 'H2D' not in fieldtypeslist:
            # Only H2D fields are stored in HDF5SAF format
            fieldslist = []
        if fieldslist == []:
            raise epygramError("no field matching '" + str(seed) + "' was found in resource " + self.container.abspath)

        if generic:
            fieldslist = [(f, {}) for f in fieldslist]

        return fieldslist

    def listfields(self, **kwargs):
        """
        Returns a list containing the netCDFSAF identifiers of all the fields of the resource.
        """
        return super(netCDFSAF, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for HDF5SAF.

        Args: \n
        - *complete*: if True method returns a list of {'netCDFSAF':netCDFSAF_fid,
                                                        'generic':generic_fid}
                      if False method return a list of netCDFSAF_fid
        """
        fieldslist = [v for v in self.nc.variables
                      if all(d in self.nc[v].dimensions for d in ('nx', 'ny'))]

        if complete:
            return [{'netCDFSAF': f, 'generic': {}} for f in fieldslist]
        else:
            return fieldslist

    def sortfields(self):
        """
        Returns a sorted list of fields with regards to their name,
        as a dict of lists.
        """

        mylist = self.listfields()

        # sort
        mylist.sort()

        outlists = {'Fields': mylist}

        return outlists

    @FileResource._openbeforedelayed
    def readfield(self, fieldidentifier, getdata=True):
        """
        Reads one field, given its identifier, and returns a Field instance.
        Interface to Fortran routines from 'ifsaux'.

        Args: \n
        - *fieldidentifier*: name.
        - *getdata*: optional, if *False*, only metadata are read, the field do not contain data.
                     Default is *True*.
        """

        if self.openmode == 'w':
            raise epygramError("cannot read fields in resource with openmode == 'w'.")
        if not isinstance(fieldidentifier, str):
            raise epygramError("fieldidentifier of a HDF5SAF field is a string.")
        if fieldidentifier not in self.listfields():
            raise epygramError("fieldidentifier is not found in the HDF5SAF file.")

        field = H2DField(fid=FPDict({self.format:fieldidentifier, 'generic':FPDict()}),
                         structure=self.geometry.structure,
                         geometry=self.geometry.deepcopy(),
                         validity=self.validity.deepcopy(),
                         processtype='observation',
                         comment="")
        if getdata:
            if len(self.nc[fieldidentifier].shape) == 3:
                # First dimension is time
                field.setdata(self.nc[fieldidentifier][0, ::-1, :])
            else:
                field.setdata(self.nc[fieldidentifier][::-1, :])

        return field

    def readfields(self, requestedfields=None, getdata=True):
        """
        Returns a :class:`epygram.base.FieldSet` containing requested fields read in the resource.

        Args: \n
        - *requestedfields*: might be \n
          - an expressions with regular expressions (e.g. '*CMS*')
          - a list of HDF5SAF fields identifiers with regular expressions (e.g. ['*TIME*', '*CMS*'])
          - if not specified, interpretated as all fields that will be found in resource
        - *getdata*: optional, if *False*, only metadata are read, the fields do not contain data.
                     Default is *True*.
        """

        requestedfields = self.find_fields_in_resource(requestedfields)
        if requestedfields == []:
            raise epygramError("unable to find requested fields in resource.")

        return super(netCDFSAF, self).readfields(requestedfields, getdata)

    def writefield(self, field):
        """
        Write a field in the resource.

        Args: \n
        - *field*: a :class:`epygram.base.Field` instance or :class:`epygram.H2DField`.
        """

        if not isinstance(field, Field):
            raise epygramError("*field* must be a Field instance.")

        raise NotImplementedError("Writing of netCDFAF file is not implemented.")

    def writefields(self, fieldset):
        """
        Write the fields of the *fieldset* in the resource.

        Args: \n
        - *fieldset*: must be a :class:`epygram.base.FieldSet` instance.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        if self.openmode == 'r':
            raise IOError("cannot write field in a netCDFAF with openmode 'r'.")

        raise NotImplementedError("Writing of netCDFAF file is not implemented.")

###########
# pre-app #
###########
    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout,
             details=False,
             sortfields=False,
             **kwargs):
        """
        Writes in file a summary of the contents of the HDF5SAF.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *sortfields*: **True** if the fields have to be sorted.
        """
        firstcolumn_width = 50
        secondcolumn_width = 16
        sepline = '{:-^{width}}'.format('', width=firstcolumn_width + secondcolumn_width + 1) + '\n'

        out.write("### FORMAT: " + self.format + "\n")
        out.write("\n")

        listfields = self.listfields()
        first_H2DField = listfields[0]
        firstfield = self.readfield(first_H2DField, getdata=False)

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
        def write_formatted_fields(dest, label):
            dest.write('{:<{width}}'.format(label, width=20) +
                       '\n')

        firstfield.what(out, vertical_geometry=False)

        out.write("######################\n")
        out.write("### LIST OF FIELDS ###\n")
        out.write("######################\n")
        if sortfields:
            sortedfields = self.sortfields()
            listoffields = []
            for k in sorted(sortedfields.keys()):
                listoffields.append(k)
                listoffields.append('--------------------')
                listoffields.extend(sortedfields[k])
                listoffields.append('--------------------')
            numfields = sum([len(v) for v in sortedfields.values()])
        else:
            listoffields = listfields
            numfields = len(listfields)
        out.write("Number: " + str(numfields) + "\n")
        if details:
            write_formatted_fields(out, "Field name")
        else:
            write_formatted_fields(out, "Field name")
        out.write(sepline)
        done = []
        for f in listoffields:
            if f not in done:
                if details:
                    write_formatted_fields(out, f)
                else:
                    write_formatted_fields(out, f)
                done.append(f)
        out.write(sepline)

##############
# the HDF5SAF WAY #
##############
    @FileResource._openbeforedelayed
    def _read_dategeom(self):
        """
        Reads netCDFSAF tags.
        """

        # Projection
        attrs = {item.split('=')[0]: item.split('=')[1]
                 for item in self.nc.cgms_projection.replace('+', '').split()}
        geo = attrs['proj'] == 'geos'

        # Date
        date = datetime.datetime.fromisoformat(self.nc.nominal_product_time.replace('Z', '+00:00'))
        self.validity = FieldValidity(basis=date, term=datetime.timedelta(hours=0))

        # Grid
        kwargs_vcoord = {'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [255]}
        vcoordinate = VGeometry(**kwargs_vcoord)
        if geo:
            # Space view
            if self.nc.satellite_identifier in ('MTGI1', 'MTI1', 'MSG3'):
                Lap = Angle(0., 'degrees')  # Certainly in degrees as others lon/lat in file
            else:
                raise NotImplementedError("Please check latitude to use for this satellite:",
                                          self.nc.satellite_identifier)
            Lop = Angle(getattr(self.nc, 'sub-satellite_longitude'), 'degrees')  # Certainly in degrees as others lon/lat in file

            rmajor, rminor = float(attrs['r_eq']) * 1000., float(attrs['r_pol']) * 1000.
            h = float(attrs['h']) * 1000. - rmajor

            # use of cfac/lfac deduced from http://www.umr-cnrm.fr/remote-sensing/IMG/pdf/lsasaf_mf_pum_al_msg_1.10.pdf
            # Here we use the sin function whereas it seems that sin(x) was approximated by x in TIFFMF ????
            # relative (absolute) difference between both formulation is about 1E9 (1E-6 meters) for resolx/resoly
            resolx = h * numpy.sin(numpy.radians(1. / (float(attrs['cfac']) * 2**-16)))
            resoly = h * numpy.sin(numpy.radians(1. / (float(attrs['lfac']) * 2**-16)))

            # In TIFFMF it was needed to shift by .5 pixel depending on the size parity
            # Here, the following input_position is found to well reproduce the lon/lat given
            # in external files, and to position correctly the field on a map
            if self.nc.region_name != 'FULL FRAME':
                raise NotImplementedError(f"Not checked on region: {self.nc.region_name}")
            nx, ny = self.nc.dimensions['nx'].size, self.nc.dimensions['ny'].size
            input_position = (float(attrs['coff']) - 1, ny - float(attrs['loff']))

            grid = {'LAMzone': None,
                    'X_resolution': resolx,
                    'Y_resolution': resoly,
                    'input_lon': Lop, 'input_lat': Lap,
                    'input_position': input_position
                    }
            projection = {'satellite_lat': Lap, 'satellite_lon': Lop, 'satellite_height': h,
                          'reference_lon': Lop,
                          'rotation': Angle(0., 'degrees')}
            dimensions = {'X': nx, 'Y': ny}
            geometryname = 'space_view'
            kwargs_geom = dict(name=geometryname,
                               grid=FPDict(grid),
                               dimensions=FPDict(dimensions),
                               projection=FPDict(projection),
                               position_on_horizontal_grid='center',
                               vcoordinate=vcoordinate,
                               geoid={'a': rmajor, 'b': rminor},
                               )
            self.geometry = ProjectedGeometry(**kwargs_geom)
        else:
            raise NotImplementedError("This projection is not implemented")
