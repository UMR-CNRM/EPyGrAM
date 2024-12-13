#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to handle HDF5 format used by SAF products.
"""

import datetime
import numpy
import sys

import footprints
from footprints import FPDict, proxy as fpx

import h5py
from epygram import epygramError, util
#from epygram import config, util
from epygram.util import Angle
from epygram.base import FieldSet, FieldValidity, Field
from epygram.geometries import VGeometry, ProjectedGeometry
from epygram.resources import FileResource
from epygram.fields import H2DField

__all__ = ['HDF5SAF']

epylog = footprints.loggers.getLogger(__name__)


class HDF5SAF(FileResource):
    """
    Class implementing all specificities for HDF5SAF resource format.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['HDF5SAF']),
                default='HDF5SAF'),
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        self.geometry = None
        self.validity = None

        super(HDF5SAF, self).__init__(*args, **kwargs)

        if not self.fmtdelayedopen:
            self.open()

    def open(self):
        """Opens a HDF5SAF, and initializes some attributes."""
        # Opening
        if self.openmode in ('r'):
            self.hdf5 = h5py.File(self.container.abspath, 'r')
            if not all([k in self.hdf5.attrs for k in ['FORECAST_STEP', 'IMAGE_ACQUISITION_TIME',
                                                       'CFAC', 'COFF', 'LFAC', 'LOFF', 'NC', 'NL',
                                                       'ORBIT_TYPE', 'PROJECTION_NAME', 'SATELLITE',
                                                       'REGION_NAME', 'SUB_SATELLITE_POINT_END_LAT',
                                                       'SUB_SATELLITE_POINT_END_LON',
                                                       'SUB_SATELLITE_POINT_START_LAT',
                                                       'SUB_SATELLITE_POINT_START_LON']]):
                self.hdf5.close()
                raise IOError('HDF5 file is not a HDF5SAF file')
            self.isopen = True
        else:
            raise NotImplementedError("HDF5SAF is only implemented for reading")

        # Reading of metadata
        if self.geometry is None:
            self._read_dategeom()

    def close(self):
        """Closes a HDF5SAF file."""
        if self.isopen:
            self.isopen = False
            # Cleanings
            self.hdf5.close()
            del self.hdf5

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
        Returns a list containing the TIFFMF identifiers of all the fields of the resource.
        """
        return super(HDF5SAF, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for HDF5SAF.

        Args: \n
        - *complete*: if True method returns a list of {'TIFFMF':TIFFMF_fid, 'generic':generic_fid}
                      if False method return a list of TIFFMF_fid
        """
        def recursive(parent, path):
            path += '' if len(path) == 0 else '/' 
            result = []
            for k, v in parent.items():
                if isinstance(v, h5py.Group):
                    result.extend(recursive(v, path + k))
                else:
                    result.append(path + k)
            return result
            
        fieldslist = recursive(self.hdf5, '')

        if complete:
            return [{'HDF5SAF':f, 'generic':{}} for f in fieldslist]
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

        outlists = {'Fields':mylist}

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
            data = self.hdf5
            for fid in fieldidentifier.split('/'):
                data = data[fid]
            field.setdata(data[...][::-1, :])

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

        return super(HDF5SAF, self).readfields(requestedfields, getdata)

    def writefield(self, field):
        """
        Write a field in the resource.

        Args: \n
        - *field*: a :class:`epygram.base.Field` instance or :class:`epygram.H2DField`.
        """

        if not isinstance(field, Field):
            raise epygramError("*field* must be a Field instance.")

        raise NotImplementedError("Writing of HDF5SAF file is not implemented.")

    def writefields(self, fieldset):
        """
        Write the fields of the *fieldset* in the resource.

        Args: \n
        - *fieldset*: must be a :class:`epygram.base.FieldSet` instance.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        if self.openmode == 'r':
            raise IOError("cannot write field in a HDF5SAF with openmode 'r'.")

        raise NotImplementedError("Writing of HDF5SAF file is not implemented.")

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
            dest.write('{:<{width}}'.format(label, width=firstcolumn_width) \
                    + ':' \
                    + '{:>{width}}'.format(str(value), width=secondcolumn_width) \
                    + '\n')
        def write_formatted_col(dest, label, value):
            dest.write('{:>{width}}'.format(label, width=firstcolumn_width) \
                    + ':' \
                    + '{:>{width}}'.format(str(value), width=secondcolumn_width) \
                    + '\n')
        def write_formatted_fields(dest, label):
            dest.write('{:<{width}}'.format(label, width=20) \
                    + '\n')

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
        Reads HDF5 tags.
        """

        #Projection
        geo = self.hdf5.attrs['ORBIT_TYPE'] == b'GEO' and \
              self.hdf5.attrs['PROJECTION_NAME'].startswith(b'GEOS(') and \
              self.hdf5.attrs['SUB_SATELLITE_POINT_END_LAT'] == \
              self.hdf5.attrs['SUB_SATELLITE_POINT_START_LAT'] and \
              self.hdf5.attrs['SUB_SATELLITE_POINT_END_LON'] == \
              self.hdf5.attrs['SUB_SATELLITE_POINT_START_LON']
        
        #Date
        if self.hdf5.attrs['FORECAST_STEP'] != 0:
            raise NotImplementedError("Does not know the forecast_step unit")
        date = self.hdf5.attrs['IMAGE_ACQUISITION_TIME'].decode('UTF8')
        if len(date) == 12:
            date = datetime.datetime.strptime(date, '%Y%m%d%H%M')
        else:
            date = datetime.datetime.strptime(date, '%Y%m%d%H%M%S')
        self.validity = FieldValidity(basis=date, term=datetime.timedelta(hours=0))

        #Grid
        cfac, coff, lfac, loff, nc, nl, satellite = [self.hdf5.attrs[k] for k in
                                                     ['CFAC', 'COFF', 'LFAC', 'LOFF', 'NC', 'NL', 'SATELLITE']]

        kwargs_vcoord = {'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [255]}
        vcoordinate = VGeometry(**kwargs_vcoord)
        if geo:
            # Space view
            Lap = Angle(self.hdf5.attrs['SUB_SATELLITE_POINT_START_LAT'], 'degrees')
            Lop = Angle(self.hdf5.attrs['SUB_SATELLITE_POINT_START_LON'], 'degrees')

            #The following values seem to be absent from the metadata in the HDF5
            #They are taken from the values obtained with a TIFFMF file (h) and ASPIC software (rmajor, rminor)
            #It could be necessary to adjust them for other satellite
            if satellite[0] not in (b'MSG4', b'MSG3'):
                raise NotImplementedError("Please check altitude and geoid to use for this satellite:", satellite[0])
            rmajor, rminor = 6378160, 6356775
            h = 35785833.88328
            
            #use of cfac/lfac deduced from http://www.umr-cnrm.fr/remote-sensing/IMG/pdf/lsasaf_mf_pum_al_msg_1.10.pdf
            #Here we use the sin function whereas it seems that sin(x) was approximated by x in TIFFMF ????
            #relative (absolute) difference between both formulation is about 1E9 (1E-6 meters) for resolx/resoly
            resolx = h * numpy.sin(numpy.radians(1. / (cfac * 2**-16)))
            resoly = h * numpy.sin(numpy.radians(1. / (lfac * 2**-16)))
            
            #In TIFFMF it was needed to shift by .5 pixel depending on the size parity
            #Here, the following input_position is found to well reproduce the lon/lat given
            #by the SAFLAND in external files, and to position correctly the field on a map
            if self.hdf5.attrs['REGION_NAME'] not in (b'Full-Disk', b'MSG-Disk'):
                raise NotImplementedError("Not checked on region: " + self.hdf5.attrs['REGION_NAME'].decode('UTF8'))
            input_position = (coff - 1, nl - loff)
            
            grid = {'LAMzone':None,
                    'X_resolution':resolx,
                    'Y_resolution':resoly,
                    'input_lon':Lop, 'input_lat':Lap,
                    'input_position':input_position
                    }
            projection = {'satellite_lat':Lap, 'satellite_lon':Lop, 'satellite_height':h,
                          'reference_lon':Lop,
                          'rotation': Angle(0., 'degrees')}
            dimensions = {'X':nc, 'Y':nl}
            geometryname = 'space_view'
            kwargs_geom = dict(name=geometryname,
                               grid=FPDict(grid),
                               dimensions=FPDict(dimensions),
                               projection=FPDict(projection),
                               position_on_horizontal_grid='center',
                               vcoordinate=vcoordinate,
                               geoid={'a':rmajor, 'b':rminor},
                               )
            self.geometry = ProjectedGeometry(**kwargs_geom)
        else:
            raise NotImplementedError("This projection is not implemented")
