#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to handle TIFFMF format.
"""

import datetime
import numpy
import sys

import footprints
from footprints import FPDict, proxy as fpx

from epygram.extra import pyexttiff

from epygram import config, epygramError, util
from epygram.util import Angle
from epygram.base import FieldSet, FieldValidity, Field
from epygram.resources import FileResource
from epygram.geometries import VGeometry, ProjectedGeometry
from epygram.fields import H2DField, make_vector_field

__all__ = ['TIFFMF']

epylog = footprints.loggers.getLogger(__name__)


class TIFFMF(FileResource):
    """
    Class implementing all specificities for TIFFMF resource format.
    """

    _footprint = dict(
        attr=dict(
            format=dict(
                values=set(['TIFFMF']),
                default='TIFFMF'),
        )
    )

    def __init__(self, *args, **kwargs):
        self.isopen = False
        self.geometry = None
        self.validity = None
        self.scan = None

        super(TIFFMF, self).__init__(*args, **kwargs)

        if not self.fmtdelayedopen:
            self.open()

    def open(self):
        """Opens a TIFFMF, and initializes some attributes."""
        # Opening
        if self.openmode in ('r'):
            self.tiff = pyexttiff.TiffFile(self.container.abspath, [(34974,)], method=2)
            try:
                self.tiff.IFDs[0].get_value(34974)
            except:
                raise epygramError('Tiff file is not a tiffmf file')
            self.isopen = True
        else:
            raise NotImplementedError("TIFFMF is only implmented for reading")

        # Reading of metadata
        if self.geometry is None:
            self._read_sections()

    def close(self):
        """Closes a TIFFMF file."""
        if self.isopen:
            self.isopen = False
            # Cleanings
            self.tiff.close()
            del self.tiff

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
            # Only H2D fields are stored in TIFFMF format
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
        return super(TIFFMF, self).listfields(**kwargs)

    @FileResource._openbeforedelayed
    def _listfields(self, complete=False):
        """
        Actual listfields() method for TIFFMF.

        Args: \n
        - *complete*: if True method returns a list of {'TIFFMF':TIFFMF_fid, 'generic':generic_fid}
                      if False method return a list of TIFFMF_fid
        """
        fieldslist = [IFD.get_value(270) for IFD in self.tiff.IFDs]

        if complete:
            return [{'TIFFMF':f, 'generic':{}} for f in fieldslist]
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
            raise epygramError("fieldidentifier of a TIFFMF field is a string.")
        if fieldidentifier not in self.listfields():
            raise epygramError("fieldidentifier is not found in th TIFFMF file.")

        if self.scan != 0:
            raise NotImplementedError("This scan mode is not implemented")
        for IFD in self.tiff.IFDs:
            if IFD.get_value(270) == fieldidentifier:
                data = IFD.get_image()[::-1, :]
                if len(data.shape) != 2:
                    if IFD.has_tag('PhotometricInterpretation') and \
                       IFD.get_entry('PhotometricInterpretation').get_value() in (2, 5, 6, 8):
                        componentNames = {2:['R', 'G', 'B'],
                                          5:['C', 'M', 'Y', 'K'],
                                          6:['Y', 'Cb', 'Cr'],
                                          8:['L*', 'a*', 'b*']}[IFD.get_entry('PhotometricInterpretation').get_value()]
                        assert len(componentNames) == data.shape[-1], "Data shape error"
                        fields = []
                        for i, cname in enumerate(componentNames):
                            field = H2DField(fid=FPDict({self.format:fieldidentifier,
                                                         'generic':FPDict(),
                                                         'color':cname}),
                                             structure=self.geometry.structure,
                                             geometry=self.geometry.deepcopy(),
                                             validity=self.validity.deepcopy(),
                                             processtype='observation',
                                             comment="")
                            if getdata:
                                field.setdata(data[..., i])
                            fields.append(field)
                        field = make_vector_field(*fields)
                    else:
                        raise NotImplementedError("This file certainly contains an image in true colors in a" + \
                                                  "color space not (yet) implemented")
                else:
                    field = H2DField(fid=FPDict({self.format:fieldidentifier, 'generic':FPDict()}),
                                     structure=self.geometry.structure,
                                     geometry=self.geometry.deepcopy(),
                                     validity=self.validity.deepcopy(),
                                     processtype='observation',
                                     comment="")
                    if getdata:
                        field.setdata(data)

        return field

    def readfields(self, requestedfields=None, getdata=True):
        """
        Returns a :class:`epygram.base.FieldSet` containing requested fields read in the resource.

        Args: \n
        - *requestedfields*: might be \n
          - an expressions with regular expressions (e.g. '*CMS*')
          - a list of TIFFMF fields identifiers with regular expressions (e.g. ['*TIME*', '*CMS*'])
          - if not specified, interpretated as all fields that will be found in resource
        - *getdata*: optional, if *False*, only metadata are read, the fields do not contain data.
                     Default is *True*.
        """

        requestedfields = self.find_fields_in_resource(requestedfields)
        if requestedfields == []:
            raise epygramError("unable to find requested fields in resource.")

        return super(TIFFMF, self).readfields(requestedfields, getdata)

    def writefield(self, field):
        """
        Write a field in the resource.

        Args: \n
        - *field*: a :class:`epygram.base.Field` instance or :class:`epygram.H2DField`.
        """

        if not isinstance(field, Field):
            raise epygramError("*field* must be a Field instance.")

        raise NotImplementedError("Writing of TIFFMF file is not implemented.")

    def writefields(self, fieldset):
        """
        Write the fields of the *fieldset* in the resource.

        Args: \n
        - *fieldset*: must be a :class:`epygram.base.FieldSet` instance.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        if self.openmode == 'r':
            raise IOError("cannot write field in a TIFFMF with openmode 'r'.")

        raise NotImplementedError("Writing of TIFFMF file is not implemented.")

###########
# pre-app #
###########
    @FileResource._openbeforedelayed
    def what(self, out=sys.stdout,
             details=False,
             sortfields=False,
             **kwargs):
        """
        Writes in file a summary of the contents of the TIFFMF.

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
# the TIFFMF WAY #
##############
    @FileResource._openbeforedelayed
    def _read_sections(self):
        """
        Reads section 1 and section 2 in the TIIF tags.
        """

        # MF tags reading
        IFD = self.tiff.IFDs[0].get_value(34974)
        typeImage = IFD.get_value(50002)
        subTypeImage = IFD.get_value(50003)
        date = IFD.get_value(50006)
        projection = IFD.get_value(50066)
        section1 = IFD.get_value(60000).view(dtype=self.tiff.dtypes.int32)
        headerSection2 = IFD.get_value(60001).view(dtype=self.tiff.dtypes.int32)
        bodySection2 = IFD.get_value(60002).view(dtype=self.tiff.dtypes.int32)

        # section decoding
        (lengthS1, version, centre, process, grid, flag,
         param, satellite, band, year, month, day, hour, minute,
         unit, P1, P2, timerange, included, missing,
         century, D) = section1
        (lengthS2, NV, PV_PL, dataRepresenationType) = headerSection2
        kwargs_vcoord = {'typeoffirstfixedsurface': 255,
                         'position_on_grid': '__unknown__',
                         'levels': [255]}
        vcoordinate = VGeometry(**kwargs_vcoord)
        if projection == 1:
            # Polar streographic
            if dataRepresenationType != 5:
                raise epygramError("Projection does not match.")
            (Nx, Ny, La1, Lo1, res, LoV, Dx, Dy, pole, scan) = bodySection2
            if typeImage == 7 and \
               subTypeImage == 13 and \
               (version, centre, process) == (1, 211, 220) and \
               (param, satellite, band) == (127, 200, 33792) and \
               (Nx, Ny) == (8192, 6144) and \
               (La1, Lo1, Dx, Dy) == (43759, -76033, 1093, 1093) and \
               LoV == 180000:
                # If tiffmf file is modified to put LoV=0. in the section 2,
                # R2 is no more able to put the image at the good location.
                #
                # So, it seems there is an error somewhere in this module, so we need to
                # reset LoV to 0. to have the good coordinates.
                # This must be corrected
                epylog.warning("LoV reset to 0. for this file.")
                LoV = 0
            else:
                raise epygramError("Do we need to reset LoV to 0.? If you encounter this " +
                                   "error, please check if we need or not to reset " +
                                   "this parameter and modify the test above.")
            if res != 64:
                raise NotImplementedError('Only res=64 is implemented')
            Nx = int(Nx)
            Ny = int(Ny)
            lat_0 = 90 if pole == 0 else -90
            m = (1. + numpy.copysign(1., lat_0) * numpy.sin(numpy.radians(lat_0))) / \
                (1. + numpy.copysign(1., lat_0) * numpy.sin(numpy.radians(60)))
            grid = {'LAMzone':None,
                    'X_resolution':float(Dx) * m,
                    'Y_resolution':float(Dy) * m,
                    'input_lon':Angle(float(Lo1) / 1000., 'degrees'),
                    'input_lat':Angle(float(La1) / 1000., 'degrees'),
                    'input_position':(0, Ny),
                    }
            if scan != 0:
                raise NotImplementedError("Polar stereographic projection with scan != 0 is not implemented.")
            projection = {'reference_lon':Angle(float(LoV) / 1000., 'degrees'),
                          'rotation': Angle(0., 'radians'),
                          'reference_lat':Angle(lat_0, 'degrees')}
            dimensions = {'X':Nx, 'Y':Ny}
            geometryname = 'polar_stereographic'
            kwargs_geom = dict(name=geometryname,
                               grid=FPDict(grid),
                               dimensions=FPDict(dimensions),
                               projection=FPDict(projection),
                               geoid=FPDict(config.default_geoid),
                               position_on_horizontal_grid='center',
                               vcoordinate=vcoordinate
                               )
            self.geometry = ProjectedGeometry(**kwargs_geom)
        elif projection == 3:
            # Mercator
            if dataRepresenationType != 1: raise epygramError("Projection does not match.")
            (Ni, Nj, La1, Lo1, res, La2, Lo2, Latin, scan, Di, Dj) = bodySection2
            raise NotImplementedError("This projection is not implemented (mercator) but could easily be done.")
        elif projection == 11:
            # Space view
            if dataRepresenationType != 90: raise epygramError("Projection does not match.")
            (Nx, Ny, Lap, Lop, res, dx, dy, Xp, Yp, scan, ort, nr, Xo, Yo) = bodySection2
            Lap = Angle(int(Lap) / 1000., 'degrees')
            Lop = Angle(int(Lop) / 1000., 'degrees')
            Nx = int(Nx)
            Ny = int(Ny)
            #rmajor, rminor = 6378169.0, 6356583.8 #old values used by an old version of the ASPIC software
            rmajor, rminor = 6378160, 6356775 #new version of ASPIC
            geoid = {'a':rmajor, 'b':rminor}
            h = int(nr) * rmajor / 1E6 - rmajor
            resolx = h * 2 * numpy.arcsin(1E6 / int(nr)) / int(dx)
            resoly = h * 2 * numpy.arcsin(1E6 / int(nr)) / int(dy)
            
            #Xp, Yp, Xo and YO are integers
            #We guess the exact position of the satellite from the parity of the apparent diameter
            #If the numbers of grid-cells is even, satellite is in between two pixels 
            grid = {'LAMzone':None,
                    'X_resolution':resolx,
                    'Y_resolution':resoly,
                    'input_lon':Lop, 'input_lat':Lap,
                    'input_position':(int(Xp) - int(Xo) - (int(dx) + 1) % 2 * 0.5,
                                      int(Yo) - int(Yp) + Ny - 1 + (int(dy) + 1) % 2 * 0.5)
                    }
            if scan != 0:
                raise NotImplementedError("Space view projection with scan != 0 is not implemented.")
            if ort != 0.:
                raise epygramError("Must be checked for rotated geometry.")
            projection = {'satellite_lat':Lap, 'satellite_lon':Lop, 'satellite_height':h,
                          'reference_lon':Lop,
                          'rotation': Angle(int(ort) / 1000., 'degrees')}
            dimensions = {'X':Nx, 'Y':Ny}
            geometryname = 'space_view'
            kwargs_geom = dict(name=geometryname,
                               grid=FPDict(grid),
                               dimensions=FPDict(dimensions),
                               projection=FPDict(projection),
                               position_on_horizontal_grid='center',
                               vcoordinate=vcoordinate,
                               geoid=geoid,
                               )
            self.geometry = ProjectedGeometry(**kwargs_geom)
        elif projection == 15:
            # Cylindric ?
            raise NotImplementedError("Cylindric projection is not implemented")
            (Ni, Nj, La1, Lo1, res, La2, Lo2, Di, Dj, scan) = bodySection2
        else:
            # In case we need it: lambert would be (Nx, Ny, La1, Lo1, res, Lon0, Lat0, K0,
            #                                       X0, Y0, Dx, Dy, pole, scan) = bodySection2
            raise NotImplementedError("This projection is not implemented")

        self.scan = scan

        # Validity
        date = datetime.datetime(int(date[:2].view(dtype=self.tiff.dtypes.uint16)),
                                 date[2], date[3], date[4], date[5])
        date2 = datetime.datetime(year + 100 * (century - 1), month, day, hour, minute)
        if date != date2:
            raise epygramError("The two dates encoded in the tiff file must be the same.")
        self.validity = FieldValidity(basis=date, term=datetime.timedelta(hours=0))
