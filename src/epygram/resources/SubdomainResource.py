#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Contains the class that handle a SubdomainResource.
This resource exposes the low level resource fields on a subdomain.
"""

from footprints import FPDict

from epygram.base import Resource, FieldSet
from epygram import epygramError
from epygram.geometries import Geometry


class SubdomainResource(Resource):
    """Class implementing a SubdomainResource."""

    _collector = ('resource_modificator', 'epyresource')
    _footprint = dict(
        attr=dict(
            resource=dict(
                type=Resource,
                info="Low level resource"),
            geometry=dict(
                type=Geometry,
                optional=True,
                default=None,
                info="Geometry on which extract fields"),
            subarray=dict(
                type=FPDict,
                optional=True,
                default=FPDict(),
                info="Dict containing the imin, imax, jmin and jmax keys delimiting the domain to extract"),
            options=dict(
                type=FPDict,
                optional=True,
                default=FPDict(),
                info="Options given to the extract_subdomain method of fields read on the low level resource"),
            name=dict(
                values=set(['Subdomain']))
        )
    )

    def __init__(self, *args, **kwargs):
        super(Resource, self).__init__(*args, **kwargs)
        if self.resource.openmode != self.openmode:
            raise epygramError("The low level resource must be opened using the same mode as this high-level resource.")

        self.format = "Subdomain"
        self.lowLevelFormat = self.resource.format
        if (self.geometry is None and len(self.subarray) == 0) or \
           (self.geometry is not None and len(self.subarray) > 0):
            raise epygramError("one and only one of geometry and subarray must be set")
        self._mode = 'subarray' if self.geometry is None else 'geometry'
        if self._mode == 'subarray':
            if set(self.subarray.keys()) != set(['imin', 'imax', 'jmin', 'jmax']):
                raise epygramError("subarray must contain exactly the imin, imax, jmin and jmax keys")

    def open(self):
        """Opens the low level resource"""
        if not self.resource.isopen:
            self.resource.open()

    @property
    def isopen(self):
        return self.resource.isopen

    def close(self):
        """Closes the low level resource."""
        try:
            self.resource.close()
        except IOError:
            pass

    def find_fields_in_resource(self, *args, **kwargs):
        """Returns a list of the fields from resource matching the request."""
        return self.resource.find_fields_in_resource(*args, **kwargs)

    def listfields(self, *args, **kwargs):
        """Lists the available fields."""
        complete = 'complete' in kwargs and kwargs['complete']
        fidlist = self.resource.listfields(*args, **kwargs)
        if complete:
            for fid in fidlist:
                fid[self.format] = fid[self.lowLevelFormat]
        return fidlist

    def sortfields(self, *args, **kwargs):
        """Returns a sorted list of fields."""
        return self.resource.sortfields(*args, **kwargs)

    def readfield(self, *args, **kwargs):
        """Read the field in the low level resource and extract subdomain."""
        result = self.readfields(*args, **kwargs)
        if len(result) != 1:
            raise epygramError(str(len(result)) + "field(s) have been found, only one expected.")
        return result[0]

    def readfields(self, *args, **kwargs):
        """Read the field in the low level resource and extract subdomain."""
        getdata = kwargs.get('getdata', True)
        fieldset = FieldSet()
        for field in self.resource.readfields(*args, **kwargs):
            if field.spectral:
                if getdata:
                    field.sp2gp()
                else:
                    field._attributes['spectral_geometry'] = None
            if self._mode == 'geometry':
                field_on_subzone = field.extract_subdomain(self.geometry,
                                                           getdata=getdata,
                                                           **self.options)
            else:
                field_on_subzone = field.extract_subarray(self.subarray['imin'],
                                                          self.subarray['imax'],
                                                          self.subarray['jmin'],
                                                          self.subarray['jmax'],
                                                          getdata=getdata,
                                                          deepcopy=False)
            field_on_subzone.fid[self.format] = field_on_subzone.fid[self.lowLevelFormat]
            fieldset.append(field_on_subzone)

        return fieldset

    def writefield(self, *args, **kwargs):
        """Write fields."""
        raise AttributeError("writefield does not exist for this resource.")

    def extractprofile(self, *args, **kwargs):
        """Extracts profiles."""
        profile = self.resource.extractprofile(*args, **kwargs)
        if self._mode == 'geometry':
            lons, lats = profile.geometry.get_lonlat_grid()
            for (lon, lat) in zip(lons.flatten(), lats.flatten()):
                if not self.geometry.point_is_inside_domain_ll(lon, lat):
                    raise epygramError("Profile is not in the subdomain geometry.")
        return profile

    def extractsection(self, *args, **kwargs):
        """Extracts sections."""
        section = self.resource.extractsection(*args, **kwargs)
        if self._mode == 'geometry':
            lons, lats = section.geometry.get_lonlat_grid()
            for (lon, lat) in zip(lons.flatten(), lats.flatten()):
                if not self.geometry.point_is_inside_domain_ll(lon, lat):
                    raise epygramError("Section is not in the subdomain geometry.")
        return section

    @property
    def spectral_geometry(self):
        """Returns the spectral_geometry."""
        return self.resource.spectral_geometry
