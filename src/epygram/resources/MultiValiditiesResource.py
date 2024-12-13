#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a MultiValiditiesResource.
"""

import copy
import numpy

from footprints import FPList, FPDict, proxy as fpx

from epygram import epygramError
from epygram.util import fmtfid
from epygram.base import Resource, FieldSet, FieldValidityList
from epygram.geometries import VGeometry
from . import open_and_close_resource


class MultiValiditiesResource(Resource):
    """Class implementing a MultiValiditiesResource."""

    _collector = ('resource_modificator', 'epyresource')
    _footprint = dict(
        attr=dict(
            resources=dict(
                type=FPList,
                info="List of resources to join."),
            name=dict(
                values=set(['MultiValidities']))
        )
    )

    def __init__(self, *args, **kwargs):
        super(Resource, self).__init__(*args, **kwargs)
        if len(self.resources) == 0:
            raise epygramError("There must be at least one resource in *resources*")
        for r in self.resources:
            if r.openmode != self.openmode:
                raise epygramError("All low level resources must be opened using the same mode as this high-level resource.")
        if len(set([r.format for r in self.resources])) != 1:
            raise epygramError("All low level resources must have the same format.")
        self.format = "MultiValidities"
        self.lowLevelFormat = self.resources[0].format
        self.isopen = True  # resource always appear to be open

    def close(self):
        """Closes all resources."""
        for r in self.resources:
            try:
                r.close()
            except Exception:
                pass

    def find_fields_in_resource(self, *args, **kwargs):
        """Call to find_fields_in_resource"""
        tmp = []
        for r in self.resources:
            with open_and_close_resource(r):
                tmp.extend(r.find_fields_in_resource(*args, **kwargs))
        result = []
        for res in tmp:
            if res not in result:
                result.append(res)
        return result

    def listfields(self, *args, **kwargs):
        """Call to listfields."""
        complete = 'complete' in kwargs and kwargs['complete']
        tmp = []
        for r in self.resources:
            with open_and_close_resource(r):
                fidlist = r.listfields(*args, **kwargs)
                if complete:
                    for fid in fidlist:
                        fid[self.format] = fid[fmtfid(r.format, fid)]
                tmp.extend(fidlist)
        result = []
        for res in tmp:
            if res not in result:
                result.append(res)
        return result

    def sortfields(self, *args, **kwargs):
        """Call to sortfields"""
        tmp = {}
        for r in self.resources:
            with open_and_close_resource(r):
                for k, v in r.sortfields(*args, **kwargs).items():
                    tmp[k] = tmp.get(k, []) + v
        result = {}
        for k, v in tmp.items():
            result[k] = list(set(v))
        return result

    def readfield(self, *args, **kwargs):
        """
        Reads the field in the different resources and join the validities.
        """
        fieldset = FieldSet()
        for r in self.resources:
            with open_and_close_resource(r):
                fieldset.append(r.readfield(*args, **kwargs))
        return self._join_validities(fieldset, **kwargs)

    def writefield(self, *args, **kwargs):
        """Write fields."""
        raise NotImplementedError("writefield is not impelemented")
        # To implement writefield we need to check that each validity is affected to only one resource

    def extractprofile(self, *args, **kwargs):
        """
        Extracts the profiles in the different resources and join the validities.
        """
        fieldset = FieldSet()
        for r in self.resources:
            with open_and_close_resource(r):
                fieldset.append(r.extractprofile(*args, **kwargs))
        return self._join_validities(fieldset)

    def extractsection(self, *args, **kwargs):
        """
        Extracts the sections in the different resources and join the validities.
        """
        fieldset = FieldSet()
        for r in self.resources:
            with open_and_close_resource(r):
                fieldset.append(r.extractsection(*args, **kwargs))
        return self._join_validities(fieldset)

    @property
    def spectral_geometry(self):
        """
        Returns the spectral_geometry
        """
        ref = self.resources[0].spectral_geometry

        if numpy.all([r.spectral_geometry == ref for r in self.resources]):
            return ref
        else:
            raise epygramError("All spectral_geometry are not identical")

    def _join_validities(self, fieldset, **kwargs):
        """
        Join the different fields
        """
        # Validities
        validities = {}
        for i in range(len(fieldset)):
            field = fieldset[i]
            if len(field.validity) != 1:
                raise NotImplementedError("Join of several multi validities field is not (yet?) implemented")
            validities[field.validity.get()] = i

        if len(validities) == 1 and len(fieldset) != 1:
            validities = {}
            for i in range(len(fieldset)):
                field = fieldset[i]
                if len(field.validity) != 1:
                    raise NotImplementedError("Join of several multi validities field is not (yet?) implemented")
                validities[field.validity.getbasis()] = i

        if len(validities) != len(fieldset):
            raise epygramError("Several resources have returned the same validity")
        sortedValidities = sorted(validities.keys())

        # Geometries
        geometry = fieldset[0].geometry.deepcopy()
        geometry.vcoordinate = VGeometry(typeoffirstfixedsurface=255, levels=[255])
        vcoordinate = fieldset[0].geometry.vcoordinate.deepcopy()
        vcoordinate.levels = []
        levels = fieldset[0].geometry.vcoordinate.levels
        joinLevels = False
        spectral = fieldset[0].spectral
        for field in fieldset[1:]:
            if spectral != field.spectral:
                raise epygramError("All fields must be gridpoint or spectral")
            geometry_field = field.geometry.deepcopy()
            geometry_field.vcoordinate = VGeometry(typeoffirstfixedsurface=255, levels=[255])
            vcoordinate_field = field.geometry.vcoordinate.deepcopy()
            vcoordinate_field.levels = []
            levels_field = field.geometry.vcoordinate.levels
            if geometry != geometry_field:
                raise epygramError("All resources must return fields with the same geometry.")
            if vcoordinate != vcoordinate_field:
                raise epygramError("All resources must return fields with the same vertical geometry.")
            if numpy.any(numpy.array(levels) != numpy.array(levels_field)):
                if len(vcoordinate.grid) > 0:
                    raise epygramError("All resources must return fields with the same vertical geometry.")
                else:
                    if len(levels) != len(levels_field):
                        raise epygramError("All resources must return fields with the same vertical geometry length.")
                    joinLevels = True

        geometry.vcoordinate = fieldset[0].geometry.vcoordinate.deepcopy()
        if joinLevels:
            if spectral:
                raise epygramError("Not sure how to merge vertical levels when spectral")

            concat = numpy.concatenate
            newLevels = []
            for i in range(len(fieldset)):
                mylevel = fieldset[validities[sortedValidities[i]]].geometry.get_levels(d4=True, nb_validities=len(field.validity))
                if isinstance(mylevel, numpy.ma.masked_array):
                    concat = numpy.ma.concatenate
                newLevels.append(mylevel)
            geometry.vcoordinate.levels = list(concat(newLevels, axis=0).swapaxes(0, 1).squeeze())

        # Other metadata
        kwargs_field = copy.deepcopy(fieldset[0].footprint_as_dict())
        kwargs_field['validity'] = FieldValidityList()
        kwargs_field['geometry'] = geometry
        for k, v in kwargs_field.items():
            if isinstance(v, FPDict):
                kwargs_field[k] = dict(v)
        field = fpx.field(**kwargs_field)
        sameProcesstype = True
        for f in fieldset[1:]:
            kwargs_field = copy.deepcopy(f.footprint_as_dict())
            kwargs_field['validity'] = FieldValidityList()
            kwargs_field['geometry'] = geometry
            sameProcesstype = sameProcesstype and kwargs_field['processtype'] == fieldset[0].processtype
            kwargs_field['processtype'] = fieldset[0].processtype  # we exclude processtype from the comparison
            kwargs_field['fid'] = field.fid
            # change FPDict to dict
            for k, v in kwargs_field.items():
                if isinstance(v, FPDict):
                    kwargs_field[k] = dict(v)
            myf = fpx.field(**kwargs_field)
            if field != myf:
                raise epygramError("All fields must be of the same kind")

        # Returned field
        fieldvaliditylist = FieldValidityList()
        fieldvaliditylist.pop()
        getdata = kwargs.get('getdata', True)
        for i in range(len(fieldset)):
            field = fieldset[validities[sortedValidities[i]]]
            fieldvaliditylist.extend(field.validity)
            if getdata:
                if i == 0:
                    data = field.getdata(d4=True)
                else:
                    if isinstance(data, numpy.ma.masked_array):
                        concat = numpy.ma.concatenate
                    else:
                        concat = numpy.concatenate
                    data = concat([data, field.getdata(d4=True)], axis=0)
        if not sameProcesstype:
            kwargs_field['processtype'] = None

        kwargs_field['validity'] = fieldvaliditylist
        kwargs_field['fid'][self.format] = kwargs_field['fid'][fmtfid(self.lowLevelFormat, kwargs_field['fid'])]
        field = fpx.field(**kwargs_field)
        if getdata:
            field.setdata(data)

        return field
