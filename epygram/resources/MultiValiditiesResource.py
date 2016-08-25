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
from epygram.base import Resource, FieldSet, FieldValidityList

class _open_and_close(object):
    def __init__(self, r):
        self.r = r
    def __enter__(self):
        self.r.open()
    def __exit__(self, t, v, tbk):
        self.r.close()

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
        """Constructor. See its footprint for arguments."""

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

#    def open(self):
#        """Opens all resources"""

#        for r in self.resources:
#            if not r.isopen: r.open()

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
            with _open_and_close(r):
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
            with _open_and_close(r):
                fidlist = r.listfields(*args, **kwargs)
                if complete:
                    for fid in fidlist:
                        fid[self.format] = fid[self.lowLevelFormat]
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
            with _open_and_close(r):
                for k, v in r.sortfields(*args, **kwargs).iteritems():
                    tmp[k] = tmp.get(k, []) + v
        result = {}
        for k, v in tmp.iteritems():
            result[k] = list(set(v))
        return result

    def readfield(self, *args, **kwargs):
        """
        Reads the field in the different resources and join the validities.
        """
        fieldset = FieldSet()
        for r in self.resources:
            with _open_and_close(r):
                fieldset.append(r.readfield(*args, **kwargs))
        return self._join_validities(fieldset)

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
            with _open_and_close(r):
                fieldset.append(r.extractprofile(*args, **kwargs))
        return self._join_validities(fieldset)

    def extractsection(self, *args, **kwargs):
        """
        Extracts the sections in the different resources and join the validities.
        """

        fieldset = FieldSet()
        for r in self.resources:
            with _open_and_close(r):
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


    def _join_validities(self, fieldset):
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
        kwargs_geom = copy.deepcopy(fieldset[0].geometry.footprint_as_dict())
        kwargs_geom['vcoordinate'] = fpx.geometry(structure='V', typeoffirstfixedsurface=255, levels=[255])
        for k, v in kwargs_geom.iteritems():
            if type(v) == type(FPDict()):
                kwargs_geom[k] = dict(v)
        geometry = fpx.geometry(**kwargs_geom)
        kwargs_vcoord = copy.deepcopy(fieldset[0].geometry.vcoordinate.footprint_as_dict())
        kwargs_vcoord['levels'] = []
        vcoordinate = fpx.geometry(**kwargs_vcoord)
        levels = fieldset[0].geometry.vcoordinate.levels
        joinLevels = False
        spectral = fieldset[0].spectral
        for field in fieldset[1:]:
            if spectral != field.spectral:
                raise epygramError("All fields must be gridpoint or spectral")
            kwargs_geom = copy.deepcopy(field.geometry.footprint_as_dict())
            kwargs_geom['vcoordinate'] = fpx.geometry(structure='V', typeoffirstfixedsurface=255, levels=[255])
            for k, v in kwargs_geom.iteritems():
                if type(v) == type(FPDict()):
                    kwargs_geom[k] = dict(v)
            geometry_field = fpx.geometry(**kwargs_geom)
            kwargs_vcoord = copy.deepcopy(field.geometry.vcoordinate.footprint_as_dict())
            kwargs_vcoord['levels'] = []
            vcoordinate_field = fpx.geometry(**kwargs_vcoord)
            levels_field = field.geometry.vcoordinate.levels
            if geometry != geometry_field:
                raise epygramError("All resources must return fields with the same geometry.")
            if vcoordinate != vcoordinate_field:
                raise epygramError("All resources must return fields with the same vertical geometry.")
            if levels != levels_field:
                if vcoordinate.grid is not None:
                    raise epygramError("All resources must return fields with the same vertical geometry.")
                else:
                    if len(levels) != len(levels_field):
                        raise epygramError("All resources must return fields with the same vertical geometry length.")
                    joinLevels = True

        kwargs_vcoord = copy.deepcopy(fieldset[0].geometry.vcoordinate.footprint_as_dict())
        if joinLevels:
            if spectral:
                raise epygramError("Not sure how to merge vertical levels when spectral")
            raise NotImplementedError("a revoir")
            # il faut revoir:
            #  mettre au clair le nb de dimensions de level
            #  utiliser get_level
            if fieldset[0].geometry.datashape['k']:
                shapeH = fieldset[0].getdata().shape[1:]
            else:
                shapeH = fieldset[0].getdata().shape
            newLevels = []
            for k in range(len(levels)):
                levelValue = numpy.ndarray(tuple([len(fieldset)] + list(shapeH)))
                for i in range(len(fieldset)):
                    field = fieldset[validities[sortedValidities[i]]]
                    levelValue[i] = field.geometry.vcoordinate.levels[k]
                newLevels.append(levelValue)
            kwargs_vcoord['levels'] = newLevels
            geometry.vcoordinate = fpx.geometry(**kwargs_vcoord)
        else:
            geometry.vcoordinate = fpx.geometry(**kwargs_vcoord)

        # Other metadata
        kwargs_field = copy.deepcopy(fieldset[0].footprint_as_dict())
        kwargs_field.pop('data')
        kwargs_field['validity'] = FieldValidityList()
        kwargs_field['geometry'] = geometry
        for k, v in kwargs_field.iteritems():
            if type(v) == type(FPDict()):
                kwargs_field[k] = dict(v)
        field = fpx.field(**kwargs_field)
        sameProcesstype = True
        for f in fieldset[1:]:
            kwargs_field = copy.deepcopy(f.footprint_as_dict())
            kwargs_field.pop('data')
            kwargs_field['validity'] = FieldValidityList()
            kwargs_field['geometry'] = geometry
            sameProcesstype = sameProcesstype and  kwargs_field['processtype'] == fieldset[0].processtype
            kwargs_field['processtype'] = fieldset[0].processtype  # we exclude processtype from the comparison
            for k, v in kwargs_field.iteritems():
                if type(v) == type(FPDict()):
                    kwargs_field[k] = dict(v)
            myf = fpx.field(**kwargs_field)
            if field != myf:
                raise epygramError("All fields must be of the same kind")

        # Returned field
        shape = tuple([len(fieldset)] + list(fieldset[0].getdata().shape))
        data = numpy.ndarray(shape)
        fieldvaliditylist = FieldValidityList()
        fieldvaliditylist.pop()
        for i in range(len(fieldset)):
            field = fieldset[validities[sortedValidities[i]]]
            fieldvaliditylist.extend(field.validity)
            data[i] = field.getdata()
        if not geometry.datashape['k']:
            shape = tuple([shape[0]] + [1] + list(shape[1:]))
            data = data.reshape(shape)
        if not sameProcesstype:
            kwargs_field['processtype'] = None

        kwargs_field['validity'] = fieldvaliditylist
        kwargs_field['fid'][self.format] = kwargs_field['fid'][self.lowLevelFormat]
        field = fpx.field(**kwargs_field)
        print data.shape
        field.setdata(data)

        return field
