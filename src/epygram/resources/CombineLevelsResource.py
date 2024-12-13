#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a CombineLevelsResource.
This resource exposes 3D fields when the low level resource only expose horizontal fields.
"""

import copy

from footprints import proxy as fpx
from footprints import FPDict

from epygram.base import Resource, FieldSet
from epygram.util import fmtfid
from epygram import epygramError


class CombineLevelsResource(Resource):
    """Class implementing a CombineLevelsResource."""

    _collector = ('resource_modificator', 'epyresource')
    _footprint = dict(
        attr=dict(
            resource=dict(
                type=Resource,
                info="Low level resource"),
            name=dict(
                values=set(['CombineLevels'])),
            virtual=dict(
                info="Must readfield return a virtual field",
                type=bool,
                optional=True,
                default=False)
        )
    )

    def __init__(self, *args, **kwargs):
        super(Resource, self).__init__(*args, **kwargs)

        if self.resource.openmode != self.openmode:
            raise epygramError("The low level resource must be opened using the same mode as this high-level resource.")

        self.format = "CombineLevels"

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
        except Exception:
            pass

    def find_fields_in_resource(self, seed=None, generic=False):
        """
        Returns a list of the fields from resource whose name match the given
        seed.

        :param seed: might be a 'handgrip', i.e. a dict where you can store all
          requested keys,
          e.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl', 'level':850},
          a list of handgrips or *None*. If *None* (default), returns the list
          of all fields in resource.
        :param generic: if True, returns a list of tuples (fid, fid) of
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

    def _create_list(self):
        """Creates the list of available fields associated with the original fids."""
        result = {}
        # Loop over the low level resource fids
        for fid in self.resource.listfields(complete=True):
            assert 'generic' in fid, \
                   "Not able to combine levels if fields do not have 'generic' fids"
            original_fid = fid[fmtfid(self.resource.format, fid)]
            generic_fid = copy.deepcopy(fid['generic'])
            level = generic_fid.pop('level', None)  # we suppress level from generic_fid
            scaledValueOfFirstFixedSurface  = generic_fid.pop('scaledValueOfFirstFixedSurface', None)
            hashable_generic_fid = tuple([(k, generic_fid[k]) for k in sorted(generic_fid.keys())])
            if hashable_generic_fid not in result:
                result[hashable_generic_fid] = {'original_fids':[], 'generic':None}
            result[hashable_generic_fid]['original_fids'].append((original_fid, level, scaledValueOfFirstFixedSurface))
            result[hashable_generic_fid]['generic'] = generic_fid
        # For fields present on only one layer, we put again the level in the generic fid
        for k, v in result.items():
            if len(v['original_fids']) == 1:
                if v['original_fids'][0][1] is not None:
                    v['generic']['level'] = v['original_fids'][0][1]
                if v['original_fids'][0][2] is not None:
                    v['generic']['scaledValueOfFirstFixedSurface'] = v['original_fids'][0][2]
        return result

    def listfields(self, onlykey=None, select=None, complete=False):
        """Lists the available fields."""
        fidlist = [v['generic'] for v in self._create_list().values()]
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

        :param sortingkey: sorting key
        :param onlykey: can be specified as a string or a tuple of strings,
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
            if category in sortedfields:
                sortedfields[category].append(field)
            else:
                sortedfields[category] = [field]
        return sortedfields

    def readfield(self, handgrip, getdata=True):
        """
        Read the field in the low level resource and join the levels.

        :param handgrip: identification of the field
        :param getdata: if False, do not read data but only metadata
        """
        result = self.readfields(handgrip=handgrip, getdata=getdata)
        if len(result) != 1:
            raise epygramError(str(len(result)) + "field(s) have been found, only one expected.")
        return result[0]

    def readfields(self, handgrip, getdata=True):
        """
        Read the field in the low level resource and join the levels.

        :param handgrip: identification of the field
        :param getdata: if False, do not read data but only metadata
        """
        fieldset = FieldSet()
        cont = self._create_list()
        for fid in self.listfields(select=handgrip):
            fid = FPDict(fid)  # footprints purpose
            found = None
            for k, v in cont.items():
                if fid == v['generic']:
                    if found is not None:
                        raise epygramError("Internal error...")
                    else:
                        found = k
            if found is None:
                raise epygramError("Internal error....")
            if 'parameterNumber' not in fid or fid['parameterNumber'] == 255:
                # We do not join levels when parameterNumber is not known
                print("readfields", handgrip)
                for original_fid in cont[found]['original_fids']:
                    field = self.resource.readfield(original_fid[0], getdata=getdata)
                    field.fid[self.format] = fid
                    fieldset.append(field)
            else:
                fidList = [original_fid[0] for original_fid in cont[found]['original_fids']]

                if len(fidList) == 1:
                    field = self.resource.readfield(fidList[0], getdata=getdata)
                    field.fid[self.format] = fid
                    fieldset.append(field)
                else:
                    structure = self.resource.readfield(fidList[0], getdata=False).structure
                    virtualField = fpx.field(fid={'generic':fid, self.format:fid},
                                             structure={'H2D':'3D',
                                                        'Point':'V1D',
                                                        'H1D':'V2D',
                                                        'V1D':'V1D',
                                                        'V2D':'V2D',
                                                        '3D':'3D'}[structure],
                                             resource=self.resource,
                                             resource_fids=fidList)
                    if self.virtual:
                        fieldset.append(virtualField)
                    else:
                        fieldset.append(virtualField.as_real_field(getdata=getdata))
        return fieldset

    def writefield(self, *args, **kwargs):
        """Write fields."""
        raise NotImplementedError("writefield is not implemented")
        # To implement writefield we need to check that each validity is affected to only one resource

    def extractprofile(self, *args, **kwargs):
        """Extracts profiles."""
        profile = self.resource.extractprofile(*args, **kwargs)
        profile.fid[self.format] = profile.fid[self.resource.format]
        return profile

    def extractsection(self, *args, **kwargs):
        """Extracts sections."""
        section = self.resource.extractsection(*args, **kwargs)
        section.fid[self.format] = section.fid[self.resource.format]
        return section

    @property
    def spectral_geometry(self):
        """Returns the spectral_geometry"""
        return self.resource.spectral_geometry
