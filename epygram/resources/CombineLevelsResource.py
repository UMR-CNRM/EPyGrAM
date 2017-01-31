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
import numpy

from footprints import FPDict, proxy as fpx

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
                values=set(['CombineLevels']))
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(Resource, self).__init__(*args, **kwargs)

        if self.resource.openmode != self.openmode:
            raise epygramError("The low level resource must be opened using the same mode as this high-level resource.")

        self.format = "CombineLevels"

    def open(self):
        """Opens the low level resource"""

        if not self.resource.isopen: self.resource.open()

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

        Args: \n
        - *seed*: might be a 'handgrip', i.e. a dict where you can store all
          requested keys,
          e.g. {'shortName':'t', 'indicatorOfTypeOfLevel':'pl', 'level':850},
          a list of handgrips or *None*. If *None* (default), returns the list
          of all fields in resource.
        - *genric*: if True, returns a list of tuples (fid, fid) of
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
            raise epygramError("no field matching '" + str(seed) + \
                               "' was found in resource ")
        if generic:
            fieldslist = [(fieldslist[i], fieldslist[i]) for i in range(len(fieldslist))]

        return fieldslist

    def _create_list(self):
        """Creates the list of available fields associated with the original fids."""

        result = {}

        # Loop over the low level resource fids
        for fid in self.resource.listfields(complete=True):
            assert 'generic' in fid.keys(), \
                   "Not able to combine levels if fields do not have 'generic' fids"
            original_fid = fid[fmtfid(self.resource.format, fid)]
            generic_fid = fid['generic']
            level = generic_fid.pop('level', None)  # we suppress level from generic_fid
            hashable_generic_fid = tuple([(k, generic_fid[k]) for k in sorted(generic_fid.keys())])
            if not hashable_generic_fid in result:
                result[hashable_generic_fid] = {'original_fids':[], 'generic':None}
            result[hashable_generic_fid]['original_fids'].append((original_fid, level))
            result[hashable_generic_fid]['generic'] = generic_fid

        # For fields present on only one layer, we put again the level in the generic fid
        for k, v in result.iteritems():
            if len(v['original_fids']) == 1 and v['original_fids'][0][1] is not None:
                v['generic']['level'] = v['original_fids'][0][1]

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
            if category in sortedfields.keys():
                sortedfields[category].append(field)
            else:
                sortedfields[category] = [field]

        return sortedfields

    def readfield(self, handgrip, getdata=True):
        """
        Read the field in the low level resource and join the levels.
        """

        result = self.readfields(handgrip=handgrip, getdata=getdata)
        if len(result) != 1:
            raise epygramError(str(len(result)) + "field(s) have been found, only one expected.")
        return result[0]

    def readfields(self, handgrip, getdata=True):
        """
        Read the field in the low level resource and join the levels.
        """

        fieldset = FieldSet()
        cont = self._create_list()
        for fid in self.listfields(select=handgrip):
            found = None
            for k, v in cont.iteritems():
                if fid == v['generic']:
                    if found is not None:
                        raise epygramError("Internal error...")
                    else:
                        found = k
            if found is None:
                raise epygramError("Internal error....")
            if not 'parameterNumber' in fid or fid['parameterNumber'] == 255:
                # We do not join levels when parameterNumber is not known
                for original_fid in cont[found]['original_fid']:
                    field = self.resource.readfield(original_fid[0])
                    field.fid[self.format] = fid
                    fieldset.append(field)
            else:
                # We join levels

                # Fields to join
                fields = self.resource.readfields([original_fid[0] for original_fid in cont[found]['original_fids']])

                if len(fields) == 1:
                    field = fields[0]
                    field.fid[self.format] = fid
                    fieldset.append(field)
                else:
                    # Geometry check and level list: part I
                    levellist = {}
                    kwargs_geom = copy.deepcopy(fields[0].geometry.footprint_as_dict())
                    kwargs_vcoord = copy.deepcopy(fields[0].geometry.vcoordinate.footprint_as_dict())
                    kwargs_vcoord['levels'] = [255]
                    kwargs_geom['vcoordinate'] = fpx.geometry(**kwargs_vcoord)
                    for k, v in kwargs_geom.iteritems():
                        if type(v) == type(FPDict()):
                            kwargs_geom[k] = dict(v)
                    ref_geometry = fpx.geometry(**kwargs_geom)
                    spectral = fields[0].spectral
                    # Other metadata check: part I
                    kwargs_field = copy.deepcopy(fields[0].footprint_as_dict())
                    kwargs_field['fid'] = {'generic':fid}
                    kwargs_field['geometry'] = ref_geometry
                    for k, v in kwargs_field.iteritems():
                        if type(v) == type(FPDict()):
                            kwargs_field[k] = dict(v)
                    ref_field = fpx.field(**kwargs_field)
                    for i, field in enumerate(fields):
                        # Geometry check and level list: part II
                        if spectral != field.spectral:
                            raise epygramError("All fields must be gridpoint or spectral")
                        kwargs_geom = copy.deepcopy(field.geometry.footprint_as_dict())
                        kwargs_vcoord = copy.deepcopy(field.geometry.vcoordinate.footprint_as_dict())
                        kwargs_vcoord['levels'] = [255]
                        kwargs_geom['vcoordinate'] = fpx.geometry(**kwargs_vcoord)
                        for k, v in kwargs_geom.iteritems():
                            if type(v) == type(FPDict()):
                                kwargs_geom[k] = dict(v)
                        geometry_field = fpx.geometry(**kwargs_geom)
                        levels_field = field.geometry.vcoordinate.levels
                        if len(levels_field) != 1:
                            raise NotImplementedError("Cannot join multi-levels fields.")
                        levellist[i] = field.geometry.vcoordinate.levels[0]
                        if ref_geometry != geometry_field:
                            raise epygramError("All resources must return fields with the same geometry.")
                        # Other metadata check: part II
                        kwargs_field = copy.deepcopy(field.footprint_as_dict())
                        kwargs_field['fid'] = {'generic':fid}
                        kwargs_field['geometry'] = geometry_field
                        for k, v in kwargs_field.iteritems():
                            if type(v) == type(FPDict()):
                                kwargs_field[k] = dict(v)
                        myf = fpx.field(**kwargs_field)
                        if ref_field != myf:
                            raise epygramError("All fields must be of the same kind")

                    # levels and data of the resulting field
                    levelsorder = sorted(levellist, key=levellist.get)
                    levels4d = ref_field.geometry.get_levels(d4=True, nb_validities=len(fields[0].validity))
                    shape = levels4d.shape  # shape of 4D with only one level
                    shape = tuple([shape[0], len(fields)] + list(shape[2:]))
                    levels = numpy.ndarray(shape, dtype=levels4d.dtype)

                    # Loop over the fields
                    first = True
                    for i in levelsorder:
                        if fields[i].spectral: fields[i].sp2gp()
                        if first:
                            first = False
                            data = fields[i].getdata(d4=True)[:, 0:1, :, :]
                        else:
                            if isinstance(data, numpy.ma.masked_array):
                                concat = numpy.ma.concatenate
                            else:
                                concat = numpy.concatenate
                            data = concat((data, fields[i].getdata(d4=True)[:, 0:1, :, :]), axis=1)
                        levels[:, i] = fields[i].geometry.get_levels(d4=True, nb_validities=len(fields[0].validity))[:, 0]
                    cst_horizontal = True
                    cst_time = True
                    for k in range(levels.shape[1]):
                        for t in range(levels.shape[0]):
                            cst_time = cst_time and numpy.all(levels[t, k] == levels[0, k])
                            cst_horizontal = cst_horizontal and \
                                             cst_time and \
                                             numpy.all(levels[t, k] == levels[t, k].flatten()[0])
                    if cst_time:
                        levels = levels[0]
                    if cst_horizontal:  # implies cst_time, so levels does not have time dimension
                        levels = [levels[k].flatten()[0] for k in range(levels.shape[0])]
                    kwargs_vcoord['levels'] = levels
                    kwargs_geom['vcoordinate'] = fpx.geometry(**kwargs_vcoord)
                    kwargs_geom['structure'] = {'H2D':'3D', 'point':'V1D', 'H1D':'V2D', 'V1D':'V1D', 'V2D':'V2D'}[kwargs_geom['structure']]
                    kwargs_field['geometry'] = fpx.geometry(**kwargs_geom)
                    kwargs_field['structure'] = kwargs_geom['structure']
                    kwargs_field['spectral_geometry'] = None
                    new_field = fpx.field(**kwargs_field)
                    new_field.setdata(data)
                    new_field.fid[self.format] = fid
                    fieldset.append(new_field)
        return fieldset

    def writefield(self, *args, **kwargs):
        """Write fields."""
        raise NotImplementedError("writefield is not implemented")
        # To implement writefield we need to check that each validity is affected to only one resource

    def extractprofile(self, *args, **kwargs):
        """
        Extracts profiles.
        """

        return self.resource.extractprofile(*args, **kwargs)

    def extractsection(self, *args, **kwargs):
        """
        Extracts sections.
        """

        return self.resource.extractsection(*args, **kwargs)

    @property
    def spectral_geometry(self):
        """
        Returns the spectral_geometry
        """

        return self.resource.spectral_geometry
