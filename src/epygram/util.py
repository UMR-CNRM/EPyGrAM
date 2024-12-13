#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some useful utilities...
"""

import math
import copy
import numpy
import sys
import datetime
import hashlib
import os
from urllib.request import urlopen

from footprints import FootprintBase
from bronx.graphics.colormapping import (add_cmap,
                                         get_norm4colorscale)
from bronx.syntax.decorators import nicedeco
from bronx.fancies import loggers

from . import config, epygramError
from .colormapping import register_colormap_from_json

epylog = loggers.getLogger(__name__)

class CheckAttribute(object):
    def add_attr(self, name, set_check, error_message, extra=None):
        """
        :param name: name of the attribute
        :param set_check: function called each time the attribute is modified to check
                          the validity of new value. This function takes two arguments:
                            - the value to test
                            - the extra parameter
                          To enable pickle to be used, functions must be named (no lambda)
                          and not local (nested functions).
        :param error_message: message to display in assertion when value is not valid
        :param extra: an argument to provide to the test function

        The method adds a control on the authorised value of an attribute
        """
        assert isinstance(name, str), "*name* must be a string"
        assert callable(set_check), "*set_check* must be callable"
        assert isinstance(error_message, str), "*error_message* must be a string"
        if not hasattr(self, '_epyattr'):
            super().__setattr__('_epyattr', {})
        self._epyattr[name] = (set_check, error_message, extra)

    def add_attr_int(self, name):
        """
        :param name: name of the attribute
        The method adds a control to check if the attribute is an integer
        """
        return self.add_attr_class(name, int)

    def add_attr_str(self, name):
        """
        :param name: name of the attribute
        The method adds a control to check if the attribute is a string
        """
        return self.add_attr_class(name, str)

    def add_attr_float(self, name):
        """
        :param name: name of the attribute
        The method adds a control to check if the attribute is a float
        """
        return self.add_attr_class(name, float)

    def add_attr_list(self, name):
        """
        :param name: name of the attribute
        The method adds a control to check if the attribute is a list
        """
        return self.add_attr_class(name, list)

    def add_attr_dict(self, name):
        """
        :param name: name of the attribute
        The method adds a control to check if the attribute is a dict
        """
        return self.add_attr_class(name, dict)

    @staticmethod
    def _inlist(x, values): return x in values
    def add_attr_inlist(self, name, values):
        """
        :param name: name of the attribute
        :param values: authorised values
        The method adds a control to check if the attribute value is allowed
        """
        assert isinstance(values, list), "*values* must be a list"
        return self.add_attr(name, self._inlist,
                             "*{name}* (={value}) must be a among {values}".format(name=name,
                                                                                   values=str(values),
                                                                                   value='{value}'),
                             values)

    def add_attr_class(self, name, cls):
        """
        :param name: name of the attribute
        :param cls: authorised class
        The method adds a control to check if the attribute is an instance of *cls*
        """
        return self.add_attr(name, isinstance,
                             "*{name}* must be a {clsc} instance (and not a {cls})".format(name=name,
                                                                                           clsc=cls.__name__,
                                                                                           cls='{cls}'),
                             cls)

    def __setattr__(self, name, value):
        if not hasattr(self, '_epyattr'):
            super().__setattr__('_epyattr', {})
        if name in self._epyattr:
            assert self._epyattr[name][0](value, self._epyattr[name][2]), \
                   self._epyattr[name][1].format(cls=type(value).__name__, value=value)
        super().__setattr__(name, value)


class RecursiveObject(CheckAttribute):
    """
    Generic abstract class implementing useful recursive properties:

    - display of object: the *__str__* method returns str(attr) for each of the
      object's attributes, with automatical indentation for readability.

    - test of (in)equality: the a == b test will be true if a and b have the
      same attributes and *a.attr == b.attr* for each attribute.
    """

    # ghost attributes are ignored when comparing 2 objects between them
    _ghost_attributes = ['_epyattr']

    def _strItem(self, item, reclevel=1):
        """Recursive display of object attributes."""

        offset = "".rjust(reclevel * 4)
        itemstring = ""
        if isinstance(item, FootprintBase):
            itemstring += item.__class__.__name__ + " containing:"
            for attr in item.__dict__.keys():
                if attr == '_attributes':
                    for i in item.__dict__[attr].keys():
                        itemstring += "\n" + offset + i + ": " + self._strItem(item.__dict__[attr][i], reclevel + 1)
                elif attr not in ('_puredict', '_observer', '_epyattr'):
                    itemstring += "\n" + offset + attr + ": " + self._strItem(item.__dict__[attr], reclevel + 1)
        elif isinstance(item, list):
            if len(item) > 0:
                if isinstance(item[0], RecursiveObject):
                    itemstring = "[" + repr(item[0]) + ", ... (" + \
                                 str(len(item)) + " objects)]"
                else:
                    itemstring = str(numpy.array(item))
            else:
                itemstring = "[]"
        elif isinstance(item, RecursiveObject):
            itemstring += item.__class__.__name__ + " containing:"
            for attr in item.__dict__.keys():
                if attr != '_epyattr':
                    itemstring += "\n" + offset + attr + ": " + self._strItem(item.__dict__[attr], reclevel + 1)
        elif isinstance(item, dict):
            for key in sorted(item.keys()):
                itemstring += "\n" + offset + str(key) + ": " + self._strItem(item[key], reclevel + 1)
        else:
            itemstring = str(item)
        return itemstring

    def __str__(self):
        """
        Recursive display of object: displays each of its attributes
        with indentation.
        """
        return self._strItem(self)

    def __eq__(self, other):
        """Test of equality by recursion on the object's attributes."""
        return self.tolerant_equal(other, tolerance=0.)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # known issue __eq__/__hash__ must be defined both or none, else inheritance is broken
        return object.__hash__(self)

    def tolerant_equal(self, other, tolerance=config.epsilon):
        """
        Test of equality by recursion on the object's attributes,
        with a **tolerance**.
        """
        ok = True
        if self.__class__ == other.__class__:
            self_attrs = set([k for k in self.__dict__.keys() if k not in self._ghost_attributes])
            other_attrs = set([k for k in other.__dict__.keys() if k not in other._ghost_attributes])
            if self_attrs != other_attrs:
                ok = False
            else:
                for attr in self_attrs:
                    if not Comparator.are_equal(self.__dict__[attr],
                                                other.__dict__[attr],
                                                tolerance):
                        ok = False
                        #print(attr)
                        #print(ok)
                        #print(Comparator.diff(self.__dict__[attr],
                        #                      other.__dict__[attr]))
                        break
        else:
            ok = False
        return ok

    def copy(self):
        """Returns a copy of the object."""
        return copy.copy(self)

    def deepcopy(self):
        """Returns a deepcopy of the object."""
        return copy.deepcopy(self)

    def recursive_diff(self, other):
        """Recursively list what differs from **other**."""
        if not self.__eq__(other):
            diff = {}
            if self.__class__ == other.__class__:
                if isinstance(self, FootprintBase):
                    self_dict = self._attributes
                    other_dict = other._attributes
                else:
                    self_dict = {k:v for k,v in self.__dict__.items() if k not in self._ghost_attributes}
                    other_dict = {k:v for k,v in other.__dict__.items() if k not in self._ghost_attributes}
                diff = Comparator.diff(self_dict, other_dict)
            else:
                diff['__class__'] = (str(self.__class__), str(other.__class__))
            return diff


class Comparator(object):
    """Helper to recursively compare objects."""

    @classmethod
    def _float_are_equal(cls, float1, float2, tolerance):
        # tolerance for floats
        return nearlyEqual(float1, float2, tolerance)

    @classmethod
    def _array_are_equal(cls, array1, array2, tolerance):
        if tolerance == 0.:
            return numpy.all(array1 == array2)
        else:
            if (array1.dtype == array2.dtype and
                    array1.dtype in [numpy.dtype(d)
                                     for d in ['float16', 'float32', 'float64']]):
                # tolerance for floats
                return numpy.all(nearlyEqualArray(array1, array2, tolerance))
            else:
                return numpy.all(array1 == array2)

    @classmethod
    def _dict_are_equal(cls, dict1, dict2, tolerance):
        if set(dict1.keys()) == set(dict2.keys()):
            ok = True
            for k in dict1.keys():
                if not cls.are_equal(dict1[k], dict2[k], tolerance):
                    # print("------------not equal", k)
                    ok = False
                    break
            return ok
        else:
            return False

    @classmethod
    def _dict_diff(cls, dict1, dict2):
        if not cls.are_equal(dict1, dict2, 1e-17):
            diff = {}
            keys1 = set(dict1.keys())
            keys2 = set(dict2.keys())
            for k in keys1.intersection(keys2):
                if not cls.are_equal(dict1[k], dict2[k], 1e-17):
                    diff[k] = cls.diff(dict1[k], dict2[k])
            for k in keys1.difference(keys2):
                diff[k] = (dict1[k], None)
            for k in keys2.difference(keys1):
                diff[k] = (None, dict2[k])
            return diff

    @classmethod
    def _list_are_equal(cls, list1, list2, tolerance):
        if len(list1) != len(list2):
            return False
        else:
            ok = True
            for i in range(len(list1)):
                if not cls.are_equal(list1[i], list2[i], tolerance):
                    ok = False
                    break
            return ok

    @classmethod
    def are_equal(cls, obj1, obj2, tolerance=0.):
        """Checks equality of objects."""
        if isinstance(obj1, float) and isinstance(obj2, float):
            return cls._float_are_equal(obj1, obj2, tolerance)
        elif isinstance(obj1, numpy.ndarray) and isinstance(obj2, numpy.ndarray):
            return cls._array_are_equal(obj1, obj2, tolerance)
        elif isinstance(obj1, dict) and isinstance(obj2, dict):
            return cls._dict_are_equal(obj1, obj2, tolerance)
        elif isinstance(obj1, (list, tuple)) and isinstance(obj2, (list, tuple)):
            return cls._list_are_equal(obj1, obj2, tolerance)
        elif isinstance(obj1, RecursiveObject) and isinstance(obj2, RecursiveObject):
            return obj1.tolerant_equal(obj2, tolerance)
        else:
            return obj1 == obj2

    @classmethod
    def diff(cls, obj1, obj2):
        """Inspect differences between objects."""
        if not cls.are_equal(obj1, obj2, 1e-17):
            if isinstance(obj1, RecursiveObject) and isinstance(obj2, RecursiveObject):
                return obj1.recursive_diff(obj2)
            elif isinstance(obj1, dict) and isinstance(obj2, dict):
                return cls._dict_diff(obj1, obj2)
            else:
                return obj1, obj2


class Angle(RecursiveObject):
    """
    This class handles an angle.
    It enables conversions of units, while saving the original unit and value
    of the angle at its construction.

    Available units: 'degrees', 'radians', 'cos_sin' (cos, sin),
                     'DMS' (Degree, Minutes, Seconds).
    """

    deg = 'degrees'
    rad = 'radians'
    trig = 'cos_sin'
    dms = 'DMS'
    units = set([deg, dms, rad, trig])

    def __init__(self, value, unit):
        """
        Constructor.
        'unit' argument must be among:
        - 'degrees',
        - 'DMS' - in which case, value is a tuple (degrees, minutes, seconds),
        - 'radians',
        - 'cos_sin' - in which case, value is a tuple (cos, sin).
        """
        if unit in Angle.units:
            self.__dict__['_' + unit] = value
            if unit in ('degrees', 'radians'):
                # get modulo around 0.
                if unit == 'degrees':
                    circle = 360.
                elif unit == 'radians':
                    circle = 2. * math.pi
                while self.__dict__['_' + unit] > circle / 2.:
                    self.__dict__['_' + unit] -= circle
                while self.__dict__['_' + unit] < -circle / 2.:
                    self.__dict__['_' + unit] += circle
        else:
            raise ValueError("this angle unit is not implemented: " + str(unit))
        self._origin_unit = unit
        self._origin_value = value

    def __hash__(self):
        return hash(self._origin_value) + hash(self._origin_unit)

    def __mul__(self, factor):
        return Angle(self.get('radians') * factor, 'radians')

    def __div__(self, factor):
        return Angle(self.get('radians') / factor, 'radians')

    def __add__(self, other):
        assert isinstance(other, Angle)
        return Angle(self.get('radians') + other.get('radians'), 'radians')

    def __sub__(self, other):
        assert isinstance(other, Angle)
        return Angle(self.get('radians') - other.get('radians'), 'radians')

    def get(self, unit=None):
        """
        Returns the angle in the requested unit.
        If no unit is supplied, the origin unit is used.
        """
        if unit is None:
            unit = self._origin_unit  # or a default one ?
        elif unit in Angle.units:
            if not '_' + unit in self.__dict__:
                self._compute(unit)
        else:
            raise ValueError("this angle unit is not implemented: " + str(unit))
        return self.__dict__['_' + unit]

    def _compute(self, unit):
        """
        Compute the angle in the requested unit, from the original one.
        See constructor for more details about units.
        """
        # conversion to degrees
        if unit == Angle.deg:
            if self._origin_unit == Angle.rad:
                self.__dict__['_' + unit] = math.degrees(self._origin_value)
            elif self._origin_unit == Angle.trig:
                self.__dict__['_' + unit] = math.degrees(
                    math.copysign(math.acos(self._origin_value[0]),
                                  self._origin_value[1]))
            elif self._origin_unit == Angle.dms:
                self.__dict__['_' + unit] = (self._origin_value[0] +
                                             self._origin_value[1] / 60. +
                                             self._origin_value[2] / 3600.)
            else:
                raise NotImplementedError("conversion from this unit (" +
                                          self._origin_unit +
                                          ") is not coded.")

        # conversion to radians
        elif unit == Angle.rad:
            if self._origin_unit == Angle.deg:
                self.__dict__['_' + unit] = math.radians(self._origin_value)
            elif self._origin_unit == Angle.trig:
                self.__dict__['_' + unit] = math.copysign(
                    math.acos(self._origin_value[0]),
                    self._origin_value[1])
            elif self._origin_unit == Angle.dms:
                self.__dict__['_' + unit] = math.radians(self._origin_value[0] +
                                                         self._origin_value[1] / 60. +
                                                         self._origin_value[2] / 3600.)
            else:
                raise NotImplementedError("conversion from this unit (" +
                                          self._origin_unit +
                                          ") is not coded.")

        # conversion to (cos, sin)
        elif unit == Angle.trig:
            if self._origin_unit == Angle.deg:
                self.__dict__['_' + unit] = (math.cos(math.radians(self._origin_value)),
                                             math.sin(math.radians(self._origin_value)))
            elif self._origin_unit == Angle.rad:
                self.__dict__['_' + unit] = (math.cos(self._origin_value),
                                             math.sin(self._origin_value))
            elif self._origin_unit == Angle.dms:
                anglerad = math.radians(self._origin_value[0] +
                                        self._origin_value[1] / 60. +
                                        self._origin_value[2] / 3600.)
                self.__dict__['_' + unit] = (math.cos(anglerad),
                                             math.sin(anglerad))
            else:
                raise NotImplementedError("conversion from this unit (" +
                                          self._origin_unit +
                                          ") is not coded.")

        # conversion to (degrees, minutes, seconds)
        elif unit == Angle.dms:
            if self._origin_unit == Angle.deg:
                decdeg = self._origin_value
            elif self._origin_unit == Angle.rad:
                decdeg = math.degrees(self._origin_value)
            elif self._origin_unit == Angle.trig:
                decdeg = math.degrees(math.copysign(math.acos(self._origin_value[0]),
                                                    self._origin_value[1]))
            else:
                raise NotImplementedError("conversion from this unit (" +
                                          self._origin_unit +
                                          ") is not coded.")
            sign = int(math.copysign(1, decdeg))
            decdeg = decdeg * sign
            degrees = int(decdeg)
            decmin = (decdeg - degrees) * 60.
            minutes = int(decmin)
            seconds = (decmin - minutes) * 60
            self.__dict__['_' + unit] = (degrees * sign, minutes, seconds)

        else:
            raise NotImplementedError("conversion to this unit (" + unit +
                                      ") is not coded.")

    def tolerant_equal(self, other, tolerance=config.epsilon):
        """
        Redefinition because of dynamism of buffering new computed values...
        """
        if not isinstance(other, Angle):
            return False
        if abs(self.get('degrees') -
               degrees_nearest_mod(other.get('degrees'),
                                   self.get('degrees'))) <= tolerance:
            ok = True
        else:
            ok = False
        return ok

    def recursive_diff(self, other):
        if self != other:
            return {'_origin_unit': Comparator.diff(self._origin_unit, other._origin_unit),
                    '_origin_value': Comparator.diff(self._origin_value, other._origin_value)}


# FUNCTIONS #
#############
def is_scalar(x):
    """
    Returns True if argument is scalar
    """
    return numpy.ndim(x) == 0


def as_numpy_array(x):
    """
    Returns a numpy (ma) array
    """
    if isinstance(x, list):
        missing = set([item.fill_value for item in x if isinstance(item, numpy.ma.masked_array)])
        if len(missing) > 0:
            return numpy.ma.atleast_1d(numpy.ma.asanyarray(x, fill_value=missing.pop() if len(missing) == 1 else None))
        else:
            return numpy.atleast_1d(numpy.asanyarray(x))
    else:
        return numpy.atleast_1d(numpy.asanyarray(x))


def find_re_in_list(regexp, a_list):
    """
    Finds all elements from a list that match a regular expression.
    The regexp and the different elements of the list must be of the same type:

    - strings
    - tuples with the same length
    - dictionnaries: all regexp keys must be keys of the list
    """

    def check_string_pattern(pattern, element):
        import re
        if not isinstance(pattern, str) or \
                not isinstance(element, str):
            raise epygramError("pattern and element must be strings in \
                                check_string_pattern function.")
        # protect '.'
        mypattern = re.subn('\.', r'\.', pattern)[0]
        # change unix '?' to python '.' (any char)
        mypattern = mypattern.replace('?', '.')
        # change unix '*' to python '.*' (several any char)
        mypattern = mypattern.replace('*', '.*')
        mypattern += '(?!.)'
        return re.match(mypattern, element.strip())

    found = []
    if isinstance(regexp, str):
        for field in a_list:
            if check_string_pattern(regexp, str(field)):
                found.append(field)
    elif isinstance(regexp, tuple):
        for field in a_list:
            if not isinstance(field, tuple):
                raise epygramError("pattern and elements of the list must be of\
                                    the same type (tuples here).")
            if len(regexp) != len(field):
                raise epygramError("pattern and elements of the list must be\
                                    tuples of the same length.")
            ok = True
            for i in range(len(regexp)):
                ok = ok and check_string_pattern(str(regexp[i]), str(field[i]))
            if ok:
                found.append(field)
    elif isinstance(regexp, dict):
        raise NotImplementedError("the dictionnary type of regexp is not yet\
                                   implemented.")
    else:
        raise NotImplementedError("this type of regexp is not (yet?)\
                                   implemented.")
    return found


def degrees_nearest_mod(d, ref):
    """Returns the angle(s) **d** in the modulo nearest to **ref**."""
    try:
        _ = len(d)
        scalar = False
    except Exception:
        d = [d]
        scalar = True
    d = numpy.array(d)
    d_inf = d - 360.
    d_sup = d + 360.
    result = d_sup
    mask = numpy.logical_and(numpy.abs(d - ref) <= numpy.abs(d_sup - ref),
                             numpy.abs(d - ref) <= numpy.abs(d_inf - ref))
    result[mask] = d[mask]
    mask = numpy.logical_and(numpy.abs(d_inf - ref) <= numpy.abs(d_sup - ref),
                             numpy.abs(d_inf - ref) <= numpy.abs(d - ref))
    result[mask] = d_inf[mask]

    if scalar:
        result = result[0]
    return result


def positive_longitudes(lons):
    """Returns numpy array *lons* forced (modulo) in [0;360[."""
    positive = lons + 360.
    mask = lons >= 0.
    positive[mask] = lons[mask]
    return positive


def longitudes_between_minus180_180(lons):
    """Returns numpy array *lons* forced (modulo) in ]-180;180]."""
    newlons = positive_longitudes(lons)
    negative_lons = newlons - 360.
    mask = newlons > 180.
    newlons[mask] = negative_lons[mask]
    return newlons


def positive_longitude(lon, unit='degrees'):
    """Returns *lon* shifted in [0;360[ or [0;2pi[ (depending on *unit*)."""
    if lon < 0.:
        if unit == 'degrees':
            lon += 360.
        elif unit == 'radians':
            lon += 2. * numpy.pi
        else:
            raise NotImplementedError()
    return lon


def get_file(url, filename, authorize_cache=True, subst=None):
    """
    Get file from url into filename.
    If authorize_cache is True and a directory is set for caching
    in user preferences, file is first searched in cache directory,
    downloaded in the cache directory if missing, then filename is
    built as a symlink to the file in cache directory.
    This way, filename must be always deleted by the caller.

    :param url: url to get
    :param filename: filename in which to put the result
    :param authorize_cache: authorize to use cache (if possible)
                            for this request
    :param subst: dictionary whose keys are searched in url to
                  be replaced by corresponding key

    example: url='https://a.tile.openstreetmap.org/${z}/${x}/${y}.png'
             subst={'${z}': 4, '${x}': 8, '${y}': 5}
             will get the file https://a.tile.openstreetmap.org/4/8/5.png
    """

    def md5(s):
        h = hashlib.md5()
        h.update(s)
        return h.hexdigest()

    # Substitution
    actual_url = url
    if subst is not None:
        for k, v in subst.items():
            actual_url = actual_url.replace(str(k), str(v))

    if authorize_cache and config.internet_cache_dir is not None:
        # Corresponding file name in cache
        url_hash = md5(url.encode('UTF8'))
        actual_url_hash = md5(actual_url.encode('UTF8'))
        if url_hash != actual_url_hash:
            directory = os.path.join(config.internet_cache_dir,
                                     url_hash)
            if not os.path.exists(directory):
                os.mkdir(directory)
            cache_filename = os.path.join(directory, actual_url_hash)
        else:
            cache_filename = os.path.join(config.internet_cache_dir,
                                          url_hash)

        # Cache filling
        if not os.path.exists(cache_filename):
            get_file(actual_url, cache_filename, authorize_cache=False, subst=None)

        # Symlink
        if os.path.exists(filename):
            os.remove(filename)
        os.symlink(cache_filename, filename)

    else:
        with open(filename, 'wb') as f:
            f.write(urlopen(actual_url).read())


def load_cmap(cmap):
    """
    Reads and registers the epygram-or-user-colormap called *cmap*,
    which must be either in config.epygram_colormaps or
    config.usercolormaps.

    Works with both old-way (.cmap) and new way (.json).
    """
    import matplotlib.pyplot as plt
    if cmap not in plt.colormaps() and cmap in config.colormaps:
        filename = config.colormaps[cmap]
        if filename.endswith('.json'):
            return register_colormap_from_json(filename)
        else:
            epylog.warning(_deprecated_cmap)
            with open(filename, 'r') as ocm:
                add_cmap(cmap, ocm)


_deprecated_cmap = ' '.join(["the use of '.cmap' user colormaps is deprecated,",
                             "(and not possible with cartoplot());",
                             "move to json format, using epygram.moves.cmap2json()"])

formatting_default_widths = (50, 20)
separation_line = '{:-^{width}}'.format('', width=sum(formatting_default_widths) + 1) + '\n'


def write_formatted(dest, label, value,
                    align='<',
                    firstcolumn_width=formatting_default_widths[0],
                    secondcolumn_width=formatting_default_widths[1]):
    dest.write(('{:' + align + '{width}}').format(label, width=firstcolumn_width) +
               ':' +
               '{:>{width}}'.format(str(value), width=secondcolumn_width) +
               '\n')


def write_formatted_fields(dest, label, value=None,
                           compression=None,
                           firstcolumn_width=formatting_default_widths[0],
                           secondcolumn_width=formatting_default_widths[1]):
    if compression is None:
        if value is None:
            dest.write('{:<{width}}'.format(label, width=20) +
                       '\n')
        else:
            dest.write('{:<{width}}'.format(label, width=20) +
                       ':' +
                       '{:^{width}}'.format(str(value), width=10) +
                       '\n')
    else:
        line = '{:<{width}}'.format(label, width=20) + \
               ':' + \
               '{:^{width}}'.format(str(value), width=10) + \
               ':' + \
               compression + \
               '\n'
        dest.write(line)


def write_formatted_dict(dest, fid, sort_function=sorted):
    name = fid.pop('name')
    dest.write('name: ' + name + '\n')
    for k in sort_function(fid.keys()):
        dest.write('  ' + str(k) + ': ' + str(fid[k]) + '\n')


def write_formatted_table(dest, table,
                          alignments=['<', '^'], precision=6, float_type='E'):
    """
    A table is meant to be :
    <str> <str> <str> ...
    <str> <num> <num> ...
     ...   ...   ...  ...
    """
    float_style = ' {: .{precision}{type}} '
    array = numpy.array(table, dtype=str)
    columns_dimension = [len(c) + 2 for c in array[0, :]]
    columns_dimension[0] = max([len(c) + 1 for c in array[:, 0]])
    float_len = len(float_style.format(2. / 3.,
                                       type=float_type,
                                       precision=precision)) + 1  # +1 for x < 1E-100
    columns_dimension[1:] = [max(float_len, c) for c in columns_dimension[1:]]

    for i in range(array.shape[0]):
        elements = []
        for elem in table[i][1:]:
            if isinstance(elem, str):
                elements.append(elem)
            elif isinstance(elem, float) or isinstance(elem, int):
                elements.append(float_style.format(elem,
                                                   precision=precision,
                                                   type=float_type))
            else:
                elements.append('-')
        line = ('{:' + alignments[0] + '{width}}').format(table[i][0], width=columns_dimension[0])
        line += ''.join([('{:' + alignments[1] + '{width}}').format(elements[j], width=columns_dimension[j + 1])
                         for j in range(len(elements))])
        dest.write(line + '\n')


def auto_meridians_parallels(geometry,
                             meridians='auto',
                             parallels='auto',
                             extent='focus'):
    """
    Compute meridians and parallels.

    *meridians* and *parallels* enable to fine-tune the choice of lines to
    plot, with either:
      - 'auto': automatic scaling to the map extents
      - 'default': range(0,360,10) and range(-90,90,10)
      - a list of values
      - a grid step, e.g. 5 to plot each 5 degree.
      - None: no one is plot

    :param extent: among 'focus' or 'global', used to determine outer limits
    """

    def Delta2delta(Delta):
        if Delta <= 10:
            delta = 1
        elif 10 < Delta <= 30:
            delta = 5
        elif 30 < Delta <= 180:
            delta = 10
        else:
            delta = 20
        return delta

    # meridians    # TODO: stereopol : all ! + parallels enough ?
    if meridians is None or meridians == 0.:
        meridians = []
    elif meridians == 'default' or (meridians == 'auto' and
                                    ('gauss' in geometry.name or
                                     geometry.name == 'polar_stereographic')):
        meridians = numpy.arange(0, 370, 10)
    elif meridians == 'auto' or isinstance(meridians, int) or isinstance(meridians, float):
        if extent == 'focus':
            minmax = geometry.minmax_ll()
        else:
            minmax = {'lonmax': 180, 'lonmin': -180, 'latmax': 90, 'latmin': -90}
        if meridians == 'auto':
            delta_lon = Delta2delta(minmax['lonmax'] - minmax['lonmin'])
        else:
            delta_lon = float(meridians)
        lonmin = minmax['lonmin'] - minmax['lonmin'] % delta_lon
        if numpy.isinf(minmax['lonmax']):
            lonmax = lonmin + 360.
        else:
            lonmax = minmax['lonmax'] - minmax['lonmax'] % delta_lon + 2 * delta_lon
        lonmax = min(lonmax, lonmin + 360.)  # space-view geometry used to return 1.E30 as lonmax
        meridians = numpy.arange(lonmin, lonmax + delta_lon, delta_lon)
        if max(meridians) > 180.:
            meridians = meridians - 180.  # FIXME: cartopy does not plot meridians > 180°
    # parallels
    if parallels is None or parallels == 0.:
        parallels = []
    elif parallels == 'default' or ('gauss' in geometry.name and parallels == 'auto'):
        parallels = numpy.arange(-90, 100, 10)
    elif parallels == 'auto' or isinstance(parallels, int) or isinstance(parallels, float):
        if extent == 'focus':
            minmax = geometry.minmax_ll()
        else:
            minmax = {'lonmax': 180, 'lonmin': -180, 'latmax': 90, 'latmin': -90}
        if parallels == 'auto':
            delta_lat = Delta2delta(minmax['latmax'] - minmax['latmin'])
        else:
            delta_lat = float(parallels)
        latmin = minmax['latmin'] - minmax['latmin'] % delta_lat
        if latmin <= -90:
            latmin += delta_lat
        if numpy.isinf(minmax['latmax']):
            latmax = 90.
        else:
            latmax = minmax['latmax'] - minmax['latmax'] % delta_lat + 2 * delta_lat
        latmax = min(latmax, 90.)  # space-view geometry used to return 1.E30 as latmax
        if latmax >= 90:
            latmax -= delta_lat
        parallels = numpy.arange(latmin, latmax + delta_lat, delta_lat)
    return list(meridians), list(parallels)


def nearlyEqual(a, b, epsilon=config.epsilon):
    """
    Function to compare floats
    http://floating-point-gui.de/errors/comparison/
    Float.MIN_NORMAL was replaced by sys.float_info.min
    Float.MAX_VALUE was replaced by sys.float_info.max
    """
    absA = numpy.abs(a)
    absB = numpy.abs(b)
    diff = numpy.abs(a - b)

    if a == b:  # shortcut, handles infinities
        return True
    elif a == 0 or b == 0 or diff < sys.float_info.min:
        # a or b is zero or both are extremely close to it
        # relative error is less meaningful here
        return diff < (epsilon * sys.float_info.min)
    else:  # use relative error
        return diff / min((absA + absB), sys.float_info.max) < epsilon


nearlyEqualArray = numpy.vectorize(nearlyEqual)
nearlyEqualArray.__doc__ = "Vector version of nearlyEqual()."


def restrain_to_index_i_of_dim_d(a, i, d, n=None):
    """
    Of an array a[d1, d2, d3, ... dn], returns the array restricted to
    index **i** of the dimension **d**.

    :param a: the input array
    :param i: index in dimension **d**
    :param d: the dimension to restrain
    :param n: specify *a priori* the number of dimensions of **a**

    A more elegant solution would have been the following, but it does
    not work when accessing netCDF variable (for which it was necessary)::

        indexes = [range(len(self._dimensions[d])) for d in variable.dimensions] # equivalent to [:, :, :, ...]
        for k in only.keys():
            indexes[variable.dimensions.index(k)] = [only[k]] # restrain to the "only" give
        return array[numpy.ix_(*indexes)]
    """
    if n is None:
        n = len(a.shape)
    if n == 1:
        ra = a[[i]]
    elif n == 2:
        if d == 0:
            ra = a[[i], :]
        else:
            ra = a[:, [i]]
    elif n == 3:
        if d == 0:
            ra = a[[i], :, :]
        elif d == 1:
            ra = a[:, [i], :]
        else:
            ra = a[:, :, [i]]
    elif n == 4:
        if d == 0:
            ra = a[[i], :, :, :]
        elif d == 1:
            ra = a[:, [i], :, :]
        elif d == 2:
            ra = a[:, :, [i], :]
        else:
            ra = a[:, :, :, [i]]
    elif n == 5:
        if d == 0:
            ra = a[[i], :, :, :, :]
        elif d == 1:
            ra = a[:, [i], :, :, :]
        elif d == 2:
            ra = a[:, :, [i], :, :]
        elif d == 3:
            ra = a[:, :, :, [i], :]
        else:
            ra = a[:, :, :, :, [i]]
    elif n == 6:
        if d == 0:
            ra = a[[i], :, :, :, :, :]
        elif d == 1:
            ra = a[:, [i], :, :, :, :]
        elif d == 2:
            ra = a[:, :, [i], :, :, :]
        elif d == 3:
            ra = a[:, :, :, [i], :, :]
        elif d == 4:
            ra = a[:, :, :, :, [i], :]
        else:
            ra = a[:, :, :, :, :, [i]]
    else:
        raise NotImplementedError("more than 5 dimensions in array.")
    return ra


def datetimes2fieldvaliditylist(datetimes, basis=None):
    """
    Return a FieldValidityList from a list of datetime.datetime instances
    (or a single datetime.datetime).

    :param basis: can be either
                  - None (default): basis = validity
                  - a single datetime.datetime
                  - a list of the same length as datetimes
    """
    from epygram.base import FieldValidityList
    if isinstance(datetimes, datetime.datetime):
        assert (isinstance(basis, datetime.datetime) or
                (isinstance(basis, list) and
                 isinstance(basis[0], datetime.datetime)) or
                basis is None)
        if isinstance(basis, list):
            fvl = FieldValidityList(date_time=[datetimes for _ in basis], basis=basis)
        else:
            fvl = FieldValidityList(date_time=datetimes, basis=basis)
    elif isinstance(datetimes, list) and isinstance(datetimes[0], datetime.datetime):
        assert (isinstance(basis, datetime.datetime) or
                basis is None or
                (isinstance(basis, list) and
                 isinstance(basis[0], datetime.datetime) and
                 len(basis) == len(datetimes)))
        if isinstance(basis, datetime.datetime) or basis is None:
            basis = [basis for _ in range(len(datetimes))]
        fvl = FieldValidityList(date_time=datetimes,
                                basis=basis)
    else:
        raise TypeError("'datetimes' must be a datetime.datetime or a list of.")
    return fvl


def ifNone_emptydict(arg):
    """
    Transforms a None into a {}.
    To be used as workaround for empty dicts in default values of methods.
    """
    if arg is None:
        arg = {}
    return arg


def fmtfid(fmt, fid):
    """
    Given a resource format name **fmt** and a (full) **fid**, returns the key
    corresponding the actual format of the field in resource.
    (Useful for distinguishing GRIB1/2)
    """
    if fmt == 'GRIB':
        if 'GRIB1' in fid:
            fmtfid = 'GRIB1'
        else:
            fmtfid = 'GRIB2'
    else:
        fmtfid = fmt
    return fmtfid


def moveaxis(a, source, destination):
    """
    Calls numpy.moveaxis(), or if numpy version is too old, emulates for simple
    cases.
    """
    return numpy.moveaxis(a, source, destination)


@nicedeco
def call_before(mtd, hook_mtd):
    """Decorator for methods: call method hook_mtd before actually calling method."""

    def hooked(self, *args, **kwargs):
        getattr(self, hook_mtd)()
        return mtd(self, *args, **kwargs)

    return hooked


def mpl_interactive_backend():
    """Return whether the matplotlib backend is an interactive one or not."""
    import matplotlib
    noninteractive_backends = ('agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template')
    if matplotlib.get_backend().lower() in noninteractive_backends:
        return False
    else:
        return True
