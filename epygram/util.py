#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some useful utilities...
"""

from __future__ import print_function, absolute_import, division , unicode_literals
import six

import math
import copy
import numpy
import sys
import datetime
import hashlib
import os
from six.moves.urllib.request import urlopen  # @UnresolvedImport
from distutils.version import LooseVersion

from footprints import FootprintBase
from bronx.graphics.colormapping import add_cmap, get_norm4colorscale
from bronx.syntax.decorators import nicedeco

from epygram import config, epygramError


class RecursiveObject(object):
    """
    Generic abstract class implementing useful recursive properties:

    - display of object: the *__str__* method returns str(attr) for each of the
      object's attributes, with automatical indentation for readability.

    - test of (in)equality: the a == b test will be true if a and b have the
      same attributes and *a.attr == b.attr* for each attribute.
    """

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
                elif attr not in ('_puredict', '_observer'):
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
        return not self == other

    def __hash__(self):
        # known issue __eq__/__hash__ must be defined both or none, else inheritance is broken
        return object.__hash__(self)

    def tolerant_equal(self, other, tolerance=config.epsilon):
        """
        Test of equality by recursion on the object's attributes,
        with a **tolerance**.
        """
        if self.__class__ == other.__class__ and \
           set(self.__dict__.keys()) == set(other.__dict__.keys()):
            ok = True
            for attr in self.__dict__.keys():
                if attr in ('_puredict', '_observer'):  # footprints special attributes
                    continue
                else:
                    if not Comparator.are_equal(self.__dict__[attr],
                                                other.__dict__[attr],
                                                tolerance):
                        ok = False
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
            if self.__class__ == other.__class__ and \
               set(self.__dict__.keys()) == set(other.__dict__.keys()):
                for attr in self.__dict__.keys():
                    if attr in ('_puredict', '_observer'):  # footprints special attributes
                        continue
                    else:
                        if not Comparator.are_equal(self.__dict__[attr],
                                                    other.__dict__[attr]):
                            if all([isinstance(obj, RecursiveObject)
                                    for obj in [self.__dict__[attr],
                                                other.__dict__[attr]]]):
                                diff[attr] = self.__dict__[attr].diff(other.__dict__[attr])
                            else:
                                diff[attr] = Comparator.diff(self.__dict__[attr],
                                                             other.__dict__[attr])
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
    def _dict_are_equal(cls, dict1, dict2):
        if set(dict1.keys()) == set(dict2.keys()):
            ok = True
            for k in dict1.keys():
                if not cls.are_equal(dict1[k], dict2[k]):
                    ok = False
                    break
            return ok
        else:
            return False

    @classmethod
    def _dict_diff(cls, dict1, dict2):
        if not cls.are_equal(dict1, dict2):
            diff = {}
            if set(dict1.keys()) == set(dict2.keys()):
                for k in dict1.keys():
                    if not cls.are_equal(dict1[k], dict2[k]):
                        diff[k] = cls.diff(dict1[k], dict2[k])
            else:
                diff = (str(dict1), str(dict2))
        else:
            diff = None
        return diff

    @classmethod
    def _list_are_equal(cls, list1, list2):
        if len(list1) != len(list2):
            return False
        else:
            ok = True
            for i in range(len(list1)):
                if not cls.are_equal(list1[i], list2[i]):
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
            return cls._dict_are_equal(obj1, obj2)
        elif isinstance(obj1, list) and isinstance(obj2, list):
            return cls._list_are_equal(obj1, obj2)
        else:
            try:
                return obj1.__eq__(obj2)
            except AttributeError:
                return obj1 == obj2

    @classmethod
    def diff(cls, obj1, obj2):
        """Inspect differences between objects."""
        if not cls.are_equal(obj1, obj2):
            if isinstance(obj1, dict) and isinstance(obj2, dict):
                return cls._dict_diff(obj1, obj2)
            else:
                return (obj1, obj2)


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

    def __eq__(self, other):
        """
        Redefinition because of dynamism of buffering new computed values...
        """
        if not isinstance(other, Angle):
            return False
        if abs(self.get('degrees') -
               degrees_nearest_mod(other.get('degrees'),
                                   self.get('degrees'))) <= config.epsilon:
            ok = True
        else:
            ok = False
        return ok

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._origin_value) + hash(self._origin_unit)

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


# FUNCTIONS #
#############
def as_numpy_array(x):
    """
    Returns a numpy array
    """
    return numpy.array(x, copy=False, ndmin=1, subok=True)


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
        if not isinstance(pattern, six.string_types) or \
           not isinstance(element, six.string_types):
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
    if isinstance(regexp, six.string_types):
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
        url_hash = md5(url)
        actual_url_hash = md5(actual_url)
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
    """
    import matplotlib.pyplot as plt
    if cmap not in plt.colormaps() and cmap in config.colormaps:
        with open(config.colormaps[cmap], 'r') as ocm:
            add_cmap(cmap, ocm)


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


def write_formatted_dict(dest, fid):
    name = fid.pop('name')
    dest.write('name: ' + name + '\n')
    for k in sorted(fid.keys()):
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
            if isinstance(elem, six.string_types):
                elements.append(elem)
            elif isinstance(elem, float) or isinstance(elem, int) :
                elements.append(float_style.format(elem,
                                                   precision=precision,
                                                   type=float_type))
            else:
                elements.append('-')
        line = ('{:' + alignments[0] + '{width}}').format(table[i][0], width=columns_dimension[0])
        line += ''.join([('{:' + alignments[1] + '{width}}').format(elements[j], width=columns_dimension[j + 1])
                         for j in range(len(elements))])
        dest.write(line + '\n')


def add_meridians_and_parallels_to(bm,
                                   meridians='auto',
                                   parallels='auto',
                                   ax=None,
                                   drawparallels_kwargs=None,
                                   drawmeridians_kwargs=None,
                                   drawequator_kwargs=None,
                                   drawgreenwich_kwargs=None):
    """
    Adds meridians and parallels to a basemap instance *bm*.

    *meridians* and *parallels* enable to fine-tune the choice of lines to
    plot, with either:
      - 'auto': automatic scaling to the basemap extents
      - 'default': range(0,360,10) and range(-90,90,10)
      - a list of values
      - a grid step, e.g. 5 to plot each 5 degree.
      - None: no one is plot
      - *meridian* == 'greenwich' // 'datechange' // 'greenwich+datechange'
        *parallel* == 'equator' // 'polarcircles' // 'tropics' or any
        combination (,) will plot only these.

    :param ax: the ax to be plotted on
    :param drawparallels_kwargs: kwargs to be passed to basemap.drawparallels()
    :param drawmeridians_kwargs: kwargs to be passed to basemap.drawgreenwich()
    :param drawequator_kwargs: draw kwargs to emphasize equator parallel
    :param drawgreenwich_kwargs: draw kwargs to emphasize greenwich meridian
    """
    try:
        parallels = float(parallels)
    except (TypeError, ValueError):
        parallels = parallels
    try:
        meridians = float(meridians)
    except (TypeError, ValueError):
        meridians = meridians
    drawparallels_kwargs = ifNone_emptydict(drawparallels_kwargs)
    drawmeridians_kwargs = ifNone_emptydict(drawmeridians_kwargs)
    drawequator_kwargs = ifNone_emptydict(drawequator_kwargs)
    drawgreenwich_kwargs = ifNone_emptydict(drawgreenwich_kwargs)

    if bm.projection == 'rotpole':
        Delta_lat = bm.ymax - bm.ymin
    else:
        Delta_lat = bm.latmax - bm.latmin
    if parallels == 'auto' or isinstance(parallels, int) or isinstance(parallels, float):
        if Delta_lat <= 10:
            delta_lat = 1
        elif 10 < Delta_lat <= 30:
            delta_lat = 5
        else:
            delta_lat = 10
        if isinstance(parallels, int) or isinstance(parallels, float):
            delta_lat = parallels
        latmin = bm.latmin - bm.latmin % delta_lat
        latmax = bm.latmax - bm.latmax % delta_lat + delta_lat
        parallels = numpy.arange(latmin, latmax, delta_lat)
    elif parallels == 'default':
        parallels = numpy.arange(-90, 90, 10)
    elif isinstance(parallels, six.string_types) or isinstance(parallels, list):
        pl = []
        if 'equator' in parallels:
            pl.append(0.)
        if 'polarcircles' in parallels:
            pl.extend([-66.5628, 66.5628])
        if 'tropics' in parallels:
            pl.extend([-23.4372, 23.4372])
        if isinstance(parallels, list):
            for p in parallels:
                try:
                    pl.append(float(p))
                except (ValueError, TypeError):
                    pass
        pl.sort()
        parallels = pl

    if bm.projection == 'rotpole':
        Delta_lon = bm.xmax - bm.xmin
    else:
        Delta_lon = bm.lonmax - bm.lonmin
    if meridians == 'auto' or isinstance(meridians, int) or isinstance(meridians, float):
        if Delta_lon <= 10:
            delta_lon = 1
        elif 10 < Delta_lon <= 30:
            delta_lon = 5
        elif 30 < Delta_lon <= 180:
            delta_lon = 10
        else:
            delta_lon = 20
        if isinstance(meridians, int) or isinstance(meridians, float):
            delta_lon = meridians
        lonmin = bm.lonmin - bm.lonmin % delta_lon
        lonmax = bm.lonmax - bm.lonmax % delta_lon + delta_lon
        meridians = numpy.arange(lonmin, lonmax, delta_lon)
    elif meridians == 'default':
        meridians = numpy.arange(0, 360, 10)
    elif isinstance(meridians, six.string_types) or isinstance(meridians, list):
        ml = []
        if 'greenwich' in meridians:
            ml.append(0.)
        if 'datechange' in meridians:
            ml.append(180.)
        if isinstance(meridians, list):
            for m in meridians:
                try:
                    ml.append(float(m))
                except (ValueError, TypeError):
                    pass
        ml.sort()
        meridians = ml

    if parallels is not None:
        if bm.projection in ('ortho', 'nsper'):
            if 'labels' not in drawparallels_kwargs.keys():
                drawparallels_kwargs['labels'] = [False, False, False, False]
        else:
            if 'labels' not in drawparallels_kwargs.keys():
                drawparallels_kwargs['labels'] = [True, False, False, False]
        bm.drawparallels(parallels, ax=ax, **drawparallels_kwargs)
        if 0. in parallels or 0 in parallels:
            if 'dashes' not in drawequator_kwargs.keys():
                drawequator_kwargs['dashes'] = [10, 1]
            drawequator_kwargs['labels'] = [False] * 4
            bm.drawparallels([0], ax=ax, **drawequator_kwargs)
    if meridians is not None:
        if bm.projection in ('spstere', 'npstere', 'stere'):
            if 'labels' not in drawmeridians_kwargs.keys():
                drawmeridians_kwargs['labels'] = [True, False, False, True]
        elif bm.projection in ('ortho', 'moll', 'nsper'):
            if 'labels' not in drawmeridians_kwargs.keys():
                drawmeridians_kwargs['labels'] = [False, False, False, False]
        else:
            if 'labels' not in drawmeridians_kwargs.keys():
                drawmeridians_kwargs['labels'] = [False, False, False, True]
        bm.drawmeridians(meridians, ax=ax, **drawmeridians_kwargs)
        if 0. in meridians or 0 in meridians:
            if 'dashes' not in drawgreenwich_kwargs.keys():
                drawgreenwich_kwargs['dashes'] = [10, 1]
            drawgreenwich_kwargs['labels'] = [False] * 4
            bm.drawmeridians([0], ax=ax, **drawgreenwich_kwargs)


def auto_meridians_parallels(geometry,
                             meridians='auto',
                             parallels='auto'):
    """
    Compute meridians and parallels.

    *meridians* and *parallels* enable to fine-tune the choice of lines to
    plot, with either:
      - 'auto': automatic scaling to the basemap extents
      - 'default': range(0,360,10) and range(-90,90,10)
      - a list of values
      - a grid step, e.g. 5 to plot each 5 degree.
      - None: no one is plot
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
    if meridians == 'default' or (meridians == 'auto' and
                                  ('gauss' in geometry.name or
                                   geometry.name == 'polar_stereographic')):
        meridians = numpy.arange(0, 370, 10)
    elif meridians == 'auto' or isinstance(meridians, int) or isinstance(meridians, float):
        minmax = geometry.minmax_ll()
        if meridians == 'auto':
            delta_lon = Delta2delta(minmax['lonmax'] - minmax['lonmin'])
        else:
            delta_lon = float(meridians)
        lonmin = minmax['lonmin'] - minmax['lonmin'] % delta_lon
        lonmax = minmax['lonmax'] - minmax['lonmax'] % delta_lon + 2 * delta_lon
        meridians = numpy.arange(lonmin, lonmax, delta_lon)
        if max(meridians > 180.):
            meridians = meridians - 180.  # FIXME: cartopy does not plot meridians > 180°
    elif meridians is None:
        meridians = []
    # parallels
    if parallels == 'default' or ('gauss' in geometry.name and parallels == 'auto'):
        parallels = numpy.arange(-90, 90, 10)
    elif parallels == 'auto' or isinstance(parallels, int) or isinstance(parallels, float):
        minmax = geometry.minmax_ll()
        if parallels == 'auto':
            delta_lat = Delta2delta(minmax['latmax'] - minmax['latmin'])
        else:
            delta_lat = float(parallels)
        latmin = minmax['latmin'] - minmax['latmin'] % delta_lat
        latmax = minmax['latmax'] - minmax['latmax'] % delta_lat + 2 * delta_lat
        parallels = numpy.arange(latmin, latmax, delta_lat)
    elif parallels is None:
        parallels = []
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


def scale_colormap(cmap, max_val=None):
    """
    Creates a matplotlib.colors.BoundaryNorm object tuned for scaled colormaps,
    i.e. discrete, irregular colorshades.

    :param cmap: name of the colormap, as found in config.colormaps_scaling
    :param max_val: if given, replaces the upper bound.

    :return: a tuple (norm, scaling), scaling being eventually modified
             according to **max_val**
    """
    bounds = copy.copy(config.colormaps_scaling.get(cmap, None))
    return get_norm4colorscale(bounds, max_val=max_val)


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


def set_map_up(bm, ax,
               drawrivers=False,
               drawcoastlines=True,
               drawcountries=True,
               meridians='auto',
               parallels='auto',
               departments=False,
               bluemarble=0.0,
               background=False,
               drawmapboundary_kwargs=None,
               fillcontinents_kwargs=None,
               drawcoastlines_kwargs=None,
               drawcountries_kwargs=None,
               drawparallels_kwargs=None,
               drawmeridians_kwargs=None,
               drawequator_kwargs=None,
               drawgreenwich_kwargs=None):
    """Cf. :meth:`H2DField.plotfield` documentation."""
    if drawmapboundary_kwargs is None:
        drawmapboundary_kwargs = dict(fill_color='lightskyblue')
    if fillcontinents_kwargs is None:
        fillcontinents_kwargs = dict(color='wheat', lake_color='skyblue',
                                     zorder=0)
    drawcoastlines_kwargs = ifNone_emptydict(drawcoastlines_kwargs)
    drawcountries_kwargs = ifNone_emptydict(drawcountries_kwargs)
    if background:
        bm.drawmapboundary(ax=ax, **drawmapboundary_kwargs)
        bm.fillcontinents(ax=ax, **fillcontinents_kwargs)
    if bluemarble:
        bm.bluemarble(alpha=bluemarble, ax=ax)
    if drawcoastlines:
        bm.drawcoastlines(ax=ax, **drawcoastlines_kwargs)
    if departments:
        if not hasattr(bm, '_epygram_departments'):
            import json
            with open(config.installdir + '/data/french_departments.json', 'r') as dp:
                depts = json.load(dp)[1]
            bm._epygram_departments = depts
        else:
            depts = bm._epygram_departments
        for d in range(len(depts)):
            for part in range(len(depts[d])):
                dlon = depts[d][part][0]
                dlat = depts[d][part][1]
                (x, y) = bm(dlon, dlat)
                bm.plot(x, y, color=drawcountries_kwargs.get('color', 'k'), ax=ax)
    elif drawcountries:
        bm.drawcountries(ax=ax, **drawcountries_kwargs)
    if drawrivers:
        bm.drawrivers(color='blue', ax=ax)
    add_meridians_and_parallels_to(bm,
                                   parallels=parallels,
                                   meridians=meridians,
                                   ax=ax,
                                   drawparallels_kwargs=drawparallels_kwargs,
                                   drawmeridians_kwargs=drawmeridians_kwargs,
                                   drawequator_kwargs=drawequator_kwargs,
                                   drawgreenwich_kwargs=drawgreenwich_kwargs)


def vtk_modify_grid(grid, grid_type, datamin=None):
    """
    Modifies the kind of grid
    Input grid must be an sgrid_point
    :param grid_type: can be:
        - sgrid_point: structured grid filled with points
        - sgrid_cell: structured grid filled with hexahedron
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
        - ugrid_point: unstructured grid filled with points
        - ugrid_cell: unstructured grid build filled with cells
                      If the field is 2D, a zero thickness is used.
                      If the field is 3D, thickness are approximately computed
    :param datamin: for an unknown reason, we need the minimum of the data
                    to transform grid into an unstructured one

    If grid_type is 'sgrid_point', the result is the grid; otherwise
    the result is the function is the last filter used.
    """
    import vtk  # @UnresolvedImport

    if grid_type not in ('sgrid_point', 'sgrid_cell', 'ugrid_point', 'ugrid_cell'):
        raise ValueError("Unknown grid type: " + grid_type)

    if grid_type in ('sgrid_cell', 'ugrid_cell'):
        interp = vtk.vtkPointDataToCellData()
        interp.SetInputData(grid)
        interp.Update()
        grid = interp
        # Values are now associated to cells

    if grid_type in ('ugrid_point', 'ugrid_cell'):
        fil = vtk.vtkThreshold()
        if grid_type in ('sgrid_cell', 'ugrid_cell'):
            fil.SetInputConnection(grid.GetOutputPort())
        else:
            fil.SetInputData(grid)
        # minScalar = grid.GetPointData().GetScalars().GetRange()[0] does not work every time
        if datamin is None:
            raise epygramError("datamin must be provided for unstructured grid types")
        minScalar = datamin
        fil.ThresholdByUpper(minScalar - 1.)
        grid = fil
        # Grid is now unstructured

    return grid


def vtk_write_png(rendering, filename, resolution_increase=1, enable_alpha=True):
    """
    Writes a png file with the present vien on vtk window
    :param rendering:  a dictionary containing, at least, the window key
    :param filename: name of the file to produce
    :param resolution_increase: vtk window resolution is multiplied
                                by this factor to get the png resolution
    """
    import vtk  # @UnresolvedImport
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(rendering['window'])
    windowToImageFilter.SetMagnification(resolution_increase)
    if enable_alpha:
        windowToImageFilter.SetInputBufferTypeToRGBA()  # also record the alpha (transparency) channel
    # windowToImageFilter.ReadFrontBufferOff() #Do not know what this means but we get bad images when uncomented
    PNGWriter = vtk.vtkPNGWriter()
    PNGWriter.SetFileName(filename)
    windowToImageFilter.Update()
    PNGWriter.SetInputConnection(windowToImageFilter.GetOutputPort())
    PNGWriter.Write()


def vtk_set_window(background_color, window_size, hide_axes=False, offscreen=False):
    """
    This function creates a simple vtk environment and returns
    a dictionary holding the different objects created
    :param background_color: must be a color name or a 3-tuple
    :param window_size: must be a 2-tuple
    :param hide_axes: True to hide the axes representation
    :param offscreen: True to hide window (useful when we only
                      want to produce png file instead of
                      interactively viewing the window)
    """
    import vtk  # @UnresolvedImport

    renderer = vtk.vtkRenderer()
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    result = dict(renderer=renderer, window=renderWin)
    if offscreen:
        renderWin.SetOffScreenRendering(True)
    else:
        renderInteractor = vtk.vtkRenderWindowInteractor()
        renderInteractor.SetRenderWindow(renderWin)
        style = vtk.vtkInteractorStyleTrackballCamera()
        renderInteractor.SetInteractorStyle(style)

        def exitCheck(obj, event):
            if obj.GetEventPending() != 0:
                obj.SetAbortRender(1)
        renderWin.AddObserver("AbortCheckEvent", exitCheck)
        renderInteractor.Initialize()
        result['interactor'] = renderInteractor

    if isinstance(background_color, tuple):
        renderer.SetBackground(*background_color)
    else:
        renderer.SetBackground(vtk.vtkNamedColors().GetColor3d(background_color))
    renderWin.SetSize(*window_size)

    axes = vtk.vtkAxesActor()
    axes.SetTotalLength([30., 30., 30.])
    # axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(colors.GetColor3d("Red"));
    # axes.SetXAxisLabelText("test");
    # axes.GetYAxisShaftProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Yellow"))
    # axes.GetYAxisTipProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Orange"))
    if hide_axes:
        axes.VisibilityOff()
    renderer.AddActor(axes)

    return result


def vtk_print_text_in_window(rendering, text, pos, fontsize=20, color='Black'):
    """
    This function write text in a vtk window
    :param rendering:  a dictionary containing, at least, the renderer key
    :param text: the text to render
    :param pos: tuple (x, y) of the position of the text (pixel unit)
    :param fontsize: fontsize to use for the text
    :param color: color of the text
    """
    import vtk  # @UnresolvedImport
    textActor = vtk.vtkTextActor()
    textActor.SetInput(text)
    textActor.SetDisplayPosition(*pos)  # SetPosition2(pos)
    textActor.GetTextProperty().SetFontSize(fontsize)
    textActor.GetTextProperty().SetColor(vtk.vtkNamedColors().GetColor3d(color))
    rendering['renderer'].AddActor2D(textActor)


def vtk_check_transform(rendering, current_typeoffirstfixedsurface, hCoord, z_factor, offset):
    """
    :param rendering:  a dictionary containing, at least, the renderer key
    :param current_typeoffirstfixedsurface: typeoffirstfixedsurface associated to the object to plot
    :param hCoord: 'll': horizontal coordinates are the lon/lat values
                   a basemap: horizontal coordinates are set according to this basemap
    :param z_factor: factor to apply on z values (to modify aspect ratio of the plot)
    :param offset: (x_offset, y_offset). Offsets are subtracted to x and y coordinates
    """

    if not hasattr(rendering['renderer'], 'epygram'):
        rendering['renderer'].epygram = dict()
    conf = rendering['renderer'].epygram

    if 'Z_axis_type' not in conf:
        conf['Z_axis_type'] = current_typeoffirstfixedsurface
    assert conf['Z_axis_type'] == current_typeoffirstfixedsurface, \
           "type of first fixed surface must be the same for all plotted objects"

    if 'hCoord' not in conf:
        conf['hCoord'] = hCoord if hCoord is not None else 'll'
    assert hCoord is None or conf['hCoord'] == hCoord, \
           "horizontal transformation of coordinate (hCoord option) must be the same for all plotted objects"

    if 'z_factor' not in conf:
        conf['z_factor'] = z_factor if z_factor is not None else 1.
    assert z_factor is None or conf['z_factor'] == z_factor, \
           "factor applied on z axis must be the same for all plotted objects"

    if 'offset' not in conf:
        conf['offset'] = offset if offset is not None else (0., 0.)
    assert offset is None or conf['offset'] == offset, \
           "offset applied on x/y axis must be the same for all plotted objects"

    return conf['hCoord'], conf['z_factor'], conf['offset']


def vtk_write_grid(grid, filename):
    """
    Writes a grid into a file
    :param grid: a vtk grid
    :param filename: filename to save in
    """
    import vtk  # @UnresolvedImport
    if isinstance(grid, vtk.vtkUnstructuredGrid):
        writer = vtk.vtkUnstructuredGridWriter()
    elif isinstance(grid, vtk.vtkStructuredGrid):
        writer = vtk.vtkStructuredGridWriter()
    else:
        epygramError('Unknown grid type')
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()


def datetimerange(*_, **__):
    """
    .. deprecated:: 1.2.11
    """
    raise DeprecationWarning("You should use function daterange/daterangex " +
                             "from bronx.stdtypes.date")


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
    # CLEANME: workaround for old installs...
    if LooseVersion(numpy.__version__) < LooseVersion('1.11.0'):
        if source == 0 and destination == -1:
            b = numpy.transpose(a, list(range(len(a.shape))[1:]) + [0])
        elif source == -1 and destination == 0:
            b = numpy.transpose(a, [-1] + list(range(len(a.shape))[:-1]))
        else:
            raise NotImplementedError('(source, destination) != (0,-1) or (-1,0) with that version of numpy')
    else:
        b = numpy.moveaxis(a, source, destination)
    return b


@nicedeco
def call_before(mtd, hook_mtd):
    """Decorator for methods: call method hook_mtd before actually calling method."""
    def hooked(self, *args, **kwargs):
        getattr(self, hook_mtd)()
        return mtd(self, *args, **kwargs)
    return hooked
