#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some useful utilities...
"""

from __future__ import print_function, absolute_import, division  # , unicode_literals
import six

import math
import copy
import numpy
import sys
import datetime

from footprints import FootprintBase
from bronx.graphics.colormapping import add_cmap, get_norm4colorscale

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
                itemstring += "\n" + offset + key + ": " + self._strItem(item[key], reclevel + 1)
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
        def comp_float(float1, float2):
            # tolerance for floats
            return abs(float1 - float2) <= config.epsilon

        def comp_array(array1, array2):
            if (array1.dtype == array2.dtype and
                array1.dtype in [numpy.dtype(d)
                                 for d in ['float16', 'float32', 'float64']]):
                # tolerance for floats
                return (abs(array1 - array2) <= config.epsilon).all()
            else:
                return numpy.all(array1 == array2)

        def comp_dict(dict1, dict2):
            if set(dict1.keys()) == set(dict2.keys()):
                ok = True
                for k in dict1.keys():
                    if not comp(dict1[k], dict2[k]):
                        ok = False
                        break
                return ok
            else:
                return False

        def comp_list(list1, list2):
            if len(list1) != len(list2):
                return False
            else:
                ok = True
                for i in range(len(list1)):
                    if not comp(list1[i], list2[i]):
                        ok = False
                        break
                return ok

        def comp(obj1, obj2):
            if isinstance(obj1, float) and isinstance(obj2, float):
                return comp_float(obj1, obj2)
            elif isinstance(obj1, numpy.ndarray) and isinstance(obj2, numpy.ndarray):
                return comp_array(obj1, obj2)
            elif isinstance(obj1, dict) and isinstance(obj2, dict):
                return comp_dict(obj1, obj2)
            elif isinstance(obj1, list) and isinstance(obj2, list):
                return comp_list(obj1, obj2)
            else:
                return obj1 == obj2

        if self.__class__ == other.__class__ and \
           set(self.__dict__.keys()) == set(other.__dict__.keys()):
            ok = True
            for attr in self.__dict__.keys():
                if attr in ('_puredict', '_observer'):
                    pass
                else:
                    if not comp(self.__dict__[attr], other.__dict__[attr]):
                        ok = False
                        break
        else:
            ok = False
        return ok

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        # known issue __eq__/must be defined both or none, else inheritance is broken
        # return super().__hash__() # CLEANME:
        return object.__hash__(self)
        # return super(RecursiveObject, self).__hash__()
        """_hash = 0
        for k, v in self.__dict__.items():
            if isinstance(v, dict) or isinstance(v, list):
                _hash += hash(k) + sum([hash(i) for i in v])
            else:
                _hash += hash(k) + hash(v)
        return _hash"""

    def copy(self):
        """Returns a copy of the object."""
        return copy.copy(self)

    def deepcopy(self):
        """Returns a deepcopy of the object."""
        return copy.deepcopy(self)


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
            raise ValueError("cannot compare instances of different classes.")
        if abs(self.get('radians') - other.get('radians')) <= config.epsilon:
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
        n = len(d)
        scalar = False
    except Exception:
        n = 1
        d = [d]
        scalar = True
    d = numpy.array(d)
    d_inf = d - 360.
    d_sup = d + 360.
    result = numpy.zeros(n)
    for i in range(n):
        if abs(d[i] - ref) <= abs(d_sup[i] - ref) and \
           abs(d[i] - ref) <= abs(d_inf[i] - ref):
            result[i] = d[i]
        elif abs(d_inf[i] - ref) <= abs(d_sup[i] - ref) and \
             abs(d_inf[i] - ref) <= abs(d[i] - ref):
            result[i] = d_inf[i]
        else:
            result[i] = d_sup[i]
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


def load_cmap(cmap):
    """
    Reads and registers the epygram-or-user-colormap called *cmap*,
    which must be either in config.epygram_colormaps or
    config.usercolormaps.
    """
    import matplotlib.pyplot as plt
    if cmap not in plt.colormaps():
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
        import json
        with open(config.installdir + '/data/departments.json', 'r') as dp:
            depts = json.load(dp)[1]
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
