#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some useful utilities... 
"""

import os
import math
import copy
import numpy
import sys
import datetime
from contextlib import contextmanager

from footprints import FootprintBase

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
            if array1.dtype == array2.dtype and array1.dtype in [numpy.dtype(d) for d in ['float16', 'float32', 'float64']]:
                # tolerance for floats
                return (abs(array1 - array2) <= config.epsilon).all()
            else:
                return numpy.all(array1 == array2)
        def comp_dict(dict1, dict2):
            if dict1.keys() == dict2.keys():
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
            if type(obj1) == type(obj2) == float:
                return comp_float(obj1, obj2)
            elif type(obj1) == type(obj2) == numpy.ndarray:
                return comp_array(obj1, obj2)
            elif type(obj1) == type(obj2) == dict:
                return comp_dict(obj1, obj2)
            elif type(obj1) == type(obj2) == list:
                return comp_list(obj1, obj2)
            else:
                return obj1 == obj2

        if self.__class__ == other.__class__ and \
           self.__dict__.keys() == other.__dict__.keys():
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

    def get(self, unit=None):
        """
        Returns the angle in the requested unit.
        If no unit is supplied, the origin unit is used.
        """

        if unit == None:
            unit = self._origin_unit  # or a default one ?
        elif unit in Angle.units:
            if not self.__dict__.has_key('_' + unit):
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
                self.__dict__['_' + unit] = math.degrees(math.copysign(math.acos(self._origin_value[0]),
                                                                       self._origin_value[1]))
            elif self._origin_unit == Angle.dms:
                self.__dict__['_' + unit] = self._origin_value[0] + self._origin_value[1] / 60. + self._origin_value[2] / 3600.
            else:
                raise NotImplementedError("conversion from this unit (" + \
                                          self._origin_unit + \
                                          ") is not coded.")

        # conversion to radians
        elif unit == Angle.rad:
            if self._origin_unit == Angle.deg:
                self.__dict__['_' + unit] = math.radians(self._origin_value)
            elif self._origin_unit == Angle.trig:
                self.__dict__['_' + unit] = math.copysign(math.acos(self._origin_value[0]),
                                                          self._origin_value[1])
            elif self._origin_unit == Angle.dms:
                self.__dict__['_' + unit] = math.radians(self._origin_value[0] + self._origin_value[1] / 60. + self._origin_value[2] / 3600.)
            else:
                raise NotImplementedError("conversion from this unit (" + \
                                          self._origin_unit + \
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
                anglerad = math.radians(self._origin_value[0] + self._origin_value[1] / 60. + self._origin_value[2] / 3600.)
                self.__dict__['_' + unit] = (math.cos(anglerad),
                                             math.sin(anglerad))
            else:
                raise NotImplementedError("conversion from this unit (" + \
                                          self._origin_unit + \
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
                raise NotImplementedError("conversion from this unit (" + \
                                          self._origin_unit + \
                                          ") is not coded.")
            sign = int(math.copysign(1, decdeg))
            decdeg = decdeg * sign
            degrees = int(decdeg)
            decmin = (decdeg - degrees) * 60.
            minutes = int(decmin)
            seconds = (decmin - minutes) * 60
            self.__dict__['_' + unit] = (degrees * sign, minutes, seconds)

        else:
            raise NotImplementedError("conversion to this unit (" + unit + \
                                      ") is not coded.")



#################
### FUNCTIONS ###
#################

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
        if (type(pattern) not in [type(""), type(u"")]) or \
           (type(element) not in [type(""), type(u"")]):
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
    if type(regexp) in [type(""), type(u"")]:
        for field in a_list:
            if check_string_pattern(regexp, str(field)):
                found.append(field)
    elif type(regexp) == type(()):
        for field in a_list:
            if type(field) != type(()):
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
    elif type(regexp) == type({}):
        raise NotImplementedError("the dictionnary type of regexp is not yet\
                                   implemented.")
    else:
        raise NotImplementedError("this type of regexp is not (yet?)\
                                   implemented.")

    return found

def nicedeco(decorator):
    """
    A decorator of decorator, for the decorated method to keep the original
    __name__, __doc__ and __dict__.
    """
    def new_decorator(f):
        g = decorator(f)
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        g.__dict__.update(f.__dict__)
        return g
    return new_decorator

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
            lon += 2.*numpy.pi
        else:
            raise NotImplementedError()

    return lon

def make_custom_cmap(filename):
    """
    Creates a custom Colormap from a set of RGB colors in a file with the
    following formating:
    
    r1,g1,b1;\n
    r2,g2,b2;\n
    ...\n
    rn,gn,bn
    
    each value being comprised between 0 and 255
    e.g. coming from http://colormap.org.
    """
    import matplotlib

    with open(filename, 'r') as f:
        colors = f.readlines()
    for i in range(len(colors)):
        colors[i] = colors[i].replace(';', '')
        colors[i] = colors[i].replace('[', '')
        colors[i] = colors[i].replace(']', '')
        colors[i] = colors[i].replace('\n', '')
        colors[i] = colors[i].split(',')
    colors = numpy.array(colors, dtype=numpy.float64)
    if colors.max() > 1.:
        cm = matplotlib.colors.ListedColormap(colors / 255.)
    else:
        cm = matplotlib.colors.ListedColormap(colors)

    return cm

def add_cmap(cmap):
    """
    Reads and registers the epygram-or-user-colormap called *cmap*,
    which must be either in config.epygram_colormaps or
    config.usercolormaps. 
    """
    import matplotlib.pyplot as plt

    if cmap not in plt.colormaps():
        if cmap in config.colormaps.keys():
            plt.register_cmap(name=cmap,
                              cmap=make_custom_cmap(config.colormaps[cmap]))
        else:
            raise epygramError("unknown colormap '" + cmap + \
                               "': must be added to userconfig, in \
                                usercolormaps.")

def printstatus(step, end, refresh_freq=1):
    """
    Print percentage of the loop it is in, with 'step' being the current
    loopstep, 'end' the final loopstep and 'refresh_freq' the frequency
    in % at which reprinting status.
    """

    status = step * 100 / end
    if status % refresh_freq == 0:
        sys.stdout.write('{:>{width}}%'.format(int(status), width=3))
        sys.stdout.flush()
        if step < end :
            sys.stdout.write('\b' * 4)
        else:
            sys.stdout.write('\n')

def read_CSV_as_dict(filename):
    """
    Reads a .csv file as a list of dict, with the assumption:
    - on first line is described the delimiter
    - on second line is described the 'priority' of the dict.
    """
    import io
    import csv

    field_dict = []
    with io.open(filename, 'r') as f:
        delimiter = str(f.readline()[0])
        file_priority = str(f.readline()[0:-1])
        field_table = csv.reader(f, delimiter=delimiter)
        for field in field_table:
            # syntax example of field description:
            # name:FIELDNAME;param:value;...
            if len(field) > 1 and field[0][0] != '#':
                fd = {}
                for kv in field:
                    try:
                        fd[kv.split(':')[0]] = int(kv.split(':')[1])
                    except ValueError:
                        fd[kv.split(':')[0]] = kv.split(':')[1]
                field_dict.append(fd)

    return field_dict, file_priority

def gfl2R(q, ql=0., qi=0., qr=0., qs=0., qg=0.):
    """
    Computes air specific gas constant R according to specific humidity,
    and hydrometeors if present.
    """

    # Constants
    Rd = config.Rd
    Rv = config.Rv

    q = numpy.array(q)
    ql = numpy.array(ql)
    qi = numpy.array(qi)
    qr = numpy.array(qr)
    qs = numpy.array(qs)
    qg = numpy.array(qg)

    R = Rd + (Rv - Rd) * q - Rd * (ql + qi + qr + qs + qg)

    return R

formatting_default_widths = (50, 20)
separation_line = '{:-^{width}}'.format('', width=sum(formatting_default_widths) + 1) + '\n'
def write_formatted(dest, label, value,
                    align='<',
                    firstcolumn_width=formatting_default_widths[0],
                    secondcolumn_width=formatting_default_widths[1]):
    dest.write(('{:' + align + '{width}}').format(label, width=firstcolumn_width) \
               + ':' \
               + '{:>{width}}'.format(str(value), width=secondcolumn_width) \
               + '\n')
def write_formatted_fields(dest, label, value=None,
                           compression=None,
                           firstcolumn_width=formatting_default_widths[0],
                           secondcolumn_width=formatting_default_widths[1]):
    if compression is None:
        if value is None:
            dest.write('{:<{width}}'.format(label, width=20)
                       + '\n')
        else:
            dest.write('{:<{width}}'.format(label, width=20) \
                      + ':' \
                      + '{:^{width}}'.format(str(value), width=10) \
                      + '\n')
    else:
        line = '{:<{width}}'.format(label, width=20) \
             + ':' \
             + '{:^{width}}'.format(str(value), width=10) \
             + ':' \
             + compression \
             + '\n'
        dest.write(line)
def write_formatted_dict(dest, fid):
    name = fid.pop('name')
    dest.write('name: ' + name + '\n')
    for k in sorted(fid.keys()):
        dest.write('  ' + str(k) + ': ' + str(fid[k]) + '\n')

def write_formatted_table(dest, table, alignments=['<', '^'], precision=6, float_type='E'):
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
    float_len = len(float_style.format(2. / 3., type=float_type, precision=precision))
    columns_dimension[1:] = [max(float_len, c) for c in columns_dimension[1:]]

    for i in range(array.shape[0]):
        elements = []
        for elem in table[i][1:]:
            if isinstance(elem, str):
                elements.append(elem)
            elif isinstance(elem, float) or isinstance(elem, int) :
                elements.append(float_style.format(elem, precision=precision, type=float_type))
            else:
                elements.append('-')
        line = ('{:' + alignments[0] + '{width}}').format(table[i][0], width=columns_dimension[0])
        line += ''.join([('{:' + alignments[1] + '{width}}').format(elements[j], width=columns_dimension[j + 1]) for j in range(len(elements))])
        dest.write(line + '\n')

def linearize(s, quotes=False):
    """
    Returns string *s* linearized, i.e. without special characters that may
    be forbidden in filenames.
    - quotes: must we also remove quotes?
    """
    replacements = [(' ', '_'), ('{', ''), ('}', ''), ("'", ''), ('*', '')]
    if quotes:
        replacements = [("'", ""), ('"', '')]
    result = s.strip()
    for repl in replacements:
        result = result.replace(*repl)

    return result

def linearize2str(o, quotes=False):
    """Returns str(*o*) linearized (cf. util.linearized)."""
    return linearize(str(o))

def str_or_int_to_datetime(dt):
    """
    Creates a datetime.datetime from a string or int YYYYMMDDHHMMSS...
    """

    dt = str(dt)
    year = int(dt[0:4])
    month = int(dt[4:6])
    day = int(dt[6:8])
    hour = minutes = seconds = 0
    try:
        hour = int(dt[8:10])
        minutes = int(dt[10:12])
        seconds = int(dt[12:14])
    except IndexError:
        pass
    dt = datetime.datetime(year, month, day, hour, minutes, seconds)

    return dt

def add_meridians_and_parallels_to(bm, meridians='auto', parallels='auto', ax=None):
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
        combination (+) will plot only these.
    """

    try:
        parallels = float(parallels)
    except (TypeError, ValueError):
        parallels = parallels
    try:
        meridians = float(meridians)
    except (TypeError, ValueError):
        meridians = meridians

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
    elif isinstance(parallels, str):
        pl = []
        if 'equator' in parallels:
            pl.append(0.)
        if 'polarcircles' in parallels:
            pl.extend([-66.5628, 66.5628])
        if 'tropics' in parallels:
            pl.extend([-23.4372, 23.4372])
        pl.sort()
        parallels = pl

    Delta_lon = bm.lonmax - bm.lonmin
    if meridians == 'auto' or isinstance(meridians, int) or isinstance(meridians, float):
        if Delta_lon <= 10:
            delta_lon = 1
        elif 10 < Delta_lon <= 30:
            delta_lon = 5
        else:
            delta_lon = 10
        if isinstance(meridians, int) or isinstance(meridians, float):
            delta_lon = meridians
        lonmin = bm.lonmin - bm.lonmin % delta_lon
        lonmax = bm.lonmax - bm.lonmax % delta_lon + delta_lon
        meridians = numpy.arange(lonmin, lonmax, delta_lon)
    elif meridians == 'default':
        meridians = numpy.arange(0, 360, 10)
    elif isinstance(meridians, str):
        ml = []
        if 'greenwich' in meridians:
            ml.append(0.)
        if 'datechange' in meridians:
            pl.append(180.)
        ml.sort()
        meridians = ml

    if parallels is not None:
        if bm.projection in ('ortho', 'nsper'):
            bm.drawparallels(parallels, labels=[False, False, False, False],
                             ax=ax)
        else:
            bm.drawparallels(parallels, labels=[True, False, False, False],
                             ax=ax)
    if meridians is not None:
        if bm.projection in ('spstere', 'npstere'):
            bm.drawmeridians(meridians, labels=[True, False, False, True],
                             ax=ax)
        elif bm.projection in ('ortho', 'moll', 'nsper'):
            bm.drawmeridians(meridians, labels=[False, False, False, False],
                             ax=ax)
        else:
            bm.drawmeridians(meridians, labels=[False, False, False, True],
                             ax=ax)
        bm.drawmeridians([0], labels=[False] * 4, linewidth=1,
                         dashes=[10, 1],
                         ax=ax)
        bm.drawparallels([0], labels=[False] * 4, linewidth=1,
                         dashes=[10, 1],
                         ax=ax)

def nearlyEqual(a, b, epsilon=config.epsilon):
    """
    Function to compare floats
    http://floating-point-gui.de/errors/comparison/
    Float.MIN_NORMAL was replaced by sys.float_info.min
    Float.MAX_VALUE was replaced by sys.float_info.max
    """
    absA = numpy.abs(a);
    absB = numpy.abs(b);
    diff = numpy.abs(a - b);

    if a == b:  # shortcut, handles infinities
        return True
    elif a == 0 or b == 0 or diff < sys.float_info.min:
        #a or b is zero or both are extremely close to it
        #relative error is less meaningful here
        return diff < (epsilon * sys.float_info.min)
    else:  # use relative error
        return diff / min((absA + absB), sys.float_info.max) < epsilon

nearlyEqualArray = numpy.vectorize(nearlyEqual)

def parse_str2dict(string, try_convert=None):
    """
    Parse a *string* (of syntax 'key1:value1,key2=value2') to a dict.
    If *try_convert* is not None, try to convert values as type *try_convert*.
    """

    d = {i.replace('=', ':').split(':')[0].strip():i.replace('=', ':').split(':')[1].strip() for i in string.split(',')}
    if try_convert is not None:
        for k, v in d.items():
            try:
                d[k] = try_convert(v)
            except ValueError:
                pass
    return d

def stretch_array(array):
    """
    Return array.flatten() or compressed(), whether the array is
    masked or not.
    """

    if isinstance(array, numpy.ma.masked_array):
        array = array.compressed()
    elif isinstance(array, numpy.ndarray):
        array = array.flatten()
    else:
        raise NotImplementedError(' '.join(['type:', type(array), 'array']))

    return array

def color_scale(cmap, max_rr=None):
    """
    Creates a matplotlib.colors.BoundaryNorm object tuned for radar colormaps.
    """

    import matplotlib.colors as colors
    if cmap == 'radar':
        bounds = [0., 0.1, 1., 3., 5., 7., 10., 15., 20., 30., 50., 70., 100., 150.]
        if max_rr <= 150.:
            max_rr = 300.
    elif cmap in ('rr1h', 'rr6h'):
        bounds = [0., 0.2, 0.5, 1, 1.5, 2., 4., 10., 25., 50., 100.]
        if max_rr <= 100.:
            max_rr = 300.
    elif cmap == 'rr24h':
        bounds = [0., 0.2, 1., 2., 4., 10., 25., 50., 100., 150., 200., 300.]
        if max_rr <= 300.:
            max_rr = 500.
    else:
        raise NotImplementedError('cmap == ' + cmap)
    bounds.append(max_rr)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds) - 1)

    return (norm, bounds)

@contextmanager
def stdout_redirected(to=os.devnull):
    '''
    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    #http://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close()  # + implicit flush()
        os.dup2(to.fileno(), fd)  # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w')  # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as f:
            _redirect_stdout(to=f)
        try:
            yield  # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout)  # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

@contextmanager
def stderr_redirected(to=os.devnull):
    '''
    import os

    with stderr_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    #Based on stdout_redirected
    fd = sys.stderr.fileno()

    def _redirect_stderr(to):
        sys.stderr.close()  # + implicit flush()
        os.dup2(to.fileno(), fd)  # fd writes to 'to' file
        sys.stderr = os.fdopen(fd, 'w')  # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stderr:
        with open(to, 'w') as f:
            _redirect_stderr(to=f)
        try:
            yield  # allow code to be run with the redirected stderr
        finally:
            _redirect_stderr(to=old_stderr)  # restore stderr.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

def restrain_to_index_i_of_dim_d(a, i, d, n=None):
    """
    Of an array a[d1, d2, d3, ... dn], returns the array restricted to index i
    of the dimension d.
    
    A more elegant solution would have been the following, except that it does
    not work when accessing netCDF variable (for which it was necessary).
    
    indexes = [range(len(self._dimensions[d])) for d in variable.dimensions] # equivalent to [:, :, :, ...]
    for k in only.keys():
        indexes[variable.dimensions.index(k)] = [only[k]] # restrain to the "only" give
    return array[numpy.ix_(*indexes)]
    """

    if n is None:
        n = a.shape
    if n == 2:
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
    else:
        raise NotImplementedError("more than 5 dimensions in array.")

    return ra

def resize_to_4D(a, indexes, ma=False):
    """
    Of an array a[d1, d2, d3, ... dn], returns the array as a 4D one restricted to index i
    of the dimension d.
    
    A more elegant solution would have been the following, except that it does
    not work when accessing netCDF variable (for which it was necessary).
    
    indexes = [range(len(self._dimensions[d])) for d in variable.dimensions] # equivalent to [:, :, :, ...]
    for k in only.keys():
        indexes[variable.dimensions.index(k)] = [only[k]] # restrain to the "only" give
    return array[numpy.ix_(*indexes)]
    """



def datetimes2fieldvaliditylist(datetimes, basis=None):
    """
    Return a FieldValidityList from a list of datetime.datetime instances
    (or a single datetime.datetime).
    
    *basis* can be either
      - None (default): basis = validity
      - a single datetime.datetime
      - a list of the same length as datetimes
    """

    from epygram.base import FieldValidityList

    if isinstance(datetimes, datetime.datetime):
        assert isinstance(basis, datetime.datetime) \
               or (isinstance(basis, list) and isinstance(basis[0], datetime.datetime)) \
               or basis is None
        if isinstance(basis, list):
            fvl = FieldValidityList(date_time=[datetimes for _ in basis], basis=basis)
        else:
            fvl = FieldValidityList(date_time=datetimes, basis=basis)
    elif isinstance(datetimes, list) and isinstance(datetimes[0], datetime.datetime):
        assert isinstance(basis, datetime.datetime) \
               or basis is None \
               or (isinstance(basis, list) \
                   and isinstance(basis[0], datetime.datetime) \
                   and len(basis) == len(datetimes))
        if isinstance(basis, datetime.datetime) or basis is None:
            basis = [basis for _ in range(len(datetimes))]
        fvl = FieldValidityList(date_time=datetimes,
                                basis=basis)

    else:
        raise TypeError("'datetimes' must be a datetime.datetime or a list of.")

    return fvl

def ifNone_emptydict(arg):
    """ Transforms a None into a {}. """
    if arg is None:
        arg = {}
    return arg

def set_DateHour_axis(axis, datetimerange,
                      showgrid=True,
                      datefmt=None,
                      xtickslabelsrotation=30.):
    """
    Set an adequate axis ticks and ticks labels for Date/Hour axis.
    
    *datetimerange* supposed to be a :class:`datetime.timedelta` instance
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    dayhourformatter = mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S')
    dayformatter = mdates.DateFormatter('%Y-%m-%d')
    if datetimerange <= datetime.timedelta(2):
        major_locator = mdates.HourLocator(interval=6)
        minor_locator = mdates.HourLocator(interval=1)
        formatter = mdates.AutoDateFormatter(major_locator)
    elif datetimerange <= datetime.timedelta(7):
        major_locator = mdates.DayLocator(interval=1)
        minor_locator = mdates.HourLocator(interval=6)
        formatter = dayhourformatter
    elif datetimerange <= datetime.timedelta(21):
        major_locator = mdates.DayLocator(interval=2)
        minor_locator = mdates.DayLocator(interval=1)
        formatter = dayhourformatter
    elif datetimerange <= datetime.timedelta(100):
        major_locator = mdates.DayLocator(interval=7)
        minor_locator = mdates.DayLocator(interval=1)
        formatter = dayformatter
    else:
        major_locator = mdates.AutoDateLocator()
        minor_locator = None
        formatter = mdates.AutoDateFormatter(major_locator)
    if datefmt is not None:
        formatter = mdates.DateFormatter(datefmt)
    axis.xaxis.set_major_locator(major_locator)
    axis.xaxis.set_major_formatter(formatter)
    axis.grid(showgrid)
    if minor_locator is not None:
        axis.xaxis.set_minor_locator(minor_locator)
        axis.grid(showgrid, which='minor', axis='x', color='grey')
    if xtickslabelsrotation != 0.:
        _ax = plt.gca()
        plt.sca(axis)
        plt.xticks(rotation=xtickslabelsrotation)
        plt.sca(_ax)

def set_figax(figure, ax, figsize=config.plotsizes):
    """
    Given existing matplotlib *figure* and an *ax* (or None),
    check consistency or generate a consistent (figure, ax) duet.
    """
    import matplotlib.pyplot as plt

    if ax is not None and figure is None:
        figure = ax.figure
    elif ax is None and figure is not None:
        if len(figure.axes) > 0:
            ax = figure.axes[0]
        else:
            ax = figure.gca()
    elif ax is not None and figure is not None:
        if ax not in figure.axes:
            raise epygramError('*over*: inconsistency between given fig and ax')
    elif figure is ax is None:
        figure, ax = plt.subplots(1, 1, figsize=figsize)

    return (figure, ax)

def set_map_up(bm, ax,
               drawrivers=False,
               drawcoastlines=True,
               drawcountries=True,
               meridians='auto',
               parallels='auto',
               departments=False,
               boundariescolor='0.25',
               bluemarble=0.0,
               background=False):
    """Cf. :meth:`H2DField.plotfield` documentation."""

    if background:
        bm.drawmapboundary(fill_color='lightskyblue', ax=ax)
        bm.fillcontinents(color='wheat', lake_color='skyblue',
                          zorder=0, ax=ax)
    if bluemarble:
        bm.bluemarble(alpha=bluemarble, ax=ax)
    if drawcoastlines:
        bm.drawcoastlines(color=boundariescolor, ax=ax)
    if departments:
        import json
        with open(config.installdir + '/data/departments.json', 'r') as dp:
            depts = json.load(dp)[1]
        for d in range(len(depts)):
            for part in range(len(depts[d])):
                dlon = depts[d][part][0]
                dlat = depts[d][part][1]
                (x, y) = bm(dlon, dlat)
                bm.plot(x, y, color=boundariescolor, ax=ax)
    elif drawcountries:
        bm.drawcountries(color=boundariescolor, ax=ax)
    if drawrivers:
        bm.drawrivers(color='blue', ax=ax)
    add_meridians_and_parallels_to(bm,
                                   parallels=parallels,
                                   meridians=meridians,
                                   ax=ax)

def datetimerange(start, stop=None, step=1, stepunit='h', tzinfo=None):
    """
    A generator of datetime.datetime objects ranging from *start* to *stop*
    (included) by *step*.
    
    Arguments syntax:\n
    - *start* and *stop* being either:\n
      - a string: 'YYYYMMDDhhmmssx', hh, mm, ss and x being optional (x = microseconds)
                  or a date/time in ISO 8601 format (cf. datetime.datetime.isoformat())
      - a tuple or list: (year, month, day[, hour[, minute[, seconde[, microsecond]]]])
      - a datetime.datetime instance
      if *stop* is None, returns [datetime(start)]
    - *step* being either an integer, which unit is specified in *stepunit*
      or a datetime.timedelta instance
    - *stepunit* among ('D', 'h', 'm', 's', 'x')
    - *tzinfo*: time zone info, cf. datetime.datetime
    """

    def parse_iterable(i):
        return datetime.datetime(*i, tzinfo=tzinfo)

    def parse_str(s):
        try:
            from dateutil.parser import parse
            dt = parse(s)
        except (ImportError, ValueError):
            hour = 0
            minute = 0
            second = 0
            microsecond = 0
            try:
                year = int(s[:4])
                month = int(s[4:6])
                day = int(s[6:8])
                if len(s) >= 10:
                    hour = int(s[8:10])
                if len(s) >= 12:
                    minute = int(s[10:12])
                if len(s) >= 12:
                    second = int(s[12:14])
                if len(s) > 12:
                    microsecond = int(s[12:])
            except ValueError:
                raise ValueError('please check syntax of date/time string.')
            dt = parse_iterable((year, month, day,
                                 hour, minute, second, microsecond))
        return dt

    if not isinstance(start, datetime.datetime):
        if isinstance(start, str):
            start = parse_str(start)
        elif isinstance(start, list) or isinstance(start, tuple):
            start = parse_iterable(start)
        else:
            raise TypeError("unknown type for *start*: " + str(type(start)))
    if not isinstance(stop, datetime.datetime):
        if stop is None:
            stop = start
        elif isinstance(stop, str):
            stop = parse_str(stop)
        elif isinstance(stop, list) or isinstance(stop, tuple):
            stop = parse_iterable(stop)
        else:
            raise TypeError("unknown type for *stop*: " + str(type(stop)))
    if not isinstance(step, datetime.timedelta):
        if isinstance(stop, str):
            step = int(step)
        assert isinstance(step, int)
        assert stepunit in ('D', 'h', 'm', 's', 'x')
        if stepunit == 'D':
            step = datetime.timedelta(step)
        elif stepunit == 'h':
            step = datetime.timedelta(0, step * 3600)
        elif stepunit == 'm':
            step = datetime.timedelta(0, step * 60)
        elif stepunit == 's':
            step = datetime.timedelta(0, step)
        elif stepunit == 'x':
            step = datetime.timedelta(0, microseconds=step)

    if start < stop:
        assert step > datetime.timedelta(0), 'step must be > 0 for start < stop'
    elif start > stop:
        assert step < datetime.timedelta(0), 'step must be < 0 for start > stop'

    rng = [start]
    dt = start + step
    if start < stop:
        while dt <= stop:
            rng.append(dt)
            dt += step
    elif start > stop:
        while dt >= stop:
            rng.append(dt)
            dt += step  # step < 0

    return rng
