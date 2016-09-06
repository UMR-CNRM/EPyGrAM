#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains some base classes of *epygram*.
"""

import numpy
import datetime
import copy
import sys

import footprints
from footprints import FootprintBase, FPDict

from epygram import epygramError, config
from epygram.util import RecursiveObject, separation_line, \
                         write_formatted

epylog = footprints.loggers.getLogger(__name__)



class Field(RecursiveObject, FootprintBase):
    """
    Generic abstract class implementing a Field, composed of an identifier and
    a data.

    The field identifier *fid* identifies a field with a set of keys.
    Each key (named after the format name) idetifies the field for a given format.
    A specific key 'generic' is a GRIB2-like description.
    E.g. *{'FA':'SURFTEMPERATURE',
           'generic':{'typeOfFirstFixedSurface':1, 'discipline':2, 'parameterCategory':3, 'parameterNumber':18}}*.
    """

    _collector = ('field',)
    _abstract = True
    _footprint = dict(
        attr=dict(
            fid=dict(
                type=FPDict,
                access='rwx'),
            data=dict(
                optional=True,
                default=None),
            comment=dict(
                type=str,
                optional=True,
                access='rwd')
        )
    )

    def getdata(self):
        """
        Returns the field data.
        Generic, default method for inheriting classes that do not overwrite it.
        """
        return self.data

    def setdata(self, data):
        """
        Sets or overwrites the field data as a numpy array.
        Mainly useful because a footprints attribute cannot be a numpy array...
        """

        if not isinstance(data, numpy.ndarray):
            data = numpy.array(data)
        self._attributes['data'] = data

    def setfid(self, fid):
        """
        Sets or overwrites the field fid given as a dict.
        """

        if not isinstance(fid, dict):
            raise epygramError("*fid* must be a dict.")
        self._attributes['fid'] = fid

    def clone(self, fid=None):
        """
        Returns a cloned field, optionally with a new *fid* given as a dict.
        """

        clone = self.deepcopy()
        if fid != None:
            clone.setfid(fid)
        return clone

###################
# pre-applicative #
###################
    def stats(self, **kwargs):
        """
        Computes some basic statistics on the field, as a dict containing:
        {'min', 'max', 'mean', 'std', 'quadmean', 'nonzero'}.

        Optional arguments can be passed, depending on the inheriting class,
        passed to getdata().
        """

        return {'min':self.min(**kwargs),
                'max':self.max(**kwargs),
                'mean':self.mean(**kwargs),
                'std':self.std(**kwargs),
                'quadmean':self.quadmean(**kwargs),
                'nonzero':self.nonzero(**kwargs)}

    def min(self, **kwargs):
        """Returns the minimum value of data."""
        data = self.getdata(**kwargs)
        return numpy.ma.masked_outside(data,
                                       - config.mask_outside,
                                       config.mask_outside).min()

    def max(self, **kwargs):
        """Returns the maximum value of data."""
        data = self.getdata(**kwargs)
        return numpy.ma.masked_outside(data,
                                       - config.mask_outside,
                                       config.mask_outside).max()

    def mean(self, **kwargs):
        """Returns the mean value of data."""
        data = self.getdata(**kwargs)
        return numpy.ma.masked_outside(data,
                                       - config.mask_outside,
                                       config.mask_outside).mean()

    def std(self, **kwargs):
        """Returns the standard deviation of data."""
        data = self.getdata(**kwargs)
        return numpy.ma.masked_outside(data,
                                       - config.mask_outside,
                                       config.mask_outside).std()

    def quadmean(self, **kwargs):
        """Returns the quadratic mean of data."""
        data = self.getdata(**kwargs)
        return numpy.sqrt((numpy.ma.masked_outside(data,
                                                   - config.mask_outside,
                                                   config.mask_outside) ** 2).mean())

    def nonzero(self, **kwargs):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        data = self.getdata(**kwargs)
        return numpy.count_nonzero(abs(numpy.ma.masked_outside(data,
                                                               - config.mask_outside,
                                                               config.mask_outside)) > config.epsilon)

#############
# OPERATORS #
#############

    def operation(self, operation, operand=None):
        """
        Makes the requested operation on the field.
        
        Implemented *operation* :
        '+', '-', '*', '/',
        'normalize', 'ceil', 'exp', 'log'... and you can try with every other
        **numpy** function.
        """

        if operand is not None:
            if isinstance(operand, self.__class__):
                self.operation_with_other(operation, operand)
            else:
                self.scalar_operation(operation, operand)
        else:
            if operation == 'normalize':
                self.setdata((self.data - self.data.min()) / \
                             (self.data.max() - self.data.min())
                             )
            else:
                try:
                    self.setdata(getattr(numpy, operation)(self.data))
                except Exception:
                    raise

    def operation_with_other(self, operation, other):
        """
        Makes a field-with-other operation, among:
        ('+', '-', '*', '/').
        """

        self._check_operands(other)

        if operation == '+':
            self.setdata(self.data + other.data)
        elif operation == '*':
            self.setdata(self.data * other.data)
        elif operation == '-':
            self.setdata(self.data - other.data)
        elif operation == '/':
            self.setdata(self.data / other.data)

    def scalar_operation(self, operation, scalar):
        """
        Makes a scalar operation on field:
        
        *operation* being one of ('+', '-', '*', '/')
        
        *scalar* being a float or int.
        """

        self._check_operands(scalar)

        if operation == '+':
            self.setdata(self.data + scalar)
        elif operation == '*':
            self.setdata(self.data * scalar)
        elif operation == '-':
            self.setdata(self.data - scalar)
        elif operation == '/':
            self.setdata(self.data / scalar)

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """

        if not isinstance(other, self.__class__):
            try:
                other = float(other)
            except Exception:
                raise ValueError("operations on " + self.__class__.__name__
                                 + " must involve either scalars \
                                    (integer/float) or "
                                 + self.__class__.__name__ + ".")
        else:
            if numpy.shape(self.data) != numpy.shape(other.data):
                raise epygramError("dimensions mismatch.")

    def _add(self, other, **kwargs):
        """
        Definition of addition, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'+'} and null validity.
        """

        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other.data
        else:
            rhs = other
        result = self.data + rhs
        newid = {'op':'+'}
        newfield = footprints.proxy.field(fid=newid,
                                          **kwargs)
        newfield.setdata(result)
        return newfield

    def _mul(self, other, **kwargs):
        """
        Definition of multiplication, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'*'} and null validity.
        """

        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other.data
        else:
            rhs = other
        result = self.data * rhs
        newid = {'op':'*'}
        newfield = footprints.proxy.field(fid=newid,
                                          **kwargs)
        newfield.setdata(result)
        return newfield

    def _sub(self, other, **kwargs):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """

        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other.data
        else:
            rhs = other
        result = self.data - rhs
        newid = {'op':'-'}
        newfield = footprints.proxy.field(fid=newid,
                                          **kwargs)
        newfield.setdata(result)
        return newfield

    def _div(self, other, **kwargs):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """

        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other.data
        else:
            rhs = other
        result = self.data / rhs
        newid = {'op':'/'}
        newfield = footprints.proxy.field(fid=newid,
                                          **kwargs)
        newfield.setdata(result)
        return newfield

    # default behaviors
    def __add__(self, other):
        return self._add(other)

    def __mul__(self, other):
        return self._mul(other)

    def __sub__(self, other):
        return self._sub(other)

    def __div__(self, other):
        return self._div(other)

    __radd__ = __add__
    __rmul__ = __mul__
    __rdiv__ = __div__
    __rsub__ = __sub__


class FieldSet(RecursiveObject, list):
    """
    Handles a set of Fields, in the manner of Python's builtin list,
    with some extra features, especially ensuring its components all are Fields.

    Constructor optional argument *fields* has to be either a :class:`Field`
    or an iterable of.
    """

    def __init__(self, fields=()):
        """
        Constructor.
        Checks that optional 'fields' argument is actually iterable and
        contains Field instances, or is a single Field.
        """

        if fields == ():
            pass
        elif isinstance(fields, Field):
            fields = (fields,)
        else:
            try:
                for item in fields:
                    if not isinstance(item, Field):
                        raise epygramError("A FieldSet can only be made out" + \
                                           " of Field instances.")
            except TypeError:
                raise epygramError("'fields' argument must be either a" + \
                                   " Field instance or an iterable of fields.")
            except Exception:
                raise

        super(FieldSet, self).__init__(fields)

    def __setitem__(self, position, field):

        if not isinstance(field, Field):
            raise epygramError("A FieldSet can contain only Field instances.")
        super(FieldSet, self).__setitem__(position, field)

    def __setslice__(self, pos1, pos2, fieldset):

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        super(FieldSet, self).__setslice__(pos1, pos2, fieldset)

    def append(self, field):
        """
        Checks that *field* is a :class:`Field` instance before appending it.
        """

        if not isinstance(field, Field):
            raise epygramError("A FieldSet can contain only Field instances.")
        super(FieldSet, self).append(field)

    def index(self, fid):
        """
        Returns the index of the first field of the FieldSet whose fid 
        matches *fid*, *fid* being a dict.
        """

        if not isinstance(fid, dict):
            raise ValueError("'fid' must be a dict.")

        idx = None
        for f in range(0, len(self)):
            if self[f].fid == fid:
                idx = f
                break
        return idx

    def extend(self, fieldset):
        """
        Checks that *fieldset* is a :class:`FieldSet` instance before extending
        with it.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        super(FieldSet, self).extend(fieldset)

    def insert(self, position, field):
        """
        Checks that *field* is a :class:`Field` instance before inserting it.
        """

        if not isinstance(field, Field):
            raise epygramError("A FieldSet can contain only Field instances.")
        super(FieldSet, self).insert(position, field)

    def remove(self, fid):
        """
        Removes from the FieldSet the first field whose fid matches *fid*,
        *fid* being a dict.
        """

        try:
            idx = self.index(fid)
            del self[idx]
        except Exception:
            pass

    def sort(self, attribute, key=None, reverse=False):
        """
        Sorts the fields of the FieldSet by the increasing criterion.

        If attribute is a string, sorting will be done according to
        *field.attribute[key]* or *field.attribute* (if *key==None*).

        If attribute is a list *[a1, a2...]*, sorting will be done according to
        *field.a1.a2[key]* or *field.a1.a2* (if *key==None*).

        If *reverse* is *True*, sorts by decreasing order.
        """

        if isinstance(attribute, str):
            if key == None:
                cmpfct = lambda x, y: cmp(x._attributes[attribute],
                                         y._attributes[attribute])
            else:
                cmpfct = lambda x, y: cmp(x._attributes[attribute][key],
                                         y._attributes[attribute][key])
        elif isinstance(attribute, list):
            a = attribute
            if isinstance(self[0]._attributes[a[0]], FootprintBase):
                if len(attribute) == 2:
                    if key == None:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]]._attributes[a[1]],
                                                 y._attributes[a[0]]._attributes[a[1]])
                    else:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]]._attributes[a[1]][key],
                                                 y._attributes[a[0]]._attributes[a[1]][key])
                elif len(attribute) == 3:
                    if key == None:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]]._attributes[a[1]].__dict__[a[2]],
                                                 y._attributes[a[0]]._attributes[a[1]].__dict__[a[2]])
                    else:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]]._attributes[a[1]].__dict__[a[2]][key],
                                                 y._attributes[a[0]]._attributes[a[1]].__dict__[a[2]][key])
                else:
                    raise NotImplementedError("len(attribute) > 3.")
            else:
                if len(attribute) == 2:
                    if key == None:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]].__dict__[a[1]],
                                                 y._attributes[a[0]].__dict__[a[1]])
                    else:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]].__dict__[a[1]][key],
                                                 y._attributes[a[0]].__dict__[a[1]][key])
                elif len(attribute) == 3:
                    if key == None:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]].__dict__[a[1]].__dict__[a[2]],
                                                 y._attributes[a[0]].__dict__[a[1]].__dict__[a[2]])
                    else:
                        cmpfct = lambda x, y: cmp(x._attributes[a[0]].__dict__[a[1]].__dict__[a[2]][key],
                                                 y._attributes[a[0]].__dict__[a[1]].__dict__[a[2]][key])
                else:
                    raise NotImplementedError("len(attribute) > 3.")

        else:
            raise TypeError("attribute must be a string or list of string.")

        super(FieldSet, self).sort(cmp=cmpfct, reverse=reverse)

    def listfields(self, fidkey=None):
        """
        Returns a list of the identifiers of the FieldSet.
        If *fidkey* is supplied, the list contains only **fid[*fidkey*]**,
        and not whole fid.
        """

        if fidkey is None:
            fieldslist = [f.fid for f in self]
        else:
            fieldslist = [f.fid[fidkey] for f in self]

        return fieldslist

    def filter(self, by, criteria):
        """
        Not Implemented Yet.

        Returns a new FieldSet filtered according to criteria specified in
        argument.

        Args:

        - by: the kind of filter; on what to filter ?
        - criteria: how to filter on that ?

        Available filters, examples:

        - by='id', criteria={typefmt:identifier} will return only fields whose
          id[typefmt] match value...
        """

        # TODO: id=, fieldtype=, spectral, (...)
        raise NotImplementedError("not yet...")



class Resource(RecursiveObject, FootprintBase):
    """Generic abstract class implementing a Resource."""

    _abstract = True
    _collector = ('epyresource',)
    _footprint = dict(
        attr=dict(
            openmode=dict(
                values=set(['r', 'read', 'w', 'write', 'a', 'append']),
                remap=dict(
                    read='r',
                    write='w',
                    append='a'),
                info="Opening mode.")
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""

        super(Resource, self).__init__(*args, **kwargs)

    def __enter__(self):
        """Context enter."""
        return self

    def __exit__(self, t, v, tbk):
        """Context exit."""
        self.close()

    def __del__(self):
        """Destructor. Closes the resource properly."""

        try:
            self.close()
        except Exception as e:
            epylog.warning(
                "Exception catched in epygram.base.Resource.__del__(): " + \
                str(e))

    def __len__(self):
        """Returns the number of fields in resource."""
        return len(self.listfields())

    def __contains__(self, fid):
        return fid in self.listfields()

    def __iter__(self):
        """Caution: iteration returns the fields one by one, read from file !"""
        for f in self.listfields():
            yield self.readfield(f)

    def readfields(self, requestedfields, getdata=True):
        """
        Returns a :class:`FieldSet` containing requested fields read in the
        resource.

        Args:

        - *requestedfields*: a field identifier of the resource format, or a
          list of.
        - *getdata*: optional, if *False*, only metadata are read, the fields
          do not contain data. Default is *True*.
        """

        fieldset = FieldSet()
        if isinstance(requestedfields, list):
            for f in requestedfields:
                fieldset.append(self.readfield(f, getdata=getdata))
        else:
            fieldset.append(self.readfield(requestedfields, getdata=getdata))

        return fieldset

    def writefields(self, fieldset):
        """
        Write the fields of the 'fieldset' in the resource;
        *fieldset* must be a :class:`FieldSet` instance.
        """

        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be a FieldSet " + \
                               "instance.")
        for field in fieldset:
            self.writefield(field)

    def listfields(self, **kwargs):
        """
        Returns a list containing the identifiers (in the resource format)
        of all the fields of the resource.

        (Generic wrapper with buffering if openmode == 'r'.)
        """

        if self.openmode == 'r' and not hasattr(self, '_bufferedlistfields'):
            self._bufferedlistfields = []
        elif self.openmode in ('w', 'a'):
            # list of fields subject to evolution; no buffering
            self._bufferedlistfields = None

        if self._bufferedlistfields and self._bufferedlistfieldsoptions != kwargs:
            self._bufferedlistfields = []

        if not self._bufferedlistfields:
            fieldslist = self._listfields(**kwargs)
            self._bufferedlistfieldsoptions = copy.deepcopy(kwargs)
            if self._bufferedlistfields != None:
                self._bufferedlistfields = fieldslist  # save
        else:
            # immutable and already read
            fieldslist = self._bufferedlistfields

        return fieldslist

    def _listfields(self):
        """
        Actual listfields() method (virtual).
        """
        raise NotImplementedError("virtual method.")

    def find_fields_in_resource_by_generic_fid(self, handgrip):
        """
        Find in resource the fields whose generic fid (if the resource is able
        to give one) matches the *handgrip*.
        """

        fieldslist = self.listfields(complete=True)
        found = []
        for f in fieldslist:
            if all([f['generic'][k] == handgrip[k] for k in handgrip.keys()]):
                found.append(f)

        return found



class FieldValidity(RecursiveObject):
    """
    This class handles a uniq temporal validity for a meteorological field:
    its date and time of validity (*date_time*), as well as the validity of
    its origin (*basis*, i.e. for a forecast field for instance, the beginning
    of the forecast) and its *term*.

    An additional optional *cumulativeduration* parameter can define
    the duration for which cumulative fields (e.g. precipitation) are valid.
    If such, optional *statistical_process_on_duration* can be supplied as a
    string or an int to describe the kind of statistical process that runs over the
    *cumulativeduration*.

    Constructor arguments: cf. *set()* method.
    """

    def __init__(self, date_time=None, basis=None, term=None,
                 cumulativeduration=None,
                 statistical_process_on_duration=None):
        """
        Constructor.
        Args:
        - date_time: has to be of type datetime.datetime;
        - basis: has to be of type datetime.datetime;
        - term: has to be of type datetime.timedelta;
        - cumulativeduration: has to be of type datetime.timedelta;
        - statistical_process_on_duration: kind of statistical process
          that runs over the cumulative duration.
        """

        self._basis = None
        self._date_time = None
        self._cumulativeduration = None
        self._statistical_process_on_duration = None

        kwargs = dict(date_time=date_time,
                      basis=basis,
                      term=term,
                      cumulativeduration=cumulativeduration,
                      statistical_process_on_duration=statistical_process_on_duration)
        if not (date_time == None and basis == None and term == None):
            self.set(**kwargs)

    def term(self, fmt=None):
        """
        This method returns the term as the difference between date and time
        of validity and basis.

        By default, it is returned as a :class:`datetime.timedelta`;
        otherwise, *fmt* argument can specify the desired return format.

        Coded versions of *fmt*: 'IntHours', 'IntSeconds', and that's all for
        now...
        """

        if fmt == None:
            out = self._date_time - self._basis
        elif fmt == 'IntHours':
            term = self._date_time - self._basis
            out = int(term.total_seconds() / 3600)
        elif fmt == 'IntSeconds':
            term = self._date_time - self._basis
            out = int(term.total_seconds())
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " + \
                                      self.__class__.__name__ + ".term().")

        return out

    def cumulativeduration(self, fmt=None):
        """
        This method returns the cumulative duration,
        i.e. the duration for which cumulative fields (e.g. precipitation) are
        valid.

        By default, it is returned as a :class:`datetime.timedelta`;
        otherwise, *fmt* argument can specify the desired return format.

        Coded versions of *fmt*: 'IntHours', 'IntSeconds', and that's all for
        now...
        """

        if fmt == None:
            out = self._cumulativeduration
        elif fmt == 'IntHours':
            out = int(self._cumulativeduration.total_seconds() / 3600)
        elif fmt == 'IntSeconds':
            out = int(self._cumulativeduration.total_seconds())
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " + \
                                      self.__class__.__name__ + \
                                      ".cumulativeduration().")

        return out

    def statistical_process_on_duration(self, asGRIB2code=False):
        """
        If the field describes a cumulative process over a cumulativeduration,
        returns the kind of statistical process that runs over the duration.
        
        If *asGRIB2code*, returned as a GRIB2 code (cf. GRIB2 table 4.10).
        """
        from .formats import grib_utilities

        if not asGRIB2code and isinstance(self._statistical_process_on_duration, int):
            out = grib_utilities.statistical_processes.get(self._statistical_process_on_duration, None)
        elif asGRIB2code and isinstance(self._statistical_process_on_duration, str):
            out = {v:k for k, v in grib_utilities.statistical_processes.items()}.get(self._statistical_process_on_duration, None)
        else:
            out = self._statistical_process_on_duration

        return out

    def get(self, fmt=None):
        """
        Returns the date and time of validity.

        By default, as a :class:`datetime.datetime`;
        otherwise, *fmt* argument can specify the desired return format.

        Coded versions of *fmt*:
        'IntStr' (e.g. '20140731104812' = 2014 july 31th at 10h, 48m, 12s).
        And that's all for now...
        """

        if fmt == None:
            out = self._date_time
        elif fmt == 'IntStr':
            out = '{:0>{width}}'.format(str(self._date_time.year), width=4) \
                + '{:0>{width}}'.format(str(self._date_time.month), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.day), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.hour), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.minute), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.second), width=2)
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " + \
                                      self.__class__.__name__ + ".get().")

        return out

    def getbasis(self, fmt=None):
        """
        Returns the date and time of origin (basis).

        By default, as a :class:`datetime.datetime`;
        otherwise, *fmt* argument can specify the desired return format.

        Coded versions of *fmt*:
        'IntStr' (e.g. '20140731104812' = 2014 july 31th at 10h, 48m, 12s).
        And that's all for now...
        """

        if fmt == None:
            out = self._basis
        elif fmt == 'IntStr':
            out = '{:^{width}}'.format(str(self._basis.year), width=4) \
                + '{:0>{width}}'.format(str(self._basis.month), width=2) \
                + '{:0>{width}}'.format(str(self._basis.day), width=2) \
                + '{:0>{width}}'.format(str(self._basis.hour), width=2) \
                + '{:0>{width}}'.format(str(self._basis.minute), width=2) \
                + '{:0>{width}}'.format(str(self._basis.second), width=2)
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " + \
                                      self.__class__.__name__ + ".getbasis().")

        return out

    def set(self, date_time=None, basis=None, term=None,
            cumulativeduration=None,
            statistical_process_on_duration=None):
        """
        Sets validity and basis according to arguments.
        A consistency check is done if the three arguments are provided
        (which is useless anyway).

        Args: \n
        - *date_time*: has to be of type :class:`datetime.datetime`;
        - *basis*: has to be of type :class:`datetime.datetime`;
        - *term*: has to be of type :class:`datetime.timedelta`;
        - *cumulativeduration*: has to be of type :class:`datetime.timedelta`;
        - *statistical_process_on_duration*: kind of statistical process
          that runs over the cumulative duration.
        """

        if isinstance(date_time, datetime.datetime):
            self._date_time = date_time
        elif date_time != None:
            raise epygramError("argument 'date_time' must be of type" + \
                               " datetime.datime")
        if isinstance(basis, datetime.datetime):
            self._basis = basis
        elif basis != None:
            raise epygramError("argument 'basis' must be of type" + \
                               " datetime.datime")
        if term != None and not isinstance(term, datetime.timedelta):
            raise epygramError("argument 'term' must be of type" + \
                               " datetime.timedelta")
        if cumulativeduration != None and\
           not isinstance(cumulativeduration, datetime.timedelta):
            raise epygramError("argument 'cumulativeduration' must be of" + \
                               " type datetime.timedelta")

        if isinstance(term, datetime.timedelta):
            if date_time != None and basis != None and term != None \
               and date_time - basis != term:
                raise epygramError("inconsistency between 'term', 'basis'" + \
                                   " and 'date_time' arguments.")

            if self._date_time == None:
                if self._basis == None:
                    raise epygramError("cannot set 'term' without 'basis'" + \
                                       " nor 'date_time'.")
                else:
                    self._date_time = self._basis + term
            else:
                if self._basis == None:
                    self._basis = self._date_time - term
                else:
                    self._date_time = self._basis + term

        if cumulativeduration != None:
            self._cumulativeduration = cumulativeduration
        if self._cumulativeduration != None and  statistical_process_on_duration != None:
            self._statistical_process_on_duration = statistical_process_on_duration


class FieldValidityList(RecursiveObject, list):
    """
    This class handles a list of temporal validity.
    """

    def __init__(self, validity_instance=None, length=1, **kwargs):
        """
        Constructor.
        
        - *validity_instance*, if given is an instance of FieldValidity
        - *length*, to build a series of validities from either the
        *validity_instance* or from an uninitialized one.
        - other kwargs: same as :class:`FieldValidity` constructor.
        
        """

        super(list, self).__init__([])

        if validity_instance is not None:
            if len(kwargs) != 0:
                raise epygramError("One can not give, at the same time, validity_instance and other argument.")
            if isinstance(validity_instance, FieldValidity):
                self.append(validity_instance)
                if length > 1:
                    for _ in range(length - 1):
                        self.append(validity_instance.deepcopy())
            elif isinstance(validity_instance, FieldValidityList):
                self.extend(validity_instance)
            else:
                raise epygramError("FieldValidityList must be built from FieldValidity or from FieldValidityList instances.")
        elif kwargs != {}:
            #Check that all lengths are equal
            length = None
            mykwargs = {}
            for k, v in kwargs.iteritems():
                mykwargs[k] = [v] if type(v) != type(list()) else v
                if length == None or length == 1:
                    length = len(mykwargs[k])
                if len(mykwargs[k]) != length:
                    raise epygramError("All the arguments must have the same length.")

            for k, v in mykwargs.iteritems():
                if len(v) == 1:
                    mykwargs[k] = mykwargs[k] * length

            #We set the different objects
            if length == None:
                length = 1
            self.extend([FieldValidity(**{key: value[i] for (key, value) in mykwargs.iteritems()}) for i in range(length)])
        elif isinstance(length, int):
            for _ in range(length):
                self.append(FieldValidity())

    def term(self, one=True, **kwargs):
        """This method returns the terms of all the validities"""

        length = len(self)
        result = [self[i].term(**kwargs) for i in range(length)]
        return result[0] if (one and length == 1) else result

    def cumulativeduration(self, one=True, **kwargs):
        """This method returns the cumulative duration of all the validities."""

        length = len(self)
        result = [self[i].cumulativeduration(**kwargs) for i in range(length)]
        return result[0] if (one and length == 1) else result

    def statistical_process_on_duration(self, one=True, **kwargs):
        """This method returns the statistical process on duration of all the validities."""

        length = len(self)
        result = [self[i].statistical_process_on_duration(**kwargs) for i in range(length)]
        return result[0] if (one and length == 1) else result

    def get(self, one=True, **kwargs):
        """Returns the date and time of all the validities."""

        length = len(self)
        result = [self[i].get(**kwargs) for i in range(length)]
        return result[0] if (one and length == 1) else result

    def getbasis(self, one=True, **kwargs):
        """Returns the date and time of origin (basis) of all the validities."""

        length = len(self)
        result = [self[i].getbasis(**kwargs) for i in range(length)]
        return result[0] if (one and length == 1) else result

    def set(self, **kwargs):
        """Sets validity objects"""

        #Check that all lengths are equal
        length = None
        mykwargs = {}
        for k, v in kwargs.iteritems():
            mykwargs[k] = [v] if type(v) != type(list()) else v
            if length == None or length == 1:
                length = len(mykwargs[k])
            if len(mykwargs[k]) != length:
                raise epygramError("All the arguments must have the same length.")

        if length == 1:
            length = len(self)

        for k, v in mykwargs.iteritems():
            if len(v) == 1:
                mykwargs[k] = mykwargs[k] * length

        #We set the different objects
        for i in range(length):
            self[i].set(**{key: value[i] for (key, value) in mykwargs.iteritems()})

    def what(self, out=sys.stdout, cumulativeduration=True):
        """
        Writes in file a summary of the validity.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *cumulativeduration*: if False, not written.
        """

        out.write("################\n")
        out.write("### VALIDITY ###\n")
        out.write("################\n")
        write_formatted(out, "Validity", self.get())
        write_formatted(out, "Basis", self.getbasis())
        try:
            write_formatted(out, "Term", self.term())
        except TypeError:
            pass
        if cumulativeduration:
            write_formatted(out, "Duration for cumulative quantities",
                            self.cumulativeduration())
            write_formatted(out, "Statistical process",
                            self.statistical_process_on_duration())
            write_formatted(out, "(GRIB2 code -- table 4.10)",
                            self.statistical_process_on_duration(asGRIB2code=True))
        out.write(separation_line)
        out.write("\n")
