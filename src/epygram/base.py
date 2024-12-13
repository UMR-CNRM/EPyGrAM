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
import hashlib

import footprints
from footprints import FootprintBase, FPDict

from epygram import epygramError, config
from epygram.util import RecursiveObject, separation_line, write_formatted

epylog = footprints.loggers.getLogger(__name__)


class Field(RecursiveObject, FootprintBase):
    """
    Generic abstract class implementing a Field, composed of an identifier and
    a data.

    The field identifier *fid* identifies a field with a set of keys.
    Each key (named after the format name) identifies the field for a given format.
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
            comment=dict(
                optional=True,
                access='rwd'),
            misc_metadata=dict(
                type=FPDict,
                optional=True,
                default=FPDict(),
                access='rwd'),
            units=dict(
                optional=True,
                access='rwd',
                default='')
        )
    )

    def __init__(self, *args, **kwargs):
        """Constructor. See its footprint for arguments."""
        super(Field, self).__init__(*args, **kwargs)
        self._data = None

    def getdata(self):
        """
        Returns the field data.
        Generic, default method for inheriting classes that do not overwrite it.
        """
        return self._data

    def setdata(self, data):
        """Sets or overwrites the field data as a numpy array."""
        if not isinstance(data, numpy.ndarray):
            data = numpy.array(data)
        self._data = data

    def deldata(self):
        """Empties the data."""
        self._data = None

    data = property(getdata, setdata, deldata, "Accessor to the field data.")

    def setfid(self, fid):
        """
        Sets or overwrites the field fid given as a dict.
        """
        if not isinstance(fid, dict):
            raise epygramError("**fid** must be a dict.")
        self._attributes['fid'] = fid

    def clone(self, fid=None):
        """
        Returns a cloned field, optionally with a new **fid** given as a dict.
        """
        clone = self.deepcopy()
        if fid is not None:
            clone.setfid(fid)
        return clone

###################
# pre-applicative #
###################
    def stats(self, **kwargs):
        """
        Computes some basic statistics on the field, as a dict containing:
        {'min', 'max', 'mean', 'std', 'quadmean', 'nonzero'}.

        See each of these methods for details.

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
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(data.min())

    def max(self, **kwargs):
        """Returns the maximum value of data."""
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(data.max())

    def mean(self, **kwargs):
        """Returns the mean value of data."""
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(data.mean())

    def std(self, **kwargs):
        """Returns the standard deviation of data."""
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(data.std())

    def quadmean(self, **kwargs):
        """Returns the quadratic mean of data."""
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(numpy.sqrt((data ** 2).mean()))

    def absmean(self, **kwargs):
        """Returns the mean of absolute value of data."""
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return float(numpy.abs(data).mean())

    def nonzero(self, **kwargs):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        data = numpy.ma.masked_outside(self.getdata(**kwargs),
                                       - config.mask_outside,
                                       config.mask_outside)
        return int(numpy.count_nonzero(abs(data) > config.epsilon))

    def sha256_checksum(self, **kwargs):
        """
        Return a SHA256 checksum of the field data.
        """
        s256 = hashlib.sha256()
        s256.update(self.getdata(**kwargs).tobytes())
        return s256.hexdigest()

#############
# OPERATORS #
#############

    def operation(self, operation, operand=None):
        """
        Makes the requested operation on the field.

        :param operation: any of '+', '-', '*', '/',
                          or 'normalize', 'ceil', 'exp', 'log'...
                          and you can try with every other **numpy** function.
        :param operand: operand for the 4 basic operations, may be a scalar or
                        another field with according geometry.
        """
        if operand is not None:
            if isinstance(operand, self.__class__):
                self.operation_with_other(operation, operand)
            else:
                self.scalar_operation(operation, operand)
        else:
            if operation == 'normalize':
                self.setdata((self._data - self._data.min()) /
                             (self._data.max() - self._data.min())
                             )
            else:
                try:
                    self.setdata(getattr(numpy, operation)(self._data))
                except Exception:
                    raise

    def operation_with_other(self, operation, other):
        """
        Makes an in-place operation with another field.

        :param operation: among ('+', '-', '*', '/')
        :param other: another field, with according dimensions
        """
        self._check_operands(other)

        if operation == '+':
            self.setdata(self._data + other._data)
        elif operation == '*':
            self.setdata(self._data * other._data)
        elif operation == '-':
            self.setdata(self._data - other._data)
        elif operation == '/':
            self.setdata(self._data / other._data)

    def scalar_operation(self, operation, scalar):
        """
        Makes an in-place scalar operation on field.

        :param operation: among ('+', '-', '*', '/')
        :param scalar: a float or int.
        """
        self._check_operands(scalar)

        if operation == '+':
            self.setdata(self._data + scalar)
        elif operation == '*':
            self.setdata(self._data * scalar)
        elif operation == '-':
            self.setdata(self._data - scalar)
        elif operation == '/':
            self.setdata(self._data / scalar)

    def compare_to(self, other):
        """
        Compare a field to another one, with several criteria:

        - bias: average of errors distribution
        - std: standard deviation of errors distribution
        - errmax: maximum absolute error
        - common_mask: indicate whether fields have the same mask or not;
          beware that above statistics apply only to commonly unmasked data

        :return ({bias, std, errmax}, common_mask)
        """
        selfdata = self._masked_data()
        otherdata = self._masked_data()
        common_mask = not numpy.any(numpy.logical_xor(selfdata.mask, otherdata.mask))
        diff = self - other
        return ({'bias':diff.mean(),
                 'std':diff.std(),
                 'errmax':max(abs(diff.min()), abs(diff.max()))},
                common_mask)

    def normalized_comparison(self, ref):
        """
        Compare field to a reference, with prior normalization (by reference
        magnitude) of both fields.
        Hence the figures of comparison can be interpretated as percentages.
        """
        refmin = ref.min()
        refmax = ref.max()
        selfmin = self.min()
        selfmax = self.max()
        if abs(refmax - refmin) <= config.epsilon:
            # ref is constant
            if refmin <= config.epsilon:
                normalizedref = ref
                if abs(selfmax - selfmin) <= config.epsilon:
                    # test is also constant
                    if selfmin <= config.epsilon:
                        # both 0. : no normalization required
                        normalizedself = self
                    else:
                        # test is not 0. vs. ref is 0.
                        normalizedself = self / selfmin  # so that normalized error = 1.
                else:
                    # test is not constant but ref is constant 0. : normalize by itself
                    normalizedself = (self - selfmin).__div__(selfmax - selfmin)  # FIXME: classical operators seem to fail ?
            else:
                # ref is constant not 0.
                normalizedref = ref.__div__(refmin)
                normalizedself = self.__div__(refmin)
        else:
            # ref is not constant
            normalizedself = (self - refmin).__div__(refmax - refmin)  # FIXME: classical operators seem to fail ?
            normalizedref = (ref - refmin).__div__(refmax - refmin)  # FIXME: classical operators seem to fail ?
        return normalizedself.compare_to(normalizedref)

    def _masked_data(self, mask_outside=config.mask_outside,
                     **kwargs):
        """
        Return self field data as a masked array.

        :param mask_outside: if None, mask is empty;
                             else, mask data outside +/- this value
        """
        if isinstance(self._data, numpy.ma.masked_array):
            mdata = self.getdata(**kwargs)
        else:
            if mask_outside is not None:
                mdata = numpy.ma.masked_outside(self.getdata(**kwargs),
                                                -mask_outside,
                                                mask_outside)
            else:
                mdata = numpy.ma.masked_array(self.getdata(**kwargs))
        return mdata

    def _masked_any(self, other, mask_outside=config.mask_outside,
                    **kwargs):
        """
        Get a copy of self data and **other** data, with masked data where any
        of them is masked.

        :param mask_outside: if None, mask is empty;
                             else, mask data outside +/- this value
        """
        data = self._masked_data(mask_outside, **kwargs)
        otherdata = other._masked_data(mask_outside, **kwargs)
        cmask = numpy.logical_or(data.mask, otherdata.mask)
        data.mask = cmask
        otherdata.mask = cmask
        return data, otherdata

    def correlation(self, other,
                    commonmask=False,
                    mask_outside=config.mask_outside):
        """
        Compute a correlation coefficient R to another field.

        :param commonmask: if True, compute distance on the subset of point that
                           are not masked for any of the two fields.
        :param mask_outside: if None, mask is empty;
                             else, mask data outside +/- this value
        """
        # TODO: treat more complicated cases, where mask is present but commonmask is False,
        # or not present but to be masked without commonmask...
        from bronx.syntax.arrays import stretch_array
        otherdata = other.data
        selfdata = self.data
        if commonmask:
            if not isinstance(selfdata, numpy.ma.masked_array):
                selfdata = numpy.ma.masked_outside(selfdata,
                                                   -mask_outside,
                                                   mask_outside)
            if not isinstance(otherdata, numpy.ma.masked_array):
                otherdata = numpy.ma.masked_outside(otherdata,
                                                    -mask_outside,
                                                    mask_outside)
            cmask = numpy.logical_or(selfdata.mask, otherdata.mask)
            selfdata.mask = cmask
            otherdata.mask = cmask
        otherdata = stretch_array(otherdata)
        selfdata = stretch_array(selfdata)
        if not commonmask and otherdata.shape != selfdata.shape:
            raise Exception('inconsistency between masks')
        r = numpy.corrcoef(selfdata, otherdata)[0,1]
        return r

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """

        if isinstance(other, self.__class__) or isinstance(self, other.__class__):
            if numpy.shape(self._data) != numpy.shape(other._data):
                raise epygramError("dimensions mismatch.")
        else:
            try:
                other = float(other)
            except Exception:
                raise ValueError("operations on " + self.__class__.__name__ +
                                 " must involve either scalars " +
                                 "(integer/float) or " +
                                 self.__class__.__name__ + ".")

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
            rhs = other._data
        else:
            rhs = other
        result = self._data + rhs
        kwargs.setdefault('fid', {'op':'+'})
        newfield = footprints.proxy.field(**kwargs)
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
            rhs = other._data
        else:
            rhs = other
        result = self._data * rhs
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
            rhs = other._data
        else:
            rhs = other
        result = self._data - rhs
        kwargs.setdefault('fid', {'op':'-'})
        newfield = footprints.proxy.field(**kwargs)
        newfield.setdata(result)
        return newfield

    def _rsub(self, other, **kwargs):
        """
        Definition of reverse substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """
        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other._data
        else:
            rhs = other
        result = rhs - self._data
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
            rhs = other._data
        else:
            rhs = other
        result = self._data / rhs
        newid = {'op':'/'}
        newfield = footprints.proxy.field(fid=newid,
                                          **kwargs)
        newfield.setdata(result)
        return newfield

    def _rdiv(self, other, **kwargs):
        """
        Definition of reverse division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """
        self._check_operands(other)
        if isinstance(other, self.__class__):
            rhs = other._data
        else:
            rhs = other
        result = rhs / self._data
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

    def __rsub__(self, other):
        return self._rsub(other)

    def __rdiv__(self, other):
        return self._rdiv(other)


class FieldSet(RecursiveObject, list):
    """
    Handles a set of Fields, in the manner of Python's builtin list,
    with some extra features, especially ensuring its components all are Fields.

    Constructor optional argument **fields** has to be either a :class:`Field`
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
                        raise epygramError("A FieldSet can only be made out" +
                                           " of Field instances.")
            except TypeError:
                raise epygramError("'fields' argument must be either a" +
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
        matches **fid**, **fid** being a dict.
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
        Checks that **fieldset** is a :class:`FieldSet` instance before extending
        with it.
        """
        if not isinstance(fieldset, FieldSet):
            raise epygramError("'fieldset' argument must be of kind FieldSet.")
        super(FieldSet, self).extend(fieldset)

    def insert(self, position, field):
        """
        Checks that **field** is a :class:`Field` instance before inserting it
        at the **position**.
        """
        if not isinstance(field, Field):
            raise epygramError("A FieldSet can contain only Field instances.")
        super(FieldSet, self).insert(position, field)

    def remove(self, fid):
        """
        Removes from the FieldSet the first field whose fid matches **fid**,
        **fid** being a dict.
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
        *field.attribute[key]* or *field.attribute* (if **key** is **None**).

        If attribute is a list *[a1, a2...]*, sorting will be done according to
        *field.a1.a2[key]* or *field.a1.a2* (if *key==None*).

        If *reverse* is *True*, sorts by decreasing order.
        """
        if isinstance(attribute, str):
            if key is None:
                def cmpfct(x, y):
                    cmp(x._attributes[attribute],
                        y._attributes[attribute])
            else:
                def cmpfct(x, y):
                    cmp(x._attributes[attribute][key],
                        y._attributes[attribute][key])
        elif isinstance(attribute, list):
            a = attribute
            if isinstance(self[0]._attributes[a[0]], FootprintBase):
                if len(attribute) == 2:
                    if key is None:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]]._attributes[a[1]],
                                y._attributes[a[0]]._attributes[a[1]])
                    else:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]]._attributes[a[1]][key],
                                y._attributes[a[0]]._attributes[a[1]][key])
                elif len(attribute) == 3:
                    if key is None:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]]._attributes[a[1]].__dict__[a[2]],
                                y._attributes[a[0]]._attributes[a[1]].__dict__[a[2]])
                    else:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]]._attributes[a[1]].__dict__[a[2]][key],
                                y._attributes[a[0]]._attributes[a[1]].__dict__[a[2]][key])
                else:
                    raise NotImplementedError("len(attribute) > 3.")
            else:
                if len(attribute) == 2:
                    if key is None:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]].__dict__[a[1]],
                                y._attributes[a[0]].__dict__[a[1]])
                    else:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]].__dict__[a[1]][key],
                                y._attributes[a[0]].__dict__[a[1]][key])
                elif len(attribute) == 3:
                    if key is None:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]].__dict__[a[1]].__dict__[a[2]],
                                y._attributes[a[0]].__dict__[a[1]].__dict__[a[2]])
                    else:
                        def cmpfct(x, y):
                            cmp(x._attributes[a[0]].__dict__[a[1]].__dict__[a[2]][key],
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

        Examples of filters:

        - by='id', criteria={typefmt:identifier} will return only fields whose
          id[typefmt] match value...
        """
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

    def __exit__(self, t, v, tbk):  # @UnusedVariables
        """Context exit."""
        self.close()

    def __del__(self):
        """Destructor. Closes the resource properly."""

        try:
            self.close()
        except Exception as e:
            epylog.warning(
                "Exception catched in epygram.base.Resource.__del__(): " +
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

        :param requestedfields: a field identifier of the resource format, or a
                                list of.
        :param getdata: optional, if *False*, only metadata are read, the fields
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
            raise epygramError("'fieldset' argument must be a FieldSet " +
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
            if self._bufferedlistfields is not None:
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
    This class handles a unique temporal validity for a meteorological field:
    its date and time of validity (*date_time*), as well as the validity of
    its origin (*basis*, i.e. for a forecast field for instance, the beginning
    of the forecast) and its *term*.

    An additional optional *cumulativeduration* parameter can define
    the duration for which cumulative fields (e.g. precipitation) are valid.

    If such, optional *statistical_process_on_duration* can be supplied as a
    string or an int (cf. GRIB2 table typeOfStatisticalProcessing) to describe
    the kind of statistical process that runs over the *cumulativeduration*.

    Constructor arguments: cf. *set()* method.
    """

    def __init__(self,
                 date_time=None,
                 basis=None,
                 term=None,
                 cumulativeduration=None,
                 statistical_process_on_duration=None,
                 statistical_time_increment=None):
        """
        Constructor.

        :param date_time: has to be of type datetime.datetime;
        :param basis: has to be of type datetime.datetime;
        :param term: has to be of type datetime.timedelta;
        :param cumulativeduration: has to be of type datetime.timedelta;
        :param statistical_process_on_duration: kind of statistical process
            that runs over the cumulative duration.
        :param statistical_time_increment: time step over used for statistical
            process.
        """
        self._basis = None
        self._date_time = None
        self._cumulativeduration = None
        self._statistical_process_on_duration = None
        self._statistical_time_increment = None

        kwargs = dict(date_time=date_time,
                      basis=basis,
                      term=term,
                      cumulativeduration=cumulativeduration,
                      statistical_process_on_duration=statistical_process_on_duration,
                      statistical_time_increment=statistical_time_increment)
        if not (date_time is None and basis is None and term is None):
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
        if fmt is None:
            out = self._date_time - self._basis
        elif fmt == 'IntHours':
            term = self._date_time - self._basis
            out = int(term.total_seconds() // 3600)
        elif fmt == 'IntSeconds':
            term = self._date_time - self._basis
            out = int(term.total_seconds())
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " +
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
        if fmt is None:
            out = self._cumulativeduration
        elif fmt == 'IntHours':
            out = int(self._cumulativeduration.total_seconds() // 3600)
        elif fmt == 'IntSeconds':
            out = int(self._cumulativeduration.total_seconds())
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " +
                                      self.__class__.__name__ +
                                      ".cumulativeduration().")
        return out

    def statistical_process_on_duration(self, asGRIB2code=False):
        """
        If the field describes a cumulative process over a cumulativeduration,
        returns the kind of statistical process that runs over the duration.

        If *asGRIB2code*, returned as a GRIB2 code (cf. GRIB2 table 4.10).
        """
        from epygram.extra import griberies
        if not asGRIB2code and isinstance(self._statistical_process_on_duration, int):
            out = griberies.tables.statistical_processes.get(self._statistical_process_on_duration, None)
        elif asGRIB2code and isinstance(self._statistical_process_on_duration, str):
            out = {v:k for k, v in griberies.tables.statistical_processes.items()}.get(self._statistical_process_on_duration, None)
        else:
            out = self._statistical_process_on_duration
        return out

    def statistical_time_increment(self, fmt=None):
        """
        This method returns the statistical_time_increment,
        i.e. the time step used for statistical process over cumulative
        duration.

        By default, it is returned as a :class:`datetime.timedelta`;
        otherwise, *fmt* argument can specify the desired return format.

        Coded versions of *fmt*: 'IntHours', 'IntSeconds', and that's all for
        now...
        """
        if fmt is None:
            out = self._statistical_time_increment
        elif fmt == 'IntHours':
            out = int(self._statistical_time_increment.total_seconds() // 3600)
        elif fmt == 'IntSeconds':
            out = int(self._statistical_time_increment.total_seconds())
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " +
                                      self.__class__.__name__ +
                                      ".statistical_time_increment().")
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
        if fmt is None:
            out = self._date_time
        elif fmt == 'IntStr':
            out = '{:0>{width}}'.format(str(self._date_time.year), width=4) \
                + '{:0>{width}}'.format(str(self._date_time.month), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.day), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.hour), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.minute), width=2) \
                + '{:0>{width}}'.format(str(self._date_time.second), width=2)
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " +
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
        if fmt is None:
            out = self._basis
        elif fmt == 'IntStr':
            out = '{:^{width}}'.format(str(self._basis.year), width=4) \
                + '{:0>{width}}'.format(str(self._basis.month), width=2) \
                + '{:0>{width}}'.format(str(self._basis.day), width=2) \
                + '{:0>{width}}'.format(str(self._basis.hour), width=2) \
                + '{:0>{width}}'.format(str(self._basis.minute), width=2) \
                + '{:0>{width}}'.format(str(self._basis.second), width=2)
        else:
            raise NotImplementedError("fmt=" + fmt + " option for " +
                                      self.__class__.__name__ + ".getbasis().")
        return out

    def set(self,
            date_time=None,
            basis=None,
            term=None,
            cumulativeduration=None,
            statistical_process_on_duration=None,
            statistical_time_increment=None):
        """
        Sets validity and basis according to arguments.
        A consistency check is done if the three arguments are provided
        (which is useless anyway).

        Args: \n
        :param date_time: has to be of type :class:`datetime.datetime`
        :param basis: has to be of type :class:`datetime.datetime`
        :param term: has to be of type :class:`datetime.timedelta`
        :param cumulativeduration: has to be of type :class:`datetime.timedelta`
        :param statistical_process_on_duration: kind of statistical process
            that runs over the cumulative duration.
            Cf. GRIB2 typeOfStatisticalProcessing
        :param statistical_time_increment: time step over used for statistical
            process.
        """
        if isinstance(date_time, datetime.datetime):
            self._date_time = date_time
        elif date_time is not None:
            raise epygramError("argument 'date_time' must be of type" +
                               " datetime.datetime")
        if isinstance(basis, datetime.datetime):
            self._basis = basis
        elif basis is not None:
            raise epygramError("argument 'basis' must be of type" +
                               " datetime.datetime")
        if term is not None and not isinstance(term, datetime.timedelta):
            raise epygramError("argument 'term' must be of type" +
                               " datetime.timedelta")
        if cumulativeduration is not None and\
           not isinstance(cumulativeduration, datetime.timedelta):
            raise epygramError("argument 'cumulativeduration' must be of" +
                               " type datetime.timedelta")
        if statistical_time_increment is not None and\
           not isinstance(statistical_time_increment, datetime.timedelta):
            raise epygramError("argument 'statistical_time_increment' must be of" +
                               " type datetime.timedelta")

        if isinstance(term, datetime.timedelta):
            if date_time is not None and basis is not None and term is not None \
               and date_time - basis != term:
                raise epygramError("inconsistency between 'term', 'basis'" +
                                   " and 'date_time' arguments.")

            if self._date_time is None:
                if self._basis is None:
                    raise epygramError("cannot set 'term' without 'basis'" +
                                       " nor 'date_time'.")
                else:
                    self._date_time = self._basis + term
            else:
                if self._basis is None:
                    self._basis = self._date_time - term
                else:
                    self._date_time = self._basis + term

        if cumulativeduration is not None:
            self._cumulativeduration = cumulativeduration
        if self._cumulativeduration is not None and \
           statistical_process_on_duration is not None:
            self._statistical_process_on_duration = statistical_process_on_duration
        if statistical_time_increment is not None:
            self._statistical_time_increment = statistical_time_increment

    def is_valid(self):
        """Check the validity is valid, i.e. not null."""
        return self.get() is not None


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
        list.__init__(self, [])

        if validity_instance is not None:
            if len(kwargs) != 0:
                raise epygramError("One can not give, at the same time, validity_instance and other argument.")
            failed = False
            if isinstance(validity_instance, FieldValidity):
                self.append(validity_instance)
                if length > 1:
                    for _ in range(length - 1):
                        self.append(validity_instance.deepcopy())
            elif isinstance(validity_instance, FieldValidityList):
                self.extend(validity_instance)
            elif isinstance(validity_instance, list):
                if all([isinstance(f, FieldValidity) for f in validity_instance]):
                    self.extend(validity_instance)
                else:
                    failed = True
            else:
                failed = True
            if failed:
                raise epygramError("FieldValidityList must be built from" +
                                   " FieldValidity, from FieldValidityList" +
                                   " instances or from a list of FieldValidity.")
        elif kwargs != {}:
            # Check that all lengths are equal
            length = None
            mykwargs = {}
            for k, v in kwargs.items():
                mykwargs[k] = [v] if not isinstance(v, list) else v
                if length is None or length == 1:
                    length = len(mykwargs[k])
                if len(mykwargs[k]) != length:
                    raise epygramError("All the arguments must have the same length.")

            for k, v in mykwargs.items():
                if len(v) == 1:
                    mykwargs[k] = mykwargs[k] * length

            # We set the different objects
            if length is None:
                length = 1
            self.extend([FieldValidity(**{key: value[i] for (key, value) in mykwargs.items()}) for i in range(length)])
        elif isinstance(length, int):
            for _ in range(length):
                self.append(FieldValidity())

    def __str__(self):
        strout = '<List of FieldValidity which date/time are:\n'
        for v in self:
            strout += str(v.get()) + '\n'
        strout += '>'
        return strout

    def __getitem__(self, key):
        result = super(FieldValidityList, self).__getitem__(key)
        if isinstance(key, slice):
            result = FieldValidityList(result)
        elif not isinstance(key, int):
            raise TypeError("*key* should be of type 'int' or 'slice'")
        return result

    def __getslice__(self, start, end):  # deprecated but not for some builtins such as list...
        result = super(FieldValidityList, self).__getitem__(slice(start, end))
        return FieldValidityList(result)

    def __eq__(self, other):
        if isinstance(other, self.__class__) and len(self) == len(other):
            return all([v == other[i] for i,v in enumerate(self)])
        else:
            return False

    def __hash__(self):
        # known issue __eq__/__hash__ must be defined both or none, else inheritance is broken
        return object.__hash__(self)

    def recursive_diff(self, other):
        """Recursively list what differs from **other**."""
        if self != other:
            if not isinstance(other, self.__class__) or len(self) != len(other):
                return (str(self), str(other))
            else:
                return [v.recursive_diff(other[i]) for i,v in enumerate(self)]

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

    def statistical_time_increment(self, one=True, **kwargs):
        """This method returns the statistical_time_increment of all the validities."""
        length = len(self)
        result = [self[i].statistical_time_increment(**kwargs) for i in range(length)]
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
        # Check that all lengths are equal
        length = None
        mykwargs = {}
        for k, v in kwargs.items():
            mykwargs[k] = [v] if not isinstance(v, list) else v
            if length is None or length == 1:
                length = len(mykwargs[k])
            if len(mykwargs[k]) != length:
                raise epygramError("All the arguments must have the same length.")

        if length == 1:
            length = len(self)

        for k, v in mykwargs.items():
            if len(v) == 1:
                mykwargs[k] = mykwargs[k] * length

        # We set the different objects
        for i in range(length):
            self[i].set(**{key: value[i] for (key, value) in mykwargs.items()})

    def is_valid(self):
        """Check the validity is valid, i.e. not null."""
        return all([self[i].is_valid() for i in range(len(self))])

    def what(self, out=sys.stdout, cumulativeduration=True):
        """
        Writes in file a summary of the validity.

        :param out: the output open file-like object.
        :param cumulativeduration: if False, not written.
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
