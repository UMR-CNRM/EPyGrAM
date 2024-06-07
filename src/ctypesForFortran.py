#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This module provides a function (ctypesForFortranFactory) which return a decorator
and another function (fortran2signature) which helps at building signature.
See these function documentations for more help on them.

The module also exposes a dlclose function to try closing lib
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import ctypes
import numpy
import subprocess
import os
import re

from _ctypes import dlclose

__all__ = []

__version__ = "1.0.0"

__license__ = 'CeCILL-C'

__authors__ = ['Sébastien Riette']


# Static values used to define input/output status of arguments
IN = 1
OUT = 2
INOUT = 3


def addReturnCode(func):
    def wrapper(*args, **kwargs):
        """
        This decorator adds an integer at the beginning of the "returned"
        signature of the Python function.
        """
        out = func(*args, **kwargs)
        out[1].insert(0, (numpy.int64, None, OUT))
        return out
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def treatReturnCode(func):
    def wrapper(*args):
        """
        This decorator raises a Python error if the integer returned by
        addReturnCode is different from 0.
        """
        result = func(*args)
        try:
            nout = len(result)
        except Exception:
            nout = 1
        if nout == 1:
            result = (result,)
        if result[0] != 0:
            raise RuntimeError("Error code " + str(result[0]) + " was raised.")
        result = result[1:]
        if len(result) == 1:
            result = result[0]
        elif len(result) == 0:
            result = None
        return result
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def get_dynamic_libs(obj):
    """Get dynamic libs from a shared object lib or executable."""
    libs = {}
    osname = str(os.uname()[0])
    if osname == 'Linux':
        _re = re.compile('((?P<libname>lib.*) => )?(?P<libpath>/.*/.*\.so(\.\d+)?)\s\\(0x.*\)')
        ldd_out = [line.decode().strip() for line in
                     subprocess.Popen(['ldd', obj],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).stdout.readlines()]
        for line in ldd_out:
            match = _re.match(line)
            if match:
                matchdict = match.groupdict()
                if matchdict.get('libname') is None:
                    matchdict['libname'] = matchdict.get('libpath').split('/')[-1]
                libs[matchdict['libname']] = matchdict.get('libpath')
    elif osname == 'Darwin':
        _re = re.compile('\s*(?P<libdir>/.*/)(?P<libname>.*(\.\d+)?\.dylib)\s+.*')
        otool_out = [line.decode().strip() for line in
                     subprocess.Popen(['otool', '-L', obj],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).stdout.readlines()]
        for line in otool_out:
            match = _re.match(line)
            if match:
                libs[match.group('libname')] = match.group('libdir') + match.group('libname')
    else:
        raise NotImplementedError("OS: " + osname)
    return libs


def ctypesForFortranFactory(solib):
    """
    solib is a shared library already opened (with ctypes) or a filename to use.

    This function returns a tuple (function, handle). The handle can be used with dlclose.
    The returned function will return a decorator used to call fortran routines or functions
    contained in <solib> using ctypes. The function take two arguments: prefix and suffix which will be
    added to the python function name to build the name of the function in the shared library.
    By default prefix is the empty string whereas suffix is '_'.

    The actual python function that must be decorated must return the signature of the fortran routine.
    The signature is a tuple. Fisrt element is the actual list of arguments which will be used to call
    the fortran routine.
    The arguments must be put in the same order as the one declared in the fortran routine.
    The second element of the signature tuple is also a list, each element of the list is
    a tuple (type, shape, in_out_status), where:
      - type must be one of str, bool, np.float64, np.int64
      - shape is - None in case of a scalar value (bool, np.float64, np.int64)
                 - a tuple which represent the shape of the array
                 - in case of str, first shape element is the string size,
                   if other elements are present, variable is a string array whose
                   shape is given by shape[1:]
      - in_out_status is one of IN, OUT or INOUT constant declared in the module
    The third, and last, element of the signature tuple is None or a tuple representing the
    output value for a function. It is a tuple with (type, shape) as described above but
    without the in_out_status element. None must be used for subroutine.

    For input arguments, type checking is done against this signature.
    All elements must be declared using numpy classes but scalars are expected to be true
    python scalars (for instance (np.float64, None, IN) is the signature for a python float variable)

    The decorated function will call the fortran routine (or function) and return:
      - a tuple of all OUT and INOUT arguments (in the declaration order,
                                                with function output in first position)
      - a single value if only one argument is OUT or INOUT

    Note on strings: scalars strings are converted from and to unicode; it is certainly not tatally wanted
                     strings arrays are declared (in signature) with str but must be
                     created with the 'S' dtype with python3 (str is OK with python2)

    Note on inout arrays: there is only one version of the array (if `a' is an array declared
                          in inout "foo(a)" will return `a' with new values (except for bool
                          arrays in some circumstances)

    Known limitations:
      - only some types have been tested, other raise an exception but this
        could normally be extended easily
      - logical arrays must be declared with KIND=1 in the fortran routine
        as bool numpy array elements are 1 byte (and not 1 bit) values
      - fortran assumed-shape arrays, assumed-rank arrays or asterisk length strings are not callable
      - there is no support for optional fortran argument. We can use python optional arguments
        but the default python value will be passed to the fortran routine (it then will appear as present)
      - unicode/ASCCI issue is more than likely... (for scalars and string arrays)
      - only scalars (integer, real and logical) can be returned by functions
      - because integer value of boolean variables vary among compilers, we need to determine the
        compiler used. For now only gfortran and ifort are supported.

    Usage:
        Fortran code:
            FUNCTION F_INT(KIN)
              IMPLICIT NONE
              INTEGER(KIND=8), INTENT(IN) :: KIN
              INTEGER(KIND=8) :: F_INT
              F_INT=KIN+1
            END FUNCTION

            FUNCTION F_REAL(PIN)
              IMPLICIT NONE
              REAL(KIND=8), INTENT(IN) :: PIN
              REAL(KIND=8) :: F_REAL
              F_REAL=PIN+1.
            END FUNCTION

            FUNCTION F_BOOL(LIN)
              IMPLICIT NONE
              LOGICAL(KIND=1), INTENT(IN) :: LIN
              LOGICAL(KIND=1) :: F_BOOL
              F_BOOL=.NOT. LIN
            END FUNCTION

            SUBROUTINE FOO(KIN, KOUT, KINOUT,       & !Integer scalars
                           KAIN, KAOUT, KAINOUT,    & !Integer arrays
                           CDIN, CDOUT, CDINOUT,    & !Strings
                           CDAIN, CDAOUT, CDAINOUT, & !String arrays
                           PIN, POUT, PINOUT,       & !Float scalars
                           PAIN, PAOUT, PAINOUT,    & !Float arrays
                           LIN, LOUT, LINOUT,       & !Logical scalars
                           LAIN, LAOUT, LAINOUT,    & !Logical arrays
                           KAIN2, KAOUT2)             !2D integer arrays

              INTEGER(KIND=8), INTENT(IN) :: KIN
              INTEGER(KIND=8), INTENT(OUT) :: KOUT
              INTEGER(KIND=8), INTENT(INOUT) :: KINOUT

              INTEGER(KIND=8), DIMENSION(KIN), INTENT(IN) :: KAIN
              INTEGER(KIND=8), DIMENSION(KIN), INTENT(OUT) :: KAOUT
              INTEGER(KIND=8), DIMENSION(KIN), INTENT(INOUT) :: KAINOUT

              CHARACTER(LEN=10), INTENT(IN) :: CDIN
              CHARACTER(LEN=20), INTENT(OUT) :: CDOUT
              CHARACTER(LEN=20), INTENT(INOUT) :: CDINOUT

              CHARACTER(LEN=10), DIMENSION(2, 3), INTENT(IN) :: CDAIN
              CHARACTER(LEN=10), DIMENSION(2, 3), INTENT(OUT) :: CDAOUT
              CHARACTER(LEN=10), DIMENSION(2, 3), INTENT(INOUT) :: CDAINOUT

              REAL(KIND=8), INTENT(IN) :: PIN
              REAL(KIND=8), INTENT(OUT) :: POUT
              REAL(KIND=8), INTENT(INOUT) :: PINOUT

              REAL(KIND=8), DIMENSION(KIN), INTENT(IN) :: PAIN
              REAL(KIND=8), DIMENSION(KIN), INTENT(OUT) :: PAOUT
              REAL(KIND=8), DIMENSION(KIN), INTENT(INOUT) :: PAINOUT

              LOGICAL(KIND=1), INTENT(IN) :: LIN
              LOGICAL(KIND=1), INTENT(OUT) :: LOUT
              LOGICAL(KIND=1), INTENT(INOUT) :: LINOUT

              LOGICAL(KIND=1), DIMENSION(40), INTENT(IN) :: LAIN
              LOGICAL(KIND=1), DIMENSION(40), INTENT(OUT) :: LAOUT
              LOGICAL(KIND=1), DIMENSION(40), INTENT(INOUT) :: LAINOUT

              INTEGER(KIND=8), DIMENSION(4, 5), INTENT(IN) :: KAIN2
              INTEGER(KIND=8), DIMENSION(4, 5), INTENT(OUT) :: KAOUT2

              KOUT=KIN+1
              KINOUT=KINOUT+1

              KAOUT(:)=KAIN(:)+1
              KAINOUT(:)=KAINOUT(:)+1

              CDOUT=CDIN // "Foo"
              CDINOUT=CDINOUT(1:5) // CDINOUT(1:15)
              CDAOUT(1,:) = CDAIN(2, 3:1:-1)
              CDAOUT(2,:) = CDAIN(1, 3:1:-1)
              CDAINOUT(1,:) = CDAINOUT(2,:)
              CDAINOUT(2,:) = CDAINOUT(1,:)

              POUT=PIN+1.
              PINOUT=PINOUT+1.

              PAOUT(:)=PAIN(:)+1.
              PAINOUT(:)=PAINOUT(:)+1.

              LOUT=.NOT.LIN
              LINOUT=.NOT.LINOUT

              LAOUT(1:10)=LAIN(1:10)
              LAOUT(11:20)=.NOT. LAIN(11:20)
              LAOUT(21:30)=LAIN(21:30)
              LAOUT(31:40)=.NOT. LAIN(31:40)
              LAINOUT(:)=.NOT. LAINOUT(:)

              KAOUT2(1,:)=KAIN2(1,:)
              KAOUT2(2,:)=-KAIN2(2,:)
              KAOUT2(3,:)=KAIN2(3,:)
              KAOUT2(4,:)=-KAIN2(4,:)
            END SUBROUTINE
          compiled to create foo.so shared library (ex with gfortran:
          "gfortran -c -fPIC foo.F90 && gfortran -shared -g -o foo.so foo.o",
          or with ifort:
          "ifort -c -fpic foo.F90 && ifort -shared -g -o foo.so foo.o"
          )

        Python code:
            import numpy
            import ctypesForFortran

            IN = ctypesForFortran.IN
            OUT = ctypesForFortran.OUT
            INOUT = ctypesForFortran.INOUT

            ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory("./foo.so")

            @ctypesFF() #With gfortran, if f_int was inside module 'toto', we would use @ctypesFF(prefix='__toto_MOD_', suffix='')
            def f_int(KIN):
                return ([KIN],
                        [(numpy.int64, None, IN)], #INTEGER(KIND=8), INTENT(IN) :: KIN
                        (numpy.int64, None))

            @ctypesFF()
            def f_real(PIN):
                return ([PIN],
                        [(numpy.float64, None, IN)], #REAL(KIND=8), INTENT(IN) :: PIN
                        (numpy.float64, None))

            @ctypesFF()
            def f_bool(LIN):
                return ([LIN],
                        [(bool, None, IN)], #LOGICAL(KIND=1), INTENT(IN) :: LIN
                        (bool, None))

            @ctypesFF()
            def foo(KIN, KINOUT,    # Integer scalars  # Only IN and INOUT variabes.
                    PIN, PINOUT,    # Float scalars    #
                    LIN, LINOUT,    # Logical scalars  # The order can be different than the
                    CDIN, CDINOUT,  # Strings          # one expected by the fortran routine.
                    CDAIN, CDAINOUT,# String arrays    # Here we put the scalars first, then
                    KAIN, KAINOUT,  # Integer arrays   # the arrays, whereas this is not the
                    PAIN, PAINOUT,  # Float arrays     # order declared in fortran code.
                    LAIN, LAINOUT,  # Logical arrays   #
                    KAIN2):         # 2D integer arrays#
                return ([KIN, KINOUT,     #
                         KAIN, KAINOUT,   # Only IN and INOUT variabes.
                         CDIN, CDINOUT,   #
                         CDAIN, CDAINOUT, #
                         PIN, PINOUT,     # Here, this *must* be the same order
                         PAIN, PAINOUT,   # as the one in the fortran declaration
                         LIN, LINOUT,     #
                         LAIN, LAINOUT,   #
                         KAIN2],          #
                        [(numpy.int64, None, IN), #INTEGER(KIND=8), INTENT(IN) :: KIN
                         (numpy.int64, None, OUT), #INTEGER(KIND=8), INTENT(OUT) :: KOUT
                         (numpy.int64, None, INOUT), #INTEGER(KIND=8), INTENT(INOUT) :: KINOUT

                         (numpy.int64, (KIN, ), IN), #INTEGER(KIND=8), DIMENSION(KIN), INTENT(IN) :: KAIN
                         (numpy.int64, (KIN, ), OUT), #INTEGER(KIND=8), DIMENSION(KIN), INTENT(OUT) :: KAOUT
                         (numpy.int64, (KIN, ), INOUT), #INTEGER(KIND=8), DIMENSION(KIN), INTENT(INOUT) :: KAINOUT

                         (str, (10, ), IN), #CHARACTER(LEN=10), INTENT(IN) :: CDIN
                         (str, (20, ), OUT), #CHARACTER(LEN=20), INTENT(OUT) :: CDOUT
                         (str, (20, ), INOUT), #CHARACTER(LEN=20), INTENT(INOUT) :: CDINOUT

                         (str, (10, 2, 3), IN), #CHARACTER(LEN=10), DIMENSION(2, 3),, INTENT(IN) :: CDAIN
                         (str, (10, 2, 3), OUT), #CHARACTER(LEN=10), DIMENSION(2, 3),, INTENT(OUT) :: CDAOUT
                         (str, (10, 2, 3), INOUT), #CHARACTER(LEN=10), DIMENSION(2, 3),, INTENT(INOUT) :: CDAINOUT

                         (numpy.float64, None, IN), #REAL(KIND=8), INTENT(IN) :: PIN
                         (numpy.float64, None, OUT), #REAL(KIND=8), INTENT(OUT) :: POUT
                         (numpy.float64, None, INOUT), #REAL(KIND=8), INTENT(INOUT) :: PINOUT

                         (numpy.float64, (KIN, ), IN), #REAL(KIND=8), DIMENSION(KIN), INTENT(IN) :: PAIN
                         (numpy.float64, (KIN, ), OUT), #REAL(KIND=8), DIMENSION(KIN), INTENT(OUT) :: PAOUT
                         (numpy.float64, (KIN, ), INOUT), #REAL(KIND=8), DIMENSION(KIN), INTENT(INOUT) :: PAINOUT

                         (bool, None, IN), #LOGICAL(KIND=1), INTENT(IN) :: LIN
                         (bool, None, OUT), #LOGICAL(KIND=1), INTENT(OUT) :: LOUT
                         (bool, None, INOUT), #LOGICAL(KIND=1), INTENT(INOUT) :: LINOUT

                         (bool, (40, ), IN), #LOGICAL(KIND=1), DIMENSION(40), INTENT(IN) :: LAIN
                         (bool, (40, ), OUT), #LOGICAL(KIND=1), DIMENSION(40), INTENT(OUT) :: LAOUT
                         (bool, (40, ), INOUT), #LOGICAL(KIND=1), DIMENSION(40), INTENT(INOUT) :: LAINOUT

                         (numpy.int64, (4, 5), IN), #INTEGER(KIND=8), DIMENSION(4, 5), INTENT(IN) :: KAIN2
                         (numpy.int64, (4, 5), OUT), #INTEGER(KIND=8), DIMENSION(4, 5), INTENT(OUT) :: KAOUT2
                        ],
                        None)

            assert f_int(5) == 6, "f_int"
            assert f_real(5.) == 6., "f_real"
            assert f_bool(True) == False and f_bool(False) == True, "f_bool"

            kin = 5
            kinout = 8
            kain = numpy.arange(kin, dtype=numpy.int64)
            kainout = numpy.arange(kin, dtype=numpy.int64) * 10
            cdin = "blabla"
            cdinout = "azertyuiop"
            cdain = numpy.ndarray((2, 3), dtype=('S', 10))
            cdainout = numpy.ndarray((2, 3), dtype=('S', 10))
            for j in range(cdain.shape[0]):
                for i in range(cdain.shape[1]):
                    cdain[j, i] = str(i) + "_" + str(j)
                    cdainout[j, i] = str(i*10) + "_" + str(j*10)
            pin = 12.
            pinout = 53.
            pain = numpy.arange(kin, dtype=numpy.float64)
            painout = numpy.arange(kin, dtype=numpy.float64) * 10.
            lin = True
            linout = False
            lain = numpy.array([True, False] * 20)
            lainout = numpy.array([True, False, False, False] * 10)
            kain2 = numpy.arange(20, dtype=numpy.int64).reshape((4, 5))

            #IN/OUT test
            args = [kin, kinout, pin, pinout, lin, linout, cdin, cdinout,
                    cdain, cdainout, kain, kainout, pain, painout]
            kwargs = dict(LAIN=lain, LAINOUT=lainout, KAIN2=kain2)

            result = foo(*args, **kwargs) #We can call the python function with keyword arguments
            (kout, kinout, kaout, kainout,
             cdout, cdinout, cdaout, cdainout,
             pout, pinout, paout, painout,
             lout, linout, laout, lainout,
             kaout2) = result

            assert kout == kin + 1, "k 1"
            assert kinout == 8 + 1, "k 2"
            assert numpy.all(kaout == kain + 1), "k 3"
            assert numpy.all(kainout == numpy.arange(kin, dtype=numpy.int64) * 10 + 1), "k 4"

            assert cdout == (cdin.ljust(10) + "Foo").ljust(20), "cd 1"
            assert cdinout == "azertyuiop".ljust(20)[0:5] + "azertyuiop".ljust(20)[0:15], "cd 2"
            assert numpy.all(cdaout[0, :] == numpy.core.defchararray.ljust(cdain, 10)[1, ::-1]) and \
                   numpy.all(cdaout[1, :] == numpy.core.defchararray.ljust(cdain, 10)[0, ::-1]), "cd 3"
            assert numpy.all(cdainout[0, :] == numpy.core.defchararray.ljust(cdainout, 10)[1, :]) and \
                   numpy.all(cdainout[1, :] == numpy.core.defchararray.ljust(cdainout, 10)[0, :]), "cd 4"

            assert pout == pin + 1., "p 1"
            assert pinout == 53. + 1., "p 2"
            assert numpy.all(paout == pain + 1.), "p 3"
            assert numpy.all(painout == numpy.arange(kin, dtype=numpy.float64) * 10. + 1.), "p 4"

            assert lout == (not lin), "l 1"
            assert linout == (not False), "l 2"
            assert numpy.all(laout[0:10]==lain[0:10]) and \
                   numpy.all(laout[10:20] == (numpy.logical_not(lain[10:20]))) and \
                   numpy.all(laout[20:30]==lain[20:30]) and \
                   numpy.all(laout[30:40] == (numpy.logical_not(lain[30:40]))), "l 3"
            assert numpy.all(lainout == numpy.logical_not(numpy.array([True, False, False, False] * 10))), "l 4"

            assert numpy.all(kaout2[0, :]==kain2[0, :]) and numpy.all(kaout2[1, :]==-kain2[1, :]) and \
                   numpy.all(kaout2[2, :]==kain2[2, :]) and numpy.all(kaout2[3, :]==-kain2[3, :]), "K 5"

            #Checks
            #Normal order args = [kin, kinout, pin, pinout, lin, linout, cdin, cdinout, cdain,
            #                     cdainout, kain, kainout, pain, painout, lain, lainout, kain2]

            for args in [[kin, kinout*1., pin, pinout, lin, linout, cdin, cdinout, cdain,
                          cdainout, kain, kainout, pain, painout, lain, lainout, kain2], # wrong type for kinout
                         [kin, kinout, pin, pinout, lin, linout, cdin.ljust(500), cdinout, cdain,
                          cdainout, kain, kainout, pain, painout, lain, lainout, kain2], # cdin string too long
                         [kin, kinout, pin, pinout, lin, linout, cdin, cdinout, cdain,
                          cdainout, kain, kainout, pain, painout, lain, lainout, kain2.reshape((5, 4))], # wrong shape for kain2
                         [kin, kinout, pin, pinout, lin, linout, cdin, cdinout, cdain,
                          cdainout, kain*1., kainout, pain, painout, lain, lainout, kain2], # wrong type for kain array
                         [kin, kinout, pin, pinout, lin, linout, cdin, cdinout, cdain,
                          cdainout, kain, kainout, pain, painout, lain.reshape((40, 1)), lainout, kain2], # wrong rank for lain
                         [kin, kinout, pin, pinout, lin, linout, cdin, cdinout, cdain,
                          cdainout, list(kain), kainout, pain, painout, lain, lainout, kain2], # type not implemented for kain
                        ]:
                try:
                    result = foo(*args)
                    raise IOError("fake IOError")
                except IOError:
                    raise RuntimeError("call has not raise error; not normal")
                except:
                    #It is normal for call to raise an error
                    pass
            ctypesForFortran.dlclose(handle)

          if OK must execute without any output
    """
    import six
    if isinstance(solib, six.string_types):
        filename = solib
        my_solib = ctypes.CDLL(solib, ctypes.RTLD_GLOBAL)
    else:
        my_solib = solib
        filename = my_solib._name
    compiler = set()
    libs = get_dynamic_libs(filename)
    for lib in libs.keys():
        if lib.startswith('libgfortran'):
            compiler.add('gfortran')
        if lib.startswith('libifport'):
            compiler.add('ifort')
    if len(compiler) == 0:
        raise IOError("Don't know which compiler was used to build the shared library")
    true, false = {'ifort': (-1, 0),
                   'gfortran': (1, 0)}[compiler.pop()]

    def ctypesFF(prefix="", suffix="_"):
        """
        This function returns the decorator to use.
        prefix (resp. suffix) is the string that we must put before (resp. after)
        the python function name to build the name of the function contained
        in the shared library.

        Please refer to ctypesForFortranFactory for a complete documentation.
        """

        def decorator(func):
            """
            This function must be used as a decorator.
            True python function is called to determine the signature of the
            underlying fortran function of same name contained in the shared library.
            Input values are checked against the signature.
            Arguments are preapred to be passed to the fortran routine and
            output arguments are processed to be retrun by the python function.

            Please refer to ctypesForFortranFactory for a complete documentation.
            """

            def wrapper(*args, **kwargs):
                sorted_args, signature, ret = func(*args, **kwargs)
                assert isinstance(signature, list), "signature must be a list"
                assert all([s is None or isinstance(s, tuple) for s in signature]), \
                    "all elements of the signature must be a tuple"
                if not all([s[0] in [str, bool, numpy.int64, numpy.float64, numpy.int32, numpy.float32] for s in signature]):
                    raise NotImplementedError("This type is not (yet?) implemented")
                assert all([s[1] is None or isinstance(s[1], tuple) for s in signature]), \
                    "second element of argument signature must be None or a tuple"
                assert all([len(s[1]) > 0 for s in signature if isinstance(s[1], tuple)]), \
                    "if second element of argument is a tuple, it must not be empty"
                assert all([all([((isinstance(item, int) or
                                   isinstance(item, numpy.int64)) and
                                  item >= 0) for item in s[1]])
                            for s in signature if isinstance(s[1], tuple)]), \
                    "if second element of argument is a tuple, it must contain " + \
                    "only positive or null integer values"
                assert all([s[2] in [IN, INOUT, OUT] for s in signature]), \
                    "third element of argument signature must be IN, INOUT or OUT"

                assert len(sorted_args) == len([s for s in signature if s[2] in [IN, INOUT]]), \
                    "Get " + str(len(sorted_args)) + " arguments, " + \
                    str(len([s for s in signature if s[2] in [IN, INOUT]])) + " expected."

                argtypes = []
                effectiveArgs = []
                additionalArgs = []
                additionalArgtypes = []
                resultArgs = []
                iarg = 0
                for s in signature:
                    if s[0] == str:
                        if not isinstance(s[1], tuple):
                            raise ValueError("Signature for string must provide a length")
                    if s[0] == str and len(s[1]) == 1:
                        argtypes.append(ctypes.c_char_p)
                        if s[2] in [IN, INOUT]:
                            argument = sorted_args[iarg].encode("utf-8")
                            iarg += 1
                            if len(argument) > s[1][0]:
                                raise ValueError("String is too long (#arg " + str(iarg) + ")")
                            argument = ctypes.create_string_buffer(argument.ljust(s[1][0]))
                        else:
                            argument = ctypes.create_string_buffer(s[1][0])
                        additionalArgs.append(ctypes.c_longlong(s[1][0]))  # | according to f2py output for gfortran and ifort
                        additionalArgtypes.append(ctypes.c_longlong)       # | not a pointer, passed by value!
                        effectiveArgs.append(argument)
                        if s[2] in [OUT, INOUT]:
                            resultArgs.append(argument)
                    else:
                        if s[1] is None:
                            # scalar value
                            if s[0] == bool:
                                if (true, false) == (1, 0):
                                    cl = ctypes.c_bool
                                else:
                                    cl = ctypes.c_int8
                            elif s[0] == numpy.int64:
                                cl = ctypes.c_longlong
                            elif s[0] == numpy.float64:
                                cl = ctypes.c_double
                            elif s[0] == numpy.int32:
                                cl = ctypes.c_long
                            elif s[0] == numpy.float32:
                                cl = ctypes.c_float
                            else:
                                raise NotImplementedError("This scalar type is not (yet?) implemented")
                            argtypes.append(ctypes.POINTER(cl))
                            if s[2] in [IN, INOUT]:
                                argument = cl(sorted_args[iarg])
                                iarg += 1
                            else:
                                argument = cl()
                            effectiveArgs.append(ctypes.byref(argument))
                            if s[2] in [OUT, INOUT]:
                                resultArgs.append(argument)
                        else:
                            # Arrays
                            if s[0] == str:
                                expected_dtype = numpy.dtype(('S', s[1][0]))
                                effective_dtype = expected_dtype
                                expected_shape = s[1][1:]
                            else:
                                if s[0] == bool and (true, false) != (1, 0):
                                    expected_dtype = s[0]
                                    effective_dtype = numpy.int8
                                else:
                                    expected_dtype = s[0]
                                    effective_dtype = expected_dtype
                                expected_shape = s[1]
                            if s[2] in [IN, INOUT]:
                                argument = sorted_args[iarg]
                                iarg += 1
                                if not isinstance(argument, numpy.ndarray):
                                    raise ValueError("Arrays must be numpy.ndarrays")
                                if argument.dtype != expected_dtype:
                                    raise ValueError("Wrong dtype for #arg " + str(iarg - 1))
                                if len(expected_shape) != len(argument.shape):
                                    raise ValueError("Wrong rank for input array (#arg " + str(iarg - 1) + ")")
                                if expected_shape != argument.shape:
                                    raise ValueError("Wrong shape for input array (#arg " +
                                                     str(iarg - 1) + "), get " + str(argument.shape) +
                                                     ", expected " + str(expected_shape))
                                if s[0] == str:
                                    argument = numpy.core.defchararray.ljust(argument, s[1][0])
                                if s[0] == bool and (true, false) != (1, 0):
                                    arr = numpy.empty_like(argument, dtype=numpy.int8, order='F')
                                    arr[argument is True] = true
                                    arr[argument is False] = false
                                    argument = arr
                                if not argument.flags['F_CONTIGUOUS']:
                                    argument = numpy.asfortranarray(argument)
                            else:
                                argument = numpy.ndarray(expected_shape, dtype=effective_dtype, order='F')
                                if s[0] == str:
                                    argument = numpy.asfortranarray(numpy.core.defchararray.ljust(argument, s[1][0]))
                            argtypes.append(numpy.ctypeslib.ndpointer(dtype=effective_dtype,
                                                                      ndim=len(argument.shape),
                                                                      flags='F_CONTIGUOUS'))
                            effectiveArgs.append(argument)
                            if s[2] in [OUT, INOUT]:
                                resultArgs.append(argument)
                f = my_solib.__getitem__(prefix + func.__name__ + suffix)
                f.argtypes = argtypes + additionalArgtypes
                if ret is not None:
                    assert len(ret) == 2, "returned value must be described by a two-values tuple"
                    if ret[1] is None or (ret[0] == str and len(ret[1]) == 1):
                        if ret[0] == bool:
                            if (true, false) == (1, 0):
                                cl = ctypes.c_bool
                            else:
                                cl = ctypes.c_int8
                        elif ret[0] == numpy.int64:
                            cl = ctypes.c_longlong
                        elif ret[0] == numpy.float64:
                            cl = ctypes.c_double
                        elif ret[0] == numpy.int32:
                            cl = ctypes.c_long
                        elif ret[0] == numpy.float32:
                            cl = ctypes.c_float
                        elif ret[0] == str:
                            cl = ctypes.c_char_p
                            raise NotImplementedError("Functions with string as output value are not working")
                        else:
                            raise NotImplementedError("This scalar type is not (yet?) implemented")
                        argument = cl
                    else:
                        raise NotImplementedError("Functions with arrays as output are not working")
                        """
                        assert isinstance(ret[1], tuple), \
                            "if second element of returned value signature is not None, it must be a tuple"
                        if ret[0] == str:
                            dtype = numpy.dtype(('S', ret[1][0]))
                            ctype = ctypes.c_char_p
                            shape = ret[1][1:]
                        elif ret[0] == bool and (true, false) != (1, 0):
                            dtype = numpy.int8
                            ctype = ctypes.c_int8
                            shape = ret[1]
                        else:
                            dtype = ret[0]
                            shape = ret[1]
                            if ret[0] == numpy.int64:
                                ctype = ctypes.c_longlong
                            elif ret[0] == numpy.float64:
                                ctype = ctypes.c_double
                        result = numpy.ndarray(shape=shape, dtype=dtype)
                        argument = result.ctypes.data_as(ctypes.POINTER(ctype))
                        argument = ctypes.POINTER(ctype)
                        argument = numpy.ctypeslib.ndpointer(dtype=dtype, shape=shape)
                        """

                    f.restype = argument

                val = f(*(effectiveArgs + additionalArgs))

                if ret is not None:
                    if ret[0] == bool and (true, false) != (1, 0):
                        result = [val == true]
                    else:
                        result = [val]
                else:
                    result = []
                iarg = 0
                for s in signature:
                    if s[2] in [OUT, INOUT]:
                        argument = resultArgs[iarg]
                        iarg += 1
                        if s[0] == str and len(s[1]) == 1:
                            argument = argument.value.decode('utf-8')
                        else:
                            if s[1] is None:
                                # scalar
                                if s[0] == bool and (true, false) != (1, 0):
                                    argument = argument.value == true
                                else:
                                    argument = argument.value
                            else:
                                # array
                                if s[0] == bool and (true, false) != (1, 0):
                                    argument = argument == true
                                pass  # If needed, we could reverse F_CONTIGOUS here (we then would need to track those changes)
                        result.append(argument)
                if len(result) > 1:
                    return tuple(result)
                elif len(result) == 1:
                    return result[0]
            wrapper.__name__ = func.__name__
            wrapper.__doc__ = func.__doc__
            return wrapper
        return decorator
    return ctypesFF, my_solib._handle


def fortran2signature(filename=None, fortran_code=None, as_string=True,
                      kind_real=None, kind_logical=None, kind_integer=None,
                      prefix="", suffix="", solib=None, only=None, **kwargs):
    """
    This functions returns the signature as a string (if as_string) or as a python object.
    In this later case, in variables must be put in **kwargs.
    The default kind (specified by compilation options) can be provided for real, logical and integer.
    prefix and suffix are used to build the symbol name to use in the shared lib (use %s in
    these strings as a placeholder for the module name)

    By default, signature for all symbols are built, you can limit to one subroutine or function
    by specifying its name in the only argument.
    When using this function with as_string, this argument becomes mandatory.

    signature is the signature as expected by ctypesForFortran.

    one of filename or fortran_code is mandatory

    Function was tested only against fortran code with relatively simple formatting,
    there are certainly issues with real fortran code.

    Example: the command line
    ctypesForFortran.py --solib=./foo.so --suffix="_" --kind_real 8 --kind_integer 8 --kind_logical 1 foo.F90
    outputs a python script that may be used in place of the signature part of the python code
    given as example in the ctypesForFortranFactory function (except that, in the example, order of argument for
    foo subroutine have been changed).

    Alternatively, the signature part of the same example (given in the ctypesForFortranFactory function) can
    be replaced by something like:
    with open('foo.F90', 'r') as f:
        fortran = f.read()

    @ctypesFF()
    def f_int(KIN):
        return ctypesForFortran.fortran2signature(fortran_code=fortran, as_string=False,
                                                  prefix="", suffix="_", only='f_int', kin=KIN)

    @ctypesFF()
    def f_real(PIN):
        return ctypesForFortran.fortran2signature(fortran_code=fortran, as_string=False,
                                                  prefix="", suffix="_", only='f_real', PIN=PIN)

    @ctypesFF()
    def f_bool(LIN):
        return ctypesForFortran.fortran2signature(fortran_code=fortran, as_string=False,
                                                  prefix="", suffix="_", only='f_bool', LIN=LIN)

    @ctypesFF()
    def foo(KIN, KINOUT,    # Integer scalars  # Only IN and INOUT variabes.
            PIN, PINOUT,    # Float scalars    #
            LIN, LINOUT,    # Logical scalars  # The order can be different than the
            CDIN, CDINOUT,  # Strings          # one expected by the fortran routine.
            CDAIN, CDAINOUT,# String arrays    # Here we put the scalars first, then
            KAIN, KAINOUT,  # Integer arrays   # the arrays, whereas this is not the
            PAIN, PAINOUT,  # Float arrays     # order declared in fortran code.
            LAIN, LAINOUT,  # Logical arrays   #
            KAIN2):         # 2D integer arrays#
        return ctypesForFortran.fortran2signature(fortran_code=fortran, as_string=False,
                                                  prefix="", suffix="_", only='foo',
                                                  kind_logical=1,
                                                  KIN=KIN, KINOUT=KINOUT,
                                                  PIN=PIN, PINOUT=PINOUT,
                                                  LIN=LIN, LINOUT=LINOUT,
                                                  CDIN=CDIN, CDINOUT=CDINOUT,
                                                  CDAIN=CDAIN, CDAINOUT=CDAINOUT,
                                                  KAIN=KAIN, KAINOUT=KAINOUT,
                                                  PAIN=PAIN, PAINOUT=PAINOUT,
                                                  LAIN=LAIN, LAINOUT=LAINOUT,
                                                  KAIN2=KAIN2)
    """
    import re

    assert filename is not None or fortran_code is not None, "one of filename or fortran_code must be provided"
    assert not (filename is not None and fortran_code is not None), "one of filename or fortran_code must be None"
    assert len(kwargs) == len(set([k.lower() for k in kwargs])), "fortran variables are case-insensitive"
    assert as_string is False or solib is not None, "solib is required if as_string"
    assert as_string or only is not None, "only must be specified when as_string is False"

    if filename is not None:
        with open(filename, 'r') as f:
            fortran_code = f.read()

    # lines will contain the source code split by instructions
    lines_tmp = fortran_code.splitlines()
    lines = []
    line = ''
    ind = 0
    in_str = False
    while len(lines_tmp) > 0:
        if line == '':
            line = lines_tmp.pop(0).strip()
            ind = 0
        # Look for first (if any) interesting character among ', ", #, &
        se = re.search('\'|"|&|!|;', line[ind:])
        if se is None:
            if line.strip() != "":
                lines.append(line.strip())
            line = ''
        else:
            c = se.group()
            if c == '!' and not in_str:
                lines.append(line.strip())
                line = ''
            elif c == '&':
                if in_str and line[ind + line[ind:].find(c) - 1] == '\\':
                    # do not count
                    ind = ind + line[ind:].find(c) + 1
                else:
                    # not in str or true continuation character
                    after = line[(ind + line[ind:].find(c) + 1):].strip()
                    if not (after == "" or after[0] == "!"):
                        raise RuntimeError("& followed by something")
                    line = line[:ind + line[ind:].find(c)].strip()
                    nextline = lines_tmp.pop(0).strip()
                    if not in_str:
                        line += " "
                    if nextline.startswith('&'):
                        line += nextline[nextline.find('&') + 1:]
                    else:
                        line += nextline
            elif c == ';' and not in_str:
                lines.append(line[:ind + line[ind:].find(c)].strip())
                line = line[ind + line[ind:].find(c) + 1:]
                ind = 0
            elif c in ['"', "'"]:
                if not in_str:
                    # Beginning of string
                    ind = ind + line[ind:].find(c) + 1
                    in_str = c
                else:
                    if c != in_str:
                        # quote or double quote in string
                        ind = ind + line[ind:].find(c) + 1
                    else:
                        if line[ind + line[ind:].find(c) - 1] == '\\':
                            # not the end of string
                            ind = ind + line[ind:].find(c) + 1
                        else:
                            # This is the end of string
                            ind = ind + line[ind:].find(c) + 1
                            in_str = False

    objs = []  # each item is a dict with keys: module, name, in_var_names (list of in var names), signatures (list of signatures), result (result signature)
    low_kwargs = {k.lower():v for k, v in kwargs.items()}
    in_module = False
    in_function = False
    in_subroutine = False
    path = []
    select = {('real', 2): numpy.float16,
              ('real', 4): numpy.float32,
              ('real', 8): numpy.float64,
              ('integer', 1): numpy.int8,
              ('integer', 2): numpy.int16,
              ('integer', 4): numpy.int32,
              ('integer', 8): numpy.int64,
              ('logical', 1): bool}
    for line in lines:
        if line.lower().startswith('module') and line[6] in [" ", "\t"]:
            if in_module:
                raise RuntimeError("Already in module")
            in_module = line[6:].lower().strip()
        elif line.lower().startswith('subroutine') and line[10] in [" ", "\t"]:
            name = line[10:]
            name = name[:name.find('(')].lower().strip()
            path.append(name)
            if len(path) == 1 and only in [None, name]:
                in_subroutine = name
                args = [arg.strip().lower() for arg in line[line.find("(") + 1:line.find(")")].split(",")]
                current_obj = {'module': in_module if in_module else "",
                               'name': in_subroutine,
                               'var_names': args,
                               'signature': {},
                               'result': None,
                               'must_be_in': set(),
                               'intents': {},
                               'result_name': None}
        elif line.lower().startswith('function') and line[8] in [" ", "\t"]:
            name = line[8:]
            name = name[:name.find('(')].lower().strip()
            path.append(name)
            if len(path) == 1 and only in [None, name]:
                in_function = name
                args = [arg.strip().lower() for arg in line[line.find("(") + 1:line.find(")")].split(",")]
                current_obj = {'module': in_module if in_module else "",
                               'name': in_function,
                               'var_names': args,
                               'signature': {},
                               'result': None,
                               'must_be_in': set(),
                               'intents': {},
                               'result_name': None}
                result_name = line[line.find(")") + 1:].strip()
                if len(result_name) > 0:
                    if not result_name.startswith("result"):
                        raise RuntimeError("Something after the function definition which is not the result?")
                    result_name = result_name[result_name.find("(") + 1:result_name.find(")")]
                    current_obj['result_name'] = result_name
                else:
                    current_obj['result_name'] = in_function
        elif '::' in line and (in_subroutine or in_function) and len(path) == 1:
            end = len(line) if '!' not in line else line.find('!')
            args = [arg.strip().lower() for arg in line[line.find('::') + 2:end].split(',')]
            options_tmp = [opt.strip().lower() for opt in line[:line.find('::')].split(',')]
            options = []
            while len(options_tmp) > 0:
                opt = options_tmp.pop(0)
                if '(' in opt:
                    while ')' not in opt:
                        opt += ", " + options_tmp.pop(0)
                options.append(opt)
            dtype = None
            shape = None
            intent = None
            for arg in args:
                if arg in current_obj['var_names'] + ([current_obj['result_name']] if in_function else []):
                    if arg in current_obj['signature']:
                        raise RuntimeError("arg already declared: " + arg)
                    if dtype is None:
                        decode_kind = False
                        for opt in options:
                            decode_kind = False
                            if opt.startswith("dimension") and opt[9] in [' ', '\t', '(']:
                                if shape is None:
                                    shape = "(" if as_string else []
                                else:
                                    if as_string:
                                        shape = shape[:-1]
                                    else:
                                        shape = list(shape)
                                dimensions = [item.strip() for item in opt[opt.find('(') + 1:opt.find(')')].split(',')]
                                for d in dimensions:
                                    if d in current_obj['var_names']:
                                        current_obj['must_be_in'].add(d)
                                        if not as_string:
                                            if d not in low_kwargs:
                                                raise ValueError(d + "must be in kwargs")
                                            d = low_kwargs[d]
                                    else:
                                        d = int(d)
                                    if as_string:
                                        shape += str(d) + ', '
                                    else:
                                        shape.append(d)
                                shape = (shape + ")") if as_string else tuple(shape)
                            elif opt.startswith("intent") and opt[6] in [' ', '\t', '(']:
                                intent = opt[opt.find('(') + 1:opt.find(')')].upper()
                                if intent not in ['IN', 'OUT', 'INOUT']:
                                    raise RuntimeError("intent must be IN, OUT or INOUT")
                                if not as_string:
                                    intent = {'IN':IN, 'OUT':OUT, 'INOUT':INOUT}[intent]
                                current_obj['intents'][arg] = intent
                            elif opt == 'optional':
                                raise RuntimeError("optional argument are not allowed")
                            elif opt.startswith("real") and (len(opt) == 4 or opt[4] in [' ', '\t', '(']):
                                dtype = "real"
                                decode_kind = True
                            elif opt.startswith("integer") and (len(opt) == 7 or opt[7] in [' ', '\t', '(']):
                                dtype = "integer"
                                decode_kind = True
                            elif opt.startswith("logical") and (len(opt) == 7 or opt[7] in [' ', '\t', '(']):
                                dtype = "logical"
                                decode_kind = True
                            elif opt.startswith("character") and (len(opt) == 9 or opt[9] in [' ', '\t', '(']):
                                if as_string:
                                    dtype = "str"
                                else:
                                    dtype = str
                                if '(' not in opt:
                                    raise RuntimeError("character type must provide a length")
                                length = opt[opt.find('(') + 1:opt.find(')')].strip()
                                if not length.startswith('len'):
                                    raise RuntimeError("character length must start with len")
                                length = length[3:].strip()
                                if length[0] != '=':
                                    raise RuntimeError("character legth must start with len=")
                                length = length[1:]
                                if length in current_obj['var_names']:
                                    current_obj['must_be_in'].add(length)
                                    if not as_string:
                                        if length not in low_kwargs:
                                            raise ValueError(length + "must be in kwargs")
                                        length = low_kwargs[length]
                                else:
                                    length = int(length)
                                if shape is None:
                                    if as_string:
                                        shape = "(" + str(length) + ", )"
                                    else:
                                        shape = (length, )
                                else:
                                    if as_string:
                                        shape = "(" + str(length) + ", " + shape[1:]
                                    else:
                                        shape = tuple([length] + list(shape))
                            if decode_kind:
                                kind = None
                                if '(' in opt:
                                    kind = opt[opt.find('(') + 1:opt.find(')')].strip()
                                    if not kind.startswith('kind'):
                                        raise RuntimeError("kind specification must start with kind")
                                    kind = kind[4:].strip()
                                    if kind[0] != "=":
                                        raise RuntimeError("kind specification must start with kind=")
                                    kind = kind[1:].strip()
                                    if kind in current_obj['var_names']:
                                        current_obj['must_be_in'].add(kind)
                                        if not as_string:
                                            if kind not in low_kwargs:
                                                raise ValueError(kind + "must be in kwargs")
                                            kind = low_kwargs[kind]
                                    else:
                                        # kind must be an int
                                        kind = int(kind)
                                else:
                                    if dtype == 'real' and kind_real is not None:
                                        kind = kind_real
                                    elif dtype == 'integer' and kind_integer is not None:
                                        kind = kind_integer
                                    elif dtype == 'logical' and kind_logical is not None:
                                        kind = kind_logical
                                if as_string:
                                    dtype = "select[('" + dtype + "', " + str(kind) + ")]"
                                else:
                                    if (dtype, kind) not in select:
                                        raise NotImplementedError("This kind is not implemented: " + str((dtype, kind)))
                                    dtype = select[(dtype, kind)]
                        if dtype is None:
                            raise RuntimeError("declaration of arg " + arg + " does not provide type information")
                        if intent is None and arg != current_obj['result_name']:
                            raise RuntimeError("declaration of arg " + arg + " does not provide intent IN/OUT information")
                    if as_string:
                        if shape is None:
                            shape = "None"
                        if arg == current_obj['result_name']:
                            current_obj['signature'][arg] = "(" + dtype + "," + shape + ")"
                        else:
                            current_obj['signature'][arg] = "(" + dtype + "," + shape + "," + intent + ")"
                    else:
                        if arg == current_obj['result_name']:
                            current_obj['signature'][arg] = (dtype, shape)
                        else:
                            current_obj['signature'][arg] = (dtype, shape, intent)
        elif len(line.lower().split()) > 1 and line.lower().split()[0:2] in [['end', 'subroutine'], ['end', 'function']]:
            if len(path) == 1 and (in_subroutine or in_function):
                current_obj['in_var_names'] = []
                current_obj['signatures'] = []
                for arg in current_obj['var_names']:
                    if arg not in current_obj['signature']:
                        raise RuntimeError("arg not found in declaration: " + arg)
                    if current_obj['intents'][arg] in ['IN', 'INOUT', IN, INOUT]:
                        current_obj['in_var_names'].append(arg)
                    current_obj['signatures'].append(current_obj['signature'][arg])
                if in_function:
                    current_obj['result'] = current_obj['signature'][current_obj['result_name']]
                    del current_obj['signature'][current_obj['result_name']]
                if current_obj['result'] is None and as_string:
                    current_obj['result'] = "None"
                for arg in current_obj['must_be_in']:
                    if arg not in current_obj['in_var_names']:
                        raise RuntimeError("An argument (" + arg + ") with INTENT(OUT) has been used for kind or dimension")
                del current_obj['signature'], current_obj['var_names'], current_obj['must_be_in'], current_obj['intents']
                objs.append(current_obj)
                in_subroutine = False
                in_function = False
            path = path[:-1]
        elif line == 'contains':
            pass

    if as_string:
        result = "import numpy\nimport ctypesForFortran\n\n"
        result += "IN = ctypesForFortran.IN\nOUT = ctypesForFortran.OUT\nINOUT = ctypesForFortran.INOUT\n\n"
        result += "select = {('real', 2): numpy.float16,\n          ('real', 4): numpy.float32,\n"
        result += "          ('real', 8): numpy.float64,\n          ('integer', 1): numpy.int8,\n"
        result += "          ('integer', 2): numpy.int16,\n          ('integer', 4): numpy.int32,\n"
        result += "          ('integer', 8): numpy.int64,\n          ('logical', 1): bool}\n\n"
        result += "pre_suf = {}\n"
        for module_name in set([obj['module'] for obj in objs]):
            result += "pre_suf['" + module_name + "'] = ('"
            result += ((prefix % module_name) if '%s' in prefix else prefix) + "', '"
            result += ((suffix % module_name) if '%s' in suffix else suffix) + "')\n"
        result += "\n\n"
        result += "ctypesFF = ctypesForFortran.ctypesForFortranFactory('" + solib + "')\n\n"
        for obj in objs:
            result += "@ctypesFF(*pre_suf['" + obj['module'] + "'])\n"
            result += "def " + obj['name'] + "(" + ', '.join(obj['in_var_names']) + "):\n"
            result += "    return ([" + ', '.join(obj['in_var_names']) + "],\n"
            result += "            [" + ',\n             '.join(obj['signatures']) + "],\n"
            result += "            " + obj['result'] + ")\n\n"
    else:
        if len(objs) != 1:
            raise ValueError("The searched symbol was not found")
        obj = objs[0]
        result = ([low_kwargs[arg] for arg in obj['in_var_names']],
                  obj['signatures'],
                  obj['result'])

    return result


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Simple fortran parser which produce signature")
    parser.add_argument('filename', metavar='filename', type=str, nargs=1,
                        help='file name of the file containing the fortran cde to parse')
    parser.add_argument('--kind_real', type=int, default=None, help='Kind to use for reals when not specified in declaration')
    parser.add_argument('--kind_integer', type=int, default=None, help='Kind to use for integers when not specified in declaration')
    parser.add_argument('--kind_logical', type=int, default=None, help='Kind to use for logicals when not specified in declaration')
    parser.add_argument('--prefix', type=str, default="", help='prefix to add to the python function name to build the symbol name as found in the shared lib')
    parser.add_argument('--suffix', type=str, default="", help='suffix to add to the python function name to build the symbol name as found in the shared lib')
    parser.add_argument('--solib', type=str, help='path the shared lib', required=True)
    options = parser.parse_args()
    print(fortran2signature(filename=options.filename[0], kind_real=options.kind_real,
                            kind_logical=options.kind_logical, kind_integer=options.kind_integer,
                            prefix=options.prefix, suffix=options.suffix, solib=options.solib))
