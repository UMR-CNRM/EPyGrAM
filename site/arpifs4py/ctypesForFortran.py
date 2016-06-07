#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
ctypesForFortran:

Contains a decorator to wrap call to ctypes functions.
"""

from ctypes import c_char_p, create_string_buffer, c_longlong, POINTER, byref
import numpy as np

IN = 1
OUT = 2
INOUT = 3

def ctypesForFortranFactory(solib):
    def ctypesForFortran(func):
        """
        Decorator to wrap call to ctypes functions.

        The real python function must return the signature of the ctypes
        function as a list whose elements are tuples (data, InOut) where data
        is a ctypes instance (this instance must contain the effective data
        when in input) and InOut must be IN, OUT or INOUT.

        In case of a numpy array, the array is directly put in first place of
        the tuple.

        If an argument is INOUT:
        - if it's a scalar or a string: argument is untouched and new value is
          returned (copy)
        - if it's an array: argument is modified ans new value is
          returned (no copy).

        Result of the call is a value or a tuple of values.

        Known limitations:
        - access to argument by name is not allowed (the use of prototypes
          could help to introduce this feature but sizes of array and strings
          will be difficult to retrieve)
        - a type check against function signature is done except for numpy
          arrays, they are pass directly to the fortran function. Type check
          for numpy arrays could be done inside signature function.
        - only numpy arrays of int64, float64 and strings are implemented
          (extension is easy)
        - fortran assumed-shape arrays are not allowed.

        Usage:
          Fortran code:
            SUBROUTINE FOO(KIN, KOUT, CDIN, CDOUT, PIN, POUT, LIN)
              INTEGER(KIND=8), INTENT(INOUT) :: KIN
              INTEGER(KIND=8), INTENT(OUT) :: KOUT
              CHARACTER(LEN=*), INTENT(IN) :: CDIN
              CHARACTER(LEN=20), INTENT(OUT) :: CDOUT
              REAL(KIND=8), DIMENSION(KIN), INTENT(IN) :: PIN
              REAL(KIND=8), DIMENSION(KIN), INTENT(OUT) :: POUT
              LOGICAL, INTENT(IN) :: LIN
              KOUT=KIN+1
              CDOUT="Foo"
              POUT(:)=PIN(:)+1
              KIN=KIN-1
            END SUBROUTINE
            compiled to create foo.so shared library
            (ex with gfortran:
            "gfortran -c -fPIC foo.F90;
             gfortran -shared -g -o foo.so foo.o")

          Python code:
            from ctypes import *
            import numpy as np
            import ctypesForFortran

            IN=ctypesForFortran.IN
            OUT=ctypesForFortran.OUT
            INOUT=ctypesForFortran.INOUT

            so=CDLL("./foo.so")

            ctypesFF=ctypesForFortran.ctypesForFortranFactory(so)

            @ctypesFF
            def foo(*args): return [(c_longlong(args[0]), INOUT),
                                    (c_longlong(), OUT),
                                    (c_char_p(args[1]), IN),
                                    (c_char_p(" "*20), OUT),
                                    (args[2], IN),
                                    (np.ndarray((args[0],), dtype=np.float64), OUT),
                                    (c_bool(args[3]), IN)]

            kin=5
            cdin="blabla"
            pin=np.arange(kin, dtype=np.float64)
            (kin2, kout, cdout, pout)=foo(kin, cdin, pin, True)
            assert kin2==kin-1, "kin2 is not kin-1" #kin remained untouched
            assert kout==kin+1, "kout is not kin+1"
            assert cdout=="Foo".ljust(20), "cdout is not 'Foo'"
            assert (pout==pin+1).all(), "pout is not pin+1"
        """
        def wrapper(*args):
            signature = func(*args)
            argtypes1 = []
            argtypes2 = []
            effectiveArgs1 = []
            effectiveArgs2 = []
            resultArgs = []
            for arg in signature:
                if arg[0].__class__.__name__ == 'c_char_p':
                    argtypes1.append(c_char_p)
                    if arg[1] == IN:
                        effectiveArgs1.append(arg[0].value)
                    else:
                        c = create_string_buffer(len(arg[0].value))
                        effectiveArgs1.append(c)
                        resultArgs.append(c)
                    argtypes2.append(c_longlong)
                    effectiveArgs2.append(c_longlong(len(arg[0].value)))
                elif arg[0].__class__.__name__ == 'ndarray':
                    argtypes1.append(np.ctypeslib.ndpointer(dtype=arg[0].dtype,
                                                            ndim=len(arg[0].shape),
                                                            flags='F_CONTIGUOUS'))
                    effectiveArgs1.append(arg[0])
                    if arg[1] != IN:
                        resultArgs.append(arg[0])
                else:
                    argtypes1.append(POINTER(arg[0].__class__))
                    effectiveArgs1.append(byref(arg[0]))
                    if arg[1] != IN:
                        resultArgs.append(arg[0])
            f = solib.__getitem__(func.__name__ + '_')
            f.argtypes = argtypes1 + argtypes2
            #f.restype=None
            f(*(effectiveArgs1 + effectiveArgs2))
            result = tuple([(arg.value if arg.__class__.__name__ != 'ndarray'
                             else arg)
                            for arg in resultArgs])
            if len(resultArgs) > 1:
                return result
            elif len(resultArgs) == 1:
                return result[0]
        return wrapper
    return ctypesForFortran
