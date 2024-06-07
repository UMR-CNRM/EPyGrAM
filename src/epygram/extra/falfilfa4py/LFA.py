#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for LFA library.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy as np

from ctypesForFortran import addReturnCode, IN, OUT
from . import ctypesFF
# Note to developers:
# Using the ctypesFF decorator, the Python function return a tuple containing:
#    tup[0]:
#        [the arguments of the Python function,
#         to be passed as in/inout arguments to the Fortran subroutine]
#    tup[1]:
#        [the python signature of all Fortran subroutine arguments]
#    tup[2]:
#        None in case of a Fortran subroutine, the output in case of a Fortran function)


def treatReturnCode_LFA(func):
    def wrapper(*args):
        """
        This decorator raises a Python error if the integer returned by
        addReturnCode is different from 0.
        The message raised is linked to the LFA library.
        """
        LFA_errors = {'-999':'Error while searching for an available logical\
                              unit.',
                      '-1':'Field not found in file.',
                      '1':'Field not found in file.',
                      '-6':'Article bigger than expected (argument kdimb).',
                      '-8':'Wrong data type (real, integer, char.).'}
        result = func(*args)
        try:
            nout = len(result)
        except Exception:
            nout = 1
        if nout == 1:
            result = (result,)
        if result[0] != 0:
            raise RuntimeError("Error code " + str(result[0]) +
                               " was raised : " +
                               LFA_errors[str(result[0])])
        result = result[1:]
        if len(result) == 1:
            result = result[0]
        elif len(result) == 0:
            result = None
        return result
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfaouv(cdnomf, cdtypo):
    """
    Open a LFA file.

    Input:
        cdnomf      file name.
        cdtypo      opening type: 'R' READ, 'W' WRITE, 'A' APPEND, 'S' SCRATCH.
    Output:
        kul         logical unit of LFA file.
    """
    return ([cdnomf, cdtypo],
            [(str, (len(cdnomf.encode('utf8')),), IN),
             (str, (len(cdtypo.encode('utf8')),), IN),
             (np.int64, None, OUT),],
            None)


@ctypesFF()
def wlfafer(kul):
    """
    Close a LFA file.

    Input:
        kul        logical unit of LFA file.
    """
    return ([kul],
            [(np.int64, None, IN)],
            None)


@ctypesFF()
def wlfaecrr(kul, cdna, preel, klong):
    """
    Write real data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             name of article to write.
        preel(klong)     real data to write.
        klong            length of article to write.
    """
    return ([kul, cdna, preel, klong],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (np.float64, (klong,), IN),
             (np.int64, None, IN)],
            None)


@ctypesFF()
def wlfaecri(kul, cdna, kentier, klong):
    """
    Write integer data of LFA file.

    Input:
        kul                  logical unit of LFA file.
        cdna                 name of article to write.
        kentier(klong)       integers to write.
        klong                length of article to write.
    """
    return ([kul, cdna, kentier, klong],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (np.float64, (klong,), IN),
             (np.int64, None, IN)],
            None)


@ctypesFF()
def wlfaecrc(kul, cdna, cdcar, klong):
    """
    Write character data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             name of article to write.
        cdcar(klong)     characters to write.
        klong            length of article to write.
    """
    return ([kul, cdna, cdcar, klong],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (str, (klong,), IN),  # FIXME: how to dimension ?
             (np.int64, None, IN)],
            None)


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfaleci(kul, cdna, kdimb):
    """
    Read integer data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array kentier.
    Output:
        kentier(klong)   integer elements read.
        klong            number of integer elements read.
    """
    return ([kul, cdna, kdimb],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (np.int64, None, IN),
             (np.int64, (kdimb,), OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfalecr(kul, cdna, kdimb):
    """
    Read real data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array preel.
    Output:
        preel(klong)     real elements read.
        klong            number of real elements read.
    """
    return ([kul, cdna, kdimb],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (np.int64, None, IN),
             (np.float64, (kdimb,), OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfalecc(kul, cdna, kdimb, klenc):
    """
    Read character data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array clcar.
        klenc            max length of strings in clcar
    Output:
        clcar            array of elements read.
        klong            number of character elements read.
    """
    return ([kul, cdna, kdimb, klenc],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (str, (klenc, kdimb), OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfacas(kul, cdna):
    """
    Get documentation about a LFA article.

    Input:
        kul               file logical unit.
        cdna:             name of required article.
    Output:
        cdtype            article type: 'R4', 'I8', 'C '.
        klong             number of elements in this article.
    """
    return ([kul, cdna],
            [(np.int64, None, IN),
             (str, (len(cdna.encode('utf8')),), IN),
             (str, (2,), OUT),
             (np.int64, None, OUT)],
            None)


@ctypesFF()
def wlfalaft(kul, kdlis, klenc):
    """
    Article list of a LFA file, on an array.

    Input:
        kul            logical unit of LFA file.
        kdlis          physical dimension of array cdlis.
        klenc          max length of strings in clcar
    Output:
        knlis          number of articles on the file. This number is also
                       the number of elements written on cclis.
        cclis          array of article names.
    """
    return ([kul, kdlis, klenc],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (str, (klenc, kdlis), OUT)],
            None)


@treatReturnCode_LFA
@ctypesFF()
@addReturnCode
def wlfatest(cdnomf):
    """
    Test if a file is a LFA one.

    Input:
        cdnomf      file name.
    Output:
        ldlfa=.true. if the file is a LFA one, .false. else case.
    """
    return ([cdnomf],
            [(str, (len(cdnomf.encode('utf8')),), IN),
             (bool, None, OUT)],
            None)
