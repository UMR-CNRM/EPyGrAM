#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for LFA library.
"""

from ctypes import c_longlong, c_char_p, c_bool
import numpy as np

from . import ctypesFF, IN, OUT, treatReturnCode_LFA, addReturnCode



@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfaouv(*args):
    """
    Open a LFA file.
    
    Input:
        cdnomf      file name.
        cdtypo      opening type: 'R' READ, 'W' WRITE, 'A' APPEND, 'S' SCRATCH.
    Output:
        kul         logical unit of LFA file.
    """
    return[(c_char_p(args[0].encode("utf-8")), IN),
           (c_char_p(args[1]), IN),
           (c_longlong(), OUT)]

@ctypesFF
def wlfafer(*args):
    """
    Close a LFA file.
    
    Input:
        kul        logical unit of LFA file.
    """
    return[(c_longlong(args[0]), IN)]

@ctypesFF
def wlfaecrr(*args):
    """
    Write real data on LFA file.
    
    Input:
        kul              logical unit of LFA file.
        cdna             name of article to write.
        preel(1,klong)   real data to write.
        klong            length of article to write.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (args[2], IN),
           (c_longlong(args[3]), IN)]

@ctypesFF
def wlfaecri(*args):
    """
    Write integer data of LFA file.

    Input:
        kul                  logical unit of LFA file.
        cdna                 name of article to write.
        kentier(1,klong)     integers to write.
        klong                length of article to write.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (args[2], IN),
           (c_longlong(args[3]), IN)]

@ctypesFF
def wlfaecrc(*args):
    """
    Write character data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             name of article to write.
        cdcar(1,klong)   characters to write.
        klong            length of article to write.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (args[2], IN),
           (c_longlong(args[3]), IN)]

@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfaleci(*args):
    """
    Read real data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array preel.
    Output:
        preel(1,klong)   real elements read.
        klong            number of real elements read.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (c_longlong(args[2]), IN),
           (np.ndarray((args[2],), dtype=np.int64), OUT),
           (c_longlong(), OUT)]

@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfalecr(*args):
    """
    Read integer data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array kentier.
    Output:
        kentier(1,klong) integer elements read.
        klong            number of integer elements read.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (c_longlong(args[2]), IN),
           (np.ndarray((args[2],), dtype=np.float64), OUT),
           (c_longlong(), OUT)]

@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfalecc(*args):
    """
    Read character data on LFA file.

    Input:
        kul              logical unit of LFA file.
        cdna             article name.
        kdimb            physical dimension of array clcar.
        klenc            max length of strings in clcar
    Output:
        kreturncode      error code
        clcar            array of elements read. 
        klong            number of character elements read.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (c_longlong(args[2]), IN),
           (c_longlong(args[3]), IN),
           (np.ndarray((args[2],), dtype=np.dtype((np.str, args[3]))), OUT),
           (c_longlong(), OUT)]

@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfacas(*args):
    """
    Get documentation about a LFA article.

    Input:
        kul               file logical unit.
        cdna:             name of required article.
    Output:
        cdtype            article type: 'R4', 'I8', 'C '.
        klong             number of elements in this article.
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].encode("utf-8")), IN),
           (c_char_p(" "*2), OUT),
           (c_longlong(), OUT)]

@ctypesFF
def wlfalaft(*args):
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
    return[(c_longlong(args[0]), IN),
           (c_longlong(args[1]), IN),
           (c_longlong(args[2]), IN),
           (c_longlong(), OUT),
           (np.ndarray((args[1],), dtype=np.dtype((np.str, args[2]))), OUT)]

@treatReturnCode_LFA
@ctypesFF
@addReturnCode
def wlfatest(*args):
    """
    Test if a file is a LFA one.

    Input:
        cdnomf      file name.
    Output:
        ldlfa=.true. if the file is a LFA one, .false. else case.
    """
    return[(c_char_p(args[0].encode("utf-8")), IN),
           (c_bool(), OUT)]
