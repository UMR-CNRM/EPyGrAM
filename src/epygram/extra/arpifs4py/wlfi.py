#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for LFI library.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy as np

from ctypesForFortran import addReturnCode, treatReturnCode, IN, OUT, INOUT
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


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfinaf(KNUMER):
    """
    Get info about number of records in file.

    Args:\n
    1) KNUMER: logical unit number associated to file

    Returns:\n
    1) KNALDO: Number of actual logical data records (holes excluded)
    2) KNTROU: Number of logical records which are holes
    3) KNARES: Number of logical records which can be written in the reserved part of index (holes included)
    4) KNAMAX: Maximum number of logical records which one can write on logical unit
    """
    return ([KNUMER],
            [(np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfipos(KNUMER):
    """
    Rewind record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    """
    return ([KNUMER],
            [(np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlficas(KNUMER, LDAVAN):
    """
    Run through records getting their names and lengths.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) LDAVAN: true if one must move forward the pointer

    Returns:\n
    1) CDNOMA: name of next record
    2) KLONG: length of next record
    3) KPOSEX: position in file of the first word of next record
    """
    return ([KNUMER, LDAVAN],
            [(np.int64, None, IN),
             (str, (16,), OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (bool, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfiouv(CDFILE, CDSTATE):
    """
    Open a LFI file.

    Args:\n
    1) CDFILE: path to file to open
    2) CDSTATE: state of file ('NEW', 'OLD', 'UNKNOWN', 'SCRATCH')

    Returns:\n
    1) KNUMER: logical unit number associated to file
    """
    return ([CDFILE, CDSTATE],
            [(str, (len(CDFILE.encode('utf8')),), IN),
             (str, (len(CDSTATE.encode('utf8')),), IN),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfifer(KNUMER, CDSTTC):
    """
    Close a LFI file.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDSTTC: close status ('KEEP', 'SCRATCH', 'DELETE')
    """
    return ([KNUMER, CDSTTC],
            [(np.int64, None, IN),
             (str, (7,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfinfo(KNUMER, CDNOMA):
    """
    Get length of a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record

    Returns:\n
    1) KLONG: length of record
    2) KPOSEX: position in file of the first word of next record
    """
    return ([KNUMER, CDNOMA],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfilec(KNUMER, CDNOMA, KLONG, LDABORT):
    """
    Read a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record
    3) KLONG: length of record
    4) LDABORT: must we raise an exception on error -21 ?

    Returns:\n
    1) KTAB: integer array read
    """
    return ([KNUMER, CDNOMA, KLONG, LDABORT],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, (KLONG,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfiecr(KNUMER, CDNOMA, KSIZE, KTAB):
    """
    Write a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of field to write
    3) KSIZE: Size of KTAB
    4) KTAB: integer array to write
    """
    return ([KNUMER, CDNOMA, KSIZE, KTAB],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (np.int64, None, IN),
             (np.int64, (KSIZE,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfiren(KNUMER, CDNOM1, CDNOM2):
    """
    Rename a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOM1: name of record to rename
    3) CDNOM2: new name of record
    """
    return ([KNUMER, CDNOM1, CDNOM2],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (str, (16,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wlfisup(KNUMER, CDNOMA):
    """
    Delete a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record to delete
    """
    return ([KNUMER, CDNOMA],
            [(np.int64, None, IN),
             (str, (16,), IN)],
            None)


@ctypesFF()
def wget_compheader(KSIZE, KDATA, KLONG):
    """
    Wrapper to GET_COMPHEADER.

    Args:\n
    1) KSIZE: Size of KDATA
    2) KDATA: (part of) integer array read from record
    3) KLONG: length of compressed data

    Returns:\n
    1) KLONU: length of uncompressed data
    2) KTYPECOMP: type of compression
    """
    return ([KSIZE, KDATA, KLONG],
            [(np.int64, None, IN),
             (np.int64, (KSIZE,), IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),],
            None)


@ctypesFF()
def wdecompress_field(KSIZE, KCOMP, KTYPECOMP, KLDECOMP):
    """
    Wrapper to DECOMPRESS_FIELD

    Args:\n
    1) KSIZE: size of KCOMP
    2) KCOMP: compressed integer array
    3) KTYPECOMP: type of compression
    4) KLDECOMP: length of decompressed data

    Returns:\n
    1) KDECOMP: decompressed data integer array
    """
    return ([KSIZE, KCOMP, KTYPECOMP, KLDECOMP],
            [(np.int64, None, IN),
             (np.int64, (KSIZE,), IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KLDECOMP,), OUT)],
            None)


@ctypesFF()
def wcompress_field(KTAB, KX, KY, KSIZEDECOMP):
    """
    Wrapper to COMPRESS_FIELD

    Args:\n
    1) KTAB: decompressed data integer array (IN)
    2,3) KX, KY: x and y dimensions
    4) KSIZEDECOMP: size of decompressed data

    Returns:\n
    1) KTAB: compressed data integer array (OUT)
    2) KSIZECOMP: size of compressed integer array
    """
    return ([KTAB, KX, KY, KSIZEDECOMP],
            [(np.int64, (KSIZEDECOMP), INOUT),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, OUT)],
            None)
