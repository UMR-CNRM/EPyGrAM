#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for LFI library.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

from ctypes import c_longlong, c_char_p, c_bool
import numpy as np

from . import ctypesFF, IN, OUT, INOUT, treatReturnCode, addReturnCode


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfinaf(*args):
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
    return [(c_longlong(args[0]), IN),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfipos(*args):
    """
    Rewind record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    """
    return [(c_longlong(args[0]), IN)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlficas(*args):
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
    return [(c_longlong(args[0]), IN),
            (c_char_p(" " * 16), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_bool(args[1]), IN)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfiouv(*args):
    """
    Open a LFI file.

    Args:\n
    1) CDFILE: path to file to open
    2) CDSTATE: state of file ('NEW', 'OLD', 'UNKNOWN', 'SCRATCH')

    Returns:\n
    1) KNUMER: logical unit number associated to file
    """
    return[(c_char_p(args[0]), IN),
           (c_char_p(args[1]), IN),
           (c_longlong(), OUT)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfifer(*args):
    """
    Close a LFI file.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDSTTC: close status ('KEEP', 'SCRATCH', 'DELETE')
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(7)), IN)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfinfo(*args):
    """
    Get length of a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record

    Returns:\n
    1) KLONG: length of record
    2) KPOSEX: position in file of the first word of next record
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(16)), IN),
           (c_longlong(), OUT),
           (c_longlong(), OUT)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfilec(*args):
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
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(16)), IN),
           (c_longlong(args[2]), IN),
           (c_bool(args[3]), IN),
           (np.ndarray((args[2],), dtype=np.int64), OUT)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfiecr(*args):
    """
    Write a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of field to write
    3) KSIZE: Size of KTAB
    4) KTAB: integer array to write
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(16)), IN),
           (c_longlong(args[2]), IN),
           (args[3], IN)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfiren(*args):
    """
    Rename a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOM1: name of record to rename
    3) CDNOM2: new name of record
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(16)), IN),
           (c_char_p(args[2].ljust(16)), IN)]


@treatReturnCode
@ctypesFF
@addReturnCode
def wlfisup(*args):
    """
    Delete a record.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record to delete
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(16)), IN)]


@ctypesFF
def wget_compheader(*args):
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
    return[(c_longlong(args[0]), IN),
           (args[1], IN),
           (c_longlong(args[2]), IN),
           (c_longlong(), OUT),
           (c_longlong(), OUT)]


@ctypesFF
def wdecompress_field(*args):
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
    return[(c_longlong(args[0]), IN),
           (args[1], IN),
           (c_longlong(args[2]), IN),
           (c_longlong(args[3]), IN),
           (np.ndarray((args[3],), dtype=np.int64), OUT)]


@ctypesFF
def wcompress_field(*args):
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
    return[(args[0], INOUT),
           (c_longlong(args[1]), IN),
           (c_longlong(args[2]), IN),
           (c_longlong(args[3]), IN),
           (c_longlong(), OUT)]
