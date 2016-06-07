#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for FA library.
"""

from ctypes import c_longlong, c_char_p, c_bool, c_double
import numpy as np

from . import ctypesFF, IN, OUT, treatReturnCode, addReturnCode



@ctypesFF
def get_facst(*args):
    """
    Export maximum sizes used for fa format.
    
    Args:\n
    1) JPXPAH
    2) JPXIND
    3) JPXGEO
    4) JPXNIV
    """
    return [(c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfaitou(*args):
    """
    Open a FA file.
    
    Args:\n
    1) CDFILE: path to file to open
    2) CDSTATE: state of file ('NEW', 'OLD')
    3) CDNOMC: name of "cadre"
    
    Returns:\n
    1) KNUMER: logical unit number associated to file
    """
    return [(c_char_p(args[0].encode("utf-8")), IN),
            (c_char_p(args[1]), IN),
            (c_longlong(), OUT),
            (c_char_p(args[2].ljust(16)), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfaveur(*args):
    """
    Get compression parameters of a field in file.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    
    Returns:\n
    1) KNGRIB: Grib encoding level (-1,0,1,2,3);
    2) KNBPDG: Number of bits by grid point value
    3) KNBCSP: Number of bits by real/imaginary part of spectral coefficient
    4) KSTRON: Unpacked under truncature
    5) KPUILA: laplacian power
    6) KDMOPL: KPUILA level of modulation
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfalsif(*args):
    """
    Get identifier of file.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    
    Returns:\n
    1) CDIDEN: identifier of file
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(" "*80), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfanion(*args):
    """
    Get the characteristics of a field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDPREF: potential prefix of record name
    3) KNIVAU: potential vertical level
    4) CDSUFF: potential suffix of record name
    
    Returns:\n
    1) LDEXIS: true if record exists
    2) LDCOSP: true if field is spectral
    3) KNGRIB: level of GRIB encoding
    4) KNBITS: numbers of encoding bits
    5) KSTRON: potential under-trocature
    6) KPUILA: potential laplacian power
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(4)), IN),
            (c_longlong(args[2]), IN),
            (c_char_p(args[3].ljust(12)), IN),
            (c_bool(), OUT), (c_bool(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT)]

@ctypesFF
def wfacies(*args):
    """
    Get parameters in header.
    
    Args:\n
    1-4) KXPAH, KXIND, KXGEO, KXNIV: maximal dimensions
    5) CDNOMC: name of "cadre"
    
    Returns:\n
    1) KTYPTR: Type of horizontal transformation
    2) PSLAPO: Sinus of latitude of pole of interest
    3) PCLOPO: Cosinus of longitude of pole of interest
    4) PSLOPO: Sinus of longitude of pole of interest
    5) PCODIL: Dilatation coefficient
    6) KTRONC: Troncature
    7) KNLATI: Nomber of latitudes (from pole to pole)
    8) KNXLON: Maximum number of longitudes by parallel
    9) KNLOPA: Number of longitudes by parallel (from north pole toward equator only)
    10) KNOZPA: Maximum zonal wave number by parallel (from north pole toward equator only)
    11) PSINLA: Sinus of latitudes of north hemisphere (from north pole toward equator only)
    12) KNIVER: Number of vertical levels
    13) PREFER: Reference pressure (multiplying factor of the first function of hybrid coordinate)
    14) PAHYBR: Values of "A" function of the hybrid coordinate at LAYERiS BOUNDARIES
    15) PBHYBR: Values of "B" function of the hybrid coordinate at LAYERiS BOUNDARIES
    16) LDGARD: True if "cadre" must be kept even after the last file attached is closed
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_char_p(args[4].ljust(16)), IN),
            (c_longlong(), OUT),
            (c_double(), OUT),
            (c_double(), OUT),
            (c_double(), OUT),
            (c_double(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (np.ndarray((args[0],), dtype=np.int64), OUT),
            (np.ndarray((args[1],), dtype=np.int64), OUT),
            (np.ndarray((args[1],), dtype=np.float64), OUT),
            (c_longlong(), OUT),
            (c_double(), OUT),
            (np.ndarray((args[3] + 1,), dtype=np.float64), OUT),
            (np.ndarray((args[3] + 1,), dtype=np.float64), OUT),
            (c_bool(), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfadies(*args):
    """
    Get the date and time of a field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    
    Returns:\n
    1) KDATEF: array of date elements
    """
    return [(c_longlong(args[0]), IN),
            (np.ndarray((11,), dtype=np.int64), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfadiex(*args):
    """
    Get the date and time of a field (precision to the second).
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    
    Returns:\n
    1) KDATEF: array of date elements (precision to the second)
    """
    return [(c_longlong(args[0]), IN),
            (np.ndarray((22,), dtype=np.int64), OUT)]

@ctypesFF
def wfacade(*args):
    """
    Set parameters in header.

    Args:\n
    1) CDNOMC: name of "cadre"
    2) KTYPTR: Type of horizontal transformation
    3) PSLAPO: Sinus of latitude of pole of interest
    4) PCLOPO: Cosinus of longitude of pole of interest
    5) PSLOPO: Sinus of longitude of pole of interest
    6) PCODIL: Dilatation coefficient
    7) KTRONC: Troncature
    8) KNLATI: Nomber of latitudes (from pole to pole)
    9) KNXLON: Maximum number of longitudes by parallel
    10) KSNLOPA: Size of KNLOPA
    11) KNLOPA: Number of longitudes by parallel (from north pole toward equator only)
    12) KSNOZPA: Size of KNOZPA
    13) KNOZPA: Maximum zonal wave number by parallel (from north pole toward equator only)
    14) KSSINLA: Size of PSINLA
    15) PSINLA: Sinus of latitudes of north hemisphere (from north pole toward equator only)
    16) KNIVER: Number of vertical levels
    17) PREFER: Reference pressure (multiplying factor of the first function of hybrid coordinate)
    18) PAHYBR: Values of "A" function of the hybrid coordinate at LAYERiS BOUNDARIES
    19) PBHYBR: Values of "B" function of the hybrid coordinate at LAYERiS BOUNDARIES
    20) LDGARD: True if "cadre" must be kept even after the last file attached is closed
    """
    return [(c_char_p(args[0].ljust(16)), IN),
            (c_longlong(args[1]), IN),
            (c_double(args[2]), IN),
            (c_double(args[3]), IN),
            (c_double(args[4]), IN),
            (c_double(args[5]), IN),
            (c_longlong(args[6]), IN),
            (c_longlong(args[7]), IN),
            (c_longlong(args[8]), IN),
            (c_longlong(args[9]), IN),
            (args[10], IN),
            (c_longlong(args[11]), IN),
            (args[12], IN),
            (c_longlong(args[13]), IN),
            (args[14], IN),
            (c_longlong(args[15]), IN),
            (c_double(args[16]), IN),
            (args[17], IN),
            (args[18], IN),
            (c_bool(args[19]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfagote(*args):
    """
    Set compression parameters of a field in file.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) KNGRIB: Grib encoding level (-1,0,1,2,3);
    3) KNBPDG: Number of bits by grid point value
    4) KNBCSP: Number of bits by real/imaginary part of spectral coefficient
    5) KSTRON: Unpacked under truncature
    6) KPUILA: laplacian power
    7) KDMOPL: KPUILA level of modulation
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_longlong(args[4]), IN),
            (c_longlong(args[5]), IN),
            (c_longlong(args[6]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfautif(*args):
    """
    Set identifier of file.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDIDEN: identifier of file
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(80)), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfandar(*args):
    """
    Set the date and time of a field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) KDATEF: array of date elements
    """
    return [(c_longlong(args[0]), IN),
            (args[1], IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfandax(*args):
    """
    Set the date and time of a field (precision to the second).
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) KDATEF: array of date elements (precision to the second)
    """
    return [(c_longlong(args[0]), IN),
            (args[1], IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfacile(*args):
    """
    Read a 2D field.
    
    Args:\n
    1) KSIZE: size of array to read
    2) KNUMER: logical unit number associated to file
    3) CDPREF: potential prefix of record name
    4) KNIVAU: potential vertical level
    5) CDSUFF: potential suffix of record name
    6) LDCOSP: true if spectral
    
    Returns:\n
    6) PCHAMP: float values read
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_char_p(args[2].ljust(4)), IN),
            (c_longlong(args[3]), IN),
            (c_char_p(args[4].ljust(12)), IN),
            (np.ndarray((args[0],), dtype=np.float64), OUT),
            (c_bool(args[5]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfacilo(*args):
    """
    Read a 2D field. With reordering of spectral fields.
    
    Args:\n
    1) KSIZE: size of array to read
    2) KNUMER: logical unit number associated to file
    3) CDPREF: potential prefix of record name
    4) KNIVAU: potential vertical level
    5) CDSUFF: potential suffix of record name
    6) LDCOSP: true if spectral
    
    Returns:\n
    6) PCHAMP: float values read
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_char_p(args[2].ljust(4)), IN),
            (c_longlong(args[3]), IN),
            (c_char_p(args[4].ljust(12)), IN),
            (np.ndarray((args[0],), dtype=np.float64), OUT),
            (c_bool(args[5]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfaienc(*args):
    """
    Write a 2D field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDPREF: potential prefix of record name
    3) KNIVAU: potential vertical level
    4) CDSUFF: potential suffix of record name
    5) KSIZE: size of PCHAMP
    6) PCHAMP: float values to write
    7) LDCOSP: true if spectral
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(4)), IN),
            (c_longlong(args[2]), IN),
            (c_char_p(args[3].ljust(12)), IN),
            (c_longlong(args[4]), IN),
            (args[5], IN),
            (c_bool(args[6]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfaieno(*args):
    """
    Write a 2D field. With reordering of spectral fields.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDPREF: potential prefix of record name
    3) KNIVAU: potential vertical level
    4) CDSUFF: potential suffix of record name
    5) KSIZE: size of PCHAMP
    6) PCHAMP: float values to write
    7) LDCOSP: true if spectral
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(4)), IN),
            (c_longlong(args[2]), IN),
            (c_char_p(args[3].ljust(12)), IN),
            (c_longlong(args[4]), IN),
            (args[5], IN),
            (c_bool(args[6]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfalais(*args):
    """
    Write a meta-field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record
    3) KLONGD: length of PDONNE
    
    Returns:\n
    1) PDONNE: data to read
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(16)), IN),
            (np.ndarray((args[2],), dtype=np.float64), OUT),
            (c_longlong(args[2]), IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfaisan(*args):
    """
    Write a meta-field.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record
    3) KSIZE: Size of PDONNE
    4) PDONNE: data to write
    """
    return [(c_longlong(args[0]), IN),
            (c_char_p(args[1].ljust(16)), IN),
            (c_longlong(args[2]), IN),
            (args[3], IN)]

@treatReturnCode
@ctypesFF
@addReturnCode
def wfairme(*args):
    """
    Close the FA file.
    
    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDSTTU: status ('KEEP', 'DELETE', 'DEFAUT')
    """
    return[(c_longlong(args[0]), IN),
           (c_char_p(args[1].ljust(7)), IN)]
