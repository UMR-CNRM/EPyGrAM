#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for FA library.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy as np

from ctypesForFortran import addReturnCode, treatReturnCode, IN, OUT
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


@ctypesFF()
def get_facst():
    """
    Export maximum sizes used for fa format.

    Returns:\n
    1) JPXPAH
    2) JPXIND
    3) JPXGEO
    4) JPXNIV
    """
    return ([],
            [(np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfaitou(CDFILE, CDSTATE, CDNOMC):
    """
    Open a FA file.

    Args:\n
    1) CDFILE: path to file to open
    2) CDSTATE: state of file ('NEW', 'OLD')
    3) CDNOMC: name of "cadre"

    Returns:\n
    1) KNUMER: logical unit number associated to file
    """
    return ([CDFILE, CDSTATE, CDNOMC],
            [(str, (len(CDFILE.encode('utf8')),), IN),
             (str, (len(CDSTATE.encode('utf8')),), IN),
             (np.int64, None, OUT),
             (str, (16,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfaveur(KNUMER):
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
    return ([KNUMER],
            [(np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfalsif(KNUMER):
    """
    Get identifier of file.

    Args:\n
    1) KNUMER: logical unit number associated to file

    Returns:\n
    1) CDIDEN: identifier of file
    """
    return ([KNUMER],
            [(np.int64, None, IN),
             (str, (80,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfanion(KNUMER, CDPREF, KNIVAU, CDSUFF):
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
    return ([KNUMER, CDPREF, KNIVAU, CDSUFF],
            [(np.int64, None, IN),
             (str, (4,), IN),
             (np.int64, None, IN),
             (str, (12,), IN),
             (bool, None, OUT),
             (bool, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),],
            None)


@ctypesFF()
def wfacies(KXPAH, KXIND, KXGEO, KXNIV, CDNOMC):
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
    return ([KXPAH, KXIND, KXGEO, KXNIV, CDNOMC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (str, (16,), IN),
             (np.int64, None, OUT),
             (np.float64, None, OUT),
             (np.float64, None, OUT),
             (np.float64, None, OUT),
             (np.float64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KXPAH,), OUT),
             (np.int64, (KXIND,), OUT),
             (np.float64, (KXIND,), OUT),
             (np.int64, None, OUT),
             (np.float64, None, OUT),
             (np.float64, (KXNIV + 1,), OUT),
             (np.float64, (KXNIV + 1,), OUT),
             (bool, None, OUT),],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfadies(KNUMER):
    """
    Get the date and time of a field.

    Args:\n
    1) KNUMER: logical unit number associated to file

    Returns:\n
    1) KDATEF: array of date elements
    """
    return ([KNUMER],
            [(np.int64, None, IN),
             (np.int64, (11), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfadiex(KNUMER):
    """
    Get the date and time of a field (precision to the second).

    Args:\n
    1) KNUMER: logical unit number associated to file

    Returns:\n
    1) KDATEF: array of date elements (precision to the second)
    """
    return ([KNUMER],
            [(np.int64, None, IN),
             (np.int64, (22,), OUT)],
            None)


@ctypesFF()
def wfacade(CDNOMC, KTYPTR,
            PSLAPO, PCLOPO, PSLOPO,
            PCODIL,
            KTRONC, KNLATI, KNXLON,
            KSNLOPA, KNLOPA,
            KSNOZPA, KNOZPA,
            KSSINLA, PSINLA,
            KNIVER, PREFER, PAHYBR, PBHYBR,
            LDGARD):
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
    return ([CDNOMC, KTYPTR,
             PSLAPO, PCLOPO, PSLOPO,
             PCODIL,
             KTRONC, KNLATI, KNXLON,
             KSNLOPA, KNLOPA,
             KSNOZPA, KNOZPA,
             KSSINLA, PSINLA,
             KNIVER, PREFER, PAHYBR, PBHYBR,
             LDGARD],
            [(str, (16,), IN),
             (np.int64, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSNLOPA,), IN),
             (np.int64, None, IN),
             (np.int64, (KSNOZPA,), IN),
             (np.int64, None, IN),
             (np.float64, (KSSINLA,), IN),
             (np.int64, None, IN),
             (np.float64, None, IN),
             (np.float64, (KNIVER + 1,), IN),
             (np.float64, (KNIVER + 1,), IN),
             (bool, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfagote(KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL):
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
    return ([KNUMER, KNGRIB, KNBPDG, KNBCSP, KSTRON, KPUILA, KDMOPL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfautif(KNUMER, CDIDEN):
    """
    Set identifier of file.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDIDEN: identifier of file
    """
    return ([KNUMER, CDIDEN],
            [(np.int64, None, IN),
             (str, (80,), IN),],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfandar(KNUMER, KDATEF):
    """
    Set the date and time of a field.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) KDATEF: array of date elements
    """
    return ([KNUMER, KDATEF],
            [(np.int64, None, IN),
             (np.int64, (11,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfandax(KNUMER, KDATEF):
    """
    Set the date and time of a field (precision to the second).

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) KDATEF: array of date elements (precision to the second)
    """
    return ([KNUMER, KDATEF],
            [(np.int64, None, IN),
             (np.int64, (22,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfacile(KSIZE, KNUMER, CDPREF, KNIVAU, CDSUFF, LDCOSP):
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
    1) PCHAMP: float values read
    """
    return ([KSIZE, KNUMER, CDPREF, KNIVAU, CDSUFF, LDCOSP],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (str, (4,), IN),
             (np.int64, None, IN),
             (str, (12,), IN),
             (np.float64, (KSIZE,), OUT),
             (bool, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfacilo(KSIZE, KNUMER, CDPREF, KNIVAU, CDSUFF, LDCOSP):
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
    1) PCHAMP: float values read
    2) LDUNDEF: true if the field has undef values
    3) PDUNDEF: undef value
    """
    return ([KSIZE, KNUMER, CDPREF, KNIVAU, CDSUFF, LDCOSP],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (str, (4,), IN),
             (np.int64, None, IN),
             (str, (12,), IN),
             (np.float64, (KSIZE,), OUT),
             (bool, None, IN),
             (bool, None, OUT),
             (np.float64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfaienc(KNUMER, CDPREF, KNIVAU, CDSUFF, KSIZE, PCHAMP, LDCOSP):
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
    return ([KNUMER, CDPREF, KNIVAU, CDSUFF, KSIZE, PCHAMP, LDCOSP],
            [(np.int64, None, IN),
             (str, (4,), IN),
             (np.int64, None, IN),
             (str, (12,), IN),
             (np.int64, None, IN),
             (np.float64, (KSIZE,), IN),
             (bool, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfaieno(KNUMER, CDPREF, KNIVAU, CDSUFF, KSIZE, PCHAMP, LDCOSP, LDUNDEF, PDUNDEF):
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
    8) LDUNDEF: true if the field has undef values
    9) PDUNDEF: undef value
    """
    return ([KNUMER, CDPREF, KNIVAU, CDSUFF, KSIZE, PCHAMP, LDCOSP, LDUNDEF, PDUNDEF],
            [(np.int64, None, IN),
             (str, (4,), IN),
             (np.int64, None, IN),
             (str, (12,), IN),
             (np.int64, None, IN),
             (np.float64, (KSIZE,), IN),
             (bool, None, IN),
             (bool, None, IN),
             (np.float64, None, IN),],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfalais(KNUMER, CDNOMA, KLONGD):
    """
    Read a meta-field.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record
    3) KLONGD: length of PDONNE

    Returns:\n
    1) PDONNE: data to read
    """
    return ([KNUMER, CDNOMA, KLONGD],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (np.float64, (KLONGD,), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfaisan(KNUMER, CDNOMA, KSIZE, PDONNE):
    """
    Write a meta-field.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDNOMA: name of record
    3) KSIZE: Size of PDONNE
    4) PDONNE: data to write
    """
    return ([KNUMER, CDNOMA, KSIZE, PDONNE],
            [(np.int64, None, IN),
             (str, (16,), IN),
             (np.int64, None, IN),
             (np.float64, (KSIZE,), IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def wfairme(KNUMER, CDSTTU):
    """
    Close the FA file.

    Args:\n
    1) KNUMER: logical unit number associated to file
    2) CDSTTU: status ('KEEP', 'DELETE', 'DEFAUT')
    """
    return ([KNUMER, CDSTTU],
            [(np.int64, None, IN),
             (str, (7,), IN),],
            None)
