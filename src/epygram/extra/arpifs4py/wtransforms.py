#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for trans/etrans library.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy as np

from . import ctypesFF, IN, OUT, treatReturnCode, addReturnCode
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
def w_etrans_inq(KSIZEI, KSIZEJ,
                 KPHYSICALSIZEI, KPHYSICALSIZEJ,
                 KTRUNCX, KTRUNCY,
                 KNUMMAXRESOL,
                 PDELATX, PDELATY):
    """
    Simplified wrapper to ETRANS_INQ.

    Args:\n
    1,2) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    3,4) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    5,6) KTRUNCX, KTRUNCY: troncatures
    7) KNUMMAXRESOL: maximum number of troncatures handled
    8,9) PDELTAX, PDELTAY: resolution along x,y axis

    Returns:\n
    1) KGPTOT: number of gridpoints
    2) KSPEC: number of spectral coefficients
    """
    return ([KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             PDELATX, PDELATY],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def w_trans_inq(KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL):
    """
    Simplified wrapper to TRANS_INQ.

    Args:\n
    1) KSIZEJ: number of latitudes in grid-point space
    2) KTRUNC: troncature
    3) KSLOEN: Size of KLOEN
    4) KLOEN: number of points on each latitude row
    5) KNUMMAXRESOL: maximum number of troncatures handled

    Returns:\n
    1) KGPTOT: number of gridpoints
    2) KSPEC: number of spectral coefficients
    3) KNMENG: cut-off zonal wavenumber
    """
    return ([KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KSLOEN,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def w_spec2gpt_lam(KSIZEI, KSIZEJ,
                   KPHYSICALSIZEI, KPHYSICALSIZEJ,
                   KTRUNCX, KTRUNCY,
                   KNUMMAXRESOL,
                   KSIZE,
                   LGRADIENT,
                   LREORDER,
                   PDELTAX, PDELTAY,
                   PSPEC):
    """
    Transform spectral coefficients into grid-point values.

    Args:\n
    1,2) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    3,4) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    5,6) KTRUNCX, KTRUNCY: troncatures
    7) KNUMMAXRESOL: maximum number of troncatures handled
    8) KSIZE: size of PSPEC
    9) LGRADIENT: gradient computation
    10) LREORDER: reorder spectral coefficients or not
    11,12) PDELTAX,PDELTAY: resolution along x,y axis
    13) PSPEC: spectral coefficient array

    Returns:\n
    1) PGPT: grid-point field
    2) PGPTM: N-S derivative if LGRADIENT
    3) PGPTL: E-W derivative if LGRADIENT
    """
    return ([KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             KSIZE,
             LGRADIENT,
             LREORDER,
             PDELTAX, PDELTAY,
             PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (bool, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.float64, (KSIZE,), IN),
             (np.float64, (KSIZEI * KSIZEJ,), OUT),
             (np.float64, (KSIZEI * KSIZEJ,), OUT),
             (np.float64, (KSIZEI * KSIZEJ,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def w_gpt2spec_lam(KSIZE,
                   KSIZEI, KSIZEJ,
                   KPHYSICALSIZEI, KPHYSICALSIZEJ,
                   KTRUNCX, KTRUNCY,
                   KNUMMAXRESOL,
                   PDELTAX, PDELTAY,
                   LREORDER,
                   PGPT):
    """
    Transform grid point values into spectral coefficients.

    Args:\n
    1) KSIZE: size of spectral field
    2,3) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    4,5) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    6,7) KTRUNCX, KTRUNCY: troncatures
    8) KNUMMAXRESOL: maximum number of troncatures handled
    9,10) PDELTAX, PDELTAY: resolution along x,y axis
    11) LREORDER: reorder spectral coefficients or not
    12) PGPT: grid-point field

    Returns:\n
    1) PSPEC: spectral coefficient array
    """
    return ([KSIZE,
             KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             PDELTAX, PDELTAY,
             LREORDER,
             PGPT],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.float64, None, IN),
             (np.float64, None, IN),
             (bool, None, IN),
             (np.float64, (KSIZEI * KSIZEJ,), IN),
             (np.float64, (KSIZE,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def w_spec2gpt_gauss(KSIZEJ,
                     KTRUNC,
                     KNUMMAXRESOL,
                     KGPTOT,
                     KSLOEN,
                     KLOEN,
                     KSIZE,
                     LGRADIENT,
                     LREORDER,
                     PSPEC):
    """
    Transform spectral coefficients into grid-point values.

    Args:\n
    1) KSIZEJ: Number of latitudes
    2) KTRUNC: troncature
    3) KNUMMAXRESOL: maximum number of troncatures handled
    4) KGPTOT: number of grid-points
    5) KSLOEN: Size of KLOEN
    6) KLOEN:
    7) KSIZE: Size of PSPEC
    8) LGRADIENT: compute derivatives
    9) LREORDER: reorder spectral coefficients or not
    10) PSPEC: spectral coefficient array

    Returns:\n
    1) PGPT: grid-point field
    2) PGPTM: N-S derivative if LGRADIENT
    3) PGPTL: E-W derivative if LGRADIENT
    """
    return ([KSIZEJ,
             KTRUNC,
             KNUMMAXRESOL,
             KGPTOT,
             KSLOEN,
             KLOEN,
             KSIZE,
             LGRADIENT,
             LREORDER,
             PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (bool, None, IN),
             (np.float64, (KSIZE,), IN),
             (np.float64, (KGPTOT,), OUT),
             (np.float64, (KGPTOT,), OUT),
             (np.float64, (KGPTOT,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def w_gpt2spec_gauss(KSPEC,
                     KSIZEJ,
                     KTRUNC,
                     KNUMMAXRESOL,
                     KSLOEN,
                     KLOEN,
                     KSIZE,
                     LREORDER,
                     PGPT):
    """
    Transform grid-point values into spectral coefficients.

    Args:\n
    1) KSPEC: size of spectral coefficients array
    2) KSIZEJ: Number of latitudes
    3) KTRUNC: troncature
    4) KNUMMAXRESOL: maximum number of troncatures handled
    5) KSLOEN: Size of KLOEN
    6) KLOEN
    7) KSIZE: Size of PGPT
    8) LREORDER: reorder spectral coefficients or not
    9) PGPT: grid-point field

    Returns:\n
    1) PSPEC: spectral coefficient array
    """
    return ([KSPEC,
             KSIZEJ,
             KTRUNC,
             KNUMMAXRESOL,
             KSLOEN,
             KLOEN,
             KSIZE,
             LREORDER,
             PGPT],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.float64, (KSIZE,), IN),
             (np.float64, (KSPEC,), OUT)],
            None)


@ctypesFF()
def w_spec2gpt_fft1d(KSIZES, KTRUNC, PSPEC, KSIZEG):
    """
    Transform spectral coefficients into grid-point values,
    for a 1D array (vertical section academic model)

    Args:\n
    1) KSIZES size of PSPEC
    2) KTRUNC: troncature
    3) PSPEC: spectral coefficient array
    4) KSIZEG: size of grid-point field (with extension zone)

    Returns:\n
    1) PGPT: grid-point field
    """
    return ([KSIZES, KTRUNC, PSPEC, KSIZEG],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.float64, (KSIZES,), IN),
             (np.int64, None, IN),
             (np.float64, (KSIZEG,), OUT)],
            None)
