#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Wrappers for trans/etrans library.
"""

from ctypes import c_longlong, c_bool, c_double
import numpy as np

from . import ctypesFF, IN, OUT, treatReturnCode, addReturnCode



@treatReturnCode
@ctypesFF
@addReturnCode
def w_etrans_inq(*args):
    """
    Simplified wrapper to ETRANS_INQ.
    
    Args:\n
    1,2) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    3,4) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    5,6) KTRUNCX, KTRUNCY: troncatures
    7) KNUMMAXRESOL: maximum number of troncatures handled
    8) PDELATX: resolution along x axis
    9) PDELATY: resolution along y axis
    
    Returns:\n
    1) KGPTOT: number of gridpoints
    2) KSPEC: number of spectral coefficients
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_longlong(args[4]), IN),
            (c_longlong(args[5]), IN),
            (c_longlong(args[6]), IN),
            (c_double(args[7]), IN),
            (c_double(args[8]), IN),
            (c_longlong(), OUT),
            (c_longlong(), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def w_trans_inq(*args):
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
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (args[3], IN),
            (c_longlong(args[4]), IN),
            (c_longlong(), OUT),
            (c_longlong(), OUT),
            (np.ndarray((args[0],), dtype=np.int64), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def w_spec2gpt_lam(*args):
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
    11) PDELTAX: resolution along x axis
    12) PDELTAY: resolution along y axis
    13) PSPEC: spectral coefficient array
    
    Returns:\n
    1) PGPT: grid-point field
    2) PGPTM: N-S derivative if LGRADIENT
    3) PGPTL: E-W derivative if LGRADIENT
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_longlong(args[4]), IN),
            (c_longlong(args[5]), IN),
            (c_longlong(args[6]), IN),
            (c_longlong(args[7]), IN),
            (c_bool(args[8]), IN),
            (c_bool(args[9]), IN),
            (c_double(args[10]), IN),
            (c_double(args[11]), IN),
            (args[12], IN),
            (np.ndarray((args[0] * args[1],), dtype=np.float64), OUT),
            (np.ndarray((args[0] * args[1],), dtype=np.float64), OUT),
            (np.ndarray((args[0] * args[1],), dtype=np.float64), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def w_gpt2spec_lam(*args):
    """
    Transform grid point values into spectral coefficients.

    Args:\n
    1) KSIZE: size of spectral field
    2,3) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    4,5) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    6,7) KTRUNCX, KTRUNCY: troncatures
    8) KNUMMAXRESOL: maximum number of troncatures handled
    9) PDELTAX: resolution along x axis
    10) PDELTAY: resolution along y axis
    11) LREORDER: reorder spectral coefficients or not
    12) PGPT: grid-point field
    
    Returns:\n
    1) PSPEC: spectral coefficient array
    """
    return[(c_longlong(args[0]), IN),
           (c_longlong(args[1]), IN),
           (c_longlong(args[2]), IN),
           (c_longlong(args[3]), IN),
           (c_longlong(args[4]), IN),
           (c_longlong(args[5]), IN),
           (c_longlong(args[6]), IN),
           (c_longlong(args[7]), IN),
           (c_double(args[8]), IN),
           (c_double(args[9]), IN),
           (c_bool(args[10]), IN),
           (args[11], IN),
           (np.ndarray((args[0],), dtype=np.float64), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def w_spec2gpt_gauss(*args):
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
    8) LREORDER: reorder spectral coefficients or not
    9) PSPEC: spectral coefficient array
    
    Returns:\n
    1) PGPT: grid-point field
    """
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_longlong(args[4]), IN),
            (args[5], IN),
            (c_longlong(args[6]), IN),
            (c_bool(args[7]), IN),
            (args[8], IN),
            (np.ndarray((args[3],), dtype=np.float64), OUT)]

@treatReturnCode
@ctypesFF
@addReturnCode
def w_gpt2spec_gauss(*args):
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
    return [(c_longlong(args[0]), IN),
            (c_longlong(args[1]), IN),
            (c_longlong(args[2]), IN),
            (c_longlong(args[3]), IN),
            (c_longlong(args[4]), IN),
            (args[5], IN),
            (c_longlong(args[6]), IN),
            (c_bool(args[7]), IN),
            (args[8], IN),
            (np.ndarray((args[0],), dtype=np.float64), OUT)]

@ctypesFF
def w_spec2gpt_fft1d(*args):
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
    return[(c_longlong(args[0]), IN),
           (c_longlong(args[1]), IN),
           (args[2], IN),
           (c_longlong(args[3]), IN),
           (np.ndarray((args[3],), dtype=np.float64), OUT)]
