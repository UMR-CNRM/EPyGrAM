#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
This module contains functions to:

- to convert vertical coordinates systems from the one to another;
- compute air specific gas constant R according to specific humidity
  and hydrometeors.
- to plot vertical profiles.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy

from . import config

# Constants
Cpd = config.Cpd
Rd = config.Rd
Rv = config.Rv


################
### FORMULAS ###
################
def hybridP2fluxpressure(A, B, Psurf):
    """
    Computes the pressure at flux levels defined by hybrid-pressure
    coefficients. Coefficients are assumed to define flux levels.

    - *A* = table of A coefficients \n
    - *B* = table of B coefficients \n
    - *Psurf* = surface pressure in Pa. \n
    A and B must not contain the first coefficient (A=B=0.)

    Returns a table of pressure.
    """

    if not len(A) == len(B):
        raise ValueError("A, B must have the same size.")
    if not isinstance(Psurf, numpy.ndarray):
        if isinstance(Psurf, int) or isinstance(Psurf, float):
            Psurf = numpy.array([Psurf])
        else:
            Psurf = numpy.array(Psurf)

    pi_tilde = numpy.zeros([len(A)] + list(Psurf.shape))
    # computation
    for k in range(len(A)):
        pi_tilde[k] = A[k] + B[k] * Psurf

    if not numpy.all(pi_tilde[1:] >= pi_tilde[:-1]):
        raise ValueError('pi_tilde array must be in ascending order.')

    return pi_tilde


def hybridP2masspressure(A, B, Psurf, vertical_mean):
    """
    Computes the pressure at mass levels defined by hybrid-pressure
    coefficients. Coefficients are assumed to define flux levels.

    - *A* = table of A coefficients \n
    - *B* = table of B coefficients \n
    - *Psurf* = surface pressure in Pa. \n
    - *vertical_mean* = kind of vertical averaging
    A and B must not contain the first coefficient (A=B=0.)

    Returns a table of pressure.
    """

    pi_tilde = hybridP2fluxpressure(A, B, Psurf)
    pi = flux2masspressures(pi_tilde, vertical_mean)

    if not numpy.all(pi[1:] >= pi[:-1]):
        raise ValueError('pi array must be in ascending order.')

    return pi


def hybridH2fluxheight(A, B, Zsurf, conv2height=False):
    """
    Computes the altitude at flux levels defined by hybrid-height coefficients.
    Coefficients are assumed to define flux levels.

    - *A* = table of A coefficients \n
    - *B* = table of B coefficients \n
    - *Zsurf* = surface altitude in m. \n
    - *conv2height* True to compute height instead of altitude
    A and B must include the underground level

    Returns a table of altitudes.
    """

    if not len(A) == len(B):
        raise ValueError("A, B must have the same size.")
    if not isinstance(Zsurf, numpy.ndarray):
        if isinstance(Zsurf, int) or isinstance(Zsurf, float):
            Zsurf = numpy.array([Zsurf])
        else:
            Zsurf = numpy.array(Zsurf)

    H_tilde = numpy.zeros([len(A)] + list(Zsurf.shape))
    # computation
    for k in range(len(A)):
        H_tilde[k] = A[k] + B[k] * Zsurf

    if conv2height:
        H_tilde -= Zsurf

    if not numpy.all(H_tilde[1:] >= H_tilde[:-1]):
        raise ValueError('H_tilde array must be in ascending order.')

    return H_tilde


def hybridH2massheight(A, B, Zsurf, conv2height=False):
    """
    Computes the altitude at mass levels defined by hybrid-height
    coefficients. Coefficients are assumed to define flux levels.

    - *A* = table of A coefficients \n
    - *B* = table of B coefficients \n
    - *Zsurf* = surface altitude in m. \n
    - *conv2height* True to compute height instead of altitude
    A and B must include the underground level

    Returns a table of altitudes.
    """

    if not len(A) == len(B):
        raise ValueError("A, B must have the same size.")

    H_tilde = hybridH2fluxheight(A, B, Zsurf, conv2height=conv2height)
    H = flux2massheights(H_tilde)

    if not numpy.all(H[1:] >= H[:-1]):
        raise ValueError('H array must be in ascending order.')

    return H


def flux2masspressures(pi_tilde, vertical_mean, Ptop=config.default_Ptop):
    """Converts pressures at flux levels to mass levels."""

    if not numpy.all(pi_tilde[1:] >= pi_tilde[:-1]):
        raise ValueError('pi_tilde array must be in ascending order.')

    L = len(pi_tilde)
    if not isinstance(pi_tilde, numpy.ndarray):
        pi_tilde = numpy.array(pi_tilde)

    pi = numpy.zeros(pi_tilde.shape)
    for k in range(1, L + 1):
        ik = k - 1  # python arranging
        if vertical_mean == 'geometric':
            if k == 1:
                pi[ik] = (pi_tilde[ik] - Ptop) / (1. + Cpd / Rd)
            else:
                pi[ik] = numpy.sqrt(pi_tilde[ik] * pi_tilde[ik - 1])
        elif vertical_mean == 'arithmetic':
            if k == 1:
                pi[ik] = (pi_tilde[ik] + Ptop) / 2.
            else:
                pi[ik] = (pi_tilde[ik] + pi_tilde[ik - 1]) / 2.
        else:
            raise NotImplementedError("vertical_mean not among" +
                                      " ('geometric', 'arithmetic').")

    return pi


def mass2fluxpressures(pi, vertical_mean, Ptop=config.default_Ptop):
    """Converts pressures at mass levels to flux levels."""

    if not numpy.all(pi[1:] >= pi[:-1]):
        raise ValueError('pi array must be in ascending order.')

    L = len(pi)
    if not isinstance(pi, numpy.ndarray):
        pi = numpy.array(pi)

    pi_tilde = numpy.zeros(pi.shape)
    for k in range(1, L + 1):
        ik = k - 1  # python arranging
        if vertical_mean == 'geometric':
            if k == 1:
                pi_tilde[ik] = (pi[ik] + Ptop) * (1. + Cpd / Rd)
            else:
                pi_tilde[ik] = pi[ik] ** 2 / pi_tilde[ik - 1]
        elif vertical_mean == 'arithmetic':
            if k == 1:
                pi_tilde[ik] = (pi[ik] + Ptop) / 2
            else:
                pi_tilde[ik] = (pi[ik] + pi[ik - 1]) / 2.
        else:
            raise NotImplementedError("vertical_mean not among" +
                                      " ('geometric', 'arithmetic').")

    return pi_tilde


def flux2massheights(H_tilde):
    """Converts altitudes at flux levels to mass levels."""

    if not numpy.all(H_tilde[1:] >= H_tilde[:-1]):
        raise ValueError('H_tilde array must be in ascending order.')

    if not isinstance(H_tilde, numpy.ndarray):
        H_tilde = numpy.array(H_tilde)

    H = numpy.zeros(H_tilde.shape)
    H[:-1] = (H_tilde[:-1] + H_tilde[1:]) / 2.
    H[-1] = 2 * H_tilde[-1] - H_tilde[-2]

    return H[1:-1]


def mass2fluxheights(H):
    """Converts altitudes at mass levels to flux levels."""

    if not numpy.all(H[1:] >= H[:-1]):
        raise ValueError('H array must be in ascending order.')

    L = len(H)
    if not isinstance(H, numpy.ndarray):
        H = numpy.array(H)

    H_tilde = numpy.zeros(H.shape)
    raise NotImplementedError('mass2fluxheights is not yet implemented.')
#    for k in range(1, L+1):
#        ik = k-1 # python arranging
#        if k == 1:
#            pi_tilde[ik] = (pi[ik] + Ptop) * (1. + Cpd/Rd)
#        else:
#            pi_tilde[ik] = pi[ik]**2 / pi_tilde[ik-1]

    return H_tilde


def pressure2altitude(R, T, vertical_mean,
                      pi=None, pi_tilde=None, Pdep=None, Phi_surf=None):
    """
    Converts a pressure profile to mass-levels altitude (= height above ground
    if *Phi_surf == 0*).

    The pressure can be given at mass levels (*pi*) or flux levels (*pi_tilde*),
    or both (assuming they are consistent), with unit: Pa.

    - *R* = table of R = specific gas constant (J/kg/K) \n
    - *T* = table of Temperature (K) \n
    - *vertical_mean* = kind of vertical averaging
    - *Pdep* = table of Pressure departure: P-Pi, where P is the hydrostatic
      pressure and Pi the Non-Hydrostatic pressure (Pa). \n
      0. if Hydrostatic (default) \n
    - *Phi_surf* = surface geopotential (J/kg)
    """

    # get full pressures
    if pi is None:
        pi = flux2masspressures(pi_tilde, vertical_mean)
    elif pi_tilde is None:
        pi_tilde = mass2fluxpressures(pi, vertical_mean)
    else:
        if not isinstance(pi, numpy.ndarray):
            pi = numpy.array(pi)
        if not isinstance(pi_tilde, numpy.ndarray):
            pi_tilde = numpy.array(pi_tilde)

    L = len(pi)
    if not len(pi) == len(pi_tilde) == len(R) == len(T):
        raise ValueError("R, T, pi, pi_tilde must have the same size.")
    if not isinstance(R, numpy.ndarray):
        R = numpy.array(R)
    if not isinstance(T, numpy.ndarray):
        T = numpy.array(T)
    if Phi_surf is None:
        myPhi_surf = numpy.zeros(T[0].shape)
    else:
        if not isinstance(Phi_surf, numpy.ndarray):
            myPhi_surf = numpy.array(Phi_surf)
        else:
            myPhi_surf = Phi_surf
    if Pdep is None:
        myPdep = numpy.zeros(R.shape)
    else:
        if not isinstance(Pdep, numpy.ndarray):
            myPdep = numpy.array(Pdep)
        else:
            myPdep = Pdep
    if not pi.shape == pi_tilde.shape:
        raise ValueError("pi and pi_tilde must have the same shape.")
    if not R.shape == T.shape == myPdep.shape:
        raise ValueError("R, T, and Pdep must have the same shape.")

    if (not T[0].shape == myPhi_surf.shape) and (not T[0].size == myPhi_surf.size == 1):
        raise ValueError("Phi_surf shape must be compatible with other shapes.")

    # Pre-computation
    delta = numpy.zeros(pi.shape)
    alpha = numpy.zeros(pi.shape)
    for k in range(1, L + 1):
        ik = k - 1  # python arranging
        if k == 1:
            delta[ik] = 1. + Cpd / Rd
            alpha[ik] = 1.
        else:
            delta[ik] = (pi_tilde[ik] - pi_tilde[ik - 1]) / pi[ik]
            alpha[ik] = 1 - pi[ik] / pi_tilde[ik]

    # Geopotential
    Phi = numpy.zeros(T.shape)
    partialsum = numpy.zeros(T.shape)
    for k in reversed(range(1, L + 1)):
        ik = k - 1  # python arranging
        if k == L:
            partialsum[ik] = 0.
        else:
            partialsum[ik] = partialsum[ik + 1] + R[ik + 1] * T[ik + 1] * \
                             delta[ik + 1] / (1. + myPdep[ik + 1] / pi[ik + 1])
        Phi[ik] = myPhi_surf + partialsum[ik] + R[ik] * T[ik] * alpha[ik] / \
                  (1. + myPdep[ik] / pi[ik])
    # Altitude
    z = Phi / config.g0

    return z


def hybridP2altitude(A, B, R, T, Psurf, vertical_mean,
                     Pdep=None, Phi_surf=None, Ptop=config.default_Ptop):
    """
    Computes the altitude of mass levels defined by hybrid-pressure
    coefficients. (= height above ground if *Phi_surf == 0*).

    - *A* = table of A coefficients \n
    - *B* = table of B coefficients \n
    - *Psurf* = surface pressure in Pa. \n
    - *R* = table of R = specific gas constant (J/kg/K) \n
    - *T* = table of Temperature (K) \n
    - *vertical_mean* = kind of vertical averaging
    - *Pdep* = table of Pressure departure: P-Pi, where P is the hydrostatic
      pressure and Pi the Non-Hydrostatic pressure (Pa). \n
      0. if Hydrostatic (default) \n
    - *Phi_surf* = surface geopotential (J/kg)
    A and B must not contain the first coefficient (A=B=0.)
    """

    if not len(A) == len(B) == len(R) == len(T):
        raise ValueError("A, B, R, T must have the same size.")

    pi_tilde = hybridP2fluxpressure(A, B, Psurf)
    pi = flux2masspressures(pi_tilde, vertical_mean, Ptop=Ptop)
    z = pressure2altitude(R, T, vertical_mean,
                          pi=pi, pi_tilde=pi_tilde, Pdep=Pdep,
                          Phi_surf=Phi_surf)

    return z
