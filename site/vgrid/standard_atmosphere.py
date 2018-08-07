#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Standard atmosphere functions."""
from __future__ import print_function, absolute_import, unicode_literals, division

import os
import numpy
import ctypesForFortran

IN = ctypesForFortran.IN
OUT = ctypesForFortran.OUT
INOUT = ctypesForFortran.INOUT

shared_objects_library = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      'vertical_discretization',
                                      'atmosta.so')
ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)


@ctypesFF()
def presta(KLEV, PSTZ):
    return ([KLEV, PSTZ],
            [(numpy.int64, None, IN),
             (numpy.float64, (KLEV,), IN),
             (numpy.float64, (KLEV,), OUT)],
            None)


@ctypesFF()
def altsta(KLEV, PREHYD):
    return ([KLEV, PREHYD],
            [(numpy.int64, None, IN),
             (numpy.float64, (KLEV,), IN),
             (numpy.float64, (KLEV,), OUT)],
            None)


def pressure_at(altitude):
    """
    Compute the pressure at a series of altitude, considering standard atmosphere.

    For more documentation, cf. arpifs' ppsta.F90
    """
    return presta(len(altitude), altitude)


def altitude_at(pressure):
    """
    Compute the altitude at a series of pressure, considering standard atmosphere.

    For more documentation, cf. arpifs' susta.F90
    """
    return altsta(len(pressure), pressure)
