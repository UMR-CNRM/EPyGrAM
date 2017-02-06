#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
This module has been developped as a temporary alternative to :mod:`pyproj`.

.. deprecated:: 1.0.0
   Use :mod:`pyproj` instead. Shall be removed from the package in a future version.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import math
import numpy

from epygram import config
from epygram.base import RecursiveObject


class Proj(RecursiveObject):
    """
    A class of homemade projections on the sphere, made to handle cases not
    handled by *pyproj* (e.g. sphere geoid for polar stereographic
    projection...).
    """

    def __init__(self, proj, geoidshape=config.myproj_default_geoid['geoidshape'],
                 geoidradius=config.myproj_default_geoid['geoidradius'],
                 x_0=0., y_0=0., **proj_params):
        """
        Initializes parameters for projection formulas.

        Args:

        - *proj*: name of the projection, among ('lambert', 'mercator', 'polar_stereographic')
        - *geoidshape*: actually, only 'sphere' is implemented.
        - *earth_radius*: can be specified, to use a specific earth radius, in m.
        - *x_0*, *y_0* : offset of Origin in (x, y) coordinates of projected map.
        - *proj_params*: a set of arguments that depends on the projection.
          (all lon/lat are in degrees):

          + lambert: lon_0 = longitude of reference point \n
                     lat_1 = first secant latitude \n
                     lat_2 = second secant latitude \n
                     if tangent, lat_1 = lat_2 = lat_0 = latitude of reference point = tangency latitude
          + mercator: lon_0 = longitude of reference point
                      lat_ts = tangency latitude (0°) or secant latitude
          + polar_stereographic: lon_0 = longitude of reference point
                                 lat_0 = +/- 90° depending on the projection pole
                                 lat_ts = secant or tangency latitude
        """

        deg2rad = math.pi / 180.

        self._args = proj_params
        self._p = {'name':proj}
        if geoidshape != 'sphere':
            raise NotImplementedError("projected geometry on ellipsoid.")
        else:
            self._p['earth_radius'] = geoidradius

        self._p['FE'] = x_0
        self._p['FN'] = y_0
        self._p['lambda0'] = proj_params['lon_0'] * deg2rad
        if proj == 'lambert':
            if abs(proj_params['lat_1'] - proj_params['lat_2']) > config.epsilon:
                self.secant = True
                m1 = math.cos(proj_params['lat_1'] * deg2rad)
                m2 = math.cos(proj_params['lat_2'] * deg2rad)
                t1 = math.tan(math.pi / 4. - proj_params['lat_1'] * deg2rad / 2.)
                t2 = math.tan(math.pi / 4. - proj_params['lat_2'] * deg2rad / 2.)
                lat_0 = (proj_params['lat_1'] + proj_params['lat_1']) / 2.
                t0 = math.tan(math.pi / 4. - lat_0 * deg2rad / 2.)
                self._p['n'] = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
                self._p['F'] = m1 / (self._p['n'] * (t1 ** self._p['n']))
                self._p['r0'] = self._p['earth_radius'] * self._p['F'] * t0 ** self._p['n']
            else:
                self.secant = False
                lat_0 = proj_params['lat_1']
                m0 = math.cos(lat_0 * deg2rad)
                t0 = math.tan(math.pi / 4. - lat_0 * deg2rad / 2.)
                self._p['n'] = math.sin(lat_0 * deg2rad)
                self._p['F'] = m0 / (self._p['n'] * t0 ** self._p['n'])
                self._p['r0'] = self._p['earth_radius'] * self._p['F'] * t0 ** self._p['n']
            self._p['sign'] = math.copysign(1., lat_0)
        elif proj == 'mercator':
            if abs(proj_params['lat_ts']) > config.epsilon:
                self.secant = True
                self._p['k0'] = math.cos(proj_params['lat_ts'] * deg2rad)
            else:
                self.secant = False
                self._p['k0'] = 1.
        elif proj == 'polar_stereographic':
            self._p['sign'] = math.copysign(1., proj_params['lat_0'])
            if abs(proj_params['lat_0'] - proj_params['lat_ts']) > config.epsilon:
                self.secant = True
                # mf = math.cos(proj_params['lat_ts'] * deg2rad)
                # tf = math.tan(math.pi / 4. - self._p['sign'] * proj_params['lat_ts'] / 2.)
                # self._p['k0'] = mf / (2.*tf)
                self._p['k0'] = (1 + self._p['sign'] * math.sin(proj_params['lat_ts'] * deg2rad)) / 2.
            else:
                self.secant = False
                self._p['k0'] = 1.

    def __call__(self, x, y, inverse=False):
        """
        Converts lon/lat coordinates (in °) to x/y coord in the projection.

        If *inverse* is *True*, makes the inverse conversion, from x/y to lon/lat (in °).
        """

        if not isinstance(x, numpy.ndarray):
            x = numpy.array(x)
        if not isinstance(y, numpy.ndarray):
            y = numpy.array(y)

        p = self._p
        if not inverse:
            lon = numpy.pi / 180. * x
            lat = numpy.pi / 180. * y
            if self._p['name'] == 'lambert':
                r = p['earth_radius'] * p['F'] * numpy.tan(numpy.pi / 4. - lat / 2.) ** p['n']
                theta = p['n'] * (lon - p['lambda0'])
                E = p['FE'] + r * numpy.sin(theta)
                N = p['FN'] + p['r0'] - r * numpy.cos(theta)
            elif self._p['name'] == 'mercator':
                E = p['FE'] + p['earth_radius'] * p['k0'] * (lon - p['lambda0'])
                N = p['FN'] + p['earth_radius'] * p['k0'] * numpy.log(numpy.tan(numpy.pi / 4. + lat / 2.))
            elif self._p['name'] == 'polar_stereographic':
                t = numpy.tan(numpy.pi / 4. - p['sign'] * lat / 2.)
                rho = 2. * p['earth_radius'] * p['k0'] * t
                E = p['FE'] + rho * numpy.sin(lon - p['lambda0'])
                N = p['FN'] - p['sign'] * rho * numpy.cos(lon - p['lambda0'])

            return (E, N)

        elif inverse:
            if self._p['name'] == 'lambert':
                rp = p['sign'] * numpy.sqrt((x - p['FE']) ** 2 + (p['r0'] - (y - p['FN'])) ** 2)
                tp = (rp / (p['earth_radius'] * p['F'])) ** (1. / p['n'])
                thetap = numpy.arctan((x - p['FE']) / (p['r0'] - (y - p['FN'])))
                lon = thetap / p['n'] + p['lambda0']
                lat = numpy.pi / 2. - 2 * numpy.arctan(tp)
            elif self._p['name'] == 'mercator':
                lon = (x - p['FE']) / (p['earth_radius'] * p['k0']) + p['lambda0']
                lat = 2 * numpy.arctan(numpy.exp((y - p['FN']) / (p['earth_radius'] * p['k0']))) - numpy.pi / 2.
            elif self._p['name'] == 'polar_stereographic':
                rhop = numpy.sqrt((x - p['FE']) ** 2 + (y - p['FN']) ** 2)
                tp = rhop / (2 * p['earth_radius'] * p['k0'])
                lon = p['lambda0'] - (0.5 + 0.5 * p['sign'] * numpy.copysign(1., y - p['FN'])) * numpy.pi \
                     + numpy.arctan(-p['sign'] * (x - p['FE']) / (y - p['FN']))
                lat = p['sign'] * (numpy.pi / 2. - 2 * numpy.arctan(tp))

            return (180. / numpy.pi * lon, 180. / numpy.pi * lat)
