#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handles spectral parameters and spectral transforms.
"""

import numpy
import os

from bronx.system import memory

from epygram import config, epygramError, epylog
from epygram.util import RecursiveObject

ectrans4py_init_env_kwargs = dict()


def nearest_greater_FFT992compliant_int(guess):
    """
    Returns the first integer *n* greater than *guess* that satisfies
    n = 2^(1+i) x 3^j x 5^k, with (i,j,k) being positive integers.
    """
    # these are defined for dimensions up to 10000 points at least
    M2 = 14
    M3 = 10
    M5 = 7
    fft_compliant_dims = numpy.zeros((M2, M3, M5))
    for i in range(M2):
        for j in range(M3):
            for k in range(M5):
                fft_compliant_dims[i, j, k] = 2 ** (1 + i) * 3 ** j * 5 ** k
    fft_compliant_dims = sorted(list(set(fft_compliant_dims.flatten())))
    for i in range(len(fft_compliant_dims)):
        if fft_compliant_dims[i] >= guess:
            result = int(fft_compliant_dims[i])
            break
    return result


def truncation_from_gridpoint_dims(dimensions, grid='linear', stretching_coef=1.):
    """
    Compute truncation from gridpoint dimensions, according to the kind of
    **grid**.

    :param dimensions: dict containing dimensions, among:
                       {'X':..., 'Y':...} for LAM grids,
                       {'lat_number':..., 'max_lon_number':...} for Gauss grids
    :param grid: how to choose the truncation, among ('linear', 'quadratic',
                 'cubic')
    :param stretching_coef: stretching or dilatation coefficient for stretched
        Gauss grids (has an impact on the gridpoint/spectral dimensions)

    Formula taken from "Spectral transforms in the cycle 45 of ARPEGE/IFS",
    http://www.umr-cnrm.fr/gmapdoc/IMG/pdf/ykts45.pdf
    """
    truncation = {}
    spfactor = {'linear':2, 'quadratic':3, 'cubic':4}[grid]
    if all([k in dimensions.keys() for k in ('X', 'Y')]):
        # LAM
        truncation['in_X'] = (dimensions['X'] - 1) // spfactor
        truncation['in_Y'] = (dimensions['Y'] - 1) // spfactor
    elif all([k in dimensions.keys() for k in ('lat_number',
                                               'max_lon_number')]):
        # Gauss
        if stretching_coef > 1.:
            truncation['max'] = min(2 * dimensions['lat_number'] - 3,
                                    dimensions['max_lon_number'] - 1) // spfactor
        else:
            truncation['max'] = (dimensions['max_lon_number'] - 1) // spfactor
    return truncation


def _guess_compliant_lonlat(nlon):
    nlon = nearest_greater_FFT992compliant_int(nlon)
    if nlon % 4 != 0:
        nlat = nlon // 2 + 1
    else:
        nlat = nlon // 2
    return nlon, nlat


def gridpoint_dims_from_truncation(truncation, grid='linear', stretching_coef=1.):
    """
    Compute truncation from gridpoint dimensions, according to the kind of
    **grid**.

    :param truncation: dict containing truncation info, among:
                       {'in_X':..., 'in_Y':...} for LAM grids,
                       {'max':...} for Gauss grids
    :param grid: how to choose the truncation, among ('linear', 'quadratic',
                 'cubic')
    :param stretching_coef: stretching or dilatation coefficient for stretched
        Gauss grids (has an impact on the gridpoint/spectral dimensions)

    Formula taken from "Spectral transforms in the cycle 45 of ARPEGE/IFS",
    http://www.umr-cnrm.fr/gmapdoc/IMG/pdf/ykts45.pdf
    """
    dimensions = {}
    spfactor = {'linear':2, 'quadratic':3, 'cubic':4}[grid]
    if all([k in truncation.keys() for k in ('in_X', 'in_Y')]):
        # LAM
        dimensions['X'] = spfactor * truncation['in_X'] + 1
        dimensions['X'] = nearest_greater_FFT992compliant_int(dimensions['X'])
        dimensions['Y'] = spfactor * truncation['in_Y'] + 1
        dimensions['Y'] = nearest_greater_FFT992compliant_int(dimensions['Y'])
    elif all([k in truncation.keys() for k in ('max',)]):
        # Gauss
        nlon = spfactor * truncation['max'] + 1
        nlon, nlat = _guess_compliant_lonlat(nlon)
        if stretching_coef > 1. and truncation['max'] > min(2*nlat-3, nlon) // spfactor:
            nlon, nlat = _guess_compliant_lonlat(nlon + 1)
            epylog.warning("Grid will be sub-{}: nlon > 2*tronc + 1 because of stretching !".format(grid))
        dimensions['max_lon_number'] = nlon
        dimensions['lat_number'] = nlat
    return dimensions


def complete_gridpoint_dimensions(longitudes, latitudes,
                                  truncation, grid, stretching_coef):
    """Complete missing gridpoint dimensions from guessing values."""
    if latitudes is None and longitudes is None:
        dims = gridpoint_dims_from_truncation({'max': truncation},
                                              grid=grid,
                                              stretching_coef=stretching_coef)
        latitudes = dims['lat_number']
        longitudes = dims['max_lon_number']
    elif longitudes is None:
        longitudes = 2 * self.latitudes
    elif latitudes is None:
        if longitudes % 4 != 0:
            latitudes = longitudes // 2 + 1
        else:
            latitudes = longitudes // 2
    return longitudes, latitudes


class SpectralGeometry(RecursiveObject):
    """Handles the spectral geometry and transforms for a H2DField."""

    # in order to avoid calling etrans_inq for large truncations
    _usual_SPdatasize_for_global_trunc = {'triangular':{1198:719400,
                                                        1798:1619100}}

    def __init__(self, space, truncation):
        """
        :param space: Name of spectral space.
                      (among ['legendre', 'bi-fourier', 'fourier'])
        :param truncation: Handles the spectral truncation parameters.
        """

        self.add_attr_inlist('space', ['legendre', 'bi-fourier', 'fourier'])
        self.add_attr_dict('truncation')

        self.space = space
        self.truncation = truncation

        super(SpectralGeometry, self).__init__()
        if os.name == 'posix':
            meminfo = memory.LinuxMemInfo()
        else:
            raise NotImplementedError('MemInfo for os.name=={}'.format(os.name))
        self._total_system_memory = meminfo.system_RAM(unit='MiB')

    @classmethod
    def ectrans4py(cls):
        """
        Accessor to the ectrans4py library, not imported before the first call.
        """
        if not hasattr(cls, '_ectrans4py'):
            import ectrans4py as _ectrans4py
            cls._ectrans4py = _ectrans4py
            if ectrans4py_init_env_kwargs.pop('trigger', False):
                _ectrans4py.init_env(**ectrans4py_init_env_kwargs)
        return cls._ectrans4py

    def _prevent_swapping(self):
        if (self.space == 'legendre' and
            (float(self.needed_memory) / (1024 ** 2.) >=
             config.prevent_swapping_legendre * self._total_system_memory)):
            needed_mem_in_mb = float(self.needed_memory) / (1024 ** 2.)
            total_mem_in_mb = float(self._total_system_memory)
            raise epygramError('Legendre spectral transforms need {:.1f} MB \
                                memory, while only {:.1f} MB is available: \
                                SWAPPING prevented !'.format(needed_mem_in_mb,
                                                             total_mem_in_mb))

    def _prevent_limited_stack(self):
        if (self.space == 'legendre' and self.truncation['max'] > 1200):
            epylog.warning('Caution: large Legendre truncation may need very large stacksize !')

    def legendre_known_spectraldata_size(self):
        """In order to avoid calling trans_inq for large truncations."""
        if self.truncation['shape'] in self._usual_SPdatasize_for_global_trunc:
            return self._usual_SPdatasize_for_global_trunc.get(
                self.truncation['shape']).get(self.truncation['max'], None)

    def trans_inq(self, gpdims):
        """
        Wrapper to IAL's TRANS_INQ.

        :param dict gpdims: gridpoints dimensions
        """
        self._prevent_swapping()
        self._prevent_limited_stack()
        return self.ectrans4py().trans_inq4py(gpdims['lat_number'],
                                              self.truncation['max'],
                                              len(gpdims['lon_number_by_lat']),
                                              numpy.array(gpdims['lon_number_by_lat']),
                                              config.KNUMMAXRESOL)

    def etrans_inq(self, gpdims):
        """
        Wrapper to IAL's ETRANS_INQ.

        :param dict gpdims: gridpoints dimensions
        """
        self._prevent_swapping()
        self._prevent_limited_stack()
        return self.ectrans4py().etrans_inq4py(gpdims['X'],
                                               gpdims['Y'],
                                               gpdims['X_CIzone'],
                                               gpdims['Y_CIzone'],
                                               self.truncation['in_X'],
                                               self.truncation['in_Y'],
                                               config.KNUMMAXRESOL,
                                               gpdims['X_resolution'],
                                               gpdims['Y_resolution'])

    @property
    def needed_memory(self):
        """Memory needed for transforms, in bytes."""
        if self.space == 'legendre':
            return self.truncation['max'] ** 3 / 2 * 8
        else:
            raise NotImplementedError('space:' + self.space)

    def sp2gp(self, data, gpdims,
              spectral_coeff_order=config.spectral_coeff_order):
        """
        Makes the transform of the spectral data contained in *data* (assumed
        this spectral geometry is that of *data*) to gridpoint space, defined
        by its dimensions contained in *gpdims*, and returns the gridpoint data.

        :param data: spectral data
        :param dict gpdims: gridpoints dimensions
        :param spectral_coeff_order: among 'model' or 'FA',
          cf. default and description in config.spectral_coeff_order

        Input and output data are both 1D.
        """
        self._prevent_swapping()
        self._prevent_limited_stack()
        if self.space == 'bi-fourier':
            gpdata = self.ectrans4py().sp2gp_lam4py(gpdims['X'],
                                                    gpdims['Y'],
                                                    gpdims['X_CIzone'],
                                                    gpdims['Y_CIzone'],
                                                    self.truncation['in_X'],
                                                    self.truncation['in_Y'],
                                                    config.KNUMMAXRESOL,
                                                    len(data),
                                                    False,  # no derivatives
                                                    spectral_coeff_order != 'model',
                                                    gpdims['X_resolution'],
                                                    gpdims['Y_resolution'],
                                                    data)[0]
        elif self.space == 'legendre':
            gpdata = self.ectrans4py().sp2gp_gauss4py(gpdims['lat_number'],
                                                      self.truncation['max'],
                                                      config.KNUMMAXRESOL,
                                                      sum(gpdims['lon_number_by_lat']),
                                                      len(gpdims['lon_number_by_lat']),
                                                      numpy.array(gpdims['lon_number_by_lat']),
                                                      len(data),
                                                      False,  # no derivatives
                                                      spectral_coeff_order != 'model',
                                                      data)[0]
        elif self.space == 'fourier':
            if self.truncation['in_Y'] > 1:
                gpdata = self.ectrans4py().sp2gp_fft1d4py(len(data),
                                                          self.truncation['in_Y'],
                                                          data,
                                                          gpdims['Y'])
            else:
                gpdata = numpy.ones(gpdims['Y']) * data[0]
        else:
            raise epygramError("unknown spectral space:" + self.space + ".")
        return gpdata

    def gp2sp(self, data, gpdims,
              spectral_coeff_order=config.spectral_coeff_order):
        """
        Makes the transform of the gridpoint data contained in *data*
        to the spectral space and truncation of this object, and returns
        the spectral data.

        :param data: gridpoint data
        :param dict gpdims: gridpoints dimensions
        :param spectral_coeff_order: among 'model' or 'FA',
          cf. default and description in config.spectral_coeff_order

        Input and output data are both 1D.
        """
        self._prevent_swapping()
        self._prevent_limited_stack()
        if self.space == 'bi-fourier':
            SPdatasize = self.etrans_inq(gpdims)[1]
            spdata = self.ectrans4py().gp2sp_lam4py(SPdatasize,
                                                    gpdims['X'],
                                                    gpdims['Y'],
                                                    gpdims['X_CIzone'],
                                                    gpdims['Y_CIzone'],
                                                    self.truncation['in_X'],
                                                    self.truncation['in_Y'],
                                                    config.KNUMMAXRESOL,
                                                    gpdims['X_resolution'],
                                                    gpdims['Y_resolution'],
                                                    spectral_coeff_order != 'model',
                                                    data)
        elif self.space == 'legendre':
            SPdatasize = self.trans_inq(gpdims)[1]
            SPdatasize *= 2  # complex coefficients
            spdata = self.ectrans4py().gp2sp_gauss4py(SPdatasize,
                                                      gpdims['lat_number'],
                                                      self.truncation['max'],
                                                      config.KNUMMAXRESOL,
                                                      len(gpdims['lon_number_by_lat']),
                                                      numpy.array(gpdims['lon_number_by_lat']),
                                                      len(data),
                                                      spectral_coeff_order != 'model',
                                                      data)
        elif self.space == 'fourier':
            # 1D case
            SPdatasize = self.etrans_inq(gpdims)[1]
            if self.truncation['in_Y'] <= 1:
                spdata = numpy.zeros(SPdatasize)
                spdata[0] = data[0]
            else:
                raise NotImplementedError("direct transform for 1D fourier" +
                                          " transform.")
        else:
            raise epygramError("unknown spectral space:" + self.space + ".")
        return spdata

    def compute_xy_spderivatives(self, data, gpdims,
                                 spectral_coeff_order=config.spectral_coeff_order):
        """
        Compute the derivatives of the spectral data contained in *data* in
        spectral space (assumed this spectral geometry is that of *data*) and
        return it in gridpoint space, defined by its dimensions contained in
        *gpdims*.

        :param data: spectral data
        :param dict gpdims: gridpoints dimensions
        :param spectral_coeff_order: among 'model' or 'FA',
          cf. default and description in config.spectral_coeff_order

        Returns: (dz/dx, dz/dy)

        Input and output data are both 1D.
        """
        self._prevent_swapping()
        self._prevent_limited_stack()
        if self.space == 'bi-fourier':
            gpdata = self.ectrans4py().sp2gp_lam4py(gpdims['X'],
                                                    gpdims['Y'],
                                                    gpdims['X_CIzone'],
                                                    gpdims['Y_CIzone'],
                                                    self.truncation['in_X'],
                                                    self.truncation['in_Y'],
                                                    config.KNUMMAXRESOL,
                                                    len(data),
                                                    True,  # derivatives
                                                    spectral_coeff_order != 'model',
                                                    gpdims['X_resolution'],
                                                    gpdims['Y_resolution'],
                                                    data)
            gpdata = (gpdata[2], gpdata[1])
        elif self.space == 'legendre':
            gpdata = self.ectrans4py().sp2gp_gauss4py(gpdims['lat_number'],
                                                      self.truncation['max'],
                                                      config.KNUMMAXRESOL,
                                                      sum(gpdims['lon_number_by_lat']),
                                                      len(gpdims['lon_number_by_lat']),
                                                      numpy.array(gpdims['lon_number_by_lat']),
                                                      len(data),
                                                      True,  # derivatives
                                                      spectral_coeff_order != 'model',
                                                      data)
            gpdata = (gpdata[2], gpdata[1])
        elif self.space == 'fourier':
            raise NotImplementedError('fourier(1D): not yet !')
        else:
            raise epygramError("unknown spectral space:" + self.space + ".")
        return gpdata
