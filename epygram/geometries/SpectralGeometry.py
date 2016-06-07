#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handles spectral parameters and spectral transforms.
"""

import numpy

from footprints import FootprintBase, FPDict

from epygram import config, epygramError
from epygram.util import RecursiveObject



class SpectralGeometry(RecursiveObject, FootprintBase):
    """Handles the spectral geometry and transforms for a H2DField."""

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            space=dict(
                access='rwx',
                info='Name of spectral space.'),
            truncation=dict(
                type=FPDict,
                access='rwx',
                info='Handles the spectral truncation parameters.'),
        )
    )

    def sp2gp(self, data, gpdims,
              spectral_coeff_order=config.spectral_coeff_order):
        """
        Makes the transform of the spectral data contained in *data* (assumed
        this spectral geometry is that of 'data') to gridpoint space, defined
        by its dimensions contained in *gpdims*, and returns the gridpoint data.

        Input and output data are both 1D.
        """
        from arpifs4py import wtransforms

        if self.space == 'bi-fourier':
            gpdata = wtransforms.w_spec2gpt_lam(gpdims['X'],
                                                gpdims['Y'],
                                                gpdims['X_CIzone'],
                                                gpdims['Y_CIzone'],
                                                self.truncation['in_X'],
                                                self.truncation['in_Y'],
                                                config.KNUMMAXRESOL,
                                                len(data),
                                                False,
                                                spectral_coeff_order != 'model',
                                                gpdims['X_resolution'],
                                                gpdims['Y_resolution'],
                                                data)[0]
        elif self.space == 'legendre':
            gpdata = wtransforms.w_spec2gpt_gauss(gpdims['lat_number'],
                                                  self.truncation['max'],
                                                  config.KNUMMAXRESOL,
                                                  sum(gpdims['lon_number_by_lat']),
                                                  len(gpdims['lon_number_by_lat']),
                                                  numpy.array(gpdims['lon_number_by_lat']),
                                                  len(data),
                                                  spectral_coeff_order != 'model',
                                                  data)
        elif self.space == 'fourier':
            if self.truncation['in_Y'] > 1:
                gpdata = wtransforms.w_spec2gpt_fft1d(len(data),
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

        Input and output data are both 1D.
        """
        from arpifs4py import wtransforms

        if self.space == 'bi-fourier':
            SPdatasize = wtransforms.w_etrans_inq(gpdims['X'],
                                                  gpdims['Y'],
                                                  gpdims['X_CIzone'],
                                                  gpdims['Y_CIzone'],
                                                  self.truncation['in_X'],
                                                  self.truncation['in_Y'],
                                                  config.KNUMMAXRESOL,
                                                  gpdims['X_resolution'],
                                                  gpdims['Y_resolution'])[1]
            spdata = wtransforms.w_gpt2spec_lam(SPdatasize,
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
            SPdatasize = wtransforms.w_trans_inq(gpdims['lat_number'],
                                                 self.truncation['max'],
                                                 len(gpdims['lon_number_by_lat']),
                                                 numpy.array(gpdims['lon_number_by_lat']),
                                                 config.KNUMMAXRESOL)[1]
            SPdatasize *= 2  # complex coefficients
            spdata = wtransforms.w_gpt2spec_gauss(SPdatasize,
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
            SPdatasize = wtransforms.w_etrans_inq(gpdims['X'], gpdims['Y'],
                                                  gpdims['X_CIzone'],
                                                  gpdims['Y_CIzone'],
                                                  self.truncation['in_X'],
                                                  self.truncation['in_Y'],
                                                  config.KNUMMAXRESOL,
                                                  gpdims['X_resolution'],
                                                  gpdims['Y_resolution'])[1]
            if self.truncation['in_Y'] <= 1:
                spdata = numpy.zeros(SPdatasize)
                spdata[0] = data[0]
            else:
                raise NotImplementedError("direct transform for 1D fourier" + \
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
        
        Returns: (dz/dx, dz/dy)

        Input and output data are both 1D.
        """
        from arpifs4py import wtransforms

        if self.space == 'bi-fourier':
            gpdata = wtransforms.w_spec2gpt_lam(gpdims['X'],
                                                gpdims['Y'],
                                                gpdims['X_CIzone'],
                                                gpdims['Y_CIzone'],
                                                self.truncation['in_X'],
                                                self.truncation['in_Y'],
                                                config.KNUMMAXRESOL,
                                                len(data),
                                                True,
                                                spectral_coeff_order != 'model',
                                                gpdims['X_resolution'],
                                                gpdims['Y_resolution'],
                                                data)
            gpdata = (gpdata[2], gpdata[1])
        elif self.space == 'legendre':
            raise NotImplementedError('legendre: not yet !')
            #gpdata = wtransforms.w_spec2gpt_gauss(gpdims['lat_number'],
            #                                    self.truncation['max'],
            #                                    config.KNUMMAXRESOL,
            #                                    sum(gpdims['lon_number_by_lat']),
            #                                    len(gpdims['lon_number_by_lat']),
            #                                    numpy.array(gpdims['lon_number_by_lat']),
            #                                    len(data),
            #                                    config.reorder_spectral_coeff,
            #                                    data)
        elif self.space == 'fourier':
            raise NotImplementedError('fourier(1D): not yet !')
            #if self.truncation['in_Y'] > 1:
            #    gpdata = wtransforms.w_spec2gpt_fft1d(len(data),
            #                                        self.truncation['in_Y'],
            #                                        data,
            #                                        gpdims['Y'])
            #else:
            #    gpdata = numpy.ones(gpdims['Y']) * data[0]
        else:
            raise epygramError("unknown spectral space:" + self.space + ".")

        return gpdata
