#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, division, unicode_literals

from unittest import main

import epygram

from . import abstract_testclasses as abtc

epygram.init_env()


class TestSpectralGaussC1(abtc.TestSpectral):
    basename = 'gaussC1.fa'


class TestSpectralGaussC2p4(abtc.TestSpectral):
    basename = 'gaussC2.4.fa'


class TestSpectralLAM(abtc.TestSpectral):
    basename = 'lambert_HN.fa'


if __name__ == '__main__':
    main(verbosity=2)
