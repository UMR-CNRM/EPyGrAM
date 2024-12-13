#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

_KNOWN_EPYGRAM_REPOS = [
    '/home/common/epygram/public/cartopy-features',  # cnrm
    '/home/gmap/mrpe/mary/public/cartopy-features',  # bullx
    '/soprano/home/marp999/epygram/cartopy-features',  # dsidev TODO:
    '/home/ms/fr/rm9/public/cartopy-features']  # ecmwf TODO:


def update_config(config):
    for p in _KNOWN_EPYGRAM_REPOS:
        if os.path.exists(p) and os.path.isdir(p):
            config['pre_existing_data_dir'] = p
            break
