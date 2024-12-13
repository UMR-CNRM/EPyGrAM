#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Some moving utilities...
"""

import json
import os

from . import config


def cmap2json(cmap_filename):
    """
    Convert a .cmap colormap to json formatting, in the spirit of
    bronx.graphics.colormapping.get_ColormapHelper_fromfile()

    :param cmap_filename: file to be processed. Must be named /path/to/{colormap}.cmap
    """
    assert cmap_filename.endswith('.cmap')
    colormap = os.path.basename(cmap_filename).rstrip('.cmap')
    path = os.path.dirname(cmap_filename)
    asdict = {}
    asdict['name'] = colormap
    with open(cmap_filename, 'r') as sourcefile:
        colors = sourcefile.readlines()
        for i in range(len(colors)):
            colors[i] = colors[i].replace(';', '')
            colors[i] = colors[i].replace('[', '')
            colors[i] = colors[i].replace(']', '')
            colors[i] = colors[i].replace('\n', '')
            colors[i] = colors[i].split(',')
            colors[i] = [float(j) for j in colors[i]]
    asdict['colors_RGB'] = colors
    if colormap in config.epygram_colormaps_scaling_labels:
        asdict['colorcenters'] = config.epygram_colormaps_scaling_labels[colormap]
        asdict['normalize'] = True
    elif colormap in config.epygram_colormaps_scaling:
        asdict['colorbounds'] = config.epygram_colormaps_scaling[colormap]
        asdict['ticks'] = 'colorbounds'
        asdict['normalize'] = True
    with open(os.path.join(path, colormap + '.json', 'w')) as out:
        json.dump(asdict, out, indent=2)
