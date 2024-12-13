#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Make colormaps illustrations for documentation purpose.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from bronx.graphics.colormapping import add_cmap

import epygram

cmaps = [
    ('Perceptually Uniform Sequential',
     ['viridis', 'plasma', 'inferno', 'magma']),
    ('Sequential',
     ['Blues', 'BuGn', 'BuPu',
      'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
      'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
      'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
    ('Sequential2',
     ['afmhot', 'autumn', 'bone', 'cool', 'copper',
      'gist_heat', 'gray', 'hot', 'pink',
      'spring', 'summer', 'winter']),
    ('Diverging',
     ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
      'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
      'seismic']),
    ('Qualitative',
     ['Accent', 'Dark2', 'Paired', 'Pastel1',
      'Pastel2', 'Set1', 'Set2', 'Set3']),
    ('Miscellaneous',
     ['gist_earth', 'terrain', 'ocean', 'gist_stern',
      'brg', 'CMRmap', 'cubehelix',
      'gnuplot', 'gnuplot2', 'gist_ncar',
      'nipy_spectral', 'jet', 'rainbow',
      'gist_rainbow', 'hsv', 'flag', 'prism']),
    ('EPyGrAM',
     epygram.config.colormaps.keys())]


def plot_cmap(cmap):
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    fig, axes = plt.subplots(figsize=(3, 0.15))
    fig.subplots_adjust(top=0.95, bottom=0.10, left=0.002, right=0.8)
    axes.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap))
    pos = list(axes.get_position().bounds)
    x_text = pos[2] + 0.02
    y_text = pos[1] + pos[3] / 2.
    fig.text(x_text, y_text, cmap, va='center', ha='left', fontsize=9)
    axes.set_axis_off()
    return fig


for (colormapstype, colormaps) in sorted(cmaps):
    for cmap in colormaps:
        if cmap not in plt.colormaps():
            with open(epygram.config.colormaps[cmap], 'r') as ocm:
                add_cmap(cmap, ocm)
        fig = plot_cmap(cmap)
        dirname = 'cmaps'
        filename = '.'.join([colormapstype.replace(' ', '_'), cmap, 'png'])
        fig.savefig(os.path.join(dirname, filename), facecolor='white')
        plt.close(fig)
