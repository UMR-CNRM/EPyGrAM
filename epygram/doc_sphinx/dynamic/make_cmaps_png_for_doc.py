#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Make EPyGrAM colormaps illustrations for documentation purpose.
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import io
import os
import numpy as np
import matplotlib.pyplot as plt
import string

from bronx.graphics.colormapping import register_colormap_from_json

import epygram

here = os.path.dirname(os.path.abspath(__file__))
doc_root = os.path.dirname(here)


def plot_cmap(cmap):
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    fig, axes = plt.subplots(figsize=(5, 0.2))
    fig.subplots_adjust(top=0.95, bottom=0.05, left=0.2, right=0.996)
    axes.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap))
    pos = list(axes.get_position().bounds)
    x_text = pos[0] - 0.01
    y_text = pos[1] + pos[3] / 2.
    fig.text(x_text, y_text, cmap, va='center', ha='right', fontsize=9)
    axes.set_axis_off()
    return fig


# read template
template_rst = os.path.join(os.path.dirname(__file__), 'cmaps.rst.tpl')
with open(template_rst, 'r') as t:
    rst = string.Template(''.join(t.readlines()))

# plot colormaps
figures = ''
for cmap in sorted(epygram.config.colormaps.keys()):
    if cmap not in plt.colormaps():
        register_colormap_from_json(epygram.config.colormaps[cmap])
    fig = plot_cmap(cmap)
    path_a = os.path.join(doc_root, 'html')
    path_b = os.path.join('_images', 'colormaps')
    dirname = os.path.join(path_a, path_b)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    filename = cmap + '.png'
    fig.savefig(os.path.join(dirname, filename), facecolor='0.66')
    plt.close(fig)
    # append RST
    figures += '.. figure:: ' + os.path.join('', path_b, filename) + '\n'
rstout = rst.substitute(epy_cmaps=figures)

# write RST
with io.open(os.path.join(doc_root, 'source', 'cmaps.rst',), 'w') as out:
    out.write(rstout)
