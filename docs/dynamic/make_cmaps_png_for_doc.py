#!/usr/bin/env python3
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
import string
import argparse

import epygram
from epygram.colormapping import register_colormap_from_json
import matplotlib.pyplot as plt

doc_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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


def main(output_rootdir=None):
    if output_rootdir is None:
        output_rootdir = os.path.join(doc_root, 'build', 'html')
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
        subdir = os.path.join('_images', 'colormaps')
        dirname = os.path.join(output_rootdir, subdir)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        filename = cmap + '.png'
        fig.savefig(os.path.join(dirname, filename), facecolor='0.66')
        plt.close(fig)
        # append RST
        figures += '.. figure:: ' + os.path.join('', subdir, filename) + '\n'
    rstout = rst.substitute(epy_cmaps=figures)

    # write RST
    with io.open(os.path.join(doc_root, 'source', 'cmaps.rst',), 'w') as out:
        out.write(rstout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Make EPyGrAM colormaps illustrations for documentation purpose.")
    parser.add_argument('-d', '--output_rootdir',
                        help="Root directory for the output doc.",
                        default=None)
    args = parser.parse_args()
    main(output_rootdir=args.output_rootdir)
