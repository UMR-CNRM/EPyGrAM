#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, unicode_literals

import argparse
import os
import webbrowser

import epygram

if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description='A command to open the local epygram documentation.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')
    args = parser.parse_args()

    # 2. Initializations
    ####################
    url = os.path.join(os.path.dirname(epygram.__file__),
                       os.path.join('doc_sphinx', 'html', 'index.html'))
    print('Open url file://' + url )

    # 3. Main
    #########
    webbrowser.open(url)
