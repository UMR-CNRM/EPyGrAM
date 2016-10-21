#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import os
import argparse
import imp

import epygram
from epygram import epylog, epygramError
from epygram.args_catalog import add_arg_to_parser, \
                                 files_management, fields_management, \
                                 misc_options, output_options, \
                                 runtime_options, graphical_options

import matplotlib.pyplot as plt



def main(open_browser=False,
         verbose=False):
    """
    Run the 'epyweb' local server.
    If *open_browser*, open a web browser tab with 'epyweb' interface.
    """

    ep_root = os.getenv('EPYGRAM_INSTALL_DIR', False)
    if True:  #not ep_root:
        _e = epygram.__file__
        ep_root = '/'.join(_e.split('/')[:-2])
    print ep_root
    os.chdir('/'.join([ep_root,
                       'site',
                       'epyweb']))
    print os.getcwd()
    epyweb = imp.load_source('epyweb', 'epyweb.py')
    from epyweb import main as epyweb_main
    epyweb_main(open_browser=open_browser,
                verbose=verbose)




# end of main() ###############################################################



if __name__ == '__main__':

    ### 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description='A web interface for plotting with epygram and vortex.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')
    parser.add_argument('-o',
                        action='store_true',
                        dest='open_browser',
                        help='opens a web browser tab with epyweb interface.',
                        default=False)
    add_arg_to_parser(parser, runtime_options['verbose'])
    args = parser.parse_args()

    ### 2. Initializations
    ######################

    ### 3. Main
    ###########
    main(open_browser=args.open_browser,
         verbose=args.verbose)

###########
### END ###
###########
