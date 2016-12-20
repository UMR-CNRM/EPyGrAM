#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import argparse

import epygram
from epygram.args_catalog import add_arg_to_parser, runtime_options



def main(open_browser=False,
         port=8080,
         vortex_mode=True,
         verbose=False):
    """
    Run the 'epyweb' local server.
    If *open_browser*, open a web browser tab with 'epyweb' interface.
    The *port* to be used for the server can be specified.
    If *vortex_mode*, describe and get resources using Vortex; else,
    as a file system.
    """
    import epyweb

    epyweb.main(open_browser=open_browser,
                port=port,
                vortex_mode=vortex_mode,
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
    parser.add_argument('-p',
                        dest='port',
                        type=int,
                        help='port to be used, default to 8080.',
                        default=8080)
    parser.add_argument('--no_vortex',
                        action='store_false',
                        dest='vortex_mode',
                        help='run the interface without Vortex, selecting resources by their local path.',
                        default=True)
    add_arg_to_parser(parser, runtime_options['verbose'])
    args = parser.parse_args()

    ### 2. Initializations
    ######################

    ### 3. Main
    ###########
    main(open_browser=args.open_browser,
         port=args.port,
         vortex_mode=args.vortex_mode,
         verbose=args.verbose)

###########
### END ###
###########
