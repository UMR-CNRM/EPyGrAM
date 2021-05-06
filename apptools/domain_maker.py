#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
from six.moves import input

import os
import sys
import argparse

# Automatically set the python path
package_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, package_path)
sys.path.insert(0, os.path.join(package_path, 'site'))
import epygram
from epygram import epylog
from epygram.geometries import domain_making as dm

from epygram.args_catalog import (add_arg_to_parser,
                                  domain_maker_options,
                                  runtime_options,
                                  graphical_options)


def main(mode,
         display=True,
         maximize_CI_in_E=False,
         french_depts=False,
         background=True):
    """
    Domain maker.

    :param mode: 'center_dims' to build domain given its center and dimensions;
                 'lonlat_included' to build domain given an included lon/lat area.
    :param display: if False, deactivates the display of domain.
    :param maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                             minimum.
    :param french_depts: draws french departments instead of countries boundaries.
    :param background: if True, set a background color to continents and oceans.
    """

    print("################")
    print("# DOMAIN MAKER #")
    print("################")
    if not epygram.util.mpl_interactive_backend():
        out = 'domain_maker.out.' + epygram.config.default_graphical_output
    else:
        out = None
    if mode == 'center_dims':
        defaults = {'Iwidth':None,
                    'tilting':0.0,
                    'resolution':'',
                    'center_lon':'',
                    'center_lat':'',
                    'Xpoints_CI':'',
                    'Ypoints_CI':''}

        retry = True
        while retry:
            # ask and build
            (geometry, defaults) = dm.ask.ask_and_build_geometry(defaults,
                                                                 maximize_CI_in_E)
            print("Compute domain...")
            print(dm.output.summary(geometry))
            if display:
                plot_lonlat_included = input("Plot a lon/lat domain over model domain ? [n]: ")
                if plot_lonlat_included in ('y', 'Y', 'yes'):
                    plot_lonlat_included = True
                else:
                    plot_lonlat_included = False
                if plot_lonlat_included:
                    proposed = dm.build.compute_lonlat_included(geometry)
                    print("Min/Max longitudes & latitudes of the lon/lat domain (defaults to that proposed above):")
                    ll_boundaries = dm.ask.ask_lonlat(proposed)
                else:
                    ll_boundaries = None
                print("Plot domain...")
                dm.output.plot_geometry(geometry,
                                        lonlat_included=ll_boundaries,
                                        out=out,
                                        departments=french_depts,
                                        background=background,
                                        plotlib='cartopy')
            # retry ?
            retry = input("Do you want to modify something ? [n] ")
            if retry in ('yes', 'y', 'Y'):
                retry = True
            else:
                retry = False

    elif mode == 'lonlat_included':
        defaults = {'Iwidth':None,
                    'lonmax':'',
                    'lonmin':'',
                    'latmax':'',
                    'latmin':'',
                    'resolution':''}

        retry = True
        while retry:
            # ask and build
            (geometry, defaults) = dm.ask.ask_lonlat_and_build_geometry(defaults, maximize_CI_in_E)
            print(dm.output.summary(geometry))
            if display:
                # plot
                print("Plot domain...")
                dm.output.plot_geometry(geometry,
                                        lonlat_included=defaults,
                                        out=out,
                                        departments=french_depts,
                                        background=background,
                                        plotlib='cartopy')
            # retry ?
            retry = input("Do you want to modify something ? [n] ")
            if retry in ('yes', 'y', 'Y'):
                retry = True
            else:
                retry = False
    else:
        raise ValueError("invalid value for 'mode' argument")

    dm.output.write_geometry_as_namelists(geometry, allinone=True)
# end of main() ###############################################################


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    epygramstr = 'EPyGrAM'
    parser = argparse.ArgumentParser(description='An interactive ' + epygramstr + " tool for defining a LAM domain, visualize it, \
                                                  and generate the needed namelist blocks.",
                                     epilog='End of help for: %(prog)s (' + epygramstr + ' v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, domain_maker_options['mode'])
    add_arg_to_parser(parser, domain_maker_options['no_display'])
    add_arg_to_parser(parser, domain_maker_options['maximize_CI_in_E'])
    add_arg_to_parser(parser, runtime_options['verbose'])
    add_arg_to_parser(parser, graphical_options['french_departments'])
    add_arg_to_parser(parser, graphical_options['background'], default=True)

    args = parser.parse_args()

    # 2. Initializations
    ####################
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 3. Main
    #########
    main(args.mode,
         display=not args.no_display,
         maximize_CI_in_E=args.maximize_CI_in_E,
         french_depts=args.depts,
         background=args.background)
