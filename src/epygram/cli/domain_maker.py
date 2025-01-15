#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An interactive EPyGrAM tool for defining a LAM domain, visualize it, and generate the needed namelist blocks."""

import argparse

import epygram
from epygram import epylog as logger
from epygram.geometries import domain_making as dm
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           domain_maker_args,
                           runtime_args,
                           graphical_args)

_description = __doc__


def main():
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    domain_maker(args.mode,
                 display=not args.no_display,
                 maximize_CI_in_E=args.maximize_CI_in_E,
                 french_depts=args.depts,
                 background=args.background,
                 truncation=args.truncation,
                 orography_subtruncation=args.orography_subtruncation)


def domain_maker(mode,
                 display=True,
                 maximize_CI_in_E=False,
                 french_depts=False,
                 background=True,
                 # for namelists only
                 truncation='linear',
                 orography_subtruncation='quadratic'):
    """
    Domain maker.

    :param mode: 'center_dims' to build domain given its center and dimensions;
                 'lonlat_included' to build domain given an included lon/lat area.
    :param display: if False, deactivates the display of domain.
    :param maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                             minimum.
    :param french_depts: draws french departments instead of countries boundaries.
    :param background: if True, set a background color to continents and oceans.
    :param truncation: the kind of truncation of spectral geometry
                       to generate, among ('linear', 'quadratic', 'cubic').
    :param orography_subtruncation: additional subtruncation for orography to
                                    be generated.
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

    dm.output.write_geometry_as_namelists(geometry,
                                          allinone=True,
                                          truncation=truncation,
                                          orography_subtruncation=orography_subtruncation)


def get_args():

    parser = argparse.ArgumentParser(description=_description, epilog=epilog)

    add_arg_to_parser(parser, domain_maker_args['mode'])
    add_arg_to_parser(parser, domain_maker_args['no_display'])
    add_arg_to_parser(parser, domain_maker_args['maximize_CI_in_E'])
    add_arg_to_parser(parser, domain_maker_args['truncation'])
    add_arg_to_parser(parser, domain_maker_args['orography_subtruncation'])
    add_arg_to_parser(parser, runtime_args['verbose'])
    add_arg_to_parser(parser, graphical_args['french_departments'])
    add_arg_to_parser(parser, graphical_args['background'], default=True)

    return parser.parse_args()

