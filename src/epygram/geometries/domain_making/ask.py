#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains functions for interactive build of a LAM domain.
"""

from .util import default_Iwidth
from .build import build_geometry, build_geometry_fromlonlat


def ask_and_build_geometry(defaults,
                           maximize_CI_in_E=False):
    """
    Ask the user for geometry params,
    then builds a proposal geometry.

    :param defaults: a dict() containing default values.
    :param maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                             minimum.
    """
    # ask for geometry
    try:
        resolution = float(input("Resolution in m [" + str(defaults['resolution']) + "]: "))
    except ValueError:
        if defaults['resolution'] != '':
            resolution = defaults['resolution']
        else:
            raise ValueError("Invalid resolution.")
    try:
        center_lon = float(input("Center of domain / longitude in degrees [" + str(defaults['center_lon']) + "]: "))
    except ValueError:
        if defaults['center_lon'] != '':
            center_lon = defaults['center_lon']
        else:
            raise ValueError("Invalid longitude.")
    try:
        center_lat = float(input("Center of domain / latitude in degrees [" + str(defaults['center_lat']) + "]: "))
    except ValueError:
        if defaults['center_lat'] != '':
            center_lat = defaults['center_lat']
        else:
            raise ValueError("Invalid latitude.")
    try:
        tilting = float(input("Optional counterclockwise tilting in degrees (lon0-lonC) [" + str(defaults['tilting']) + "]: "))
    except ValueError:
        tilting = defaults['tilting']

    # and dimensions
    try:
        Xpoints_CI = int(input("C+I zonal (X) dimension in pts [" + str(defaults['Xpoints_CI']) + "]: "))
    except ValueError:
        if defaults['Xpoints_CI'] != '':
            Xpoints_CI = defaults['Xpoints_CI']
        else:
            raise ValueError("Invalid dimension.")
    try:
        Ypoints_CI = int(input("C+I meridian (Y) dimension in pts [" + str(defaults['Ypoints_CI']) + "]: "))
    except ValueError:
        if defaults['Ypoints_CI'] != '':
            Ypoints_CI = defaults['Ypoints_CI']
        else:
            raise ValueError("Invalid dimension.")
    if defaults['Iwidth'] is None:
        defaults['Iwidth'] = default_Iwidth(resolution)
    try:
        Iwidth = int(input("I-zone width in pts [" + str(defaults['Iwidth']) + "]: "))
    except Exception:
        Iwidth = defaults['Iwidth']

    geometry = build_geometry(center_lon, center_lat,
                              Xpoints_CI, Ypoints_CI,
                              resolution,
                              Iwidth=Iwidth,
                              tilting=tilting,
                              maximize_CI_in_E=maximize_CI_in_E,
                              interactive=True)

    print("--------------------------------------------------")
    defaults = {'Iwidth':geometry.dimensions['X_Iwidth'],
                'tilting':geometry.projection['reference_lon'].get('degrees') -
                          geometry.grid['input_lon'].get('degrees'),
                'resolution':geometry.grid['X_resolution'],
                'center_lon':geometry.grid['input_lon'].get('degrees'),
                'center_lat':geometry.grid['input_lat'].get('degrees'),
                'Xpoints_CI':geometry.dimensions['X_CIzone'],
                'Ypoints_CI':geometry.dimensions['Y_CIzone']}

    return geometry, defaults


def ask_lonlat(defaults):
    """Ask a lon/lat geometry."""
    try:
        lonmin = float(input("Minimum (Western) longitude in degrees [" + str(defaults['lonmin']) + "]: "))
    except ValueError:
        if str(defaults['lonmin']) != '':
            lonmin = defaults['lonmin']
        else:
            raise ValueError("Invalid longitude.")
    try:
        lonmax = float(input("Maximum (Eastern) longitude in degrees [" + str(defaults['lonmax']) + "]: "))
    except ValueError:
        if str(defaults['lonmax']) != '':
            lonmax = defaults['lonmax']
        else:
            raise ValueError("Invalid longitude.")
    try:
        latmin = float(input("Minimum (Southern) latitude in degrees [" + str(defaults['latmin']) + "]: "))
    except ValueError:
        if str(defaults['latmin']) != '':
            latmin = defaults['latmin']
        else:
            raise ValueError("Invalid latitude.")
    try:
        latmax = float(input("Maximum (Northern) latitude in degrees [" + str(defaults['latmax']) + "]: "))
    except ValueError:
        if str(defaults['latmax']) != '':
            latmax = defaults['latmax']
        else:
            raise ValueError("Invalid latitude.")

    return {'lonmin':lonmin, 'lonmax':lonmax, 'latmin':latmin, 'latmax':latmax}


def ask_lonlat_and_build_geometry(defaults,
                                  maximize_CI_in_E=False):
    """
    Ask the user for lonlat-included geometry params,
    then builds a proposal geometry.

    :param defaults: a dict() containing default values.
    :param maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                             minimum.
    """
    # ask for geometry
    try:
        resolution = float(input("Model resolution in m [" + str(defaults['resolution']) + "]: "))
    except ValueError:
        if defaults['resolution'] != '':
            resolution = defaults['resolution']
        else:
            raise ValueError("Invalid resolution.")
    print("Min/Max longitudes & latitudes that must be included in domain:")
    ll_boundaries = ask_lonlat(defaults)
    lonmin = ll_boundaries['lonmin']
    lonmax = ll_boundaries['lonmax']
    latmin = ll_boundaries['latmin']
    latmax = ll_boundaries['latmax']
    if defaults['Iwidth'] is None:
        defaults['Iwidth'] = default_Iwidth(resolution)
    try:
        Iwidth = int(input("I-zone width in pts [" + str(defaults['Iwidth']) + "]: "))
    except Exception:
        Iwidth = defaults['Iwidth']

    geometry = build_geometry_fromlonlat(lonmin, lonmax,
                                         latmin, latmax,
                                         resolution,
                                         Iwidth=Iwidth,
                                         force_projection=None,
                                         maximize_CI_in_E=maximize_CI_in_E,
                                         interactive=True)

    print("--------------------------------------------------")
    defaults = {'Iwidth':Iwidth,
                'lonmax':lonmax,
                'lonmin':lonmin,
                'latmax':latmax,
                'latmin':latmin,
                'resolution':resolution}

    return geometry, defaults
