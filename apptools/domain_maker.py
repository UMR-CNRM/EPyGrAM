#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import argparse
import numpy
import math

import footprints

import epygram
from epygram import epylog

import matplotlib.pyplot as plt
from epygram.args_catalog import add_arg_to_parser, \
                                 domain_maker_options, runtime_options, \
                                 graphical_options


"""
Note: I tried to compute a lon/lat optimal post-processing subdomain automatically from
a given domain; I failed. Is it possible to do so, considering that the domain
can be tilted ? Open question...
"""

# parameters
Ezone_minimum_width = 11
threshold_mercator_lambert = 1.0
maxdims_security_barrier = 10000
vkw = {'structure': 'V',
       'typeoffirstfixedsurface': 1,
       'levels': [1]}
vgeom = footprints.proxy.geometry(**vkw)



def nearest_greater_FFTcompliant_int(guess):
    """
    Returns the first integer *n* greater than *guess* that satisfies
    n = 2^(1+i) x 3^j x 5^k, with (i,j,k) being positive integers. 
    """

    # these are defined for dimensions up to 10000 points at least
    M2 = 14
    M3 = 10
    M5 = 7
    fft_compliant_dims = numpy.zeros((M2, M3, M5))
    for i in range(M2):
        for j in range(M3):
            for k in range(M5):
                fft_compliant_dims[i, j, k] = 2 ** (1 + i) * 3 ** j * 5 ** k
    fft_compliant_dims = sorted(list(set(fft_compliant_dims.flatten())))
    for i in range(len(fft_compliant_dims)):
        if fft_compliant_dims[i] >= guess:
            result = int(fft_compliant_dims[i])
            break
    return result



def ask_and_build_geometry(defaults,
                           maximize_CI_in_E=False):
    """
    Ask the user for geometry params,
    then builds a proposal geometry.
    
    Args:
        defaults: a dict() containing default values.
        maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                          minimum.
    """

    projections = {'L':'lambert', 'M':'mercator', 'PS':'polar_stereographic'}
    print_projections = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}

    # ask for geometry
    try:
        resolution = float(raw_input("Resolution in m [" + str(defaults['resolution']) + "]: "))
    except ValueError:
        if defaults['resolution'] != '':
            resolution = defaults['resolution']
        else:
            raise ValueError("Invalid resolution.")
    try:
        center_lon = float(raw_input("Center of domain / longitude in degrees [" + str(defaults['center_lon']) + "]: "))
    except ValueError:
        if defaults['center_lon'] != '':
            center_lon = defaults['center_lon']
        else:
            raise ValueError("Invalid longitude.")
    try: center_lat = float(raw_input("Center of domain / latitude in degrees [" + str(defaults['center_lat']) + "]: "))
    except ValueError:
        if defaults['center_lat'] != '':
            center_lat = defaults['center_lat']
        else:
            raise ValueError("Invalid latitude.")
    try:
        tilting = float(raw_input("Optional counterclockwise tilting in degrees [" + str(defaults['tilting']) + "]: "))
    except ValueError:
        tilting = defaults['tilting']

    # and dimensions
    try:
        Xpoints_CI = int(raw_input("C+I zonal (X) dimension in pts [" + str(defaults['Xpoints_CI']) + "]: "))
    except ValueError:
        if defaults['Xpoints_CI'] != '':
            Xpoints_CI = defaults['Xpoints_CI']
        else:
            raise ValueError("Invalid dimension.")
    try:
        Ypoints_CI = int(raw_input("C+I meridian (Y) dimension in pts [" + str(defaults['Ypoints_CI']) + "]: "))
    except ValueError:
        if defaults['Ypoints_CI'] != '':
            Ypoints_CI = defaults['Ypoints_CI']
        else:
            raise ValueError("Invalid dimension.")
    try:
        Iwidth = int(raw_input("I-zone width in pts [" + str(defaults['Iwidth']) + "]: "))
    except Exception:
        Iwidth = defaults['Iwidth']

    # begin to build a horizontal geometry
    Xpoints_CIE = nearest_greater_FFTcompliant_int(Xpoints_CI + Ezone_minimum_width)
    Ypoints_CIE = nearest_greater_FFTcompliant_int(Ypoints_CI + Ezone_minimum_width)
    if maximize_CI_in_E:
            Xpoints_CI = Xpoints_CIE - Ezone_minimum_width
            Ypoints_CI = Ypoints_CIE - Ezone_minimum_width
    # dimensions
    dimensions = {'X':Xpoints_CIE,
                  'Y':Ypoints_CIE,
                  'X_CIzone':Xpoints_CI,
                  'Y_CIzone':Ypoints_CI,
                  'X_Czone':Xpoints_CI - 2 * Iwidth,
                  'Y_Czone':Ypoints_CI - 2 * Iwidth,
                  'X_CIoffset':0,
                  'Y_CIoffset':0,
                  'X_Iwidth':Iwidth,
                  'Y_Iwidth':Iwidth
                      }
    # coordinates
    projection = {'reference_lon':epygram.util.Angle(center_lon + tilting, 'degrees'),
                  'reference_lat':epygram.util.Angle(center_lat, 'degrees'),
                  'rotation':epygram.util.Angle(0.0, 'degrees'),
                  }
    # grid
    grid = {'X_resolution':resolution,
            'Y_resolution':resolution,
            'LAMzone':'CIE',
            'input_lon':epygram.util.Angle(center_lon, 'degrees'),
            'input_lat':epygram.util.Angle(center_lat, 'degrees'),
            'input_position':(float(dimensions['X_CIzone'] - 1) / 2.,
                              float(dimensions['Y_CIzone'] - 1) / 2.)
            }

    # try to guess best projection
    if abs(center_lat) >= threshold_mercator_lambert:  # lambert or polar stereographic
        # we make a "first guess" with Lambert Geometry, just to check the pole is not inside the domain
        geometryname = 'lambert'
        geometry = footprints.proxy.geometry(structure='H2D',
                                             name=geometryname,
                                             grid=footprints.FPDict(grid),
                                             dimensions=footprints.FPDict(dimensions),
                                             projection=footprints.FPDict(projection),
                                             vcoordinate=vgeom,
                                             position_on_horizontal_grid='center')

        # test if pole in lambert grid
        pole_in_domain = geometry.point_is_inside_domain_ll(0, 90) or geometry.point_is_inside_domain_ll(0, -90)
        if pole_in_domain:
            projname = 'PS'
        else:
            projname = 'L'
    else:
        projname = 'M'
        pole_in_domain = False

    # proposal
    if not pole_in_domain:
        print "Advised projection for this center and these dimensions is:"
        print "=>", print_projections[projname]
        print "Other choices:", print_projections
        accepted_projection = raw_input("Chosen projection [" + projname + "]: ")
        if accepted_projection == '':
            accepted_projection = projname
    else:
        print "Pole is inside domain: only available projection for this center and these dimensions is:"
        print print_projections[projname]
        accepted_projection = projname

    if accepted_projection == 'PS':
        reference_lat = math.copysign(90.0, center_lat)
    elif accepted_projection == 'M':
        reference_lat = 0.0
        projection['reference_lon'] = epygram.util.Angle(center_lon, 'degrees')
        if tilting != 0.0:
            epylog.warning("! Tilting ignored: not available for Mercator projection.")
    elif accepted_projection == 'L':
        print "Advised reference latitude for Lambert domain is center latitude:"
        accepted_lat = raw_input("Reference latitude [" + str(center_lat) + "]: ")
        if accepted_lat == '':
            reference_lat = center_lat
        else:
            reference_lat = float(accepted_lat)
    projection['reference_lat'] = epygram.util.Angle(reference_lat, 'degrees')
    geometryname = projections[accepted_projection]
    geometry = footprints.proxy.geometry(structure='H2D',
                                         name=geometryname,
                                         grid=footprints.FPDict(grid),
                                         dimensions=footprints.FPDict(dimensions),
                                         projection=footprints.FPDict(projection),
                                         vcoordinate=vgeom,
                                         position_on_horizontal_grid='center',)

    print "--------------------------------------------------"
    defaults = {'Iwidth':Iwidth,
                 'tilting':tilting,
                 'resolution':resolution,
                 'center_lon':center_lon,
                 'center_lat':center_lat,
                 'Xpoints_CI':Xpoints_CI,
                 'Ypoints_CI':Ypoints_CI}
    return geometry, defaults



def ask_lonlat_and_build_geometry(defaults,
                                  maximize_CI_in_E=False):
    """
    Ask the user for lonlat-included geometry params,
    then builds a proposal geometry.
    
    Args:
        defaults: a dict() containing default values.
        maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                          minimum.
    """

    projections = {'L':'lambert', 'M':'mercator', 'PS':'polar_stereographic'}
    print_projections = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}

    # ask for geometry
    try:
        resolution = float(raw_input("Model resolution in m [" + str(defaults['resolution']) + "]: "))
    except ValueError:
        if defaults['resolution'] != '':
            resolution = defaults['resolution']
        else:
            raise ValueError("Invalid resolution.")
    print "Min/Max longitudes & latitudes that must be included in domain:"
    ll_boundaries = ask_lonlat(defaults)
    lonmin = ll_boundaries['lonmin']
    lonmax = ll_boundaries['lonmax']
    latmin = ll_boundaries['latmin']
    latmax = ll_boundaries['latmax']
    try:
        Iwidth = int(raw_input("I-zone width in pts [" + str(defaults['Iwidth']) + "]: "))
    except Exception:
        Iwidth = defaults['Iwidth']

    # begin to build a horizontal geometry
    if lonmin > lonmax:
        lonmax += 360.
    center_lon = (lonmax + lonmin) / 2.
    if center_lon > 180.:
        lonmin -= 360.
        lonmax -= 360.
    elif center_lon < -180.:
        lonmin += 360.
        lonmax += 360.
    center_lon = (lonmax + lonmin) / 2.
    center_lat = (latmax + latmin) / 2.
    if abs(center_lat) >= threshold_mercator_lambert:
        projname = 'L'
    else:
        projname = 'M'

    # proposal
    print "Advised projection with regards to domain center latitude is:"
    print "=>", print_projections[projname]
    print "Other choices:", print_projections
    accepted_projection = raw_input("Chosen projection [" + projname + "]: ")
    if accepted_projection == '':
        accepted_projection = projname
    if accepted_projection == 'PS':
        reference_lat = math.copysign(90.0, center_lat)
    elif accepted_projection == 'M':
        reference_lat = 0.0
    elif accepted_projection == 'L':
        print "Advised center latitude for Lambert domain is mean(Northern, Southern):"
        accepted_lat = raw_input("Center latitude [" + str(center_lat) + "]: ")
        if accepted_lat != '':
            center_lat = float(accepted_lat)
        print "Advised reference latitude for Lambert domain is center latitude:"
        accepted_lat = raw_input("Reference latitude [" + str(center_lat) + "]: ")
        if accepted_lat != '':
            reference_lat = float(accepted_lat)
        else:
            reference_lat = center_lat

    Xpoints_CI = 2 * Iwidth + 1
    Ypoints_CI = 2 * Iwidth + 1
    lonlat_included = False
    while not lonlat_included:
        Xpoints_CIE = nearest_greater_FFTcompliant_int(Xpoints_CI + Ezone_minimum_width)
        Ypoints_CIE = nearest_greater_FFTcompliant_int(Ypoints_CI + Ezone_minimum_width)
        if maximize_CI_in_E:
            Xpoints_CI = Xpoints_CIE - Ezone_minimum_width
            Ypoints_CI = Ypoints_CIE - Ezone_minimum_width
        # dimensions
        dimensions = {'X':Xpoints_CIE,
                      'Y':Ypoints_CIE,
                      'X_CIzone':Xpoints_CI,
                      'Y_CIzone':Ypoints_CI,
                      'X_Czone':Xpoints_CI - 2 * Iwidth,
                      'Y_Czone':Ypoints_CI - 2 * Iwidth,
                      'X_CIoffset':0,
                      'Y_CIoffset':0,
                      'X_Iwidth':Iwidth,
                      'Y_Iwidth':Iwidth
                      }
        # coordinates
        projection = {'reference_lon':epygram.util.Angle(center_lon, 'degrees'),
                      'reference_lat':epygram.util.Angle(reference_lat, 'degrees'),
                      'rotation':epygram.util.Angle(0.0, 'degrees'),
                      }
        # grid
        grid = {'X_resolution':resolution,
                'Y_resolution':resolution,
                'LAMzone':'CIE',
                'input_lon':epygram.util.Angle(center_lon, 'degrees'),
                'input_lat':epygram.util.Angle(center_lat, 'degrees'),
                'input_position':(float(dimensions['X_CIzone'] - 1) / 2.,
                                  float(dimensions['Y_CIzone'] - 1) / 2.)
                }
        geometryname = projections[accepted_projection]
        # first guess for lambert, to check that pole is not in domain
        if geometryname == 'lambert':
            geometry = footprints.proxy.geometry(structure='H2D',
                                                 name=geometryname,
                                                 grid=footprints.FPDict(grid),
                                                 dimensions=footprints.FPDict(dimensions),
                                                 projection=footprints.FPDict(projection),
                                                 vcoordinate=vgeom,
                                                 position_on_horizontal_grid='center')
            pole_in_domain = (geometry.point_is_inside_domain_ll(0, 90) or geometry.point_is_inside_domain_ll(0, -90))
            if pole_in_domain:
                epylog.warning("Pole is inside Lambert domain => shifted to Polar Stereographic projection !")
                geometryname = 'polar_stereographic'
                projection['reference_lat'] = epygram.util.Angle(math.copysign(90.0, center_lat), 'degrees')

        # guess
        geometry = footprints.proxy.geometry(structure='H2D',
                                             name=geometryname,
                                             grid=footprints.FPDict(grid),
                                             dimensions=footprints.FPDict(dimensions),
                                             projection=footprints.FPDict(projection),
                                             vcoordinate=vgeom,
                                             position_on_horizontal_grid='center')
        # test whether the lonlat corners are inside domain
        points_to_test = [(lonmax, latmax), (lonmax, latmin),
                          (lonmin, latmax), (lonmin, latmin),
                          (lonmax, (latmin + latmax) / 2.), (lonmin, (latmin + latmax) / 2.),
                          ((lonmin + lonmax) / 2., latmax), ((lonmin + lonmax) / 2., latmin)]
        xlonlat_included = all([1.0 + dimensions['X_Iwidth'] < geometry.ll2ij(*c)[0] < dimensions['X_Czone'] - 1.0 \
                                for c in points_to_test])
        ylonlat_included = all([1.0 + dimensions['Y_Iwidth'] < geometry.ll2ij(*c)[1] < dimensions['Y_Czone'] - 1.0 \
                                for c in points_to_test])
        lonlat_included = xlonlat_included and ylonlat_included
        if not lonlat_included:
            if not xlonlat_included:
                Xpoints_CI = Xpoints_CIE - Ezone_minimum_width + 1
            if not ylonlat_included:
                Ypoints_CI = Ypoints_CIE - Ezone_minimum_width + 1
            if Xpoints_CI > maxdims_security_barrier \
            or Ypoints_CI > maxdims_security_barrier:
                raise epygram.epygramError("Domain is too large, > " + str(maxdims_security_barrier) + " points.")

    print "--------------------------------------------------"
    defaults = {'Iwidth':Iwidth,
                'lonmax':lonmax,
                'lonmin':lonmin,
                'latmax':latmax,
                'latmin':latmin,
                'resolution':resolution}

    return geometry, defaults



def show_geometry(geometry):
    """
    Returns a summary of geometry as a character string.
    
    Args:
        geometry: a H2DGeometry instance
    """

    invprojections = {'lambert':'L', 'mercator':'M', 'polar_stereographic':'PS'}
    print_projections = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}

    disp = ""
    disp += "# Geometry Summary #" + '\n'
    disp += "  ================  " + '\n'
    disp += "Center Longitude: " + str(geometry.getcenter()[0].get('degrees')) + '\n'
    disp += "Center Latitude:  " + str(geometry.getcenter()[1].get('degrees')) + '\n'
    disp += "Tilting:          " + \
        str(geometry.projection['reference_lon'].get('degrees') - geometry.getcenter()[0].get('degrees')) + '\n'
    disp += "  => Reference longitude: " + str(geometry.projection['reference_lon'].get('degrees')) + '\n'
    disp += "Projection: " + print_projections[invprojections[geometry.name]] + '\n'
    disp += "  Reference latitude: " + str(geometry.projection['reference_lat'].get('degrees')) + '\n'
    disp += "Resolution: " + str(geometry.grid['X_resolution']) + '\n'
    mapfactor = geometry.map_factor_field().getdata(subzone='CI')
    mapfactor_range = [mapfactor.min(), mapfactor.max()]
    disp += "Map factor range on C+I domain: [" + \
          '{:.{precision}{type}}'.format(mapfactor_range[0], type='E', precision=3) + " -> " + \
          '{:.{precision}{type}}'.format(mapfactor_range[1], type='E', precision=3) + "]" + '\n'
    disp += "---" + '\n'
    disp += "Dimensions    " + '{:^{width}}'.format("C+I", width=10) + \
          '{:^{width}}'.format("C+I+E", width=10) + \
          '{:^{width}}'.format("E-zone", width=10) + '\n'
    Ezone_Xwidth = geometry.dimensions['X'] - geometry.dimensions['X_CIzone']
    Ezone_Ywidth = geometry.dimensions['Y'] - geometry.dimensions['Y_CIzone']
    if Ezone_Xwidth == Ezone_minimum_width:
        disp += "X:            " + '{:^{width}}'.format(str(geometry.dimensions['X_CIzone']), width=10) + \
              '{:^{width}}'.format(str(geometry.dimensions['X']), width=10) + \
              '{:^{width}}'.format(str(Ezone_Xwidth), width=10) + \
              ' (optimal)' + '\n'
    else:
        disp += "X:            " + '{:^{width}}'.format(str(geometry.dimensions['X_CIzone']), width=10) + \
              '{:^{width}}'.format(str(geometry.dimensions['X']), width=10) + \
              '{:^{width}}'.format(str(Ezone_Xwidth), width=10) + \
              '/ ' + str(Ezone_minimum_width) + ' (optimal) => advised C+I zonal (X) dimension =' + \
              '{:^{width}}'.format(str(geometry.dimensions['X'] - Ezone_minimum_width), width=10) + '\n'
    if Ezone_Ywidth == Ezone_minimum_width:
        disp += "Y:            " + '{:^{width}}'.format(str(geometry.dimensions['Y_CIzone']), width=10) + \
              '{:^{width}}'.format(str(geometry.dimensions['Y']), width=10) + \
              '{:^{width}}'.format(str(Ezone_Ywidth), width=10) + \
              ' (optimal)' + '\n'
    else:
        disp += "Y:            " + '{:^{width}}'.format(str(geometry.dimensions['Y_CIzone']), width=10) + \
              '{:^{width}}'.format(str(geometry.dimensions['Y']), width=10) + \
              '{:^{width}}'.format(str(Ezone_Ywidth), width=10) + \
              '/ ' + str(Ezone_minimum_width) + ' (optimal) => advised C+I meridian (Y) dimension =' + \
              '{:^{width}}'.format(str(geometry.dimensions['Y'] - Ezone_minimum_width), width=10) + '\n'
    disp += "I zone width: " + str(geometry.dimensions['X_Iwidth']) + '\n'
    if geometry.projection['reference_lon'].get('degrees') - geometry.getcenter()[0].get('degrees') <= epygram.config.epsilon:
        disp += "---" + '\n'
        disp += "The domain contains (at least) the following lon/lat regular area:" + '\n'
        ll_included = compute_lonlat_included(geometry)
        disp += "Longitudes: " + '{:.{precision}{type}}'.format(ll_included['lonmin'], type='F', precision=4) + \
                " <--> " + '{:.{precision}{type}}'.format(ll_included['lonmax'], type='F', precision=4) + '\n'
        disp += "Latitudes:  " + '{:.{precision}{type}}'.format(ll_included['latmax'], type='F', precision=4) + '\n'
        disp += "            " + "  ^" + '\n'
        disp += "            " + "  |" + '\n'
        disp += "            " + "  v" + '\n'
        disp += "            " + '{:.{precision}{type}}'.format(ll_included['latmin'], type='F', precision=4) + '\n'
    disp += "--------------------------------------------------"

    return disp

def compute_lonlat_included(geometry):
    """Computes a lon/lat domain included in the C zone of the model domain."""

    (longrid, latgrid) = geometry.get_lonlat_grid(subzone='C')
    lonmin = max(longrid[:, 1])
    lonmax = min(longrid[:, -2])
    latmax = min(latgrid[-2, :])
    latmin = max(latgrid[1, :])

    return {'lonmin':lonmin, 'lonmax':lonmax, 'latmin':latmin, 'latmax':latmax}

def ask_lonlat(defaults):
    """Ask a lon/lat geometry."""

    try:
        lonmin = float(raw_input("Minimum (Western) longitude in degrees [" + str(defaults['lonmin']) + "]: "))
    except ValueError:
        if defaults['lonmin'] != '':
            lonmin = defaults['lonmin']
        else:
            raise ValueError("Invalid longitude.")
    try:
        lonmax = float(raw_input("Maximum (Eastern) longitude in degrees [" + str(defaults['lonmax']) + "]: "))
    except ValueError:
        if defaults['lonmax'] != '':
            lonmax = defaults['lonmax']
        else:
            raise ValueError("Invalid longitude.")
    try:
        latmin = float(raw_input("Minimum (Southern) latitude in degrees [" + str(defaults['latmin']) + "]: "))
    except ValueError:
        if defaults['latmin'] != '':
            latmin = defaults['latmin']
        else:
            raise ValueError("Invalid latitude.")
    try:
        latmax = float(raw_input("Maximum (Northern) latitude in degrees [" + str(defaults['latmax']) + "]: "))
    except ValueError:
        if defaults['latmax'] != '':
            latmax = defaults['latmax']
        else:
            raise ValueError("Invalid latitude.")

    return {'lonmin':lonmin, 'lonmax':lonmax, 'latmin':latmin, 'latmax':latmax}

def build_lonlat_field(ll_boundaries, fid={'lon/lat':'template'}):
    """
    Build a lonlat field empty except on the border, given lon/lat boundaries.
    """

    llwidth = 1000
    llgrid = {'input_lon':epygram.util.Angle(ll_boundaries['lonmin'], 'degrees'),
              'input_lat':epygram.util.Angle(ll_boundaries['latmin'], 'degrees'),
              'input_position':(0, 0),
              'X_resolution':epygram.util.Angle((ll_boundaries['lonmax'] - ll_boundaries['lonmin']) / llwidth, 'degrees'),
              'Y_resolution':epygram.util.Angle((ll_boundaries['latmax'] - ll_boundaries['latmin']) / llwidth, 'degrees')
              }
    lldims = {'X':llwidth + 1, 'Y':llwidth + 1}
    llgeometry = footprints.proxy.geometry(structure='H2D',
                                           name='regular_lonlat',
                                           grid=footprints.FPDict(llgrid),
                                           dimensions=footprints.FPDict(lldims),
                                           vcoordinate=vgeom,
                                           position_on_horizontal_grid='center')
    lldomain = footprints.proxy.field(structure='H2D',
                                      geometry=llgeometry,
                                      fid=footprints.FPDict(fid))
    data = numpy.zeros((lldims['Y'], lldims['X']))
    data[1:-1, 1:-1] = 1.
    lldomain.setdata(data)

    return lldomain



### MAIN #######################################################################

def main(mode,
         no_display=False,
         maximize_CI_in_E=False,
         gisquality='i'):
    """
    Args:
        mode: 'center_dims' to build domain given its center and dimensions;
              'lonlat_included' to build domain given an included lon/lat area.
        no_display: if True, activates the display of domain.
        maximize_CI_in_E: boolean deciding to force the E-zone to be at its
                          minimum.
        gisquality: quality of coastlines and countries boundaries.
    """

    print "################"
    print "# DOMAIN MAKER #"
    print "################"

    if mode == 'center_dims':
        defaults = {'Iwidth':8,
                    'tilting':0.0,
                    'resolution':'',
                    'center_lon':'',
                    'center_lat':'',
                    'Xpoints_CI':'',
                    'Ypoints_CI':''}

        retry = True
        while retry:
            # ask and build
            (geometry, defaults) = ask_and_build_geometry(defaults, maximize_CI_in_E)
            print "Compute domain..."
            print show_geometry(geometry)
            if not no_display:
                # plot
                CIEdomain = footprints.proxy.field(structure='H2D',
                                                   geometry=geometry,
                                                   fid=footprints.FPDict({'zone':'C+I+E'}))
                data = numpy.ones((geometry.dimensions['Y'], geometry.dimensions['X'])) * 2.0
                data[0:geometry.dimensions['Y_CIzone'], 0:geometry.dimensions['X_CIzone']] = 1.0
                data[geometry.dimensions['Y_Iwidth']:geometry.dimensions['Y_CIzone'] - geometry.dimensions['Y_Iwidth'],
                     geometry.dimensions['X_Iwidth']:geometry.dimensions['X_CIzone'] - geometry.dimensions['X_Iwidth']] = 0.0
                CIEdomain.setdata(data)
                print "Plot domain..."
                bm = CIEdomain.geometry.make_basemap(specificproj=('nsper', {'sat_height':5000}))
                fig = CIEdomain.plotfield(existingbasemap=bm,
                                          levelsnumber=6,
                                          minmax=[-1.0, 3.0],
                                          colorbar=False,
                                          title='Domain: C+I+E',
                                          meridians=None,
                                          parallels=None,
                                          gisquality=gisquality,
                                          bluemarble=args.bluemarble)
                plot_lonlat_included = raw_input("Plot a lon/lat domain over model domain ? [n]: ")
                if plot_lonlat_included in ('y', 'Y', 'yes'):
                    plot_lonlat_included = True
                else:
                    plot_lonlat_included = False
                if plot_lonlat_included:
                    proposed = compute_lonlat_included(geometry)
                    print "Min/Max longitudes & latitudes of the lon/lat domain (defaults to that proposed above):"
                    ll_boundaries = ask_lonlat(proposed)
                    ll_domain = build_lonlat_field(ll_boundaries)
                    ll_domain.plotfield(existingfigure=fig,
                                        existingbasemap=bm,
                                        graphicmode='contourlines',
                                        title='Domain: C+I+E \n Red contour: required lon/lat',
                                        levelsnumber=2,
                                        contourcolor='red',
                                        contourwidth=2,
                                        contourlabel=False,
                                        gisquality=gisquality)
                plt.show()

            # retry ?
            retry = raw_input("Do you want to modify something ? [n] ")
            if retry in ('yes', 'y', 'Y'):
                retry = True
            else:
                retry = False

    elif mode == 'lonlat_included':
        defaults = {'Iwidth':8,
                    'lonmax':'',
                    'lonmin':'',
                    'latmax':'',
                    'latmin':'',
                    'resolution':''}

        retry = True
        while retry:
            # ask and build
            (geometry, defaults) = ask_lonlat_and_build_geometry(defaults, maximize_CI_in_E)
            print show_geometry(geometry)
            if not no_display:
                # plot
                CIEdomain = footprints.proxy.field(structure='H2D',
                                                   geometry=geometry,
                                                   fid=footprints.FPDict({'zone':'C+I+E'}))
                data = numpy.ones((geometry.dimensions['Y'], geometry.dimensions['X'])) * 2.0
                data[0:geometry.dimensions['Y_CIzone'], 0:geometry.dimensions['X_CIzone']] = 1.0
                data[geometry.dimensions['Y_Iwidth']:geometry.dimensions['Y_CIzone'] - geometry.dimensions['Y_Iwidth'],
                     geometry.dimensions['X_Iwidth']:geometry.dimensions['X_CIzone'] - geometry.dimensions['X_Iwidth']] = 0.0
                CIEdomain.setdata(data)
                print "Plot domain..."
                bm = CIEdomain.geometry.make_basemap(specificproj=('nsper', {'sat_height':5000}))
                fig = CIEdomain.plotfield(existingbasemap=bm,
                                          levelsnumber=6,
                                          minmax=[-1.0, 3.0],
                                          colorbar=False,
                                          gisquality=gisquality,
                                          bluemarble=args.bluemarble)
                lldomain = build_lonlat_field(defaults, fid={'lon/lat':'included'})
                lldomain.plotfield(existingfigure=fig,
                                   existingbasemap=bm,
                                   graphicmode='contourlines',
                                   title='Domain: C+I+E \n Red contour: required lon/lat',
                                   levelsnumber=2,
                                   contourcolor='red',
                                   contourwidth=2,
                                   contourlabel=False,
                                   gisquality=gisquality)
                plt.show()

            # retry ?
            retry = raw_input("Do you want to modify something ? [n] ")
            if retry in ('yes', 'y', 'Y'):
                retry = True
            else:
                retry = False

    # output routines
    def write_blocks(out, blocks):
        for b in blocks.keys():
            out.write("&" + b + "\n")
            for v in sorted(blocks[b].keys()):
                out.write("  " + v + "=" + str(blocks[b][v]) + "," + "\n")
            out.write("/\n")

    def write_namelist(out, name, blocks):
        out.write("------------------------------" + "\n")
        out.write(name + "\n")
        out.write("------------------------------" + "\n")
        write_blocks(out, blocks)

    # compute additionnal parameters
    Xtruncation_lin = int(numpy.floor((geometry.dimensions['X'] - 1) / 2))
    Ytruncation_lin = int(numpy.floor((geometry.dimensions['Y'] - 1) / 2))
    Xtruncation_quad = int(numpy.floor((geometry.dimensions['X'] - 1) / 3))
    Ytruncation_quad = int(numpy.floor((geometry.dimensions['Y'] - 1) / 3))

    # open output
    outputfilename = "new_domain.namelists_blocks"
    out = open(outputfilename, 'w')
    out.write(show_geometry(geometry) + '\n')
    out.write("# Namelists blocks #\n")
    out.write("  ================\n")

    # PGD namelist
    namelist_name = 'namel_pre_pgd'
    blocks = {'NAM_CONF_PROJ':{}, 'NAM_CONF_PROJ_GRID':{}}
    blocks['NAM_CONF_PROJ']['XLON0'] = geometry.projection['reference_lon'].get('degrees')
    blocks['NAM_CONF_PROJ']['XLAT0'] = geometry.projection['reference_lat'].get('degrees')
    blocks['NAM_CONF_PROJ']['XRPK'] = geometry.projection['reference_lat'].get('cos_sin')[1]
    blocks['NAM_CONF_PROJ']['XBETA'] = geometry.projection['reference_lon'].get('degrees') - geometry.getcenter()[0].get('degrees')

    blocks['NAM_CONF_PROJ_GRID']['XLONCEN'] = geometry.getcenter()[0].get('degrees')
    blocks['NAM_CONF_PROJ_GRID']['XLATCEN'] = geometry.getcenter()[1].get('degrees')
    blocks['NAM_CONF_PROJ_GRID']['NIMAX'] = geometry.dimensions['X_CIzone']
    blocks['NAM_CONF_PROJ_GRID']['NJMAX'] = geometry.dimensions['Y_CIzone']
    blocks['NAM_CONF_PROJ_GRID']['XDX'] = geometry.grid['X_resolution']
    blocks['NAM_CONF_PROJ_GRID']['XDY'] = geometry.grid['Y_resolution']

    write_namelist(out, namelist_name, blocks)

    # quadclim namelist
    namelist_name = 'namel_mens_quad'
    blocks = {'NAMDIM':{}, 'NEMGEO':{}}
    blocks['NAMDIM']['NDLON'] = geometry.dimensions['X']
    blocks['NAMDIM']['NDLUXG'] = geometry.dimensions['X_CIzone']
    blocks['NAMDIM']['NDGLG'] = geometry.dimensions['Y']
    blocks['NAMDIM']['NDGUXG'] = geometry.dimensions['Y_CIzone']
    blocks['NAMDIM']['NMSMAX'] = Xtruncation_quad
    blocks['NAMDIM']['NSMAX'] = Ytruncation_quad

    blocks['NEMGEO']['ELON0'] = geometry.projection['reference_lon'].get('degrees')
    blocks['NEMGEO']['ELAT0'] = geometry.projection['reference_lat'].get('degrees')
    blocks['NEMGEO']['ELONC'] = geometry.getcenter()[0].get('degrees')
    blocks['NEMGEO']['ELATC'] = geometry.getcenter()[1].get('degrees')
    blocks['NEMGEO']['EDELX'] = geometry.grid['X_resolution']
    blocks['NEMGEO']['EDELY'] = geometry.grid['Y_resolution']

    write_namelist(out, namelist_name, blocks)

    # linclim namelist
    namelist_name = 'namel_mens_lin'
    blocks['NAMDIM']['NMSMAX'] = Xtruncation_lin
    blocks['NAMDIM']['NSMAX'] = Ytruncation_lin

    write_namelist(out, namelist_name, blocks)

    # couplingsurf namelist
    namelist_name = 'namel_e927_surf'
    blocks = {'NAMFPD':{}, 'NAMFPG':{}}

    blocks['NAMFPD']['NLON'] = geometry.dimensions['X']
    blocks['NAMFPD']['NFPLUX'] = geometry.dimensions['X_CIzone']
    blocks['NAMFPD']['NLAT'] = geometry.dimensions['Y']
    blocks['NAMFPD']['NFPGUX'] = geometry.dimensions['Y_CIzone']
    blocks['NAMFPD']['RLONC'] = geometry.getcenter()[0].get('degrees')
    blocks['NAMFPD']['RLATC'] = geometry.getcenter()[1].get('degrees')
    blocks['NAMFPD']['RDELX'] = geometry.grid['X_resolution']
    blocks['NAMFPD']['RDELY'] = geometry.grid['Y_resolution']

    blocks['NAMFPG']['FPLON0'] = geometry.projection['reference_lon'].get('degrees')
    blocks['NAMFPG']['FPLAT0'] = geometry.projection['reference_lat'].get('degrees')
    blocks['NAMFPG']['NMFPMAX'] = Xtruncation_lin
    blocks['NAMFPG']['NFPMAX'] = Ytruncation_lin

    write_namelist(out, namelist_name, blocks)

    out.close()

# end of main() ###############################################################



if __name__ == '__main__':

    ### 1. Parse arguments
    ######################
    epygramstr = 'EPyGrAM'
    parser = argparse.ArgumentParser(description='An interactive ' + epygramstr + " tool for defining a LAM domain, visualize it, \
                                                  and generate the needed namelist blocks.",
                                     epilog='End of help for: %(prog)s (' + epygramstr + ' v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, domain_maker_options['mode'])
    add_arg_to_parser(parser, domain_maker_options['no_display'])
    add_arg_to_parser(parser, domain_maker_options['maximize_CI_in_E'])
    add_arg_to_parser(parser, runtime_options['verbose'])
    add_arg_to_parser(parser, graphical_options['gis_quality'])
    add_arg_to_parser(parser, graphical_options['bluemarble'], default=1.0)

    args = parser.parse_args()

    ### 2. Initializations
    ######################
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    ### 3. Main
    ###########
    main(args.mode,
         no_display=args.no_display,
         maximize_CI_in_E=args.maximize_CI_in_E,
         gisquality=args.gisquality)

###########
### END ###
###########
