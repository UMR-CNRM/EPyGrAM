#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains functions for building a LAM domain.
"""

import math
import numpy

from footprints import FPDict
from footprints import proxy as fpx

from .util import (Ezone_minimum_width, maxdims_security_barrier,
                   threshold_mercator_lambert, threshold_pole_distance_lambert,
                   vgeom,
                   default_Iwidth,
                   projections_g2p, projections_g2s,
                   projections_s2g, projections_s2p)
from epygram import epygramError, epylog
from epygram.config import epsilon, margin_points_within_Czone
from epygram.util import Angle
from epygram.geometries.SpectralGeometry import SpectralGeometry, nearest_greater_FFT992compliant_int
from epygram.geometries import ProjectedGeometry, RegLLGeometry


def build_geometry(center_lon, center_lat,
                   Xpoints_CI, Ypoints_CI,
                   resolution,
                   Iwidth=None,
                   tilting=0.,
                   reference_lat=None,
                   force_projection=None,
                   maximize_CI_in_E=False,
                   interactive=False):
    """
    Build an *ad hoc* geometry from the input given parameters.
    Beware only secant projection are available here.

    :param center_lon: longitude of the domain center
    :param center_lat: latitude of the domain center
    :param Xpoints_CI: number of gridpoints in C+I zone, zonal dimension X
    :param Ypoints_CI: number of gridpoints in C+I zone, meridian dimension Y
    :param resolution: resolution of the grid in m
    :param Iwidth: width of the I-zone
    :param tilting: optional inclination of the grid, in degrees
                    (only for 'polar_stereographic' and 'lambert' projections)
    :param reference_lat: (disadvised use) reference latitude of the projection;
                          beware of consistency with *force_projection*
    :param force_projection: force projection among ('polar_stereographic',
                             'lambert', 'mercator')
    :param maximize_CI_in_E: extend the C+I zone inside the C+I+E zone, in order
                             to have a E-zone width of 11 points
    :param interactive: interactive mode, to fine-tune the projection
    """
    Xpoints_CI = int(Xpoints_CI)
    Ypoints_CI = int(Ypoints_CI)
    if Iwidth is not None:
        Iwidth = int(Iwidth)
    # begin to build a horizontal geometry
    Xpoints_CIE = nearest_greater_FFT992compliant_int(Xpoints_CI + Ezone_minimum_width)
    Ypoints_CIE = nearest_greater_FFT992compliant_int(Ypoints_CI + Ezone_minimum_width)
    if maximize_CI_in_E:
        Xpoints_CI = Xpoints_CIE - Ezone_minimum_width
        Ypoints_CI = Ypoints_CIE - Ezone_minimum_width
    # dimensions
    if Iwidth is None:
        Iwidth = default_Iwidth(resolution)
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
    if reference_lat in (90., -90.):
        assert force_projection in (None, 'polar_stereographic')
    elif reference_lat == 0.:
        assert force_projection in (None, 'mercator')
    elif reference_lat is not None:
        assert force_projection in (None, 'lambert')
    elif reference_lat is None:
        reference_lat = center_lat
    projection = {'reference_lon':Angle(center_lon + tilting, 'degrees'),
                  'reference_lat':Angle(reference_lat, 'degrees'),
                  'rotation':Angle(0.0, 'degrees'),
                  }
    # grid
    grid = {'X_resolution':resolution,
            'Y_resolution':resolution,
            'LAMzone':'CIE',
            'input_lon':Angle(center_lon, 'degrees'),
            'input_lat':Angle(center_lat, 'degrees'),
            'input_position':(float(dimensions['X_CIzone'] - 1) / 2.,
                              float(dimensions['Y_CIzone'] - 1) / 2.)
            }

    # try to guess best projection
    if abs(center_lat) >= threshold_mercator_lambert:  # lambert or polar stereographic
        # we make a "first guess" with Lambert Geometry, just to check the pole is not inside the domain
        geometryname = 'lambert'
        geometry = ProjectedGeometry(name=geometryname,
                                     grid=grid,
                                     dimensions=dimensions,
                                     projection=projection,
                                     vcoordinate=vgeom,
                                     position_on_horizontal_grid='center')

        # test if pole in lambert grid
        pole_in_domain = (geometry.point_is_inside_domain_ll(0, 90) or
                          geometry.point_is_inside_domain_ll(0, -90))
        if pole_in_domain:
            projname = 'polar_stereographic'
        else:
            projname = 'lambert'
    else:
        projname = 'mercator'
        pole_in_domain = False
    if force_projection is not None:
        projname = force_projection
    if interactive:
        if force_projection is not None:
            print("User-requested projection for this center and these dimensions is:")
        else:
            print("Advised projection for this center and these dimensions is:")
        print("=>", projections_g2p[projname], '(', projections_g2s[projname], ')')
        print("Other choices:", projections_s2p)
        force_projection = input("Chosen projection [" + projections_g2s[projname] + "]: ")
        if force_projection != '':
            projname = projections_s2g[force_projection]
    assert (not (pole_in_domain and projname != 'polar_stereographic')), \
        "requested projection [" + projname + "] is not possible: " + \
        "a pole is inside domain, the only possible projection is 'polar_stereographic'."
    assert projname in projections_g2s.keys(), \
        "unknown projection [" + projname + "]"

    if projname == 'polar_stereographic':
        reference_lat = math.copysign(90.0, center_lat)
    elif projname == 'mercator':
        reference_lat = 0.0
        projection['reference_lon'] = Angle(center_lon, 'degrees')
        if tilting != 0.0:
            epylog.warning("! Tilting ignored: not available for Mercator projection.")
    elif projname == 'lambert':
        if interactive:
            print("Advised reference latitude for Lambert domain is center latitude:")
            accepted_lat = input("Reference latitude [" + str(center_lat) + "]: ")
            if accepted_lat != '':
                reference_lat = float(accepted_lat)
            else:
                reference_lat = float(center_lat)
    projection['reference_lat'] = Angle(reference_lat, 'degrees')
    geometry = ProjectedGeometry(name=projname,
                                 grid=grid,
                                 dimensions=dimensions,
                                 projection=projection,
                                 vcoordinate=vgeom,
                                 position_on_horizontal_grid='center',)
    return geometry


def build_CIE_field(geometry):
    """
    Build a field according to geometry, with constant values for each of the
    C,I,E zones.
    """
    CIEdomain = fpx.field(structure='H2D',
                          geometry=geometry,
                          fid=FPDict({'zone':'C+I+E'}))
    data = numpy.ones((geometry.dimensions['Y'],
                       geometry.dimensions['X'])) * 2.0
    data[0:geometry.dimensions['Y_CIzone'],
         0:geometry.dimensions['X_CIzone']] = 1.0
    data[geometry.dimensions['Y_Iwidth']:geometry.dimensions['Y_CIzone'] - geometry.dimensions['Y_Iwidth'],
         geometry.dimensions['X_Iwidth']:geometry.dimensions['X_CIzone'] - geometry.dimensions['X_Iwidth']] = 0.0
    CIEdomain.setdata(data)
    return CIEdomain


def build_lonlat_geometry(ll_boundaries, resolution=None):
    """
    Build a lonlat geometry given lon/lat boundaries and optionally
    resolution in degrees.
    """
    if resolution is not None:
        if isinstance(resolution, (list, tuple)):
            xres, yres = resolution
        else:
            xres = resolution
            yres = resolution
        xwidth = (ll_boundaries['lonmax'] - ll_boundaries['lonmin']) / xres
        ywidth = (ll_boundaries['latmax'] - ll_boundaries['latmin']) / yres
        if (xwidth - round(xwidth)) > 1e-6:
            raise ValueError('resolution:{} cannot divide span [lonmin:lonmax]'.
                             format(resolution))
        else:
            xwidth = int(round(xwidth))
        if (ywidth - round(ywidth)) > 1e-6:
            raise ValueError('resolution:{} cannot divide span [latmin:latmax]'.
                             format(resolution))
        else:
            ywidth = int(round(ywidth))
    elif resolution is None:
        xwidth = 1000
        ywidth = 1000
    xres = (ll_boundaries['lonmax'] - ll_boundaries['lonmin']) / xwidth
    yres = (ll_boundaries['latmax'] - ll_boundaries['latmin']) / ywidth
    llgrid = {'input_lon':Angle(ll_boundaries['lonmin'], 'degrees'),
              'input_lat':Angle(ll_boundaries['latmin'], 'degrees'),
              'input_position':(0, 0),
              'X_resolution':Angle(xres, 'degrees'),
              'Y_resolution':Angle(yres, 'degrees')
              }
    lldims = {'X':xwidth + 1, 'Y':ywidth + 1}
    llgeometry = RegLLGeometry(name='regular_lonlat',
                               grid=llgrid,
                               dimensions=lldims,
                               vcoordinate=vgeom,
                               position_on_horizontal_grid='center')
    return llgeometry


def build_lonlat_field(ll_boundaries,
                       fid={'lon/lat':'template'},
                       resolution=None):
    """
    Build a lonlat field empty except on the border, given lon/lat boundaries
    and optionally resolution in degrees.
    """
    ll_geometry = build_lonlat_geometry(ll_boundaries, resolution=resolution)
    ll_domain = ll_geometry.make_field(fid=fid)
    return ll_domain


def build_geometry_fromlonlat(lonmin, lonmax,
                              latmin, latmax,
                              resolution,
                              Iwidth=None,
                              force_projection=None,
                              maximize_CI_in_E=False,
                              interactive=False):
    """
    Build an *ad hoc* geometry from the input given parameters.
    Beware only secant projection are available here.

    :param lonmin: minimum longitude of the domain
    :param lonmax: maximum longitude of the domain
    :param latmin: minimum latitude of the domain
    :param latmax: maximum latitude of the domain
    :param resolution: resolution of the grid in m
    :param Iwidth: width of the I-zone
    :param force_projection: force projection among ('polar_stereographic',
                             'lambert', 'mercator')
    :param maximize_CI_in_E: extend the C+I zone inside the C+I+E zone, in order
                             to have a E-zone width of 11 points
    :param interactive: interactive mode, to fine-tune the projection
    """
    if Iwidth is not None:
        Iwidth = int(Iwidth)
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
        if latmax < 90. - threshold_pole_distance_lambert and \
           latmin > -90. + threshold_pole_distance_lambert:
            projname = 'lambert'
        else:
            projname = 'polar_stereographic'
    else:
        projname = 'mercator'

    if interactive:
        print("Advised projection with regards to domain center latitude is:")
        print("=>", projections_g2p[projname], '(', projections_g2s[projname], ')')
        print("Other choices:", projections_s2p)
        accepted_projection = input("Chosen projection [" + projections_g2s[projname] + "]: ")
        if accepted_projection != '':
            projname = projections_s2g[accepted_projection]
    if force_projection is not None:
        projname = force_projection
        print("Projection forced by dummy argument to:", force_projection)
    if projname == 'polar_stereographic':
        reference_lat = math.copysign(90.0, center_lat)
    elif projname == 'mercator':
        reference_lat = 0.0
    elif projname == 'lambert':
        if interactive:
            print("Advised center latitude for Lambert domain is mean(Northern, Southern):")
            accepted_lat = input("Center latitude [" + str(center_lat) + "]: ")
            if accepted_lat != '':
                center_lat = float(accepted_lat)
        reference_lat = center_lat
        if interactive:
            print("Advised reference latitude for Lambert domain is center latitude:")
            accepted_lat = input("Reference latitude [" + str(reference_lat) + "]: ")
            if accepted_lat != '':
                reference_lat = float(accepted_lat)

    if Iwidth is None:
        Iwidth = default_Iwidth(resolution)
    Xpoints_CI = 2 * Iwidth + 1
    Ypoints_CI = 2 * Iwidth + 1
    lonlat_included = False
    while not lonlat_included:
        Xpoints_CIE = nearest_greater_FFT992compliant_int(Xpoints_CI + Ezone_minimum_width)
        Ypoints_CIE = nearest_greater_FFT992compliant_int(Ypoints_CI + Ezone_minimum_width)
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
        projection = {'reference_lon':Angle(center_lon, 'degrees'),
                      'reference_lat':Angle(reference_lat, 'degrees'),
                      'rotation':Angle(0.0, 'degrees'),
                      }
        # grid
        grid = {'X_resolution':resolution,
                'Y_resolution':resolution,
                'LAMzone':'CIE',
                'input_lon':Angle(center_lon, 'degrees'),
                'input_lat':Angle(center_lat, 'degrees'),
                'input_position':(float(dimensions['X_CIzone'] - 1) / 2.,
                                  float(dimensions['Y_CIzone'] - 1) / 2.)
                }
        # first guess for lambert, to check that pole is not in domain
        if projname == 'lambert':
            geometry = ProjectedGeometry(name=projname,
                                         grid=grid,
                                         dimensions=dimensions,
                                         projection=projection,
                                         vcoordinate=vgeom,
                                         position_on_horizontal_grid='center')
            pole_in_domain = (geometry.point_is_inside_domain_ll(0, 90) or geometry.point_is_inside_domain_ll(0, -90))
            if pole_in_domain:
                epylog.warning("Pole is inside Lambert domain => shifted to Polar Stereographic projection !")
                projname = 'polar_stereographic'
                projection['reference_lat'] = Angle(math.copysign(90.0, center_lat), 'degrees')

        # guess
        geometry = ProjectedGeometry(name=projname,
                                     grid=grid,
                                     dimensions=dimensions,
                                     projection=projection,
                                     vcoordinate=vgeom,
                                     position_on_horizontal_grid='center')
        # test whether the lonlat corners are inside domain
        points_to_test = [(lonmax, latmax), (lonmax, latmin),
                          (lonmin, latmax), (lonmin, latmin),
                          (lonmax, (latmin + latmax) / 2.), (lonmin, (latmin + latmax) / 2.),
                          ((lonmin + lonmax) / 2., latmax), ((lonmin + lonmax) / 2., latmin)]
        IminC, JminC = geometry.gimme_corners_ij('C')['ll']
        ImaxC, JmaxC = geometry.gimme_corners_ij('C')['ur']
        # we cannot use geometry.point_is_inside_domain_ll() because we need to
        # know in which direction we need to extend
        margin = float(margin_points_within_Czone)
        xlonlat_included = all(
            [margin + IminC < geometry.ll2ij(*c)[0] < ImaxC - margin
             for c in points_to_test])
        ylonlat_included = all(
            [margin + JminC < geometry.ll2ij(*c)[1] < JmaxC - margin
             for c in points_to_test])
        lonlat_included = xlonlat_included and ylonlat_included
        if not lonlat_included:
            if not xlonlat_included:
                Xpoints_CI = Xpoints_CIE - Ezone_minimum_width + 1
            if not ylonlat_included:
                Ypoints_CI = Ypoints_CIE - Ezone_minimum_width + 1
            if Xpoints_CI > maxdims_security_barrier or \
               Ypoints_CI > maxdims_security_barrier:
                raise epygramError("Domain is too large, > " + str(maxdims_security_barrier) + " points.")

    return geometry


def build_geom_from_e923nam(nam):
    """
    Build geometry and spectral geometry objects, given e923-like namelist
    blocks.
    """
    if nam['NAMCT0']['LRPLANE']:
        geometryclass = ProjectedGeometry
        if nam['NEMGEO']['ELAT0'] <= epsilon:
            geometryname = 'mercator'
        elif 90. - abs(nam['NEMGEO']['ELAT0']) <= epsilon:
            geometryname = 'polar_stereographic'
        elif epsilon < abs(nam['NEMGEO']['ELAT0']) < 90. - epsilon:
            geometryname = 'lambert'
        kwargs = dict(
            projection=dict(reference_lat=Angle(nam['NEMGEO']['ELAT0'], 'degrees'),
                            reference_lon=Angle(nam['NEMGEO']['ELON0'], 'degrees'),
                            rotation=Angle(0.,'degrees')),
            grid=dict(input_lat=Angle(nam['NEMGEO']['ELATC'], 'degrees'),
                      input_lon=Angle(nam['NEMGEO']['ELONC'], 'degrees'),
                      input_position=((nam['NAMDIM']['NDLUXG'] - 1) / 2,
                                      (nam['NAMDIM']['NDGUXG'] - 1) / 2),
                      X_resolution=nam['NEMGEO']['EDELX'],
                      Y_resolution=nam['NEMGEO']['EDELY'],
                      LAMzone='CI'),
            dimensions=dict(X=nam['NAMDIM']['NDLUXG'],
                            Y=nam['NAMDIM']['NDGUXG'],
                            X_CIzone=nam['NAMDIM']['NDLUXG'],
                            Y_CIzone=nam['NAMDIM']['NDGUXG'],
                            X_Czone=nam['NAMDIM']['NDLUXG'] - 2 * nam['NEMDIM']['NBZONL'],
                            Y_Czone=nam['NAMDIM']['NDGUXG'] - 2 * nam['NEMDIM']['NBZONG'],
                            X_Iwidth=nam['NEMDIM']['NBZONL'],
                            Y_Iwidth=nam['NEMDIM']['NBZONG']))
    else:
        geometryclass = RegLLGeometry
        geometryname = 'regular_lonlat'
        kwargs = dict(
            grid=dict(input_lat=Angle(nam['NEMGEO']['ELATC'], 'degrees'),
                      input_lon=Angle(nam['NEMGEO']['ELONC'], 'degrees'),
                      input_position=((nam['NAMDIM']['NDLON'] - 1) / 2,
                                      (nam['NAMDIM']['NDGLG'] - 1) / 2),
                      X_resolution=Angle(nam['NEMGEO']['EDELX'], 'degrees'),
                      Y_resolution=Angle(nam['NEMGEO']['EDELY'], 'degrees')),
            dimensions=dict(X=nam['NAMDIM']['NDLON'],
                            Y=nam['NAMDIM']['NDGLG']))
    geom = geometryclass(name=geometryname,
                         vcoordinate=vgeom,
                         position_on_horizontal_grid='center',
                         **kwargs
                         )
    if 'NMSMAX' in nam['NAMDIM'].keys() and 'NSMAX' in nam['NAMDIM'].keys():
        spgeom = SpectralGeometry(space='bi-fourier',
                                  truncation=dict(in_X=nam['NAMDIM']['NMSMAX'],
                                                  in_Y=nam['NAMDIM']['NSMAX']))
    else:
        spgeom = None
    return (geom, spgeom)


def compute_lonlat_included(geometry):
    """
    Computes a lon/lat domain included in the C zone of the model domain,
    with a 1 gridpoint margin.
    """
    (longrid, latgrid) = geometry.get_lonlat_grid(subzone='C')
    lonmin = max(longrid[:, 1])
    lonmax = min(longrid[:, -2])
    latmax = min(latgrid[-2, :])
    latmin = max(latgrid[1, :])

    return {'lonmin':lonmin, 'lonmax':lonmax, 'latmin':latmin, 'latmax':latmax}
