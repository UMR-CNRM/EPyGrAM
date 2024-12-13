#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend Geometry with plotting methods using cartopy.
"""

import numpy

from cartopy import crs as ccrs

import footprints

from epygram import epygramError

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([default_cartopy_CRS, get_cartopy_extent],
                               plugin_name)
    from epygram.geometries import ProjectedGeometry
    # defaults arguments for cartopy plots
    ProjectedGeometry.default_cartopy_CRS = default_cartopy_CRS
    ProjectedGeometry.get_cartopy_extent = get_cartopy_extent


def default_cartopy_CRS(self):
    """
    Create a cartopy.crs appropriate to the Geometry.
    """
    globe = ccrs.Globe(semimajor_axis=self.geoid['a'], semiminor_axis=self.geoid['b'])
    if self.name == 'lambert':
        if self.secant_projection:
            lat_0 = (self.projection['secant_lat1'].get('degrees') +
                     self.projection['secant_lat2'].get('degrees')) / 2.
            secant_lats = (self.projection['secant_lat1'].get('degrees'),
                           self.projection['secant_lat2'].get('degrees'))
        else:
            lat_0 = self.projection['reference_lat'].get('degrees')
            secant_lats = (lat_0, lat_0)
        crs = ccrs.LambertConformal(
            central_longitude=self.projection['reference_lon'].get('degrees'),
            central_latitude=lat_0,
            standard_parallels=secant_lats,
            globe=globe)
    elif self.name == 'mercator':
        if self.secant_projection:
            lat = 'secant_lat'
        else:
            lat = 'reference_lat'
        crs = ccrs.Mercator(
            central_longitude=self._center_lon.get('degrees'),
            latitude_true_scale=self.projection[lat].get('degrees'),
            globe=globe)
    elif self.name == 'polar_stereographic':
        if self.secant_projection:
            lat = 'secant_lat'
        else:
            lat = 'reference_lat'
        crs = ccrs.Stereographic(
            central_latitude=numpy.copysign(90., self.projection[lat].get('degrees')),
            central_longitude=self.projection['reference_lon'].get('degrees'),
            true_scale_latitude=self.projection[lat].get('degrees'),
            globe=globe)
    elif self.name == 'space_view':
        if self.projection['satellite_lat'].get('degrees') == 0.:
            crs = ccrs.Geostationary(
                central_longitude=self.projection['satellite_lon'].get('degrees'),
                satellite_height=self.projection['satellite_height'],
                globe=globe)
        else:
            raise NotImplementedError("Implementation to be tested")
            #NearsidePerspective projection does not handle elliptical globes
            #Thus, we do not set the globe attribute from the geometry geoid
            crs = ccrs.NearsidePerspective(
                central_longitude=self.projection['satellite_lon'].get('degrees'),
                central_latitude=self.projection['satellite_lat'].get('degrees'),
                satellite_height=self.projection['satellite_height'])
    else:
        raise epygramError("Projection name unknown.")
    return crs


def get_cartopy_extent(self, subzone=None):
    """
    Gets the extension of the geometry in the default crs
    :param subzone: defines the LAM subzone to be included, in LAM case,
                    among: 'C', 'CI'.
    :return: (xmin, xmax, ymin, ymax) as expected by matplotlib's ax.set_extent
    """
    # Actually, this method transforms projected coordinates from pyproj to
    # projected coordinates from cartopy's CRS
    # Both make use of proj4 but the transformation is not exactly the same.
    # Easiest method is to transform pyproj coordinates into lat/lon then
    # into crs coordinates. But this method does not work when the point
    # does not have lat/lon. This is the case in space_view where corners
    # can be outside the sphere.
    crs = self.default_cartopy_CRS()
    if self.name != 'space_view':
        (llcrnrlon, llcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ll'])
        (urcrnrlon, urcrnrlat) = self.ij2ll(*self.gimme_corners_ij(subzone)['ur'])
        ll = crs.transform_points(ccrs.PlateCarree(),
                                  numpy.array(llcrnrlon),
                                  numpy.array(llcrnrlat))[0]
        ur = crs.transform_points(ccrs.PlateCarree(),
                                  numpy.array(urcrnrlon),
                                  numpy.array(urcrnrlat))[0]
    else:
        # one point that can be viewed from satellite
        lon0 = self.projection['satellite_lon'].get('degrees')
        lat0 = self.projection['satellite_lat'].get('degrees')
        # projected coordinates in both systems of this point
        pos_crs0 = crs.transform_points(ccrs.PlateCarree(),
                                        numpy.array(lon0),
                                        numpy.array(lat0))[0]
        pos_epy0 = self.ll2xy(lon0, lat0)
        # Corners
        ll = self.ij2xy(*self.gimme_corners_ij(subzone)['ll'])
        ur = self.ij2xy(*self.gimme_corners_ij(subzone)['ur'])
        ll = (ll[0] - pos_epy0[0] + pos_crs0[0],
              ll[1] - pos_epy0[1] + pos_crs0[1])
        ur = (ur[0] - pos_epy0[0] + pos_crs0[0],
              ur[1] - pos_epy0[1] + pos_crs0[1])
    return (ll[0], ur[0], ll[1], ur[1])
