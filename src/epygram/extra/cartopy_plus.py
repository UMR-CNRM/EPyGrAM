#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bronx.fancies import loggers
import cartopy.crs as ccrs

logger = loggers.getLogger(__name__)

natural_earth_countries_departments = [dict(category='cultural',
                                            name='admin_0_countries',
                                            facecolor='none'),
                                       dict(category='cultural',
                                            name='admin_1_states_provinces',
                                            facecolor='none',
                                            linestyle=':')]


def formatted_meridian(value):
    return '{:g}°{}'.format(abs(value), 'E' if value >= 0 else 'W')


def formatted_parallel(value):
    return '{:g}°{} '.format(abs(value), 'N' if value >= 0 else 'S')


def lambert_parallels_meridians_labels(ax, geometry, projection,
                                       meridians, parallels,
                                       subzone=None):
    """
    An Epygram workaround to cartopy's inability to write labels to parallels
    and meridians in Lambert projection.
    """
    pc = ccrs.PlateCarree()
    corners = geometry.gimme_corners_ll(subzone=subzone)
    lons, lats = geometry.get_lonlat_grid(subzone=subzone)
    left_border = (lons[:, 0], lats[:, 0])
    bottom_border = (lons[0, :], lats[0, :])
    corners = {k:projection.transform_point(*v, src_crs=pc)
               for k, v in corners.items()}
    (x0, y0) = corners['ll']
    # MERIDIANS
    # filter: only keep the meridians that cross the bottom border
    meridians = [m for m in meridians
                 if bottom_border[0].min() <= m <= bottom_border[0].max()]
    for m in meridians:
        # find closest gridpoint on the bottom border, take its latitude
        dist = bottom_border[0] - m
        dmin = 360
        for i, d in enumerate(dist):
            if abs(d) <= dmin:
                dmin = abs(d)
            else:
                break
        lat_crossborder = bottom_border[1][i]
        x, _ = projection.transform_point(m, lat_crossborder, src_crs=pc)
        ax.text(x, y0 , '\n\n' + formatted_meridian(m),
                horizontalalignment='center',
                verticalalignment='center')
    # PARALLELS
    # filter: only keep the parallels that cross the left border
    parallels = [p for p in parallels
                 if left_border[1].min() <= p <= left_border[1].max()]
    for p in parallels:
        # find closest gridpoint on the left border, take its longitude
        dist = left_border[1] - p
        dmin = 90
        for i, d in enumerate(dist):
            if abs(d) <= dmin:
                dmin = abs(d)
            else:
                break
        lon_crossborder = left_border[0][i]
        # compute y coordinate and write label
        _, y = projection.transform_point(lon_crossborder, p, src_crs=pc)
        ax.text(x0, y, formatted_parallel(p),
                horizontalalignment='right',
                verticalalignment='center')
