#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend H2DField with plotting methods using cartopy.
"""

import numpy
import copy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import footprints

from epygram import epygramError, config, util
from epygram import colormapping

epylog = footprints.loggers.getLogger(__name__)


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([cartoplot, cartoplot_fig_init,
                                cartoplot_background],
                               plugin_name)
    from epygram.fields import H2DField
    # defaults arguments for cartopy plots
    H2DField.default_NEfeatures = [dict(category='cultural',
                                        name='admin_0_countries',
                                        facecolor='none',
                                        edgecolor='k'),]
    H2DField.default_scatter_kw = {'s':20,
                                   'marker':',',
                                   'linewidths':0}
    H2DField.default_contour_kw = {'linewidths':1}
    H2DField.default_contourf_kw = {}
    H2DField.default_pcolormesh_kw = {}
    H2DField.default_clabel_kw = {'fmt':'%0i'}
    H2DField.default_gridlines_kw = {'draw_labels':True,
                                     'linewidth':1,
                                     'color':'k',
                                     'linestyle':'--'}
    # methods
    H2DField.cartoplot_fig_init = cartoplot_fig_init
    H2DField.cartoplot_background = cartoplot_background
    H2DField._cartoplot_treat_minmax = classmethod(_cartoplot_treat_minmax)
    H2DField._cartoplot_mesh_coords = _cartoplot_mesh_coords
    H2DField._cartoplot_get_coords = _cartoplot_get_coords
    H2DField._cartoplot_shape = _cartoplot_shape
    H2DField._cartoplot_get_ColormapHelper = classmethod(_cartoplot_get_ColormapHelper)
    H2DField._cartoplot_colormap_kwargs = classmethod(_cartoplot_colormap_kwargs)
    H2DField._cartoplot_actualplot = _cartoplot_actualplot
    H2DField._cartoplot_colorbar = classmethod(_cartoplot_colorbar)
    H2DField._cartoplot_text = _cartoplot_text
    H2DField.cartoplot = cartoplot


def meshable_geometry(geometry):
    """Returns True is the geometry is meshable (rectangularly)."""
    return geometry.rectangular_grid and geometry.name != 'unstructured'


def cartoplot_fig_init(self,
                       fig=None,
                       ax=None,
                       projection=None,
                       figsize=config.plotsizes,
                       rcparams=config.default_rcparams,
                       set_global=False):
    """Consistently set figure, ax and projection."""
    import matplotlib.pyplot as plt
    for args, kwargs in rcparams:
        plt.rc(*args, **kwargs)
    # fig, ax, proj
    if fig is None:
        assert ax is None, "Cannot specify an **ax** without a **fig**."
        if projection is None:
            projection = self.geometry.default_cartopy_CRS()
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1,1,1, projection=projection)
    else:
        if ax is None:
            if projection is None:
                projection = self.geometry.default_cartopy_CRS()
            ax = fig.add_subplot(1,1,1, projection=projection)
        else:
            assert hasattr(ax, 'projection'), "Cannot provide an **ax** which does not hold a projection."
            if projection is not None:
                assert ax.projection == projection
            else:
                projection = ax.projection
    # enlarge ?
    if set_global:
        ax.set_global()
    return fig, ax, projection


def cartoplot_background(self,
                         ax,
                         projection,
                         cartopy_features=[],
                         natural_earth_features='__default__',
                         meridians='auto',
                         parallels='auto',
                         gridlines_kw=None,
                         epygram_departments=False,
                         subzone=None,
                         extent='focus'):
    """Set cartography features, such as borders, coastlines, meridians and parallels..."""
    import cartopy
    from cartopy.mpl.gridliner import (LATITUDE_FORMATTER,
                                       LONGITUDE_FORMATTER)
    import matplotlib.ticker as mticker
    if natural_earth_features == '__default__':
        natural_earth_features = self.default_NEfeatures,
    if 'gauss' in self.geometry.name:
        default_scale = '110m'
    else:
        default_scale = '50m'
    for f in cartopy_features:
        if isinstance(f, str):
            f = getattr(cfeature, f)
        ax.add_feature(f.with_scale(default_scale))
    for f in natural_earth_features:
        f = copy.copy(f)
        if 'scale' not in f:
            f['scale'] = default_scale
        ax.add_feature(cfeature.NaturalEarthFeature(**f))
    # meridians and parallels
    meridians, parallels = util.auto_meridians_parallels(self.geometry,
                                                         meridians,
                                                         parallels,
                                                         extent=extent)
    if gridlines_kw is None:
        gridlines_kw = copy.copy(self.default_gridlines_kw)
    if (meridians, parallels) != ([], []):
        if not isinstance(projection, (ccrs.PlateCarree, ccrs.Mercator)):
            # only these projections have labels available in cartopy
            gridlines_kw.pop('draw_labels', False)
            # home-made workaround for Lambert
            if isinstance(projection, ccrs.LambertConformal) and not 'gauss' in self.geometry.name:
                from epygram.extra.cartopy_plus import lambert_parallels_meridians_labels
                lambert_parallels_meridians_labels(ax, self.geometry, projection,
                                                   meridians, parallels,
                                                   subzone=subzone)
            else:
                pass  # TODO: workarounds for other geometries ?
        gl = ax.gridlines(xlocs=meridians,
                          ylocs=parallels,
                          **gridlines_kw)
        if cartopy.__version__ <= '0.17.0':
            gl.xlabels_top = False
            gl.ylabels_right = False
        else:
            gl.top_labels = False
            gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        # stereo-polar, get more round circles
        if isinstance(projection, (ccrs.Stereographic,
                                   ccrs.NorthPolarStereo,
                                   ccrs.SouthPolarStereo,
                                   ccrs.NearsidePerspective)):
            gl.n_steps = 360
    if epygram_departments:
        if not isinstance(epygram_departments, dict):
            epygram_departments = dict(color='k')
        import json
        with open(config.installdir + '/data/french_departments.json', 'r') as dp:
            depts = json.load(dp)[1]
        for d in list(range(len(depts)))[:]:
            for part in range(len(depts[d])):
                dlon = numpy.array(depts[d][part][0])
                dlat = numpy.array(depts[d][part][1])
                xyz = projection.transform_points(ccrs.PlateCarree(), dlon, dlat)
                x = xyz[..., 0]
                y = xyz[..., 1]
                ax.plot(x, y, **epygram_departments)


def _cartoplot_treat_minmax(cls,
                            data,
                            mask_threshold,
                            minmax,
                            colormap_helper,
                            minmax_along_colorbar,
                            mask_outside_minmax):
    """Mask according to threshold, get min/max, mask according to colormap."""
    # 1/ mask values or explicit thresholding
    mask_outside = {'min':-config.mask_outside,
                    'max':config.mask_outside}
    if mask_threshold is not None:
        mask_outside.update(mask_threshold)
    data = numpy.ma.masked_outside(data,
                                   mask_outside['min'],
                                   mask_outside['max'])
    # 2/ get min/max for plot
    pmin = data.min()
    pmax = data.max()
    if minmax_along_colorbar and (minmax is not None or
                                  colormap_helper.explicit_colorbounds is not None):
        minmax_along_colorbar = '(min: {: .{precision}{type}} // max: {: .{precision}{type}})'.format(
            pmin, pmax, type='E', precision=3)
    else:
        minmax_along_colorbar = ''
    if minmax is not None:
        try:
            pmin = float(minmax[0])
        except ValueError:
            pass
        try:
            pmax = float(minmax[1])
        except ValueError:
            pass
    if colormap_helper.explicit_colorbounds is not None:
        if minmax is not None:
            epylog.warning('**minmax** overwritten by **colormap_helper.explicit_colorbounds**')
        pmin = min(colormap_helper.explicit_colorbounds)
        pmax = max(colormap_helper.explicit_colorbounds)
    # 3/ mask outside colorbounds/minmax
    if mask_outside_minmax:
        data = numpy.ma.masked_outside(data, pmin, pmax)
    return data, pmin, pmax, minmax_along_colorbar


def _cartoplot_mesh_coords(self, subzone):
    """Get coordinates of mesh corners."""
    if meshable_geometry(self.geometry):
        # In this case we need the coordinates of mesh borders
        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone,
                                                   position='lower-left')
        corners_ij = self.geometry.gimme_corners_ij(subzone=subzone)
        Imax, Jmax = corners_ij['ur']
        Imin, Jmin = corners_ij['ll']
        ii = numpy.arange(Imin, Imax + 1)
        jj = numpy.arange(Jmin, Jmax + 1)
        # compute upper line, upper border
        up_lon, up_lat = self.geometry.ij2ll(ii,
                                             numpy.full((len(ii),), Jmax, dtype=int),
                                             position='upper-right')
        ul = self.geometry.ij2ll(Imin, Jmax,
                                 position='upper-left')
        up_lon = numpy.ma.hstack((ul[0], up_lon))
        up_lat = numpy.ma.hstack((ul[1], up_lat))
        up_lon = up_lon.reshape((1, len(up_lon)))
        up_lat = up_lat.reshape((1, len(up_lat)))
        # compute right line, right border
        right_lon, right_lat = self.geometry.ij2ll(numpy.full((len(jj),), Imax, dtype=int),
                                                   jj,
                                                   position='upper-right')
        lr = self.geometry.ij2ll(Imax, Jmin,
                                 position='lower-right')
        right_lon = numpy.ma.hstack((lr[0], right_lon))
        right_lat = numpy.ma.hstack((lr[1], right_lat))
        right_lon = right_lon.reshape((len(right_lon), 1))
        right_lat = right_lat.reshape((len(right_lat), 1))
        # augment numpy array
        lons = numpy.ma.vstack((lons, up_lon[:, :-1]))
        lats = numpy.ma.vstack((lats, up_lat[:, :-1]))
        lons = numpy.ma.hstack((lons, right_lon))
        lats = numpy.ma.hstack((lats, right_lat))
        mask = numpy.ma.array(lats > 90.).filled(False)  # don't replace value where masked
        lats[mask] = 90.
        mask = numpy.ma.array(lats < -90.).filled(False)  # don't replace value where masked
        lats[mask] = -90.
    else:
        raise NotImplementedError("plot_method='pcolormesh' for Gauss or unstructured grids. Any idea ?")  # TODO: ?
    return lons, lats


def _cartoplot_get_coords(self,
                          plot_method,
                          subzone,
                          projection):
    # 1/ coordinates
    if plot_method == 'pcolormesh':
        # In this case we need the coordinates of mesh borders
        lons, lats = self._cartoplot_mesh_coords(subzone)
    else:
        # else the center
        lons, lats = self.geometry.get_lonlat_grid(subzone=subzone)
    if self.geometry.name == 'academic':
        x = (lons - lons.min()) * self.geometry.grid['X_resolution']
        y = (lats - lats.min()) * self.geometry.grid['Y_resolution']
    else:
        if isinstance(lons, numpy.ma.masked_array):
            lons.data[lons.mask] = 0.
            lats.data[lats.mask] = 0.
        xyz = projection.transform_points(ccrs.PlateCarree(), lons, lats)
        x = numpy.ma.array(xyz[..., 0])
        y = numpy.ma.array(xyz[..., 1])
        masked = numpy.logical_or(x == numpy.inf, y == numpy.inf)
        x[masked] = numpy.ma.masked
        y[masked] = numpy.ma.masked
        if isinstance(lons, numpy.ma.masked_array):
            x = numpy.ma.masked_where(lons.mask, x)
            y = numpy.ma.masked_where(lons.mask, y)
    return x, y


def _cartoplot_shape(self,
                     x, y,
                     data,
                     plot_method):
    """Shape according to plot_method."""
    xf = x
    yf = y
    zf = data
    if plot_method in ('contourf', 'contour'):
        if not self.geometry.rectangular_grid:
            # gauss grid
            # There are three kinds of masked data:
            #  - due to the shape of the gaussian grid (already masked in data, x and y):
            #      => must be suppressed from the output arrays
            #  - due to the projection (already masked in x and y):
            #      => must be suppressed from the output arrays
            #  - due to undefined values (masked in data, eg: SST over continents):
            #      => must be kept as missing in the output arrays (to actually have a hole on the plot)
            _fill_value = 1e20
            _shp = [1, 1] + list(data.shape)
            zf = self.geometry.fill_maskedvalues(data.reshape(_shp),  # protect actually masked values
-                                                _fill_value).reshape(data.shape)
            zf = numpy.ma.masked_where(xf.mask, numpy.ma.masked_where(yf.mask, zf)) # zf masked if xf or yf is masked
            xf = numpy.ma.array(numpy.ma.masked_where(zf.mask, x).compressed()) # xf masked where zf is masked, then compressed
            yf = numpy.ma.array(numpy.ma.masked_where(zf.mask, y).compressed()) # yf masked where zf is masked, then compressed
            zf = numpy.ma.masked_equal(zf.compressed(), _fill_value)  # flatten and re-masked
        elif self.geometry.dimensions['Y'] == 1:
            # unstructured grid
            xf = x.flatten()
            yf = y.flatten()
            zf = data.flatten()
    elif plot_method == 'scatter':
        # flatten data for scatter
        xf = x.flatten()
        yf = y.flatten()
        zf = data.flatten()
    return xf, yf, zf


def _cartoplot_get_ColormapHelper(cls,
                                  colormap,
                                  colorbounds):
    """Build ColormapHelper object."""
    if colormap in config.colormaps:
        cmapfile = config.colormaps[colormap]
        if cmapfile.endswith('.json'):
            colormaphelper = colormapping.get_ColormapHelper_fromfile(cmapfile)  # potentially already loaded
        else:
            raise NotImplementedError("colormaps taken from config must be .json files (and be suffixed so)")
    else:
        colormaphelper = colormapping.ColormapHelper(colormap,
                                                     explicit_colorbounds=colorbounds)
    return colormaphelper


def _cartoplot_colormap_kwargs(cls,
                               colormaphelper,
                               m, M,
                               plot_method,
                               colorsnumber,
                               colorstep,
                               center_cmap_on_0,
                               contourcolor):
    """Get colors kwargs."""
    from matplotlib.colors import cnames
    plot_kwargs = colormaphelper.kwargs_for_plot(plot_method,
                                                 (m, M),
                                                 center_cmap_on_0,
                                                 number=colorsnumber,
                                                 step=colorstep)
    if plot_method == 'contour':
        if contourcolor in cnames or contourcolor[0] == '#':
            colormap = None
        else:
            util.load_cmap(contourcolor)
            colormap = contourcolor
            contourcolor = None
        plot_kwargs.update(colors=contourcolor,
                           cmap=colormap)
    return plot_kwargs


def _cartoplot_actualplot(self,
                          ax,
                          x, y,
                          data,
                          plot_method,
                          plot_kwargs,
                          uniform,
                          colormap_helper,
                          scatter_kw,
                          contour_kw,
                          contourf_kw,
                          pcolormesh_kw,
                          contourlabel,
                          clabel_kw):
    """Actual call to matplotlib plotting functions."""
    from matplotlib.colors import cnames
    from matplotlib import tri
    if plot_method == 'contourf':
        if contourf_kw is None:
            contourf_kw = self.default_contourf_kw
        plot_kwargs.update(contourf_kw)
        if not self.geometry.rectangular_grid or self.geometry.dimensions['Y'] == 1:
            # triangulate for plotting
            triangulation = tri.Triangulation(x, y)
            if isinstance(data.mask, numpy.ndarray):  # masked values
                triangles_summits_mask = numpy.where(data.mask[triangulation.triangles], True, False)
                mask = numpy.any(triangles_summits_mask, axis=1)  # any masked summit masks triangle
                triangulation.set_mask(mask)
            pf = ax.tricontourf(triangulation, data,
                                **plot_kwargs)
        else:
            pf = ax.contourf(x, y, data,
                             **plot_kwargs)
    elif plot_method == 'scatter':
        # treat uniform field case
        if uniform:
            if colormap_helper.colormap in cnames or len(colormap_helper.colormap) == 1:
                colors = colormap_helper.colormap
            else:
                colors = 'silver'
        else:
            colors = data
        # plot
        if scatter_kw is None:
            scatter_kw = self.default_scatter_kw
        plot_kwargs.update(scatter_kw)
        pf = ax.scatter(x, y, c=colors,
                        **plot_kwargs)
    elif plot_method == 'pcolormesh':
        if pcolormesh_kw is None:
            pcolormesh_kw = self.default_pcolormesh_kw
        plot_kwargs.update(pcolormesh_kw)
        if any([isinstance(arr, numpy.ma.masked_array) for arr in [x, y]]):
            # pcolormesh cannot plot if x or y contain masked values
            # We mask data values instead of x or y
            data = numpy.ma.array(data)
            m = numpy.logical_or(numpy.ma.getmaskarray(x),
                                 numpy.ma.getmaskarray(y))
            m = numpy.logical_or(numpy.logical_or(m[:-1, :-1],
                                                  m[1:, :-1]),
                                 numpy.logical_or(m[:-1, 1:],
                                                  m[1:, 1:]))
            data[m] = numpy.ma.masked
            x = x.filled()
            y = y.filled()
        pf = ax.pcolormesh(x, y, data,
                           **plot_kwargs)
    elif plot_method == 'contour':
        if contour_kw is None:
            contour_kw = self.default_contour_kw
        plot_kwargs.update(contour_kw)
        if not self.geometry.rectangular_grid or self.geometry.dimensions['Y'] == 1:
            pf = ax.tricontour(x, y, data,
                               **plot_kwargs)
        else:
            pf = ax.contour(x, y, data,
                            **plot_kwargs)
        if contourlabel:
            if clabel_kw is None:
                clabel_kw = self.default_clabel_kw
            ax.clabel(pf,
                      **clabel_kw)
    else:
        raise NotImplementedError('plot_method=="{}"'.format(plot_method))
    return pf


def _cartoplot_colorbar(cls,
                        ax,
                        pf,
                        colorbar,
                        colorbar_over,
                        colorbar_ax_kw,
                        colormap_helper,
                        colorbar_legend,
                        minmax_along_colorbar,
                        m, M,
                        colorsnumber,
                        colorstep,
                        colorbar_kwargs):
    """Add colorbar."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # position of colorbar
    if colorbar_over is None:
        if colorbar_ax_kw is None:
            colorbar_ax_kw = dict(size="5%", pad=0.2)
        cax = make_axes_locatable(ax).append_axes(colorbar,
                                                  axes_class=plt.Axes,
                                                  **colorbar_ax_kw)
    else:
        cax = colorbar_over
    orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
    # ticks
    ticks_label = colormap_helper.ticks_label((m, M),
                                              number=colorsnumber,
                                              step=colorstep)
    ticks_position = colormap_helper.ticks_position((m, M),
                                                    number=colorsnumber,
                                                    step=colorstep)
    cb = plt.colorbar(pf,
                      orientation=orientation,
                      ticks=ticks_position,
                      cax=cax,
                      **colorbar_kwargs)
    if not numpy.all(ticks_label == ticks_position):
        cax.set_yticklabels(ticks_label)
    if minmax_along_colorbar and not colorbar_legend:
        cb.set_label(minmax_along_colorbar)
    elif colorbar_legend and not minmax_along_colorbar:
        cb.set_label(colorbar_legend)
    elif colorbar_legend and minmax_along_colorbar:
        raise ValueError("Arguments 'minmax_along_colorbar' and 'colorbar_legend' are mutually exclusive.")
    return cb

def _cartoplot_text(self, ax,
                    title,
                    uniform,
                    uniformvalue):
    """Set legend text: title and others."""
    if title is None:
        fid = self.fid.get('short',
                           self.fid[sorted(self.fid.keys())[0]])
        title = "\n".join([str(fid),
                           str(self.validity.get())])
    ax.set_title(title)
    if uniform:
        ax.text(1.02, 0.5, '(uniform field: value={})'.format(uniformvalue),
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                rotation=90.)


def cartoplot(self,
              # figure
              fig=None,
              ax=None,
              figsize=config.plotsizes,
              rcparams=config.default_rcparams,
              title=None,
              # geometry
              projection=None,
              subzone=None,
              extent='__default__',
              # graphical settings
              plot_method='__default__',
              minmax=None,
              mask_threshold=None,
              scatter_kw=None,
              contour_kw=None,
              contourf_kw=None,
              pcolormesh_kw=None,
              contourlabel=False,
              clabel_kw=None,
              # cartography
              meridians='auto',
              parallels='auto',
              gridlines_kw=None,
              cartopy_features=[],
              natural_earth_features='__default__',
              epygram_departments=False,
              # colormapping
              colormap_helper=None,
              colormap='plasma',
              colorbounds=None,
              colorsnumber=None,
              colorstep=None,
              center_cmap_on_0=False,
              contourcolor='blue',
              # colorbar
              colorbar='right',
              colorbar_over=None,
              colorbar_ax_kw=None,
              colorbar_kw={},
              colorbar_legend=None,
              minmax_along_colorbar=True,
              # takeover
              takeover=False
              ):
    """
    Plot field with **cartopy**. Returns (figure, axis).
    
    Figure settings:

    :param fig: any existing figure to be used for the plot.
        A matplotlib *figure* is a frame containing eventually several
        subplots (axes).
    :param ax: any existing axis to be used to plot on.
        A matplotlib *axis* is the subplot on which the drawing is done.
        If given, **fig** and **ax** must be consistent, i.e. **ax** being
        one of the **fig** axes.
    :param figsize: figure sizes in inches, e.g. (5, 8.5).
        Default figsize is config.plotsizes.
    :param rcparams: list of (*args, **kwargs) to be passed to pyplot.rc()
        defaults to [(('font',), dict(family='serif')),]
    :param title: title for the plot. Default is field identifier and validity.

    Geometry settings:

    :param projection: a cartopy.crs projection to be used for plot.
        Defaults to the field.geometry.default_cartopy_CRS()
    :param subzone: [LAM fields only] among ('C', 'CI'), plots the data
        resp. on the C or C+I zone.
        Default is no subzone, i.e. the whole field.
    :param extent: tune the extent of the map.
        Among ('focus', 'global', '__default__').
        'focus' will focus on the field geometry extent.
        'global' will de-zoom to have the whole globe, if possible
        (call cartopy GeoAxes.set_global()).
        '__default__' will choose one of these depending on the geometry.

    Graphical settings:

    :param plot_method: choice of the matplotlib plotting function to be
        used, among ('contourf', 'contour', 'scatter', 'pcolormesh', None).
        If None, prepare the blank figure, but skip actual plot.
        Default is 'pcolormesh' if the geometry is "meshable",
        else 'contourf'.
    :param minmax: defines the min and max values for the plot colorbar.
        Syntax: [min, max]. Strings 'min' and 'max' (default)
        will take the min and max of the field.
    :param mask_threshold: dict with min and/or max value(s) to mask outside.
    :param scatter_kw: kwargs to be passed to matplotlib's ax.scatter().
        Only for plot_method = 'scatter'.
    :param contour_kw: kwargs to be passed to matplotlib's ax.contour().
        Only for plot_method = 'contour'.
    :param contourf_kw: kwargs to be passed to matplotlib's ax.contourf().
        Only for plot_method = 'contourf'.
    :param pcolormesh_kw: kwargs to be passed to matplotlib's ax.pcolormesh().
        Only for plot_method = 'pcolormesh'.
    :param contourlabel: displays labels on contours.
        Only for plot_method = 'contour'.
    :param clabel_kw: kwargs to be passed to matplotlib's ax.clabel().
        Only for plot_method = 'contour'.

    Cartography settings:

    :param meridians: enable to fine-tune the choice of lines to
        plot, with either:
          - 'auto': automatic scaling to the map extents
          - 'default': every 10 degree
          - a list of values
          - a grid step, e.g. 5 to plot each 5 degree.
          - None: no one is plot
    :param parallels: cf. **meridians**
    :param gridlines_kw: graphical characteristics of meridians/parallels,
        arguments to be passed to cartopy's ax.gridlines(...)
    :param cartopy_features: list of cartopy.feature.??? features.
    :param natural_earth_features: list of dicts, each of them containing
        arguments to instanciate a cartopy.feature.NaturalEarthFeature(...).
        E.g. [dict(category='cultural', name='admin_1_states_provinces', facecolor='none', linestyle=':'),]
        will add states/departments/provinces.
        Cf. https://scitools.org.uk/cartopy/docs/latest/matplotlib/feature_interface.html#cartopy.feature.NaturalEarthFeature
        for details.
    :param epygram_departments: add high-resolution french departments
        limits stored in epygram. May be a dict containing lines plotting
        arguments, such as linewidth etc...
        Warning: not consistent with natural_earth_features !

    Colormap settings:

    :param colormap_helper: an instance of the class
        epygram.colormapping.ColormapHelper.
        Has priority on the following arguments.
    :param colormap: name of the ``matplotlib`` colormap to use (or an
        ``epygram`` one, or a user-defined one, cf.
        config.usercolormaps).
    :param colorbounds: levels on which to shift color.
    :param colorsnumber: number of levels for contours and colorbar.
    :param colorstep: step in value between two color shift.
    :param center_cmap_on_0: aligns the colormap center on the value 0.
    :param contourcolor: color or colormap to be used for 'contourlines'
        plot_method. It can be either a legal html color name, or a
        colormap name.

    Colorbar settings:

    :param colorbar: if *False*, hide colorbar the plot; else, defines the
        colorbar position, among ('bottom', 'right'). Defaults to 'right'.
    :param colorbar_over: an optional existing ax to plot the colorbar on.
    :param colorbar_ax_kw: kwargs to be passed to
        make_axes_locatable(ax).append_axes(colorbar, **kwargs)
    :param colorbar_kw: kwargs to be passed to matpltolib's colorbar()
    :param colorbar_legend: legend for the colorbar; mutually exclusive with arg 'minmax_along_colorbar'
    :param minmax_along_colorbar: if True and minmax is not None,
        adds min and max values along colorbar.
        
    Takeover:
    
    :param takeover: give the user more access to the objects used in the
        plot, by returning a dict containing them instead of only fig/ax
    """
    # 0/ defaults pre-sets
    if plot_method == '__default__':
        if meshable_geometry(self.geometry):
            plot_method = 'pcolormesh'
        else:
            plot_method = 'contourf'
    if natural_earth_features == '__default__':
        natural_earth_features = self.default_NEfeatures
    if not self.geometry.grid.get('LAMzone', False):
        subzone = None
    if extent == '__default__':
        if self.geometry.isglobal:
            extent = 'global'
        else:
            extent = 'focus'
    # 1/ geometry and figure
    fig, ax, projection = self.cartoplot_fig_init(fig,
                                                  ax,
                                                  projection,
                                                  figsize,
                                                  rcparams,
                                                  set_global=(extent == 'global'))
    result = dict(fig=fig, ax=ax)
    # 2/ background
    self.cartoplot_background(ax,
                              projection,
                              cartopy_features,
                              natural_earth_features,
                              meridians,
                              parallels,
                              gridlines_kw,
                              epygram_departments,
                              subzone=subzone,
                              extent=extent)
    # 3/ get data to plot
    if self.spectral:
        self.sp2gp()
    data = self.getdata(subzone=subzone)
    # 4/ handle min/max values
    if colormap_helper is None:
        colormap_helper = self._cartoplot_get_ColormapHelper(colormap,
                                                             colorbounds)
    data, pmin, pmax, minmax_along_colorbar = self._cartoplot_treat_minmax(data,
                                                                           mask_threshold,
                                                                           minmax,
                                                                           colormap_helper,
                                                                           minmax_along_colorbar,
                                                                           plot_method != 'contourf')
    if abs(float(pmin) - float(pmax)) < config.epsilon and plot_method is not None:
        epylog.warning("uniform field: plot as 'points'.")
        plot_method = 'scatter'
        uniform = True
        uniformvalue = pmin
    else:
        uniform = False
        uniformvalue = None
    # 5/ get coordinates
    x, y = self._cartoplot_get_coords(plot_method,
                                      subzone,
                                      projection)
    # FIXME: (CLEANME) fixed by Sebastien ?
    assert numpy.inf not in [x.min(), x.max(), y.min(), y.max()], \
        "Domain is too large to be plotted on space view"
    assert -numpy.inf not in [x.min(), x.max(), y.min(), y.max()], \
        "Domain is too large to be plotted on space view"
    # 6/ shape data/lons/lats
    x, y, data = self._cartoplot_shape(x, y,
                                       data,
                                       plot_method)
    if extent == 'focus':
        xyz = ax.projection.transform_points(projection, x, y)
        xyz = numpy.ma.masked_where(numpy.abs(xyz) == numpy.inf, xyz)
        ax.set_xlim((xyz[..., 0].min(), xyz[..., 0].max()))
        ax.set_ylim((xyz[..., 1].min(), xyz[..., 1].max()))
        del xyz
    # 7/ colormapping
    plot_kwargs = self._cartoplot_colormap_kwargs(colormap_helper,
                                                  pmin, pmax,
                                                  plot_method,
                                                  colorsnumber,
                                                  colorstep,
                                                  center_cmap_on_0,
                                                  contourcolor)
    if plot_method is not None:
        # 8/ plot
        elements = self._cartoplot_actualplot(ax,
                                              x, y,
                                              data,
                                              plot_method,
                                              plot_kwargs,
                                              uniform,
                                              colormap_helper,
                                              scatter_kw,
                                              contour_kw,
                                              contourf_kw,
                                              pcolormesh_kw,
                                              contourlabel,
                                              clabel_kw)
        result['plot_elements'] = elements
        # 9/ colorbar
        if colorbar and plot_method != 'contour':
            cb = self._cartoplot_colorbar(ax,
                                          elements,
                                          colorbar,
                                          colorbar_over,
                                          colorbar_ax_kw,
                                          colormap_helper,
                                          colorbar_legend,
                                          minmax_along_colorbar,
                                          pmin, pmax,
                                          colorsnumber,
                                          colorstep,
                                          colorbar_kw)
            result['colorbar'] = cb
    # 10/ texts
    self._cartoplot_text(ax, title, uniform, uniformvalue)
    return result if takeover else (fig, ax)
