#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy
import copy

import footprints
from bronx.graphics.axes import set_figax

from epygram import config, util, epygramError
from epygram.geometries import H2DGeometry
from . import gimme_one_point
from .D3Field import D3Field
from .PointField import PointField

epylog = footprints.loggers.getLogger(__name__)


class _H2DCartopyPlot(object):
    """
    Plugin for H2DField for plotting with Cartopy.
    """
    # defaults arguments for cartopy plots
    default_scatter_kw = {'s':20,
                          'marker':',',
                          'linewidths':0}
    default_contour_kw = {'linewidths':1}
    default_clabel_kw = {'fmt':'%0i'}
    default_gridlines_kw = {'draw_labels':True,
                            'linewidth':1,
                            'linestyle':'--'}

    def _cartoplot_fig_init(self,
                            fig,
                            ax,
                            projection,
                            figsize,
                            rcparams,
                            set_global):
        """Consistently set figure, ax and projection."""
        import matplotlib.pyplot as plt
        if rcparams is None:
            rcparams = [(('font',), dict(family='serif')), ]
        for args, kwargs in rcparams:
            plt.rc(*args, **kwargs)
        # fig, ax, proj
        if figsize is None:
            figsize = config.plotsizes
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

    def _cartoplot_background(self,
                              # fig,
                              ax,
                              projection,
                              natural_earth_features,
                              meridians,
                              parallels,
                              gridlines_kw,
                              epygram_departments):
        """Set cartography features, such as borders, coastlines, meridians and parallels..."""
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from cartopy.mpl.gridliner import (LATITUDE_FORMATTER,
                                           LONGITUDE_FORMATTER)
        for f in natural_earth_features:
            f = copy.copy(f)
            if 'scale' not in f:
                if 'gauss' in self.geometry.name:
                    f['scale'] = '110m'
                else:
                    f['scale'] = '10m'
            ax.add_feature(cfeature.NaturalEarthFeature(**f))
        # meridians and parallels
        meridians, parallels = util.auto_meridians_parallels(self.geometry,
                                                             meridians,
                                                             parallels)
        if gridlines_kw is None:
            gridlines_kw = copy.copy(self.default_gridlines_kw)
        if (meridians, parallels) != ([], []):
            if not isinstance(projection, (ccrs.PlateCarree, ccrs.Mercator)):
                # only these projections have labels available in cartopy
                gridlines_kw.pop('draw_labels', False)
                # from cartopy_plus import lambert_ticks_workaround
                # draw_labels = gridlines_kw.pop('draw_labels', False)
                # FIXME: canvas drawing can be over-consuming !!!
                # lambert_ticks_workaround(fig, ax, projection,
                #                          draw_labels, parallels, meridians)
            gl = ax.gridlines(xlocs=meridians,
                              ylocs=parallels,
                              **gridlines_kw)
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            # stereo-polar, get more round circles
            if isinstance(projection, (ccrs.Stereographic,
                                       ccrs.NorthPolarStereo,
                                       ccrs.SouthPolarStereo)):
                gl.n_steps = 90
        if epygram_departments:
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
                    ax.plot(x, y, color='k')

    @classmethod
    def _cartoplot_treat_minmax(cls,
                                data,
                                mask_threshold,
                                minmax,
                                colorbounds,
                                colormap,
                                minmax_along_colorbar):
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
                                      colormap in config.colormaps_scaling or
                                      colorbounds is not None):
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
        if colorbounds is not None:
            if minmax is not None:
                epylog.warning('**minmax** overwritten by **colorbounds**')
            pmin = min(colorbounds)
            pmax = max(colorbounds)
        # 3/ mask outside colorbounds/minmax
        data = numpy.ma.masked_outside(data, pmin, pmax)
        return data, pmin, pmax, minmax_along_colorbar

    def _cartoplot_mesh_coords(self, subzone):
        """Get coordinates of mesh corners."""
        if self.geometry.rectangular_grid:
            # In this case we need the coordinates of mesh borders
            lons, lats = self.geometry.get_lonlat_grid(subzone=subzone,
                                                       position='lower-left')
            corners_ij = self.geometry.gimme_corners_ij(subzone=subzone)
            Imax, Jmax = corners_ij['ur']
            Imin, Jmin = corners_ij['ll']
            ii = list(range(Imin, Imax + 1))
            jj = list(range(Jmin, Jmax + 1))
            # compute upper line, upper border
            up_lon, up_lat = self.geometry.ij2ll(ii,
                                                 [Jmax for _ in ii],
                                                 position='upper-right')
            ul = self.geometry.ij2ll(Imin, Jmax,
                                     position='upper-left')
            up_lon = numpy.hstack((numpy.array([ul[0]]), up_lon))
            up_lat = numpy.hstack((numpy.array([ul[1]]), up_lat))
            up_lon = up_lon.reshape((1, len(up_lon)))
            up_lat = up_lat.reshape((1, len(up_lat)))
            # compute right line, right border
            right_lon, right_lat = self.geometry.ij2ll([Imax for _ in jj],
                                                       jj,
                                                       position='upper-right')
            lr = self.geometry.ij2ll(Imax, Jmin,
                                     position='lower-right')
            right_lon = numpy.hstack((numpy.array([lr[0]]), right_lon))
            right_lat = numpy.hstack((numpy.array([lr[1]]), right_lat))
            right_lon = right_lon.reshape((len(right_lon), 1))
            right_lat = right_lat.reshape((len(right_lat), 1))
            # augment numpy array
            lons = numpy.vstack((lons, up_lon[:, :-1]))
            lats = numpy.vstack((lats, up_lat[:, :-1]))
            lons = numpy.hstack((lons, right_lon))
            lats = numpy.hstack((lats, right_lat))
            lats[lats > 90.] = 90.
            lats[lats < -90.] = -90.
        else:
            raise NotImplementedError("plot_method='pcolormesh' for Gauss grids. Any idea ?")  # TODO: ?
        return lons, lats

    def _cartoplot_get_coords(self,
                              plot_method,
                              subzone,
                              projection):
        import cartopy.crs as ccrs
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
            xyz = projection.transform_points(ccrs.PlateCarree(), lons, lats)
            x = xyz[..., 0]
            y = xyz[..., 1]
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
                _fill_value = 1e20
                _shp = [1,1] + list(data.shape)
                zf = self.geometry.fill_maskedvalues(data.reshape(_shp),  # protect actually masked values
                                                     _fill_value)
                zf = numpy.ma.masked_equal(zf.compressed(), _fill_value)  # flatten and re-masked
                xf = numpy.ma.masked_where(zf.mask, x.compressed())
                yf = numpy.ma.masked_where(zf.mask, y.compressed())
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

    @classmethod
    def _cartoplot_colormap(cls,
                            m, M,
                            colormap,
                            colorbounds,
                            plot_method,
                            colorsnumber,
                            colorstep,
                            center_cmap_on_0,
                            contourcolor):
        """Set colors settings."""
        from matplotlib.colors import cnames
        from bronx.graphics import colormapping
        util.load_cmap(colormap)
        if colormap in config.epygram_colormaps_scaling_labels:
            colormaphelper = colormapping.CenteredColormapHelper(
                colormap,
                fixed_colorcenters=config.epygram_colormaps_scaling_labels[colormap])
        elif colormap in config.epygram_colormaps_scaling:
            colormaphelper = colormapping.ColormapHelper(
                colormap,
                fixed_colorbounds=config.epygram_colormaps_scaling[colormap])
        else:
            colormaphelper = colormapping.ColormapHelper(
                colormap,
                fixed_colorbounds=colorbounds)
        plot_kwargs = colormaphelper.kwargs_for_plot((m, M), center_cmap_on_0)
        if plot_method not in ('scatter', 'pcolormesh'):
            plot_kwargs['levels'] = colormaphelper.colorbounds(minmax=(m, M),
                                                               number=colorsnumber,
                                                               step=colorstep)
        ticks = colormaphelper.ticks((m, M))
        if plot_method == 'contour':
            if contourcolor in cnames or contourcolor[0] == '#':
                colormap = None
            else:
                util.load_cmap(contourcolor)
                colormap = contourcolor
                contourcolor = None
            plot_kwargs.update(colors=contourcolor,
                               cmap=colormap)
        return plot_kwargs, ticks

    def _cartoplot_actualplot(self,
                              ax,
                              x, y, data,
                              plot_method,
                              plot_kwargs,
                              uniform,
                              colormap,
                              scatter_kw,
                              contour_kw,
                              contourlabel,
                              clabel_kw):
        """Actual call to matplotlib plotting functions."""
        from matplotlib.colors import cnames
        if plot_method == 'contourf':
            if not self.geometry.rectangular_grid or self.geometry.dimensions['Y'] == 1:
                # triangulate for plotting
                pf = ax.tricontourf(x, y, data,
                                    **plot_kwargs)
            else:
                pf = ax.contourf(x, y, data,
                                 **plot_kwargs)
        elif plot_method == 'scatter':
            # treat uniform field case
            if uniform:
                if colormap in cnames or len(colormap) == 1:
                    colors = colormap
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
            pf = ax.pcolormesh(x, y, data,
                               **plot_kwargs)
        elif plot_method == 'contour':
            if contour_kw is None:
                contour_kw = self.default_contour_kw
            plot_kwargs.update(contour_kw)
            if not self.geometry.rectangular_grid or self.geometry.dimensions['Y'] == 1:
                pf = ax.tricontour(x, y, data,
                                   **plot_kwargs)
                epylog.warning(" ".join(["There is a remaining bug in",
                                         "plotting data from non-rectangular",
                                         "grids in basemap, cf.",
                                         "https://github.com/matplotlib/basemap/issues/265"]))
            else:
                pf = ax.contour(x, y, data,
                                **plot_kwargs)
            if contourlabel:
                if clabel_kw is None:
                    clabel_kw = self.default_clabel_kw
                ax.clabel(pf,
                          colors=plot_kwargs['colors'],
                          cmap=plot_kwargs['cmap'],
                          **clabel_kw)
        else:
            raise NotImplementedError('plot_method=="{}"'.format(plot_method))
        return pf

    @classmethod
    def _cartoplot_colorbar(cls,
                            ax,
                            pf,
                            colorbar,
                            colorbar_over,
                            colorbar_ax_kw,
                            colorbar_ticks,
                            minmax_along_colorbar):
        """Add colorbar."""
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        if colorbar_over is None:
            if colorbar_ax_kw is None:
                colorbar_ax_kw = dict(size="5%", pad=0.2)
            cax = make_axes_locatable(ax).append_axes(colorbar,
                                                      axes_class=plt.Axes,
                                                      **colorbar_ax_kw)
        else:
            cax = colorbar_over
        orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
        cb = plt.colorbar(pf,
                          orientation=orientation,
                          ticks=colorbar_ticks,
                          cax=cax)
        if minmax_along_colorbar:
            cb.set_label(minmax_along_colorbar)

    def _cartoplot_text(self, ax,
                        title,
                        uniform,
                        uniformvalue):
        """Set legend text: title and others."""
        if title is None:
            ax.set_title("\n".join([str(self.fid[sorted(self.fid.keys())[0]]),
                                    str(self.validity.get())]))
        else:
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
                  figsize=None,
                  rcparams=None,
                  title=None,
                  # geometry
                  projection=None,
                  subzone=None,
                  set_global=False,
                  focus_extent=False,
                  # graphical settings
                  plot_method='contourf',
                  minmax=None,
                  mask_threshold=None,
                  scatter_kw=None,
                  contour_kw=None,
                  contourlabel=False,
                  clabel_kw=None,
                  # cartography
                  meridians='auto',
                  parallels='auto',
                  gridlines_kw=None,
                  natural_earth_features=[dict(category='cultural',
                                               name='admin_0_countries',
                                               facecolor='none'),],
                  epygram_departments=False,
                  # colormapping
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
                  minmax_along_colorbar=True,
                  ):
        """
        Plot field with **cartopy**.

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
        :param title: title for the plot. Default is field identifier.

        Geometry settings:
        :param projection: a cartopy.crs projection to be used for plot.
            Defaults to the field.geometry.default_cartopy_CRS()
        :param subzone: [LAM fields only] among ('C', 'CI'), plots the data
            resp. on the C or C+I zone.
            Default is no subzone, i.e. the whole field.
        :param set_global: call cartopy GeoAxes.set_global()
        :param focus_extent: force to focus the map boundaries to the field
            extent. Overrides **set_global**.

        Graphical settings
        :param plot_method: choice of the matplotlib plotting function to be
            used, among ('contourf', 'contour', 'scatter', 'pcolormesh').
        :param minmax: defines the min and max values for the plot colorbar.
            Syntax: [min, max]. Strings 'min' and 'max' (default)
            will take the min and max of the field.
        :param mask_threshold: dict with min and/or max value(s) to mask outside.
        :param scatter_kw: kwargs to be passed to matplotlib's ax.scatter().
            Only for plot_method = 'scatter'.
        :param contour_kw: kwargs to be passed to matplotlib's ax.contour().
            Only for plot_method = 'contour'.
        :param contourlabel: displays labels on contours.
            Only for plot_method = 'contour'.
        :param clabel_kw: kwargs to be passed to matplotlib's ax.clabel().
            Only for plot_method = 'contour'.

        Cartography settings
        :param meridians: enable to fine-tune the choice of lines to
            plot, with either:
              - 'auto': automatic scaling to the basemap extents
              - 'default': every 10 degree
              - a list of values
              - a grid step, e.g. 5 to plot each 5 degree.
              - None: no one is plot
        :param parallels: cf. **meridians**
        :param gridlines_kw: graphical characteristics of meridians/parallels,
            arguments to be passed to cartopy's ax.gridlines(...)
        :param natural_earth_features: list of dicts, each of them containing
            arguments to instanciate a cartopy.feature.NaturalEarthFeature(...).
            E.g. [dict(category='cultural', name='admin_1_states_provinces', facecolor='none', linestyle=':'),]
            will add states/departments/provinces.
            Cf. https://scitools.org.uk/cartopy/docs/latest/matplotlib/feature_interface.html#cartopy.feature.NaturalEarthFeature
            for details.
        :param epygram_departments: add high-resolution french departments
            limits stored in epygram.
            Warning: not consistent with natural_earth_features !

        Colormap settings:
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
        :param minmax_along_colorbar: if True and minmax is not None,
            adds min and max values along colorbar.
        """
        # 1/ geometry and figure
        fig, ax, projection = self._cartoplot_fig_init(fig,
                                                       ax,
                                                       projection,
                                                       figsize,
                                                       rcparams,
                                                       set_global)
        # 2/ background
        self._cartoplot_background(ax,
                                   projection,
                                   natural_earth_features,
                                   meridians,
                                   parallels,
                                   gridlines_kw,
                                   epygram_departments)
        # 3/ get data to plot
        if self.spectral:
            self.sp2gp()
        data = self.getdata(subzone=subzone)
        # 4/ handle min/max values
        data, pmin, pmax, minmax_along_colorbar = self._cartoplot_treat_minmax(data,
                                                                               mask_threshold,
                                                                               minmax,
                                                                               colorbounds,
                                                                               colormap,
                                                                               minmax_along_colorbar)
        if abs(float(pmin) - float(pmax)) < config.epsilon:
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
        # 6/ shape data/lons/lats
        x, y, data = self._cartoplot_shape(x, y,
                                           data,
                                           plot_method)
        if focus_extent:
            ax.set_extent([x.min(), x.max(), y.min(), y.max()], crs=projection)
        # 7/ colormapping
        plot_kwargs, colorbar_ticks = self._cartoplot_colormap(pmin, pmax,
                                                               colormap,
                                                               colorbounds,
                                                               plot_method,
                                                               colorsnumber,
                                                               colorstep,
                                                               center_cmap_on_0,
                                                               contourcolor)
        # 8/ plot
        plot = self._cartoplot_actualplot(ax,
                                          x, y, data,
                                          plot_method,
                                          plot_kwargs,
                                          uniform,
                                          colormap,
                                          scatter_kw,
                                          contour_kw,
                                          contourlabel,
                                          clabel_kw)
        # 9/ colorbar
        if colorbar and plot_method != 'contour':
            self._cartoplot_colorbar(ax,
                                     plot,
                                     colorbar,
                                     colorbar_over,
                                     colorbar_ax_kw,
                                     colorbar_ticks,
                                     minmax_along_colorbar)
        # 10/ texts
        self._cartoplot_text(ax, title, uniform, uniformvalue)
        return fig, ax

    def cartoplot_animation(self, *args, **kwargs):
        return self.animate_plot('cartopy', *args, **kwargs)


class _H2DBasemapPlot(object):
    """
    Plugin for H2DField for plotting with Basemap.

    .. deprecated:: 1.3.9
    """

    def basemap_plot(self,
                     subzone=None,
                     title=None,
                     gisquality='i',
                     specificproj=None,
                     zoom=None,
                     over=(None, None),
                     colorbar_over=None,
                     use_basemap=None,
                     minmax=None,
                     graphicmode='colorshades',
                     levelsnumber=21,
                     colormap='jet',
                     center_cmap_on_0=False,
                     drawrivers=False,
                     drawcoastlines=True,
                     drawcountries=True,
                     meridians='auto',
                     parallels='auto',
                     colorbar='right',
                     minmax_in_title=True,
                     departments=False,
                     boundariescolor='0.25',
                     pointsize=20,
                     contourcolor='blue',
                     contourwidth=1,
                     contourlabel=True,
                     bluemarble=0.0,
                     background=False,
                     mask_threshold=None,
                     contourlabelfmt='%0i',
                     pointsmarker=',',
                     figsize=None,
                     drawmapboundary_kwargs=None,
                     fillcontinents_kwargs=None,
                     drawcoastlines_kwargs=None,
                     drawcountries_kwargs=None,
                     drawparallels_kwargs=None,
                     drawmeridians_kwargs=None,
                     drawequator_kwargs=None,
                     drawgreenwich_kwargs=None,
                     rcparams=None,
                     colorbar_ax_kwargs=None,
                     force_colorbar_ticks_positions=None,
                     force_colorbar_ticks_labels=None):
        """
        Makes a simple plot of the field, with a number of options.

        .. deprecated:: 1.3.9

        :param subzone: among ('C', 'CI'), for LAM fields only, plots the
                        data resp. on the C or C+I zone.
                        Default is no subzone, i.e. the whole field.
        :param gisquality: among ('c', 'l', 'i', 'h', 'f') -- by increasing
                           quality. Defines the quality for GIS elements
                           (coastlines, countries boundaries...).
        :param specificproj: enables to make basemap on the specified projection,
                         among: 'kav7', 'cyl', 'ortho', ('nsper', {...})
                         (cf. Basemap doc). \n
                         In 'nsper' case, the {} may contain: ('sat_height' =
                         satellite height in km; 'lon' = longitude of nadir
                         in degrees; 'lat' = latitude of nadir in degrees.
        :param zoom: specifies the lon/lat borders of the map, implying
                         hereby a 'cyl' projection. Must be a dict(lonmin=,
                         lonmax=, latmin=, latmax=).\n
                         Overwrites *specificproj*.
        :param over: any existing figure and/or ax to be used for the
                         plot, given as a tuple (fig, ax), with None for
                         missing objects. *fig* is the frame of the
                         matplotlib figure, containing eventually several
                         subplots (axes); *ax* is the matplotlib axes on
                         which the drawing is done. When given (is not None),
                         these objects must be coherent, i.e. ax being one of
                         the fig axes.
        :param colorbar_over: an optional existing ax to plot the colorbar on.
        :param use_basemap: a basemap.Basemap object used to handle the
                         projection of the map. If given, the map projection
                         options (*specificproj*, *zoom*, *gisquality* ...)
                         are ignored, keeping the properties of the
                         *use_basemap* object.
        :param title: title for the plot. Default is field identifier.
        :param minmax: defines the min and max values for the plot
                         colorbar. \n
                         Syntax: [min, max]. [0.0, max] also works. Default
                         is min/max of the field.
        :param graphicmode: among ('colorshades', 'contourlines', 'points').
        :param levelsnumber: number of levels for contours and colorbar.
                             (A list of levels also works).
        :param colormap: name of the ``matplotlib`` colormap to use (or an
                         ``epygram`` one, or a user-defined one, cf.
                         config.usercolormaps).
        :param center_cmap_on_0: aligns the colormap center on the value 0.
        :param drawrivers: to add rivers on map.
        :param drawcoastlines: to add coast lines on map.
        :param drawcountries: to add countries on map.
        :param colorbar: if *False*, hide colorbar the plot; else, defines the
                         colorbar position, among ('bottom', 'right').
                         Defaults to 'right'.
        :param meridians: enable to fine-tune the choice of lines to
                         plot, with either:\n
                         - 'auto': automatic scaling to the basemap extents
                         - 'default': range(0,360,10)
                         - a list of values
                         - a grid step, e.g. 5 to plot each 5 degree.
                         - None: no one is plot
                         - *meridians* == 'greenwich' // 'datechange' //
                           'greenwich+datechange' or any combination (,)
                           will plot only these.
        :param parallels: enable to fine-tune the choice of lines to
                         plot, with either:\n
                         - 'auto': automatic scaling to the basemap extents
                         - 'default': range(-90,90,10)
                         - a list of values
                         - a grid step, e.g. 5 to plot each 5 degree.
                         - None: no one is plot
                         - 'equator' // 'polarcircles' // 'tropics'
                           or any combination (,) will plot only these.
        :param minmax_in_title: if True and minmax is not None, adds min and max
                         values in title.
        :param departments: if True, adds the french departments on map (instead
                         of countries).
        :param boundariescolor: color of lines for boundaries (countries,
                         departments, coastlines)
        :param pointsize: size of points for *graphicmode* == 'points'.
        :param contourcolor: color or colormap to be used for 'contourlines'
                         graphicmode. It can be either a legal html color
                         name, or a colormap name.
        :param contourwidth: width of contours for 'contourlines' graphicmode.
        :param contourlabel: displays labels on contours.
        :param bluemarble: if > 0.0 (and <=1.0), displays NASA's "blue marble"
                         as background. The numerical value sets its
                         transparency.
        :param background: if True, set a background color to
                           continents and oceans.
        :param mask_threshold: dict with min and/or max value(s) to mask outside.
        :param contourlabelfmt: format of the contour labels: e.g. 273.15 will
                         appear: '%0i' => 273, '%0f' => 273.150000,
                         '%0.2f' => 273.15, '%04i' => 0273,
                         '%0.5e' => 2.731500e+02
        :param pointsmarker: shape of the points if graphicmode='points'.
                         Cf. matplotlib.scatter() for possible markers.
        :param figsize: figure sizes in inches, e.g. (5, 8.5).
                        Default figsize is config.plotsizes.
        :param drawmapboundary_kwargs: kwargs to be passed to basemap.drawmapboundary()
        :param fillcontinents_kwargs: kwargs to be passed to basemap.fillcontinents()
        :param drawcoastlines_kwargs: kwargs to be passed to basemap.drawcoastlines()
        :param drawcountries_kwargs: kwargs to be passed to basemap.drawcountries()
        :param drawparallels_kwargs: kwargs to be passed to basemap.drawparallels()
        :param drawmeridians_kwargs: kwargs to be passed to basemap.drawgreenwich()
        :param drawequator_kwargs: draw kwargs to emphasize equator parallel
        :param drawgreenwich_kwargs: draw kwargs to emphasize greenwich meridian
        :param rcparams: list of (*args, **kwargs) to be passed to pyplot.rc()
                         defaults to [(('font',), dict(family='serif')),]
        :param colorbar_ax_kwargs: kwargs to be passed to
                                   make_axes_locatable(ax).append_axes(colorbar,
                                                                       **kwargs)
        :param force_colorbar_ticks_position: as a list, or 'center' to center
            it between the color shifting levels of the colormap
        :param force_colorbar_ticks_labels: as a list

        This method uses (hence requires) 'matplotlib' and 'basemap' libraries.
        """
        epylog.info("The 'plotfield' method is deprecated, please use 'cartoplot'")
        # 0. Initializations
        #####################
        # 0.1 matplotlib initializations
        import matplotlib.pyplot as plt
        from matplotlib.colors import cnames
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        if rcparams is None:
            rcparams = [(('font',), dict(family='serif')), ]
        for args, kwargs in rcparams:
            plt.rc(*args, **kwargs)
        if figsize is None:
            figsize = config.plotsizes
        drawcoastlines_kwargs = util.ifNone_emptydict(drawcoastlines_kwargs)
        drawcountries_kwargs = util.ifNone_emptydict(drawcountries_kwargs)
        if colorbar_ax_kwargs is None:
            colorbar_ax_kwargs = dict(size="5%", pad=0.2)

        # 0.2 checkings
        if self.spectral:
            raise epygramError("please convert to gridpoint with sp2gp()" +
                               " method before plotting.")
        if len(self.validity) > 1:
            raise epygramError('to plot H2DField with time dimension, use plotanimation().')

        # 0.3 add custom colormap if necessary
        if colormap not in plt.colormaps():
            if colormap in config.colormaps:
                util.load_cmap(colormap)
            else:
                from mpl_toolkits.basemap import cm
                if colormap in cm.datad:
                    plt.register_cmap(name=colormap, cmap=cm.__dict__[colormap])
        if graphicmode == 'contourlines':  # color of contours
            if contourcolor in cnames or contourcolor[0] == '#':
                colormap = None
            else:
                if contourcolor not in plt.colormaps():
                    util.load_cmap(colormap)
                colormap = contourcolor
                contourcolor = None

        if zoom in (None, {}):  # actual build of figure
            # 1. Figure, ax
            ###############
            fig, ax = set_figax(*over, figsize=figsize)

            # 2. Set up the map
            ###################
            academic = self.geometry.name == 'academic'
            if (academic and use_basemap is not None):
                epylog.warning('*use_basemap* is ignored for academic geometries fields')
            if use_basemap is None and not academic:
                bm = self.geometry.make_basemap(gisquality=gisquality,
                                                subzone=subzone,
                                                specificproj=specificproj,
                                                zoom=zoom)
            elif use_basemap is not None:
                bm = use_basemap
            if not academic:
                if 'color' not in drawcoastlines_kwargs.keys():
                    drawcoastlines_kwargs['color'] = boundariescolor
                if 'color' not in drawcountries_kwargs.keys():
                    drawcountries_kwargs['color'] = boundariescolor
                util.set_map_up(bm, ax,
                                drawrivers=drawrivers,
                                drawcoastlines=drawcoastlines,
                                drawcountries=drawcountries,
                                meridians=meridians,
                                parallels=parallels,
                                departments=departments,
                                bluemarble=bluemarble,
                                background=background,
                                drawmapboundary_kwargs=drawmapboundary_kwargs,
                                fillcontinents_kwargs=fillcontinents_kwargs,
                                drawcoastlines_kwargs=drawcoastlines_kwargs,
                                drawcountries_kwargs=drawcountries_kwargs,
                                drawparallels_kwargs=drawparallels_kwargs,
                                drawmeridians_kwargs=drawmeridians_kwargs,
                                drawequator_kwargs=drawequator_kwargs,
                                drawgreenwich_kwargs=drawgreenwich_kwargs)

            # 3. Prepare data
            #################
            if academic:
                (lons, lats) = self.geometry.get_lonlat_grid(subzone=subzone)
                x = (lons - lons.min()) * self.geometry.grid['X_resolution']
                y = (lats - lats.min()) * self.geometry.grid['Y_resolution']
            else:
                (lons, lats) = self.geometry.get_lonlat_grid(subzone=subzone)
                x, y = bm(lons, lats)
            # mask values
            mask_outside = {'min':-config.mask_outside,
                            'max':config.mask_outside}
            if mask_threshold is not None:
                mask_outside.update(mask_threshold)
                if not self.geometry.rectangular_grid and graphicmode != 'points':
                    epylog.warning("shift to *graphicmode*='points' for mask_threshold to be accounted for in unrectangular_grid.")
                    graphicmode = 'points'
            data = numpy.ma.masked_outside(self.getdata(subzone=subzone),
                                           mask_outside['min'],
                                           mask_outside['max'])

            # 4. Plot configuration
            #######################
            # handle min/max values
            m = data.min()
            M = data.max()
            if minmax_in_title and (minmax is not None or
                                    colormap in config.colormaps_scaling or
                                    isinstance(levelsnumber, list)):
                minmax_in_title = '(min: ' + \
                                  '{: .{precision}{type}}'.format(m, type='E', precision=3) + \
                                  ' // max: ' + \
                                  '{: .{precision}{type}}'.format(M, type='E', precision=3) + ')'
            else:
                minmax_in_title = ''
            if minmax is not None:
                try:
                    m = float(minmax[0])
                except ValueError:
                    m = data.min()
                try:
                    M = float(minmax[1])
                except ValueError:
                    M = data.max()
            if isinstance(levelsnumber, list):
                if minmax is not None:
                    epylog.warning('**minmax** overwritten by **levelsnumber**')
                m = min(levelsnumber)
                M = max(levelsnumber)

            if abs(float(m) - float(M)) < config.epsilon:
                epylog.warning("uniform field: plot as 'points'.")
                graphicmode = 'points'
                uniform = True
            else:
                uniform = False
            if center_cmap_on_0:
                vmax = max(abs(m), M)
                vmin = -vmax
            else:
                vmin = m
                vmax = M
            # set levels and ticks levels
            if isinstance(levelsnumber, list):
                levels = levelsnumber
                ticks_positions = levelsnumber
            else:
                levels = numpy.linspace(m, M, levelsnumber)
                L = int((levelsnumber - 1) // 15) + 1
                ticks_positions = [levels[l]
                                   for l in range(len(levels) - (L // 3 + 1))
                                   if l % L == 0] + [levels[-1]]
            if colormap in config.colormaps_scaling:
                (norm, levels) = util.scale_colormap(colormap)
                if not isinstance(levelsnumber, list):
                    ticks_positions = levels
                vmin = vmax = None
            else:
                norm = None
            ticks_labels = None
            if force_colorbar_ticks_labels is not None:
                ticks_labels = [str(l) for l in force_colorbar_ticks_labels]
            if force_colorbar_ticks_positions is not None:
                if isinstance(force_colorbar_ticks_positions, list):
                    ticks_positions = force_colorbar_ticks_positions
                elif force_colorbar_ticks_positions == 'center':
                    ticks_positions = [(levels[i] + levels[i + 1]) / 2. for i in range(len(levels) - 1)]
            if ticks_labels is not None:
                assert len(ticks_labels) == len(ticks_positions), \
                    str(len(ticks_labels)) + '!=' + str(len(ticks_positions))
            # 5. Plot
            #########
            if graphicmode == 'colorshades':
                if not self.geometry.rectangular_grid:
                    xf = numpy.ma.masked_where(data.mask, x).compressed()
                    yf = numpy.ma.masked_where(data.mask, y).compressed()
                    zf = data.compressed()
                    tri = True
                elif self.geometry.dimensions['Y'] == 1:
                    xf = x.flatten()
                    yf = y.flatten()
                    zf = data.flatten()
                    tri = True
                else:
                    xf = x
                    yf = y
                    zf = data
                    tri = False
                plot_kwargs = dict(cmap=colormap,
                                   vmin=vmin, vmax=vmax,
                                   norm=norm,
                                   tri=tri)
                if academic:
                    pf = ax.contourf(xf, yf, zf, levels,
                                     **plot_kwargs)
                else:
                    pf = bm.contourf(xf, yf, zf, levels, ax=ax,
                                     **plot_kwargs)
                if colorbar:
                    if colorbar_over is None:
                        cax = make_axes_locatable(ax).append_axes(colorbar,
                                                                  **colorbar_ax_kwargs)
                    else:
                        cax = colorbar_over
                    orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
                    cb = plt.colorbar(pf,
                                      orientation=orientation,
                                      ticks=ticks_positions,
                                      cax=cax)
                    if ticks_labels is not None:
                        cax.set_yticklabels(ticks_labels)
                    if minmax_in_title:
                        cb.set_label(minmax_in_title)
            elif graphicmode == 'contourlines':
                if not self.geometry.rectangular_grid:
                    xf = x.compressed()
                    yf = y.compressed()
                    zf = data.compressed()
                    tri = True
                    epylog.warning(" ".join(["There is a remaining bug in",
                                             "plotting data from non-rectangular",
                                             "grids in basemap, cf.",
                                             "https://github.com/matplotlib/basemap/issues/265"]))
                elif self.geometry.dimensions['Y'] == 1:
                    xf = x.flatten()
                    yf = y.flatten()
                    zf = data.flatten()
                    tri = True
                    epylog.warning(" ".join(["There is a remaining bug in",
                                             "plotting data from non-rectangular",
                                             "grids in basemap, cf.",
                                             "https://github.com/matplotlib/basemap/issues/265"]))
                else:
                    xf = x
                    yf = y
                    zf = data
                    tri = False
                plot_kwargs = dict(levels=levels,
                                   colors=contourcolor,
                                   cmap=colormap,
                                   linewidths=contourwidth,
                                   tri=tri)
                if academic:
                    pf = ax.contour(xf, yf, zf,
                                    **plot_kwargs)
                else:
                    pf = bm.contour(xf, yf, zf, ax=ax,
                                    **plot_kwargs)
                if contourlabel:
                    ax.clabel(pf, colors=contourcolor, cmap=colormap,
                              fmt=contourlabelfmt)
            elif graphicmode == 'points':
                xf = x.flatten()
                yf = y.flatten()
                zf = numpy.ma.masked_outside(data.flatten(), m, M)
                if uniform:
                    if colormap in cnames or len(colormap) == 1:
                        zf = colormap
                    else:
                        zf = 'seagreen'
                plot_kwargs = dict(s=pointsize,
                                   norm=norm,
                                   marker=pointsmarker,
                                   linewidths=0,
                                   cmap=colormap,
                                   vmin=vmin, vmax=vmax)
                if academic:
                    pf = ax.scatter(xf, yf, c=zf,
                                    **plot_kwargs)
                else:
                    pf = bm.scatter(xf, yf, c=zf, ax=ax,
                                    **plot_kwargs)
                if colorbar and not uniform:
                    if colorbar_over is None:
                        cax = make_axes_locatable(ax).append_axes(colorbar,
                                                                  **colorbar_ax_kwargs)
                    else:
                        cax = colorbar_over
                    orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
                    cb = plt.colorbar(pf,
                                      orientation=orientation,
                                      ticks=ticks_positions,
                                      cax=cax)
                    if ticks_labels is not None:
                        cax.set_yticklabels(ticks_labels)
                    if minmax_in_title != '':
                        cb.set_label(minmax_in_title)
                elif uniform:
                    ax.text(1.02, 0.5, '(uniform field)',
                            horizontalalignment='center',
                            verticalalignment='center',
                            transform=ax.transAxes,
                            rotation=90.)
            if title is None:
                ax.set_title("\n".join([str(self.fid[sorted(self.fid.keys())[0]]),
                                        str(self.validity.get())]))
            else:
                ax.set_title(title)
        else:
            # zoom: create zoom_field and plot it
            zoom_field = self.extract_zoom(zoom,
                                           # in regLL case, to be sure gridpoints of the border are included
                                           extra_10th=use_basemap is not None)
            # get args
            import inspect
            frame = inspect.currentframe()
            args, _, _, values = inspect.getargvalues(frame)
            kwargs = {k:values[k] for k in args if k != 'self'}
            kwargs.pop('zoom')
            kwargs.pop('subzone')
            # forward plot
            fig, ax = zoom_field.plotfield(**kwargs)

        return (fig, ax)

    def plotanimation(self, *args, **kwargs):
        return self.animate_plot('basemap', *args, **kwargs)


class H2DField(_H2DBasemapPlot, _H2DCartopyPlot, D3Field):
    """
    Horizontal 2-Dimensions field class.
    A field is defined by its identifier 'fid',
    its data, its geometry (gridpoint and optionally spectral),
    and its validity.

    The natural being of a field is gridpoint, so that:
    a field always has a gridpoint geometry, but it has a spectral geometry only
    in case it is spectral.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['H2D'])),
            geometry=dict(
                type=H2DGeometry),
        )
    )

##############
# ABOUT DATA #
##############
    def getlevel(self, level=None, k=None):
        """Returns self. Useless but for compatibility reasons."""
        if k is None and level is None:
            raise epygramError("You must give k or level.")
        if k is not None and level is not None:
            raise epygramError("You cannot give, at the same time, k and level")
        if level is not None:
            if level not in self.geometry.vcoordinate.levels:
                raise epygramError("The requested level does not exist.")
        return self

    def extract_point(self, lon, lat,
                      interpolation='nearest',
                      external_distance=None):
        """
        Extract a point as a PointField.

        Cf. getvalue_ll() doc for other arguments.
        """
        value = [self.getvalue_ll(lon, lat, validity=self.validity[i],
                                  interpolation=interpolation,
                                  external_distance=external_distance)
                 for i in range(len(self.validity))]

        pt = gimme_one_point(lon, lat,
                             field_args={'validity':self.validity.deepcopy(),
                                         'fid':self.fid},
                             geometry_args={'vcoordinate':self.geometry.vcoordinate},
                             vertical_geometry_args=None)
        pt.setdata(value)
        return pt

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def plotfield(self, *args, **kwargs):
        return self.basemap_plot(*args, **kwargs)

    def animate_plot(self,
                     maplib,
                     title='__auto__',
                     repeat=False,
                     interval=1000,
                     **kwargs):
        """
        Plot the field with animation with regards to time dimension.
        Returns a :class:`matplotlib.animation.FuncAnimation`.

        In addition to those specified below, all plot method (cartoplot/basemap_plot)
        arguments can be provided.

        :param maplib: map library to be used: cartopy/basemap
        :param title: title for the plot. '__auto__' (default) will print
          the current validity of the time frame.
        :param repeat: to repeat animation
        :param interval: number of milliseconds between two validities
        """
        import matplotlib.animation as animation
        if len(self.validity) == 1:
            raise epygramError("plotanimation can handle only field with several validities.")
        # title
        if title is not None:
            if title == '__auto__':
                title_prefix = ''
            else:
                title_prefix = title
            title = title_prefix + '\n' + self.validity[0].get().isoformat(sep=b' ')
        else:
            title_prefix = None
        # initialisation
        field0 = self.deepcopy()
        field0.validity = self.validity[0]
        field0.setdata(self.getdata()[0, ...])
        mindata = self.getdata(subzone=kwargs.get('subzone')).min()
        maxdata = self.getdata(subzone=kwargs.get('subzone')).max()
        kwargs.setdefault(minmax=(mindata, maxdata))
        if maplib == 'basemap':
            if self.geometry.name != 'academic':
                bm = kwargs.get('use_basemap')
                if bm is None:
                    bm = self.geometry.make_basemap(gisquality=kwargs.get('gisquality', 'i'),
                                                    subzone=kwargs.get('subzone'),
                                                    specificproj=kwargs.get('specificproj'),
                                                    zoom=kwargs.get('zoom'))
                kwargs['use_basemap'] = bm
            fig, ax = field0.basemap_plot(title=title,
                                          **kwargs)
        elif maplib == 'cartopy':
            kwargs.setdefault(projection=self.geometry.default_cartopy_CRS())
            fig, ax = field0.cartoplot(title=title,
                                       **kwargs)
        if kwargs.get('colorbar_over') is None:
            kwargs['colorbar_over'] = fig.axes[-1]  # the last being created, in plotfield()
        kwargs['over'] = (fig, ax)

        # update function, step after step
        def update(i, ax, myself, fieldi, title_prefix, kwargs):
            if i < len(myself.validity):
                ax.clear()
                fieldi.validity = myself.validity[i]
                fieldi.setdata(myself.getdata()[i, ...])
                if title_prefix is not None:
                    title = title_prefix + '\n' + fieldi.validity.get().isoformat(sep=b' ')
                if maplib == 'basemap':
                    fieldi.basemap_plot(title=title,
                                        **kwargs)
                elif maplib == 'cartopy':
                    fieldi.cartoplot(title=title,
                                     **kwargs)
        # build animation
        anim = animation.FuncAnimation(fig, update,
                                       fargs=[ax, self, field0, title_prefix, kwargs],
                                       frames=list(range(len(self.validity) + 1)),  # AM: don't really understand why but needed for the last frame to be shown
                                       interval=interval,
                                       repeat=repeat)
        return anim

    def morph_with_points(self, points,
                          alpha=1.,
                          morphing='nearest',
                          increment=False,
                          **kwargs):
        """
        Perturb the field values with the values of a set of points.

        :param points: a list/fieldset of PointField.
        :param alpha: is the blending value, ranging from 0. to 1.:
          e.g. with 'nearest' morphing, the final value of the point is
          alpha*point.value + (1-alpha)*original_field.value
        :param morphing: is the way the point modify the field:\n
          - 'nearest': only the nearest point is modified
          - 'exp_decay': modifies the surrounding points with an isotrop
            exponential decay weighting
          - 'gaussian': modifies the surrounding points with an isotrop
            gaussian weighting. Its standard deviation *sigma* must then
            be passed as argument, in meters. Weightingassumed to be 0. from
            3*sigma.
        :param increment: if True, the final value of the point is
          original_field.value + point.value
        """
        assert isinstance(points, list)
        assert all([isinstance(p, PointField) for p in points])
        if len(self.validity) != 1:
            raise NotImplementedError('morph field with a time dimension.')

        for p in points:
            if morphing == 'nearest':
                lon = p.geometry.grid['longitudes'][0]
                lat = p.geometry.grid['latitudes'][0]
                i, j = self.geometry.ll2ij(lon, lat)
                i = int(i)
                j = int(j)
                if increment:
                    self._data[:, :, j, i] = self._data[:, :, j, i] + p.getdata()
                else:
                    self._data[:, :, j, i] = alpha * p.getdata() + (1. - alpha) * self._data[:, :, j, i]
            elif morphing == 'gaussian':
                sigma = kwargs['sigma']
                lon = p.geometry.grid['longitudes'][0]
                lat = p.geometry.grid['latitudes'][0]
                zero_radius = 3. * sigma
                selection_points_ij = self.geometry.nearest_points(lon, lat, interpolation=('square:radius', zero_radius))
                selection_points_ll = self.geometry.ij2ll([sp[0] for sp in selection_points_ij], [sp[1] for sp in selection_points_ij])
                selection_points_ll = [(selection_points_ll[0][k], selection_points_ll[1][k]) for k in range(len(selection_points_ll[0]))]
                for sp in selection_points_ll:
                    distance = self.geometry.distance((lon, lat), sp)
                    if distance <= zero_radius:
                        (i, j) = selection_points_ij[selection_points_ll.index(sp)]
                        if increment:
                            local_alpha = numpy.exp(-distance ** 2. / (2. * sigma ** 2))
                            self._data[:, :, j, i] = self._data[:, :, j, i] + local_alpha * p.getdata()
                        else:
                            local_alpha = alpha * numpy.exp(-distance ** 2. / (2. * sigma ** 2))
                            self._data[:, :, j, i] = local_alpha * p.getdata() + (1. - local_alpha) * self._data[:, :, j, i]
            else:
                raise NotImplementedError("morphing: " + morphing + " :not yet.")

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """
        if isinstance(other, self.__class__):
            if self.spectral != other.spectral:
                raise epygramError("cannot operate a spectral field with a" +
                                   " non-spectral field.")
            if self.geometry.rectangular_grid:
                if self.geometry.dimensions['X'] != other.geometry.dimensions['X'] or\
                   self.geometry.dimensions['Y'] != other.geometry.dimensions['Y']:
                    raise epygramError("operations on fields cannot be done if" +
                                       " fields do not share their gridpoint" +
                                       " dimensions.")
            else:
                if self.geometry.dimensions != other.geometry.dimensions:
                    raise epygramError("operations on fields cannot be done if" +
                                       " fields do not share their gridpoint" +
                                       " dimensions.")
            if self.spectral_geometry != other.spectral_geometry:
                raise epygramError("operations on fields cannot be done if" +
                                   " fields do not share their spectral" +
                                   " geometry.")
        else:
            super(D3Field, self)._check_operands(other)
