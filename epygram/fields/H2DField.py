#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import numpy
import copy

import footprints

from epygram import config, util, epygramError
from epygram.geometries import H2DGeometry
from . import gimme_one_point
from .D3Field import D3Field
from .PointField import PointField
from ._BasemapFieldsPlot import _H2DBasemapPlot

epylog = footprints.loggers.getLogger(__name__)


class _H2DCartopyPlot(object):
    """
    Plugin for H2DField for plotting with Cartopy.
    """
    # defaults arguments for cartopy plots
    default_NEfeatures = [dict(category='cultural',
                               name='admin_0_countries',
                               facecolor='none',
                               edgecolor='k'),]
    default_scatter_kw = {'s':20,
                          'marker':',',
                          'linewidths':0}
    default_contour_kw = {'linewidths':1}
    default_clabel_kw = {'fmt':'%0i'}
    default_gridlines_kw = {'draw_labels':True,
                            'linewidth':1,
                            'color':'k',
                            'linestyle':'--'}

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

    @classmethod
    def cartoplot_background(cls,
                             geometry,
                             ax,
                             projection,
                             cartopy_features=[],
                             natural_earth_features=default_NEfeatures,
                             meridians='auto',
                             parallels='auto',
                             gridlines_kw=None,
                             epygram_departments=False):
        """Set cartography features, such as borders, coastlines, meridians and parallels..."""
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from cartopy.mpl.gridliner import (LATITUDE_FORMATTER,
                                           LONGITUDE_FORMATTER)
        if 'gauss' in geometry.name:
            default_scale = '110m'
        else:
            default_scale = '50m'
        for f in cartopy_features:
            if isinstance(f, six.string_types):
                f = getattr(cfeature, f)
            ax.add_feature(f.with_scale(default_scale))
        for f in natural_earth_features:
            f = copy.copy(f)
            if 'scale' not in f:
                f['scale'] = default_scale
            ax.add_feature(cfeature.NaturalEarthFeature(**f))
        # meridians and parallels
        meridians, parallels = util.auto_meridians_parallels(geometry,
                                                             meridians,
                                                             parallels)
        if gridlines_kw is None:
            gridlines_kw = copy.copy(cls.default_gridlines_kw)
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

    @property
    def _meshable(self):
        return self.geometry.rectangular_grid and self.geometry.name != 'unstructured'

    def _cartoplot_mesh_coords(self, subzone):
        """Get coordinates of mesh corners."""
        if self._meshable:
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
            raise NotImplementedError("plot_method='pcolormesh' for Gauss or unstructured grids. Any idea ?")  # TODO: ?
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
        import matplotlib.pyplot as plt
        from matplotlib.colors import cnames
        from epygram import colormapping
        if colormap not in plt.colormaps() and colormap in config.colormaps:
            cmapfile = config.colormaps[colormap]
            if cmapfile.endswith('.json'):
                colormaphelper = colormapping.get_ColormapHelper_fromfile(cmapfile)
            elif cmapfile.endswith('.cmap'):
                # deprecated
                raise epygramError(util._deprecated_cmap)
        else:
            colormaphelper = colormapping.ColormapHelper(colormap,
                                                         explicit_colorbounds=colorbounds)
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
        return plot_kwargs, colormaphelper

    def _cartoplot_actualplot(self,
                              ax,
                              x, y,
                              data,
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
                            colormap_helper,
                            minmax_along_colorbar,
                            m, M,
                            colorsnumber,
                            colorstep):
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
                          cax=cax)
        if ticks_label != ticks_position:
            cax.set_yticklabels(ticks_label)
        if minmax_along_colorbar:
            cb.set_label(minmax_along_colorbar)

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
                  set_global=False,
                  focus_extent=True,
                  # graphical settings
                  plot_method='__default__',
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
                  cartopy_features=[],
                  natural_earth_features=default_NEfeatures,
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
        :param title: title for the plot. Default is field identifier and validity.

        Geometry settings:

        :param projection: a cartopy.crs projection to be used for plot.
            Defaults to the field.geometry.default_cartopy_CRS()
        :param subzone: [LAM fields only] among ('C', 'CI'), plots the data
            resp. on the C or C+I zone.
            Default is no subzone, i.e. the whole field.
        :param set_global: call cartopy GeoAxes.set_global().
            Overrides **focus_extent**.
        :param focus_extent: force to focus the map boundaries to the field
            extent. Otherwise, let matplotlib decide of the map boundaries.

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
        :param contourlabel: displays labels on contours.
            Only for plot_method = 'contour'.
        :param clabel_kw: kwargs to be passed to matplotlib's ax.clabel().
            Only for plot_method = 'contour'.

        Cartography settings:

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
        # 0/ defaults pre-sets
        if plot_method == '__default__':
            if self._meshable:
                plot_method = 'pcolormesh'
            else:
                plot_method = 'contourf'
        # 1/ geometry and figure
        fig, ax, projection = self.cartoplot_fig_init(fig,
                                                      ax,
                                                      projection,
                                                      figsize,
                                                      rcparams,
                                                      set_global)
        # 2/ background
        self.cartoplot_background(self.geometry,
                                  ax,
                                  projection,
                                  cartopy_features,
                                  natural_earth_features,
                                  meridians,
                                  parallels,
                                  gridlines_kw,
                                  epygram_departments)
        # 3/ get data to plot
        if self.spectral:
            self.sp2gp()
        if not self.geometry.grid.get('LAMzone', False):
            subzone = None
        data = self.getdata(subzone=subzone)
        # 4/ handle min/max values
        data, pmin, pmax, minmax_along_colorbar = self._cartoplot_treat_minmax(data,
                                                                               mask_threshold,
                                                                               minmax,
                                                                               colorbounds,
                                                                               colormap,
                                                                               minmax_along_colorbar)
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
        assert numpy.inf not in [x.min(), x.max(), y.min(), y.max()], \
            "Domain is too large to be plotted on space view"
        assert -numpy.inf not in [x.min(), x.max(), y.min(), y.max()], \
            "Domain is too large to be plotted on space view"
        # 6/ shape data/lons/lats
        x, y, data = self._cartoplot_shape(x, y,
                                           data,
                                           plot_method)
        if focus_extent and not set_global and 'gauss' not in self.geometry.name:
            ax.set_extent([x.min(), x.max(), y.min(), y.max()],
                          crs=projection)
        # 7/ colormapping
        plot_kwargs, colormap_helper = self._cartoplot_colormap(pmin, pmax,
                                                                colormap,
                                                                colorbounds,
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
                                                  colormap,
                                                  scatter_kw,
                                                  contour_kw,
                                                  contourlabel,
                                                  clabel_kw)
            # 9/ colorbar
            if colorbar and plot_method != 'contour':
                self._cartoplot_colorbar(ax,
                                         elements,
                                         colorbar,
                                         colorbar_over,
                                         colorbar_ax_kw,
                                         colormap_helper,
                                         minmax_along_colorbar,
                                         pmin, pmax,
                                         colorsnumber,
                                         colorstep)
        # 10/ texts
        self._cartoplot_text(ax, title, uniform, uniformvalue)
        return fig, ax


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
        #return self.cartoplot(*args, **kwargs)
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

        .. deprecated:: 1.3.11
        Animations shall be made externally. Will not be maintained.

        In addition to those specified below, all plot method (cartoplot/basemap_plot)
        arguments can be provided.

        :param maplib: map library to be used: 'cartopy'/'basemap' (basemap is deprecated!)
        :param title: title for the plot. '__auto__' (default) will print
          the current validity of the time frame.
        :param repeat: to repeat animation
        :param interval: number of milliseconds between two validities
        """
        import matplotlib.animation as animation
        epylog.info('deprecated:: 1.3.11 - Animations shall be made externally. Will not be maintained.')
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
