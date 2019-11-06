#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the Basemap-interfaced facilities to plot field.

.. deprecated:: 1.3.11
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import numpy

import footprints
from bronx.graphics.axes import set_figax
from bronx.syntax.arrays import stretch_array

from epygram import config, util, epygramError

epylog = footprints.loggers.getLogger(__name__)


class _H2DBasemapPlot(object):
    """
    .. deprecated:: 1.3.11

    Plugin for H2DField for plotting with Basemap.
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
        epylog.warning("The 'basemap_plot' method is deprecated, and will be removed in future versions of epygram; please use 'cartoplot'")
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


class _H2DVectorBasemapPlot(object):

    def basemap_plot(self,
                     over=(None, None),
                     subzone=None,
                     title=None,
                     gisquality='i',
                     specificproj=None,
                     zoom=None,
                     use_basemap=None,
                     drawcoastlines=True,
                     drawcountries=True,
                     drawrivers=False,
                     departments=False,
                     boundariescolor='0.25',
                     parallels='auto',
                     meridians='auto',
                     subsampling=1,
                     symbol='barbs',
                     symbol_options={'color':'k', },
                     plot_module=True,
                     plot_module_options=None,
                     bluemarble=0.,
                     background=False,
                     quiverkey=None,
                     quiver_options=None,
                     components_are_projected_on='grid',
                     map_factor_correction=True,
                     mask_threshold=None,
                     figsize=None,
                     drawmapboundary_kwargs=None,
                     fillcontinents_kwargs=None,
                     drawcoastlines_kwargs=None,
                     drawcountries_kwargs=None,
                     drawparallels_kwargs=None,
                     drawmeridians_kwargs=None,
                     drawequator_kwargs=None,
                     drawgreenwich_kwargs=None,
                     rcparams=None):
        """
        Makes a simple plot of the field, with a number of options.

        .. deprecated:: 1.3.11

        :param over: to plot the vectors over an existing figure
          (e.g. colorshades).
          Any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several
          subplots (axes); *ax* is the matplotlib axes on
          which the drawing is done. When given (is not None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        :param subzone: among ('C', 'CI'), for LAM fields only, plots the data
          resp. on the C or C+I zone. \n
          Default is no subzone, i.e. the whole field.
        :param gisquality: among ('c', 'l', 'i', 'h', 'f') -- by increasing
          quality. Defines the quality for GIS elements (coastlines, countries
          boundaries...). Default is 'i'. Cf. 'basemap' doc for more details.
        :param specificproj: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        :param zoom: specifies the lon/lat borders of the map, implying hereby
          a 'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=).\n
          Overwrites *specificproj*.
        :param use_basemap: a basemap.Basemap object used to handle the
          projection of the map. If given, the map projection
          options (*specificproj*, *zoom*, *gisquality* ...)
          are ignored, keeping the properties of the
          *use_basemap* object. (because making Basemap is the most
          time-consuming step).
        :param drawrivers: to add rivers on map.
        :param departments: if True, adds the french departments on map (instead
          of countries).
        :param boundariescolor: color of lines for boundaries (countries,
          departments, coastlines)
        :param drawcoastlines: to add coast lines on map.
        :param drawcountries: to add countries on map.
        :param title: title for the plot. Default is field identifier.
        :param meridians: enable to fine-tune the choice of lines to
          plot, with either:\n
          - 'auto': automatic scaling to the basemap extents
          - 'default': range(0,360,10)
          - a list of values
          - a grid step, e.g. 5 to plot each 5 degree.
          - None: no one is plot
          - *meridian* == 'greenwich' // 'datechange' // 'greenwich+datechange'
            combination (,) will plot only these.
        :param parallels: enable to fine-tune the choice of lines to
          plot, with either:\n
          - 'auto': automatic scaling to the basemap extents
          - 'default': range(-90,90,10)
          - a list of values
          - a grid step, e.g. 5 to plot each 5 degree.
          - None: no one is plot
          - 'equator' // 'polarcircles' // 'tropics' or any
            combination (,) will plot only these.
        :param subsampling: to subsample the number of gridpoints to plot.
          Ex: *subsampling* = 10 will only plot one gridpoint upon 10.
        :param symbol: among ('barbs', 'arrows', 'stream')
        :param symbol_options: a dict of options to be passed to **barbs** or
          **quiver** method.
        :param plot_module: to plot module as colorshades behind vectors.
        :param plot_module_options: options (dict) to be passed to module.plotfield().
        :param bluemarble: if > 0.0 (and <=1.0), displays NASA's "blue marble"
          as background. The numerical value sets its transparency.
        :param background: if True, set a background color to
          continents and oceans.
        :param quiverkey: to activate quiverkey; must contain arguments to be
          passed to pyplot.quiverkey(), as a dict.
        :param components_are_projected_on: inform the plot on which axes the
          vector components are projected on ('grid' or 'lonlat').
        :param map_factor_correction: if True, applies a correction of magnitude
          to vector due to map factor.
        :param mask_threshold: dict with min and/or max value(s) to mask outside.
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

        This method uses (hence requires) 'matplotlib' and 'basemap' libraries.
        """
        epylog.info('basemap is deprecated ! use cartoplot() method.')
        import matplotlib.pyplot as plt
        if rcparams is None:
            rcparams = [(('font',), {'family':'serif'}),]
        for args, kwargs in rcparams:
            plt.rc(*args, **kwargs)
        if figsize is None:
            figsize = config.plotsizes
        drawcoastlines_kwargs = util.ifNone_emptydict(drawcoastlines_kwargs)
        drawcountries_kwargs = util.ifNone_emptydict(drawcountries_kwargs)

        plot_module_options = util.ifNone_emptydict(plot_module_options)
        quiver_options = util.ifNone_emptydict(quiver_options)

        if self.spectral:
            raise epygramError("please convert to gridpoint with sp2gp()" +
                               " method before plotting.")

        # 1. Figure, ax
        if not plot_module:
            fig, ax = set_figax(*over, figsize=figsize)

        # 2. Set up the map
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
        elif academic:
            raise NotImplementedError('plot VectorField in academic geometry')
            bm = None
        if not academic:
            if plot_module:
                module = self.to_module()
                if 'gauss' in self.geometry.name and self.geometry.grid['dilatation_coef'] != 1.:
                    if map_factor_correction:
                        module.operation_with_other('*', self.geometry.map_factor_field())
                    else:
                        epylog.warning('check carefully *map_factor_correction* w.r.t. dilatation_coef')
                fig, ax = module.plotfield(use_basemap=bm,
                                           over=over,
                                           subzone=subzone,
                                           specificproj=specificproj,
                                           title=title,
                                           drawrivers=drawrivers,
                                           drawcoastlines=drawcoastlines,
                                           drawcountries=drawcountries,
                                           meridians=meridians,
                                           parallels=parallels,
                                           departments=departments,
                                           boundariescolor=boundariescolor,
                                           bluemarble=bluemarble,
                                           background=background,
                                           **plot_module_options)
            else:
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
        # mask values
        mask_outside = {'min':-config.mask_outside,
                        'max':config.mask_outside}
        if mask_threshold is not None:
            mask_outside.update(mask_threshold)
        data = [numpy.ma.masked_outside(data,
                                        mask_outside['min'],
                                        mask_outside['max']) for data in
                self.getdata(subzone=subzone)]
        if data[0].ndim == 1:  # self.geometry.dimensions['Y'] == 1:
            u = data[0][::subsampling]
            v = data[1][::subsampling]
        else:
            u = data[0][::subsampling, ::subsampling]
            v = data[1][::subsampling, ::subsampling]
        (lons, lats) = self.geometry.get_lonlat_grid(subzone=subzone)
        if lons.ndim == 1:  # self.geometry.dimensions['Y'] == 1:
            lons = lons[::subsampling]
            lats = lats[::subsampling]
        else:
            lons = lons[::subsampling, ::subsampling]
            lats = lats[::subsampling, ::subsampling]
        if isinstance(u, numpy.ma.masked_array) \
        or isinstance(v, numpy.ma.masked_array):
            assert isinstance(u, numpy.ma.masked_array) == isinstance(u, numpy.ma.masked_array)
            common_mask = u.mask + v.mask
            u.mask = common_mask
            v.mask = common_mask
            lons = numpy.ma.masked_where(common_mask, lons)
            lats = numpy.ma.masked_where(common_mask, lats)
        x, y = bm(lons, lats)

        # Calculate the orientation of the vectors
        assert components_are_projected_on in ('grid', 'lonlat')
        if components_are_projected_on == 'grid' and 'gauss' not in self.geometry.name \
           and (specificproj is None and zoom is None):
            # map has same projection than components: no rotation necessary
            u_map = u
            v_map = v
        else:
            # (1or2) rotation(s) is(are) necessary
            if components_are_projected_on == 'lonlat' or self.geometry.name == 'regular_lonlat':
                (u_ll, v_ll) = (stretch_array(u), stretch_array(v))
            else:
                # wind is projected on a grid that is not lonlat: rotate to lonlat
                (u_ll, v_ll) = self.geometry.reproject_wind_on_lonlat(stretch_array(u), stretch_array(v),
                                                                      stretch_array(lons), stretch_array(lats),
                                                                      map_factor_correction=map_factor_correction)
            # rotate from lonlat to map projection
            (u_map, v_map) = bm.rotate_vector(u_ll,
                                              v_ll,
                                              stretch_array(lons),
                                              stretch_array(lats))
            # go back to 2D if necessary
            if symbol == 'stream':
                u_map = u_map.reshape(u.shape)
                v_map = v_map.reshape(v.shape)

        if symbol == 'stream':
            if self.geometry.rectangular_grid:
                xf = x[0, :]  # in basemap space, x is constant on a column
                yf = y[:, 0]  # in basemap space, y is constant on a row
                u = u_map
                v = v_map
                speed_width = 2 * numpy.sqrt(u ** 2 + v ** 2) / min(u.max(), v.max())
            else:
                raise NotImplementedError("matplotlib's streamplot need an evenly spaced grid.")
        else:
                xf = stretch_array(x)
                yf = stretch_array(y)
                u = stretch_array(u_map)
                v = stretch_array(v_map)
        if symbol == 'barbs':
            bm.barbs(xf, yf, u, v, ax=ax, **symbol_options)
        elif symbol == 'arrows':
            q = bm.quiver(xf, yf, u, v, ax=ax, **symbol_options)
            if quiverkey:
                ax.quiverkey(q, **quiverkey)
        elif symbol == 'stream':
            bm.streamplot(xf, yf, u, v, ax=ax, linewidth=speed_width, **symbol_options)
        if title is None:
            ax.set_title(str(self.fid) + "\n" + str(self.validity.get()))
        else:
            ax.set_title(title)

        return (fig, ax)

    def plotanimation(self,
                      title='__auto__',
                      repeat=False,
                      interval=1000,
                      **kwargs):
        """
        Plot the field with animation with regards to time dimension.
        Returns a :class:`matplotlib.animation.FuncAnimation`.

        .. deprecated:: 1.3.11

        In addition to those specified below, all :meth:`plotfield` method
        arguments can be provided.

        :param title: title for the plot. '__auto__' (default) will print
          the current validity of the time frame.
        :param repeat: to repeat animation
        :param interval: number of milliseconds between two validities
        """
        import matplotlib.animation as animation

        if len(self.validity) == 1:
            raise epygramError("plotanimation can handle only field with several validities.")

        if title is not None:
            if title == '__auto__':
                title_prefix = ''
            else:
                title_prefix = title
            title = title_prefix + '\n' + self.validity[0].get().isoformat(sep=b' ')
        else:
            title_prefix = None
        field0 = self.deepcopy()
        for c in field0.components:
            c.validity = self.validity[0]
        field0.validity = field0.components[0].validity
        field0.setdata([d[0, ...] for d in self.getdata()])
        if kwargs.get('plot_module', True):
            module = self.to_module()
            mindata = module.getdata(subzone=kwargs.get('subzone')).min()
            maxdata = module.getdata(subzone=kwargs.get('subzone')).max()
            plot_module_options = kwargs.get('plot_module_options', {})
            if plot_module_options == {}:
                kwargs['plot_module_options'] = {}
            minmax = plot_module_options.get('minmax')
            if minmax is None:
                minmax = (mindata, maxdata)
            kwargs['plot_module_options']['minmax'] = minmax
        academic = self.geometry.name == 'academic'
        if not academic:
            bm = kwargs.get('use_basemap')
            if bm is None:
                bm = self.geometry.make_basemap(gisquality=kwargs.get('gisquality', 'i'),
                                                subzone=kwargs.get('subzone'),
                                                specificproj=kwargs.get('specificproj'),
                                                zoom=kwargs.get('zoom'))
            kwargs['use_basemap'] = bm
        fig, ax = field0.plotfield(title=title,
                                   **kwargs)
        if kwargs.get('plot_module', True):
            if kwargs['plot_module_options'].get('colorbar_over') is None:
                kwargs['plot_module_options']['colorbar_over'] = fig.axes[-1]  # the last being created, in plotfield()
        kwargs['over'] = (fig, ax)

        def update(i, ax, myself, fieldi, title_prefix, kwargs):
            if i < len(myself.validity):
                ax.clear()
                for c in fieldi.components:
                    c.validity = myself.validity[i]
                fieldi.validity = fieldi.components[0].validity
                fieldi.setdata([d[i, ...] for d in myself.getdata()])
                if title_prefix is not None:
                    title = title_prefix + '\n' + fieldi.validity.get().isoformat(sep=b' ')
                fieldi.plotfield(title=title,
                                 **kwargs)

        anim = animation.FuncAnimation(fig, update,
                                       fargs=[ax, self, field0, title_prefix, kwargs],
                                       frames=list(range(len(self.validity) + 1)),  # AM: don't really understand why but needed for the last frame to be shown
                                       interval=interval,
                                       repeat=repeat)

        return anim
