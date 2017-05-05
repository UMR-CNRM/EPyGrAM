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

import footprints

from .D3Field import D3Field
from .PointField import PointField
from epygram import config, util, epygramError
from epygram.geometries import H2DGeometry
from . import gimme_one_point

epylog = footprints.loggers.getLogger(__name__)


class H2DField(D3Field):
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
                             field_args={'validity':self.validity,
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
    def plotfield(self,
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
                  pointsmarker=','):
        """
        Makes a simple plot of the field, with a number of options.

        Requires :mod:`matplotlib`

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        subzone          among ('C', 'CI'), for LAM fields only, plots the
                         data resp. on the C or C+I zone. \n
                         Default is no subzone, i.e. the whole field.
        gisquality       among ('c', 'l', 'i', 'h', 'f') -- by increasing
                         quality. Defines the quality for GIS elements
                         (coastlines, countries boundaries...).
        specificproj     enables to make basemap on the specified projection,
                         among: 'kav7', 'cyl', 'ortho', ('nsper', {...})
                         (cf. Basemap doc). \n
                         In 'nsper' case, the {} may contain: ('sat_height' =
                         satellite height in km; 'lon' = longitude of nadir
                         in degrees; 'lat' = latitude of nadir in degrees.
        zoom             specifies the lon/lat borders of the map, implying
                         hereby a 'cyl' projection. Must be a dict(lonmin=,
                         lonmax=, latmin=, latmax=).\n
                         Overwrites *specificproj*.
        over             any existing figure and/or ax to be used for the
                         plot, given as a tuple (fig, ax), with None for
                         missing objects. *fig* is the frame of the
                         matplotlib figure, containing eventually several
                         subplots (axes); *ax* is the matplotlib axes on
                         which the drawing is done. When given (is not None),
                         these objects must be coherent, i.e. ax being one of
                         the fig axes.
        colorbar_over    an optional existing ax to plot the colorbar on.
        use_basemap      a basemap.Basemap object used to handle the
                         projection of the map. If given, the map projection
                         options (*specificproj*, *zoom*, *gisquality* ...)
                         are ignored, keeping the properties of the
                         *use_basemap* object.
        title            title for the plot. Default is field identifier.
        minmax           defines the min and max values for the plot
                         colorbar. \n
                         Syntax: [min, max]. [0.0, max] also works. Default
                         is min/max of the field.
        graphicmode      among ('colorshades', 'contourlines', 'points').
        levelsnumber     number of levels for contours and colorbar.
        colormap         name of the ``matplotlib`` colormap to use (or an
                         ``epygram`` one, or a user-defined one, cf.
                         config.usercolormaps).
        center_cmap_on_0 aligns the colormap center on the value 0.
        drawrivers       to add rivers on map.
        drawcoastlines   to add coast lines on map.
        drawcountries    to add countries on map.
        colorbar         if *False*, hide colorbar the plot; else, defines the
                         colorbar position, among ('bottom', 'right').
                         Defaults to 'right'.
        meridians        enable to fine-tune the choice of lines to
                         plot, with either:\n
                         - 'auto': automatic scaling to the basemap extents
                         - 'default': range(0,360,10) and range(-90,90,10)
                         - a list of values
                         - a grid step, e.g. 5 to plot each 5 degree.
                         - None: no one is plot
                         - *meridians* == 'greenwich' // 'datechange' //
                           'greenwich+datechange' or any combination (+)
                           will plot only these.
        parallels        cf. meridians doc, with specific ones being
                         *parallels* == 'equator' // 'polarcircles' //
                         'tropics'
        minmax_in_title  if True and minmax is not None, adds min and max
                         values in title.
        departments      if True, adds the french departments on map (instead
                         of countries).
        boundariescolor  color of lines for boundaries (countries,
                         departments, coastlines)
        pointsize        size of points for *graphicmode* == 'points'.
        contourcolor     color or colormap to be used for 'contourlines'
                         graphicmode. It can be either a legal html color
                         name, or a colormap name.
        contourwidth     width of contours for 'contourlines' graphicmode.
        contourlabel     displays labels on contours.
        bluemarble       if > 0.0 (and <=1.0), displays NASA's "blue marble"
                         as background. The numerical value sets its
                         transparency.
        background       if True, set a background color to
                         continents and oceans.
        mask_threshold   dict with min and/or max value(s) to mask outside.
        contourlabelfmt  format of the contour labels: e.g. 273.15 will
                         appear: '%0i' => 273, '%0f' => 273.150000,
                         '%0.2f' => 273.15, '%04i' => 0273,
                         '%0.5e' => 2.731500e+02
        pointsmarker     shape of the points if graphicmode='points'.
                         Cf. matplotlib.scatter() for possible markers.

        ==============   ====================================================

        This method uses (hence requires) 'matplotlib' and 'basemap' libraries.
        """

        # 0. Initializations
        #####################
        # 0.1 matplotlib initializations
        import matplotlib.pyplot as plt
        from matplotlib.colors import cnames
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)

        # 0.2 checkings
        if self.spectral:
            raise epygramError("please convert to gridpoint with sp2gp()" +
                               " method before plotting.")
        if len(self.validity) > 1:
            raise epygramError('to plot H2DField with time dimension, use plotanimation().')

        # 0.3 add custom colormap if necessary
        if colormap not in plt.colormaps():
            if colormap in config.colormaps:
                util.add_cmap(colormap)
            else:
                from mpl_toolkits.basemap import cm
                if colormap in cm.datad:
                    plt.register_cmap(name=colormap, cmap=cm.__dict__[colormap])
        if graphicmode == 'contourlines':  # color of contours
            if contourcolor in cnames or contourcolor[0] == '#':
                colormap = None
            else:
                if contourcolor not in plt.colormaps():
                    util.add_cmap(contourcolor)
                colormap = contourcolor
                contourcolor = None

        if zoom in (None, {}):  # actual build of figure
            # 1. Figure, ax
            ###############
            fig, ax = util.set_figax(*over)

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
                util.set_map_up(bm, ax,
                                drawrivers=drawrivers,
                                drawcoastlines=drawcoastlines,
                                drawcountries=drawcountries,
                                meridians=meridians,
                                parallels=parallels,
                                departments=departments,
                                boundariescolor=boundariescolor,
                                bluemarble=bluemarble,
                                background=background)

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
            if minmax_in_title and (minmax is not None or colormap in config.colormaps_scaling):
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

            if abs(float(m) - float(M)) < config.epsilon:
                raise epygramError("cannot plot uniform field.")
            if center_cmap_on_0:
                vmax = max(abs(m), M)
                vmin = -vmax
            else:
                vmin = m
                vmax = M
            # set levels and ticks levels
            levels = numpy.linspace(m, M, levelsnumber)
            L = int((levelsnumber - 1) // 15) + 1
            tick_levels = [levels[l]
                           for l in range(len(levels) - (L // 3 + 1))
                           if l % L == 0] + [levels[-1]]
            if colormap in config.colormaps_scaling:
                (norm, levels) = util.color_scale(colormap)
                tick_levels = levels
                vmin = vmax = None
            else:
                norm = None

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
                                                                  size="5%",
                                                                  pad=0.2)
                    else:
                        cax = colorbar_over
                    orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
                    cb = plt.colorbar(pf,
                                      orientation=orientation,
                                      ticks=tick_levels,
                                      cax=cax)
                    if minmax_in_title:
                        cb.set_label(minmax_in_title)
            elif graphicmode == 'contourlines':
                if not self.geometry.rectangular_grid:
                    xf = x.compressed()
                    yf = y.compressed()
                    zf = data.compressed()
                    tri = True
                    # FIXME: problem of duplicate points with arpege grid
                elif self.geometry.dimensions['Y'] == 1:
                    xf = x.flatten()
                    yf = y.flatten()
                    zf = data.flatten()
                    tri = True
                    # FIXME: problem of triangulation for contourlines
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
                if colorbar:
                    if colorbar_over is None:
                        cax = make_axes_locatable(ax).append_axes(colorbar,
                                                                  size="5%",
                                                                  pad=0.2)
                    else:
                        cax = colorbar_over
                    orientation = 'vertical' if colorbar in ('right', 'left') else 'horizontal'
                    cb = plt.colorbar(pf,
                                      orientation=orientation,
                                      ticks=tick_levels,
                                      cax=cax)
                    if minmax_in_title != '':
                        cb.set_label(minmax_in_title)
            if title is None:
                ax.set_title("\n".join([str(self.fid[sorted(self.fid.keys())[0]]),
                                        str(self.validity.get())]))
            else:
                ax.set_title(title)
        else:
            # zoom: create zoom_field and plot it
            zoom_field = self.extract_zoom(zoom)
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

    def plotanimation(self, title='__auto__', repeat=False, interval=1000, **kwargs):
        """
        Plot the field with animation with regards to time dimension.
        Returns a :class:`matplotlib.animation.FuncAnimation`.

        In addition to those specified below, all :meth:`plotfield` method
        arguments can be provided.

        Args:\n
        - *title* = title for the plot. '__auto__' (default) will print
          the current validity of the time frame.
        - *repeat*: to repeat animation
        - *interval*: number of milliseconds between two validities
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
        field0.validity = self.validity[0]
        field0.setdata(self.getdata()[0, ...])
        mindata = self.getdata(subzone=kwargs.get('subzone')).min()
        maxdata = self.getdata(subzone=kwargs.get('subzone')).max()

        minmax = kwargs.get('minmax')
        if minmax is None:
            minmax = (mindata, maxdata)
        kwargs['minmax'] = minmax
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
        if kwargs.get('colorbar_over') is None:
            kwargs['colorbar_over'] = fig.axes[-1]  # the last being created, in plotfield()
        kwargs['over'] = (fig, ax)

        def update(i, ax, myself, fieldi, title_prefix, kwargs):
            if i < len(myself.validity):
                ax.clear()
                fieldi.validity = myself.validity[i]
                fieldi.setdata(myself.getdata()[i, ...])
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

    def morph_with_points(self, points, alpha=1., morphing='nearest',
                          increment=False,
                          **kwargs):
        """
        Perturb the field values with the values of a set of points.

        *points* meant to be a list/fieldset of PointField.
        *alpha* is the blending value, ranging from 0. to 1.:
          e.g. with 'nearest' morphing, the final value of the point is
          alpha*point.value + (1-alpha)*original_field.value
        *morphing* is the way the point modify the field:
          - 'nearest': only the nearest point is modified
          - 'exp_decay': modifies the surrounding points with an isotrop
            exponential decay weighting
          - 'gaussian': modifies the surrounding points with an isotrop
            gaussian weighting. Its standard deviation *sigma* must then
            be passed as argument, in meters. Weightingassumed to be 0. from
            3*sigma.
        *increment*: if True, the final value of the point is
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
