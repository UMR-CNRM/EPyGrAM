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


def activate():
    """Activate extension."""
    from epygram.fields import H2DVectorField
    H2DVectorField.basemap_plot = basemap_plot
    H2DVectorField.plotfield = plotfield
    H2DVectorField.plotanimation = plotanimation


def plotfield(self, *args, **kwargs):
        return self.basemap_plot(*args, **kwargs)


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
