#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

import numpy
import datetime

import footprints
from bronx.graphics.axes import set_figax, set_nice_time_axis

from epygram import config, util, epygramError
from epygram.base import FieldSet
from epygram.geometries import Geometry
from .D3Field import D3Field

epylog = footprints.loggers.getLogger(__name__)


class H1DField(D3Field):
    """
    Horizontal 1-Dimensions field class.
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
                values=set(['H1D'])),
            geometry=dict(
                type=Geometry),
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

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]
    def plotfield(self, *args, **kwargs):
        """
        Interface method to methods plotprofiles() and plotverticalhovmoller(),
        depending on time dimension (or not).

        Cf. these functions for arguments.
        """
        if len(self.validity) == 1:
            return self.plottransects(*args, **kwargs)
        else:
            return self.plothorizontalhovmoller(*args, **kwargs)

    def plottransects(self, *args, **kwargs):
        """Cf. eponymous function of module for arguments."""
        return plottransects(self, *args, **kwargs)

    def plothorizontalhovmoller(self, *args, **kwargs):
        """Cf. eponymous function of module for arguments."""
        return plothorizontalhovmoller(self, *args, **kwargs)

    def plotanimation(self, *args, **kwargs):
        """Cf. eponymous function of module for arguments."""
        return plotanimation(self, *args, **kwargs)


# FUNCTIONS #
#############
def plothorizontalhovmoller(transect,
                            over=(None, None),
                            fidkey=None,
                            title=None,
                            logscale=False,
                            zoom=None,
                            colorbar='vertical',
                            graphicmode='colorshades',
                            minmax=None,
                            levelsnumber=21,
                            center_cmap_on_0=False,
                            colormap='jet',
                            minmax_in_title=True,
                            contourcolor='k',
                            contourwidth=1,
                            contourlabel=True,
                            datefmt=None,
                            showgrid=True,
                            x_is='distance',
                            figsize=None,
                            rcparams=None):
    """
    Makes a simple vertical Hovmöller plot of the field.

    :param transect: being a :class:`epygram.fields.H1DField`
    :param over: any existing figure and/or ax to be used for the
      plot, given as a tuple (fig, ax), with None for
      missing objects. *fig* is the frame of the
      matplotlib figure, containing eventually several
      subplots (axes); *ax* is the matplotlib axes on
      which the drawing is done. When given (!= None),
      these objects must be coherent, i.e. ax being one of
      the fig axes.
    :param fidkey: type of fid for entitling the plot with *fid[fidkey]*,
                 if title is *None*;
                 if *None*, labels with raw fid.
    :param title: title for the plot.
    :param logscale: to set Y logarithmic scale
    :param zoom: a dict containing optional limits to zoom on the plot. \n
      Syntax: e.g. {'ymax':500, ...}.
    :param colorbar: if *False*, hide colorbar the plot; else, befines the
      colorbar orientation, among ('horizontal', 'vertical').
      Defaults to 'vertical'.
    :param graphicmode: among ('colorshades', 'contourlines').
    :param minmax: defines the min and max values for the plot colorbar. \n
      Syntax: [min, max]. [0.0, max] also works. Default is min/max of the
      field.
    :param levelsnumber: number of levels for contours and colorbar.
    :param center_cmap_on_0: aligns the colormap center on the value 0.
    :param colormap: name of the **matplotlib** colormap to use.
    :param minmax_in_title: if True and minmax != None, adds min and max
      values in title
    :param contourcolor: color or colormap to be used for 'contourlines'
      graphicmode. It can be either a legal html color name, or a colormap
      name.
    :param contourwidth: width of contours for 'contourlines' graphicmode.
    :param contourlabel: displays labels on contours.
    :param datefmt: date format to use, e.g. "%Y-%m-%d %H:%M:%S %Z"
    :param showgrid: True/False to show grid or not
    :param x_is: abscissa to be among:
      - 'distance': distance from first point of transect
      - 'lon': longitude of points
      - 'lat': latitude of points
    :param figsize: figure sizes in inches, e.g. (5, 8.5).
                    Default figsize is config.plotsizes.
    :param rcparams: list of (*args, **kwargs) to be passed to pyplot.rc()
                     defaults to [(('font',), dict(family='serif')),]

    Warning: requires **matplotlib**.
    """

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if rcparams is None:
        rcparams = [(('font',), dict(family='serif')),]
    for args, kwargs in rcparams:
        plt.rc(*args, **kwargs)
    if figsize is None:
        figsize = config.plotsizes

    # User colormaps
    if colormap not in plt.colormaps():
        util.load_cmap(colormap)

    # Figure, ax
    fig, ax = set_figax(*over, figsize=figsize)

    # coords
    y = numpy.zeros((len(transect.validity),
                     transect.geometry.dimensions['X']))

    validities = {transect.validity[i].get():i for i in range(len(transect.validity))}
    yaxis_label = 'Validity'
    if len(validities) == 1 and len(transect.validity) != 1:
        yaxis_label = 'Basis'
        validities = {transect.validity[i].getbasis():i for i in len(transect.validity)}
    epoch = datetime.datetime(1970, 1, 1)
    for i in range(len(transect.validity)):
        d = transect.validity[i].get() if yaxis_label == 'Validity' else transect.validity[i].getbasis()
        timedelta = d - epoch
        p = (timedelta.microseconds + (timedelta.seconds + timedelta.days * 24 * 3600) * 10 ** 6) / 1e6
        if hasattr(mdates, 'epoch2num'):
            y[i, :] = mdates.epoch2num(p)
        else:
            y[i, :] = mdates.date2num(mdates.datetime.datetime.utcfromtimestamp(p))
    x = numpy.zeros((len(transect.validity),
                     transect.geometry.dimensions['X']))
    lonlat = zip(*transect.geometry.get_lonlat_grid())
    p0 = lonlat[0]
    plast = p0
    distance = 0
    for i in range(transect.geometry.dimensions['X']):
        if x_is == 'distance':
            p = lonlat[i]
            distance += transect.geometry.distance((plast[0], plast[1]),
                                                   (p[0], p[1]))
            x[:, i] = distance
            plast = p
        elif x_is == 'lon':
            x[:, i] = lonlat[i][0]
        elif x_is == 'lat':
            x[:, i] = lonlat[i][1]
    plast = lonlat[-1]
    data = transect.getdata()
    # min/max
    m = data.min()
    M = data.max()
    if minmax is not None:
        if minmax_in_title:
            minmax_in_title = '(min: ' + \
                              '{: .{precision}{type}}'.format(m, type='E', precision=3) + \
                              ' // max: ' + \
                              '{: .{precision}{type}}'.format(M, type='E', precision=3) + ')'
        try:
            m = float(minmax[0])
        except Exception:
            m = data.min()
        try:
            M = float(minmax[1])
        except Exception:
            M = data.max()
    else:
        minmax_in_title = ''
    if abs(m - M) > config.epsilon:
        levels = numpy.linspace(m, M, levelsnumber)
        vmin = vmax = None
        if center_cmap_on_0:
            vmax = max(abs(m), M)
            vmin = -vmax
    else:
        raise epygramError("cannot plot uniform field.")
    L = int((levelsnumber - 1) // 15) + 1
    hlevels = [levels[l] for l in range(len(levels) - L / 3) if
               l % L == 0] + [levels[-1]]
    # plot
    if logscale:
        ax.set_yscale('log')
    ax.grid()
    if graphicmode == 'colorshades':
        pf = ax.contourf(x, y, data, levels, cmap=colormap,
                          vmin=vmin, vmax=vmax)
        if colorbar:
            position = 'right' if colorbar == 'vertical' else 'bottom'
            cax = make_axes_locatable(ax).append_axes(position,
                                                      size="5%",
                                                      pad=0.1)
            cb = plt.colorbar(pf,
                              orientation=colorbar,
                              ticks=hlevels,
                              cax=cax)
            if minmax_in_title != '':
                cb.set_label(minmax_in_title)
    elif graphicmode == 'contourlines':
        pf = ax.contour(x, y, data, levels=levels, colors=contourcolor,
                        linewidths=contourwidth)
        if contourlabel:
            ax.clabel(pf, colors=contourcolor)
    # time
    set_nice_time_axis(ax, 'y',
                       showgrid=showgrid, datefmt=datefmt)
    # decoration
    if x_is == 'distance':
        ax.set_xlabel('Distance from left-end point (m).')
    elif x_is == 'lon':
        ax.set_xlabel(u'Longitude (\u00B0).')
    elif x_is == 'lat':
        ax.set_xlabel(u'Latitude (\u00B0).')
    ax.set_ylabel(yaxis_label)
    if zoom is not None:
        ykw = {}
        xkw = {}
        for pair in (('bottom', 'ymin'), ('top', 'ymax')):
            try:
                ykw[pair[0]] = zoom[pair[1]]
            except Exception:
                pass
        for pair in (('left', 'xmin'), ('right', 'xmax')):
            try:
                xkw[pair[0]] = zoom[pair[1]]
            except Exception:
                pass
        ax.set_ylim(**ykw)
        ax.set_xlim(**xkw)
    if title is None:
        if fidkey is None:
            fid = transect.fid[sorted(transect.fid.keys())[0]]
        else:
            fid = transect.fid[fidkey]
        title = u'Horizontal Hovmöller of ' + str(fid) + ' between \n' + \
                '<- (' + str(p0[0]) + ', ' + \
                str(p0[1]) + ')' + ' and ' + \
                '(' + str(plast[0]) + ', ' + \
                str(plast[1]) + ') ->'
    ax.set_title(title)

    return (fig, ax)


def plottransects(transects,
                  over=(None, None),
                  labels=None,
                  fidkey=None,
                  unit='SI',
                  title=None,
                  logscale=False,
                  zoom=None,
                  x_is='distance',
                  figsize=None,
                  rcparams=None):
    """
    To plot a series of transects. Returns a tuple of :mod:`matplotlib`
    (*Figure*, *ax*).

    :param transects: being a :class:`epygram.base.FieldSet` of
                      :class:`epygram.fields.H1DField`, or a single
                      :class:`epygram.fields.H1DField`. \n
                      All transects are supposed to have the same unit,
                      and the same horizontal coordinate.
    :param over: any existing figure and/or ax to be used for the
      plot, given as a tuple (fig, ax), with None for
      missing objects. *fig* is the frame of the
      matplotlib figure, containing eventually several
      subplots (axes); *ax* is the matplotlib axes on
      which the drawing is done. When given (!= None),
      these objects must be coherent, e.g. ax being one of
      the fig axes.
    :param labels: a list of labels for the profiles (same length and same order).
    :param fidkey: key of fid for labelling the curve with *fid[fidkey]*;
                  if *None*, labels with raw fid.
    :param unit: label for X coordinate.
    :param title: title for the plot.
    :param logscale: to set Y logarithmic scale
    :param zoom: a dict containing optional limits to zoom on the plot. \n
      Syntax: e.g. {'ymax':500, ...}.
    :param x_is: abscissa to be among:
      - 'distance': distance from first point of transect
      - 'lon': longitude of points
      - 'lat': latitude of points
    :param figsize: figure sizes in inches, e.g. (5, 8.5).
                    Default figsize is config.plotsizes.
    :param rcparams: list of (*args, **kwargs) to be passed to pyplot.rc()
                     defaults to [(('font',), dict(family='serif')),
                                  (('figure',), dict(autolayout=True))]
    """
    import matplotlib.pyplot as plt
    if rcparams is None:
        rcparams = [(('font',), dict(family='serif')),
                    (('figure',), dict(autolayout=True))]
    for args, kwargs in rcparams:
        plt.rc(*args, **kwargs)
    plt.rc('figure', autolayout=True)
    if figsize is None:
        figsize = config.plotsizes

    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'darkolivegreen',
              'yellow', 'salmon', 'black']
    linestyles = ['-', '--', '-.', ':']

    if not isinstance(transects, FieldSet):
        p = transects.deepcopy()
        transects = FieldSet()
        transects.append(p)
    t0 = transects[0]
    x = numpy.zeros((t0.geometry.dimensions['X'],))
    lonlat = list(zip(*t0.geometry.get_lonlat_grid()))
    p0 = lonlat[0]
    plast = p0
    distance = 0
    for i in range(t0.geometry.dimensions['X']):
        if x_is == 'distance':
            p = lonlat[i]
            distance += t0.geometry.distance((plast[0], plast[1]),
                                               (p[0], p[1]))
            x[i] = distance
            plast = p
        elif x_is == 'lon':
            x[i] = lonlat[i][0]
        elif x_is == 'lat':
            x[i] = lonlat[i][1]
    plast = lonlat[-1]

    # Figure
    fig, ax = set_figax(*over, figsize=figsize)
    if logscale:
        ax.set_yscale('log')
    i = 0
    for t in transects:
        if len(t.validity) != 1:
            raise epygramError("plottransects can handle only profiles with one validity.")
        if labels is not None:
            label = labels[i]
        else:
            if fidkey is not None:
                label = t.fid.get(fidkey, t.fid)
            else:
                label = str(t.fid)
        data = t.getdata()
        plot_kwargs = {}
        if len(transects) > 1:
            plot_kwargs['color'] = colors[i % len(colors)]
            plot_kwargs['linestyle'] = linestyles[i // len(colors)]
        ax.plot(x, data.flatten(), label=label,
                **plot_kwargs)
        i += 1

    # Decoration
    if zoom is not None:
        if 'ymin' in zoom.keys():
            ax.set_ylim(bottom=zoom['ymin'])
        if 'ymax' in zoom.keys():
            ax.set_ylim(top=zoom['ymax'])
        if 'xmin' in zoom.keys():
            ax.set_xlim(left=zoom['xmin'])
        if 'xmax' in zoom.keys():
            ax.set_xlim(right=zoom['xmax'])
    if title is None:
        title = u'Transect between \n' + \
                '<- (' + str(p0[0]) + ', ' + \
                str(p0[1]) + ')' + ' and ' + \
                '(' + str(plast[0]) + ', ' + \
                str(plast[1]) + ') ->'
    ax.set_title(title)
    legend = ax.legend(loc='upper right', shadow=True)
    for label in legend.get_texts():
        label.set_fontsize('medium')
    ax.set_ylabel(r'$' + unit + '$')
    if x_is == 'distance':
        ax.set_xlabel('Distance from left-end point (m).')
    elif x_is == 'lon':
        ax.set_xlabel(u'Longitude (\u00B0).')
    elif x_is == 'lat':
        ax.set_xlabel(u'Latitude (\u00B0).')
    ax.grid()

    return (fig, ax)


def plotanimation(transect,
                  title='__auto__',
                  repeat=False,
                  interval=1000,
                  **kwargs):
    """
    To plot a time-dependent transect as an animation.
    Returns a :class:`matplotlib.animation.FuncAnimation`.

    :param transect: being a :class:`epygram.fields.H1DField`.
    :param title: title for the plot. '__auto__' (default) will print
      the current validity of the time frame.
    :param repeat: to repeat animation
    :param interval: number of milliseconds between two validities

    Other kwargs passed to plottransects().
    """
    import matplotlib.animation as animation

    if len(transect.validity) == 1:
        raise epygramError("plotanimation can handle only transect with several validities.")

    if title is not None:
        if title == '__auto__':
            title_prefix = ''
        else:
            title_prefix = title
        title = title_prefix + '\n' + transect.validity[0].get().isoformat(sep=' ')
    else:
        title_prefix = None
    transect0 = transect.getvalidity(0)
    mindata = transect.getdata().min()
    maxdata = transect.getdata().max()
    mindata -= (maxdata - mindata) / 10
    maxdata += (maxdata - mindata) / 10

    zoom = kwargs.get('zoom')
    zoom = util.ifNone_emptydict(zoom)
    if 'ymax' not in zoom.keys():
        zoom.update(ymax=maxdata)
    if 'ymin' not in zoom.keys():
        zoom.update(ymin=mindata)
    kwargs['zoom'] = zoom

    fig, ax = plottransects(transect0,
                            title=title,
                            **kwargs)
    if kwargs.get('colorbar_over') is None:
        kwargs['colorbar_over'] = fig.axes[-1]  # the last being created, in plotfield()
    kwargs['over'] = (fig, ax)

    def update(i, ax, myself, transecti, title_prefix, kwargs):
        print("updating")
        if i < len(myself.validity):
            ax.clear()
            transecti = myself.getvalidity(i)
            if title_prefix is not None:
                title = title_prefix + '\n' + transecti.validity.get().isoformat(sep=' ')
            transecti.plotfield(title=title, **kwargs)

    anim = animation.FuncAnimation(fig, update,
                                   fargs=[ax, transect, transect0, title_prefix, kwargs],
                                   frames=range(len(transect.validity) + 1),  # AM: don't really understand why but needed for the last frame to be shown
                                   interval=interval,
                                   repeat=repeat)

    return anim
