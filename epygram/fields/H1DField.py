#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Horizontal 2D field.
"""

import numpy

import footprints

from .D3Field import D3Field
from epygram.base import FieldSet
from epygram import config, util, epygramError
from epygram.geometries import H1DGeometry
import datetime
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
                type=H1DGeometry),
        )
    )

##############
# ABOUT DATA #
##############
    def getlevel(self, level=None, k=None):
        """Returns self. Useless but for compatibility reasons."""

        if k == None and level == None:
            raise epygramError("You must give k or level.")
        if k != None and level != None:
            raise epygramError("You cannot give, at the same time, k and level")
        if level != None:
            if level not in self.geometry.vcoordinate.levels:
                raise epygramError("The requested level does not exist.")
        return self

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """

        if isinstance(other, self.__class__):
            if self.spectral != other.spectral:
                raise epygramError("cannot operate a spectral field with a" + \
                                   " non-spectral field.")
            if self.geometry.rectangular_grid:
                if self.geometry.dimensions['X'] != other.geometry.dimensions['X'] or\
                   self.geometry.dimensions['Y'] != other.geometry.dimensions['Y']:
                    raise epygramError("operations on fields cannot be done if" + \
                                       " fields do not share their gridpoint" + \
                                       " dimensions.")
            else:
                if self.geometry.dimensions != other.geometry.dimensions:
                    raise epygramError("operations on fields cannot be done if" + \
                                       " fields do not share their gridpoint" + \
                                       " dimensions.")
            if self.spectral_geometry != other.spectral_geometry:
                raise epygramError("operations on fields cannot be done if" + \
                                   " fields do not share their spectral" + \
                                   " geometry.")
        else:
            super(D3Field, self)._check_operands(other)

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def plotfield(self,
                  over=(None, None),
                  labels=None,
                  fidkey=None,
                  unit='SI',
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
                  linecolor='black',
                  linestyle='-',
                  force_mode=False,
                  repeat=False,
                  interval=1000,
                  x_is='distance',
                  **ignored_kwargs):
        """
        Makes a simple (transect) plot of the field.
        Help on arguments can be found in actual plot functions docstrings.
        
        Args: \n
        - *force_mode* = *False* (default) to plot an hovmoller plot
                         if several validities are available,
                         a simple profile plot otherwise.
                         *hovmoller* to force an hovmoller plot
                         *transect* to force a simple transect plot
                         *animation* to force an animation plot
        """

        if force_mode != False:
            mode = force_mode
        else:
            if len(self.validity) == 1:
                mode = "transect"
            else:
                mode = "hovmoller"

        useless_args = []
        if mode == 'transect':
            if colorbar != 'vertical': useless_args.append('colorbar')
            if graphicmode != 'colorshades': useless_args.append('graphicmode')
            if minmax != None: useless_args.append('minmax')
            if levelsnumber != 21: useless_args.append('levelsnumber')
            if center_cmap_on_0 != False: useless_args.append('center_cmap_on_0')
            if colormap != 'jet': useless_args.append('colormap')
            if minmax_in_title != True: useless_args.append('minmax_in_title')
            if contourcolor != 'k': useless_args.append('contourcolor')
            if contourwidth != 1: useless_args.append('contourwidth')
            if contourlabel != True: useless_args.append('contourlabel')
            if datefmt != "%Y%m%d %H:%M:%S %Z": useless_args.append('datefmt')
            if linecolor != 'black': useless_args.append('linecolor')
            if linestyle != '-': useless_args.append('linestyle')
            if repeat != False: useless_args.append('repeat')
            if interval != 1000: useless_args.append('interval')
            return plottransects(self,
                                 over=over,
                                 labels=labels,
                                 fidkey=fidkey,
                                 unit=unit,
                                 title=title,
                                 logscale=logscale,
                                 zoom=zoom,
                                 x_is=x_is)
        elif mode == 'hovmoller':
            if labels != None: useless_args.append('labels')
            if unit != 'SI': useless_args.append('unit')
            if linecolor != 'black':useless_args.append('linecolor')
            if linestyle != '-':useless_args.append('linestyle')
            if repeat != False: useless_args.append('repeat')
            if interval != 1000: useless_args.append('interval')
            return plothorizontalhovmoller(self,
                                           over=over,
                                           fidkey=fidkey,
                                           title=title,
                                           logscale=logscale,
                                           zoom=zoom,
                                           colorbar=colorbar,
                                           graphicmode=graphicmode,
                                           minmax=minmax,
                                           levelsnumber=levelsnumber,
                                           center_cmap_on_0=center_cmap_on_0,
                                           colormap=colormap,
                                           minmax_in_title=minmax_in_title,
                                           contourcolor=contourcolor,
                                           contourwidth=contourwidth,
                                           contourlabel=contourlabel,
                                           datefmt=datefmt,
                                           x_is=x_is)
        elif mode == 'animation':
            if labels != None: useless_args.append('labels')
            if colorbar != 'vertical': useless_args.append('colorbar')
            if graphicmode != 'colorshades': useless_args.append('graphicmode')
            if minmax != None: useless_args.append('minmax')
            if levelsnumber != 21: useless_args.append('levelsnumber')
            if center_cmap_on_0 != False: useless_args.append('center_cmap_on_0')
            if colormap != 'jet': useless_args.append('colormap')
            if minmax_in_title != True: useless_args.append('minmax_in_title')
            if contourcolor != 'k': useless_args.append('contourcolor')
            if contourwidth != 1: useless_args.append('contourwidth')
            if contourlabel != True: useless_args.append('contourlabel')
            if datefmt != "%Y%m%d %H:%M:%S %Z": useless_args.append('datefmt')
            return plotanimation(self,
                                 over=over,
                                 fidkey=fidkey,
                                 unit=unit,
                                 title='__auto__' if title is None else title,
                                 logscale=logscale,
                                 zoom=zoom,
                                 repeat=repeat,
                                 interval=interval,
                                 x_is=x_is)
        else:
            raise NotImplementedError("This graphic mode is not implemented")

        if len(useless_args) != 0:
            epylog.warning("Some arguments to plotfield are useless and will not be used: " + str(useless_args))

    def plottransects(self, *args, **kwargs):
        return plottransects(self, *args, **kwargs)

    def plothorizontalhovmoller(self, *args, **kwargs):
        return plothorizontalhovmoller(self, *args, **kwargs)

    def plotanimation(self, *args, **kwargs):
        return plotanimation(self, *args, **kwargs)

#################
### FUNCTIONS ###
#################

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
                            x_is='distance'):
        """
        Makes a simple vertical Hovmöller plot of the field.

        Args: \n
        - *transect* being a :class:`epygram.fields.H1DField`
        - *over* = any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several 
          subplots (axes); *ax* is the matplotlib axes on 
          which the drawing is done. When given (!= None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        - *fidkey* = type of fid for entitling the plot with *fid[fidkey]*,
                     if title is *None*;
                     if *None*, labels with raw fid.
        - *title* = title for the plot.
        - *logscale* = to set Y logarithmic scale
        - *zoom*: a dict containing optional limits to zoom on the plot. \n
          Syntax: e.g. {'ymax':500, ...}.
        - *colorbar*: if *False*, hide colorbar the plot; else, befines the
          colorbar orientation, among ('horizontal', 'vertical').
          Defaults to 'vertical'.
        - *graphicmode*: among ('colorshades', 'contourlines').
        - *minmax*: defines the min and max values for the plot colorbar. \n
          Syntax: [min, max]. [0.0, max] also works. Default is min/max of the
          field.
        - *levelsnumber*: number of levels for contours and colorbar.
        - *center_cmap_on_0*: aligns the colormap center on the value 0.
        - *colormap*: name of the **matplotlib** colormap to use.
        - *minmax_in_title*: if True and minmax != None, adds min and max
          values in title
        - *contourcolor*: color or colormap to be used for 'contourlines'
          graphicmode. It can be either a legal html color name, or a colormap
          name.
        - *contourwidth*: width of contours for 'contourlines' graphicmode.
        - *contourlabel*: displays labels on contours.
        - *datefmt*: date format to use, e.g. "%Y-%m-%d %H:%M:%S %Z"
        - *showgrid*: True/False to show grid or not
        - *x_is*: abscissa to be among:
          - 'distance': distance from first point of transect
          - 'lon': longitude of points
          - 'lat': latitude of points

        Warning: requires **matplotlib**.
        """

        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)

        # User colormaps
        if colormap not in plt.colormaps():
            util.add_cmap(colormap)

        # Figure, ax
        fig, ax = util.set_figax(*over)

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
            y[i, :] = matplotlib.dates.epoch2num(p)
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
        if minmax != None:
            if minmax_in_title:
                minmax_in_title = '(min: ' + \
                                  '{: .{precision}{type}}'.format(m, type='E', precision=3) + \
                                  ' // max: ' + \
                                  '{: .{precision}{type}}'.format(M, type='E', precision=3) + ')'
            try: m = float(minmax[0])
            except Exception: m = data.min()
            try: M = float(minmax[1])
            except Exception: M = data.max()
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
        ymin = mdates.num2date(ax.axis()[2]).replace(tzinfo=None)
        ymax = mdates.num2date(ax.axis()[3]).replace(tzinfo=None)
        util.set_DateHour_axis(ax, ymax - ymin, 'y',
                               showgrid=showgrid, datefmt=datefmt)
        # decoration
        if x_is == 'distance':
            ax.set_xlabel('Distance from left-end point (m).')
        elif x_is == 'lon':
            ax.set_xlabel(u'Longitude (\u00B0).')
        elif x_is == 'lat':
            ax.set_xlabel(u'Latitude (\u00B0).')
        ax.set_ylabel(yaxis_label)
        if zoom != None:
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
                 x_is='distance'):
    """
    To plot a series of transects. Returns a tuple of :mod:`matplotlib`
    (*Figure*, *ax*).

    Args: \n
    - *transects* being a :class:`epygram.base.FieldSet` of
      :class:`epygram.fields.H1DField`, or a single
      :class:`epygram.fields.H1DField`. \n
      All transects are supposed to have the same unit, and the same horizontal
      coordinate.
    - *over* = any existing figure and/or ax to be used for the
      plot, given as a tuple (fig, ax), with None for
      missing objects. *fig* is the frame of the
      matplotlib figure, containing eventually several 
      subplots (axes); *ax* is the matplotlib axes on 
      which the drawing is done. When given (!= None),
      these objects must be coherent, e.g. ax being one of
      the fig axes.
    - *labels* = a list of labels for the profiles (same length and same order).
    - *fidkey* = key of fid for labelling the curve with *fid[fidkey]*;
                  if *None*, labels with raw fid.
    - *unit* = label for X coordinate.
    - *title* = title for the plot.
    - *logscale* = to set Y logarithmic scale
    - *zoom*: a dict containing optional limits to zoom on the plot. \n
      Syntax: e.g. {'ymax':500, ...}.
    - *x_is*: abscissa to be among:
      - 'distance': distance from first point of transect
      - 'lon': longitude of points
      - 'lat': latitude of points

    """
    import matplotlib.pyplot as plt

    plt.rc('font', family='serif')
    plt.rc('figure', autolayout=True)

    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'darkolivegreen',
              'yellow', 'salmon', 'black']
    linestyles = ['-', '--', '-.', ':']

    if not isinstance(transects, FieldSet):
        p = transects.deepcopy()
        transects = FieldSet()
        transects.append(p)
    t0 = transects[0]
    x = numpy.zeros((t0.geometry.dimensions['X'],))
    lonlat = zip(*t0.geometry.get_lonlat_grid())
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
    fig, ax = util.set_figax(*over, figsize=(6., 9.))
    if logscale:
        ax.set_yscale('log')
    i = 0
    for t in transects:
        if len(t.validity) != 1:
            raise epygramError("plottransects can handle only profiles with one validity.")
        if labels != None:
            label = labels[i]
        else:
            if fidkey != None:
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

    Args: \n
    - *transect* being a :class:`epygram.fields.H1DField`.
    - *title* = title for the plot. '__auto__' (default) will print
      the current validity of the time frame.
    - *repeat*: to repeat animation
    - *interval*: number of milliseconds between two validities
    
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
    if not 'ymax' in zoom.keys():
        zoom.update(ymax=maxdata)
    if not 'ymin' in zoom.keys():
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
                title = title_prefix + '\n' + transecti.validity.get().isoformat(sep=b' ')
            transecti.plotfield(title=title,
                             **kwargs)

    anim = animation.FuncAnimation(fig, update,
                                   fargs=[ax, transect, transect0, title_prefix, kwargs],
                                   frames=range(len(transect.validity) + 1),  # AM: don't really understand why but needed for the last frame to be shown
                                   interval=interval,
                                   repeat=repeat)

    return anim