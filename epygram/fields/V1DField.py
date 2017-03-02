#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Vertical 1D field.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import datetime
import numpy
import six

import footprints

from .D3Field import D3Field
from epygram import epygramError, config, util
from epygram.geometries import V1DGeometry

epylog = footprints.loggers.getLogger(__name__)


class V1DField(D3Field):
    """
    Vertical 1-Dimension (column) field class.
    A field is defined by its identifier 'fid',
    its data, its geometry, and its validity.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V1D'])),
            geometry=dict(
                type=V1DGeometry,
                access='rwx'),
        )
    )


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
                  Ycoordinate=None,
                  unit='SI',
                  title=None,
                  logscale=False,
                  ema=False,
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
                  **ignored_kwargs):
        """
        Makes a simple (profile) plot of the field.
        Help on arguments can be found in actual plot functions docstrings.

        Args: \n
        - *force_mode* = *False* (default) to plot an hovmoller plot
                         if several validities are available,
                         a simple profile plot otherwise.
                         *hovmoller* to force an hovmoller plot
                         *profile* to force a simple profile plot
                         *animation* to force an animation plot
        """

        if force_mode is not False:
            mode = force_mode
        else:
            if len(self.validity) == 1:
                mode = "profile"
            else:
                mode = "hovmoller"

        useless_args = []
        if mode == 'profile':
            if colorbar != 'vertical':
                useless_args.append('colorbar')
            if graphicmode != 'colorshades':
                useless_args.append('graphicmode')
            if minmax is not None:
                useless_args.append('minmax')
            if levelsnumber != 21:
                useless_args.append('levelsnumber')
            if center_cmap_on_0 is not False:
                useless_args.append('center_cmap_on_0')
            if colormap != 'jet':
                useless_args.append('colormap')
            if minmax_in_title is not True:
                useless_args.append('minmax_in_title')
            if contourcolor != 'k':
                useless_args.append('contourcolor')
            if contourwidth != 1:
                useless_args.append('contourwidth')
            if contourlabel is not True:
                useless_args.append('contourlabel')
            if datefmt != "%Y%m%d %H:%M:%S %Z":
                useless_args.append('datefmt')
            if linecolor != 'black':
                useless_args.append('linecolor')
            if linestyle != '-':
                useless_args.append('linestyle')
            if repeat is not False:
                useless_args.append('repeat')
            if interval != 1000:
                useless_args.append('interval')
            return plotprofiles(self,
                                over=over,
                                labels=labels,
                                fidkey=fidkey,
                                Ycoordinate=Ycoordinate,
                                unit=unit,
                                title=title,
                                logscale=logscale,
                                ema=ema,
                                zoom=zoom)
        elif mode == 'hovmoller':
            if labels is not None:
                useless_args.append('labels')
            if unit != 'SI':
                useless_args.append('unit')
            if ema is not False:
                useless_args.append('ema')
            if linecolor != 'black':
                useless_args.append('linecolor')
            if linestyle != '-':
                useless_args.append('linestyle')
            if repeat is not False:
                useless_args.append('repeat')
            if interval != 1000:
                useless_args.append('interval')
            return plotverticalhovmoller(self,
                                         over=over,
                                         fidkey=fidkey,
                                         Ycoordinate=Ycoordinate,
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
                                         datefmt=datefmt)
        elif mode == 'animation':
            if labels is not None:
                useless_args.append('labels')
            if colorbar != 'vertical':
                useless_args.append('colorbar')
            if graphicmode != 'colorshades':
                useless_args.append('graphicmode')
            if minmax is not None:
                useless_args.append('minmax')
            if levelsnumber != 21:
                useless_args.append('levelsnumber')
            if center_cmap_on_0 is not False:
                useless_args.append('center_cmap_on_0')
            if colormap != 'jet':
                useless_args.append('colormap')
            if minmax_in_title is not True:
                useless_args.append('minmax_in_title')
            if contourcolor != 'k':
                useless_args.append('contourcolor')
            if contourwidth != 1:
                useless_args.append('contourwidth')
            if contourlabel is not True:
                useless_args.append('contourlabel')
            if datefmt != "%Y%m%d %H:%M:%S %Z":
                useless_args.append('datefmt')
            return plotanimation(self,
                                 over=over,
                                 fidkey=fidkey,
                                 Ycoordinate=Ycoordinate,
                                 unit=unit,
                                 title=title,
                                 logscale=logscale,
                                 ema=ema,
                                 zoom=zoom)
        else:
            raise NotImplementedError("This graphic mode is not implemented")

        if len(useless_args) != 0:
            epylog.warning("Some arguments to plotfield are useless and will not be used: " + str(useless_args))


#############
# FUNCTIONS #
#############
def plotverticalhovmoller(profile,
                          over=(None, None),
                          fidkey=None,
                          Ycoordinate=None,
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
                          showgrid=True):
        """
        Makes a simple vertical Hovmöller plot of the field.

        Args: \n
        - *profile* being a :class:`epygram.fields.V1DField`
        - *over* = any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several
          subplots (axes); *ax* is the matplotlib axes on
          which the drawing is done. When given (is not None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        - *fidkey* = type of fid for entitling the plot with *fid[fidkey]*,
                     if title is *None*;
                     if *None*, labels with raw fid.
        - *Ycoordinate* = label for the Y coordinate.
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
        - *minmax_in_title*: if True and minmax is not None, adds min and max
          values in title
        - *contourcolor*: color or colormap to be used for 'contourlines'
          graphicmode. It can be either a legal html color name, or a colormap
          name.
        - *contourwidth*: width of contours for 'contourlines' graphicmode.
        - *contourlabel*: displays labels on contours.
        - *datefmt*: date format to use, e.g. "%Y-%m-%d %H:%M:%S %Z"
        - *showgrid*: True/False to show grid or not

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
        z = numpy.zeros((len(profile.validity),
                         len(profile.geometry.vcoordinate.levels)))
        for k in range(len(profile.geometry.vcoordinate.levels)):
            z[:, k] = profile.geometry.vcoordinate.levels[k]
        x = numpy.zeros((len(profile.validity),
                         len(profile.geometry.vcoordinate.levels)))
        validities = {profile.validity[i].get():i for i in range(len(profile.validity))}
        xaxis_label = 'Validity'
        if len(validities) == 1 and len(profile.validity) != 1:
            xaxis_label = 'Basis'
            validities = {profile.validity[i].getbasis():i for i in len(profile.validity)}
        epoch = datetime.datetime(1970, 1, 1)
        for i in range(len(profile.validity)):
            d = profile.validity[i].get() if xaxis_label == 'Validity' else profile.validity[i].getbasis()
            timedelta = d - epoch
            p = (timedelta.microseconds + (timedelta.seconds + timedelta.days * 24 * 3600) * 1e6) / 1e6
            x[i, :] = matplotlib.dates.epoch2num(p)
        data = profile.getdata()
        if profile.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100):
            reverseY = True
        else:
            reverseY = False
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
        else:
            raise epygramError("cannot plot uniform field.")
        if center_cmap_on_0:
            vmax = max(abs(m), M)
            vmin = -vmax
        else:
            vmin = m
            vmax = M
        L = int((levelsnumber - 1) // 15) + 1
        hlevels = [levels[l]
                   for l in range(len(levels) - L // 3)
                   if l % L == 0] + [levels[-1]]
        # plot
        if reverseY and not ax.yaxis_inverted():
            ax.invert_yaxis()
        if logscale:
            ax.set_yscale('log')
        ax.grid()
        if graphicmode == 'colorshades':
            pf = ax.contourf(x, z, data, levels, cmap=colormap,
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
            pf = ax.contour(x, z, data, levels=levels, colors=contourcolor,
                            linewidths=contourwidth)
            if contourlabel:
                ax.clabel(pf, colors=contourcolor)
        # time
        xmin = mdates.num2date(ax.axis()[0]).replace(tzinfo=None)
        xmax = mdates.num2date(ax.axis()[1]).replace(tzinfo=None)
        util.set_DateHour_axis(ax, xmax - xmin,
                               showgrid=showgrid, datefmt=datefmt)
        # decoration
        surf = z[-1, :]
        bottom = max(surf) if reverseY else min(surf)
        ax.fill_between(x[-1, :], surf, numpy.ones(len(surf)) * bottom,
                        color='k')
        if Ycoordinate is None:
            if profile.geometry.vcoordinate.typeoffirstfixedsurface == 119:
                Ycoordinate = 'Level \nHybrid-Pressure \ncoordinate'
            elif profile.geometry.vcoordinate.typeoffirstfixedsurface == 100:
                Ycoordinate = 'Pressure (hPa)'
            elif profile.geometry.vcoordinate.typeoffirstfixedsurface == 102:
                Ycoordinate = 'Altitude (m)'
            elif profile.geometry.vcoordinate.typeoffirstfixedsurface == 103:
                Ycoordinate = 'Height (m)'
            elif profile.geometry.vcoordinate.typeoffirstfixedsurface == 118:
                Ycoordinate = 'Level \nHybrid-Height \ncoordinate'
            elif profile.geometry.vcoordinate.typeoffirstfixedsurface == 109:
                Ycoordinate = 'Potential \nvortex \n(PVU)'
            else:
                Ycoordinate = 'unknown \ncoordinate'
        ax.set_xlabel(xaxis_label)
        ax.set_ylabel(Ycoordinate)
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
                fid = profile.fid[sorted(profile.fid.keys())[0]]
            else:
                fid = profile.fid[fidkey]
            title = u'Vertical Hovmöller of ' + str(fid)
        ax.set_title(title)

        return (fig, ax)


def plotprofiles(profiles,
                 over=(None, None),
                 labels=None,
                 fidkey=None,
                 Ycoordinate=None,
                 unit='SI',
                 title=None,
                 logscale=False,
                 ema=False,
                 zoom=None):
    """
    To plot a series of profiles. Returns a tuple of :mod:`matplotlib`
    (*Figure*, *ax*).

    Args: \n
    - *profiles* being a :class:`epygram.base.FieldSet` of
      :class:`epygram.fields.V1DField`, or a single
      :class:`epygram.fields.V1DField`. \n
      All profiles are supposed to have the same unit, and the same vertical
      coordinate.
    - *over* = any existing figure and/or ax to be used for the
      plot, given as a tuple (fig, ax), with None for
      missing objects. *fig* is the frame of the
      matplotlib figure, containing eventually several
      subplots (axes); *ax* is the matplotlib axes on
      which the drawing is done. When given (is not None),
      these objects must be coherent, e.g. ax being one of
      the fig axes.
    - *labels* = a list of labels for the profiles (same length and same order).
    - *fidkey* = key of fid for labelling the curve with *fid[fidkey]*;
                  if *None*, labels with raw fid.
    - *Ycoordinate* = label for the Y coordinate.
    - *unit* = label for X coordinate.
    - *title* = title for the plot.
    - *logscale* = to set Y logarithmic scale
    - *ema* = to make emagram-like plots of Temperature
    - *zoom*: a dict containing optional limits to zoom on the plot. \n
      Syntax: e.g. {'ymax':500, ...}.
    """
    import matplotlib.pyplot as plt

    plt.rc('font', family='serif')
    plt.rc('figure', autolayout=True)

    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'darkolivegreen',
              'yellow', 'salmon', 'black']
    linestyles = ['-', '--', '-.', ':']

    if isinstance(profiles, V1DField):
        profiles = [profiles]
    if isinstance(labels, six.string_types):
        labels = [labels]
    p0 = profiles[0]
    if p0.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100):
        reverseY = True
    else:
        reverseY = False
    if p0.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
        Y = p0.geometry.vcoordinate.levels
    if Ycoordinate is None:
        if p0.geometry.vcoordinate.typeoffirstfixedsurface == 119:
            Ycoordinate = 'Level \nHybrid-Pressure \ncoordinate'
        elif p0.geometry.vcoordinate.typeoffirstfixedsurface == 100:
            Ycoordinate = 'Pressure (hPa)'
        elif p0.geometry.vcoordinate.typeoffirstfixedsurface == 102:
            Ycoordinate = 'Altitude (m)'
        elif p0.geometry.vcoordinate.typeoffirstfixedsurface == 103:
            Ycoordinate = 'Height (m)'
        elif p0.geometry.vcoordinate.typeoffirstfixedsurface == 118:
            Ycoordinate = 'Level \nHybrid-Height \ncoordinate'
        elif p0.geometry.vcoordinate.typeoffirstfixedsurface == 109:
            Ycoordinate = 'Potential \nvortex \n(PVU)'
        else:
            Ycoordinate = 'unknown \ncoordinate'

    # Figure
    fig, ax = util.set_figax(*over, figsize=(6., 9.))
    if logscale:
        ax.set_yscale('log')
    if reverseY and not ax.yaxis_inverted():
        ax.invert_yaxis()
    for i, p in enumerate(profiles):
        if len(p.validity) != 1:
            raise epygramError("plotprofiles can handle only profiles with one validity.")
        Y = numpy.array(p.geometry.vcoordinate.levels).flatten()
        if labels is not None:
            label = labels[i]
        else:
            if fidkey is not None:
                label = p.fid.get(fidkey, p.fid)
            else:
                label = str(p.fid)
        data = p.getdata()
        if ema:
            mindata = numpy.inf
            maxdata = -numpy.inf
            alpha = 0.75
            templines = numpy.arange(round(min(data), -1) - 10,
                                     round(max(data), -1) + 10 + 50,
                                     10)
            for t in templines:
                ax.plot([t, t + (max(data) - min(data)) * alpha],
                        [max(Y), min(Y)],
                        color='grey', linestyle=':')
            ax.set_yticks(numpy.linspace(0, 1000, 11))
            data = data + abs(Y - Y[-1]) * (data.max() - data.min()) / \
                                           (Y.max() - Y.min()) * alpha
            unit = 'K'
            mindata = min(mindata, data.min())
            maxdata = max(maxdata, data.max())
        plot_kwargs = {}
        if len(profiles) > 1:
            plot_kwargs['color'] = colors[i % len(colors)]
            plot_kwargs['linestyle'] = linestyles[i // len(colors)]
        ax.plot(data.flatten(), Y.flatten(), label=label,
                **plot_kwargs)

    # Decoration
    if reverseY:
        ax.set_ylim(bottom=numpy.array(Y).max())
    else:
        ax.set_ylim(bottom=numpy.array(Y).min())
    if zoom is not None:
        if 'ymin' in zoom:
            ax.set_ylim(bottom=zoom['ymin'])
        if 'ymax' in zoom:
            ax.set_ylim(top=zoom['ymax'])
        if 'xmin' in zoom:
            ax.set_xlim(left=zoom['xmin'])
        if 'xmax' in zoom:
            ax.set_xlim(right=zoom['xmax'])
    if title is not None:
        ax.set_title(title)
    legend = ax.legend(loc='upper right', shadow=True)
    for label in legend.get_texts():
        label.set_fontsize('medium')
    ax.set_xlabel(r'$' + unit + '$')
    ax.set_ylabel(Ycoordinate)
    if ema:
        ax.grid(axis='y')
        ax.set_xlim(mindata - 10, maxdata + 10)
    else:
        ax.grid()

    return (fig, ax)


def plotanimation(profile,
                  title='__auto__',
                  repeat=False,
                  interval=1000,
                  **kwargs):
    """
    To plot a time-dependent profile as an animation.
    Returns a :class:`matplotlib.animation.FuncAnimation`.

    Args: \n
    - *profile* being a :class:`epygram.fields.V1DField`.
    - *title* = title for the plot. '__auto__' (default) will print
      the current validity of the time frame.
    - *repeat*: to repeat animation
    - *interval*: number of milliseconds between two validities

    Other kwargs passed to plotprofiles().
    """
    import matplotlib.animation as animation

    if len(profile.validity) == 1:
        raise epygramError("plotanimation can handle only profile with several validities.")

    if title is not None:
        if title == '__auto__':
            title_prefix = ''
        else:
            title_prefix = title
        title = title_prefix + '\n' + profile.validity[0].get().isoformat(sep=b' ')
    else:
        title_prefix = None
    profile0 = profile.deepcopy()
    profile0.validity = profile.validity[0]
    profile0.setdata(profile.getdata()[0, ...])
    mindata = profile.getdata().min()
    maxdata = profile.getdata().max()
    mindata -= (maxdata - mindata) / 10.
    maxdata += (maxdata - mindata) / 10.

    if kwargs.get('ema', False):
        epylog.warning("'ema' option not fully tested in animation: min/max may not be optimised.")
    zoom = kwargs.get('zoom')
    zoom = util.ifNone_emptydict(zoom)
    if 'xmax' not in zoom:
        zoom.update(xmax=maxdata)
    if 'xmin' not in zoom:
        zoom.update(xmin=mindata)
    kwargs['zoom'] = zoom

    fig, ax = plotprofiles(profile0,
                           title=title,
                           **kwargs)

    def update(i, ax=None, profile=None, title_prefix=None):
        if i < len(profile.validity):
            ax.lines[0].set_xdata(profile.getdata()[i, ...])
            if title_prefix is not None:
                ax.set_title(title_prefix + '\n' + profile.validity[i].get().isoformat(sep=b' '))

        return ax.lines[0],

    anim = animation.FuncAnimation(fig, update,
                                   fargs=[ax, profile, title_prefix],
                                   frames=list(range(len(profile.validity) + 1)),  # AM: don't really understand why but needed for the last frame to be shown
                                   interval=interval,
                                   repeat=repeat)

    return anim
