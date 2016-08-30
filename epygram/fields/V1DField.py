#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Vertical 1D field.
"""

import datetime
import numpy

from .D3Field import D3Field
from epygram import epygramError, config, util
from epygram.base import FieldSet
from epygram.geometries import V1DGeometry
import footprints
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
                  datefmt="%Y%m%d %H:%M:%S %Z",
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

        if force_mode != False:
            mode = force_mode
        else:
            if len(self.validity) == 1:
                mode = "profile"
            else:
                mode = "hovmoller"

        useless_args = []
        if mode == 'profile':
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
            return plotprofiles(self,
                                labels=labels,
                                fidkey=fidkey,
                                Ycoordinate=Ycoordinate,
                                unit=unit,
                                title=title,
                                logscale=logscale,
                                ema=ema,
                                zoom=zoom)
        elif mode == 'hovmoller':
            if labels != None: useless_args.append('labels')
            if unit != 'SI': useless_args.append('unit')
            if ema != False: useless_args.append('ema')
            if linecolor != 'black':useless_args.append('linecolor')
            if linestyle != '-':useless_args.append('linestyle')
            if repeat != False: useless_args.append('repeat')
            if interval != 1000: useless_args.append('interval')
            return plotverticalhovmoller(self,
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

#################
### FUNCTIONS ###
#################

def plotverticalhovmoller(profile,
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
                          datefmt="%Y%m%d %H:%M:%S %Z"):
        """
        Makes a simple vertical Hovmöller plot of the field.

        Args: \n
        - *profile* being a :class:`epygram.fields.V1DField`
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
        - *minmax_in_title*: if True and minmax != None, adds min and max
          values in title
        - *contourcolor*: color or colormap to be used for 'contourlines'
          graphicmode. It can be either a legal html color name, or a colormap
          name.
        - *contourwidth*: width of contours for 'contourlines' graphicmode.
        - *contourlabel*: displays labels on contours.
        - *datefmt*: date format to use

        Warning: requires **matplotlib**.
        """

        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)
        # User colormaps
        if colormap not in plt.colormaps():
            util.add_cmap(colormap)

        f = plt.figure()
        # coords
        z = numpy.zeros((len(profile.validity),
                         len(profile.geometry.vcoordinate.levels)))
        for k in range(len(profile.geometry.vcoordinate.levels)):
            z[:, k] = profile.geometry.vcoordinate.levels[k]
        x = numpy.zeros((len(profile.validity),
                         len(profile.geometry.vcoordinate.levels)))
        validities = {profile.validity[i].get():i for i in range(len(profile.validity))}
        date = 'Validity'
        if len(validities) == 1 and len(profile.validity) != 1:
            date = 'Basis'
            validities = {profile.validity[i].getbasis():i for i in len(profile.validity)}
        epoch = datetime.datetime(1970, 1, 1)
        for i in range(len(profile.validity)):
            d = profile.validity[i].get() if date == 'Validity' else profile.validity[i].getbasis()
            timedelta = d - epoch
            p = (timedelta.microseconds + (timedelta.seconds + timedelta.days * 24 * 3600) * 10 ** 6) / 1e6
            x[i, :] = matplotlib.dates.epoch2num(p)
        data = profile.data
        if profile.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100):
            reverseY = True
        else:
            reverseY = False
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
        if reverseY:
            plt.gca().invert_yaxis()
        if logscale:
            f.axes[0].set_yscale('log')
        plt.grid()
        if graphicmode == 'colorshades':
            pf = plt.contourf(x, z, data, levels, cmap=colormap,
                              vmin=vmin, vmax=vmax)
            if colorbar:
                cb = plt.colorbar(pf, orientation=colorbar, ticks=hlevels)
                if minmax_in_title != '':
                    cb.set_label(minmax_in_title)
        elif graphicmode == 'contourlines':
            pf = plt.contour(x, z, data, levels=levels, colors=contourcolor,
                             linewidths=contourwidth)
            if contourlabel:
                f.axes[0].clabel(pf, colors=contourcolor)
        f.axes[0].xaxis.set_major_formatter(mdates.DateFormatter(datefmt))
        matplotlib.pyplot.xticks(rotation='vertical')
        # decoration
        surf = z[-1, :]
        bottom = max(surf) if reverseY else min(surf)
        plt.fill_between(x[-1, :], surf, numpy.ones(len(surf)) * bottom,
                         color='k')
        if Ycoordinate == None:
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
        f.axes[0].set_xlabel(date + ' date.')
        f.axes[0].set_ylabel(Ycoordinate)
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
            f.axes[0].set_ylim(**ykw)
            f.axes[0].set_xlim(**xkw)
        if title is None:
            if fidkey is None:
                fid = profile.fid
            else:
                fid = profile.fid[fidkey]
            title = u'Vertical Hovmöller of ' + str(fid)
        f.axes[0].set_title(title)

        return f

def plotprofiles(profiles,
                 labels=None,
                 fidkey=None,
                 Ycoordinate=None,
                 unit='SI',
                 title=None,
                 logscale=False,
                 ema=False,
                 zoom=None):
    """
    To plot a series of profiles. Returns a :mod:`matplotlib` *Figure*.

    Args: \n
    - *profiles* being a :class:`epygram.base.FieldSet` of
      :class:`epygram.fields.V1DField`, or a single
      :class:`epygram.fields.V1DField`. \n
      All profiles are supposed to have the same unit, and the same vertical
      coordinate.
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

    if not isinstance(profiles, FieldSet):
        p = profiles.deepcopy()
        profiles = FieldSet()
        profiles.append(p)
    p0 = profiles[0]
    if p0.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100):
        reverseY = True
    else:
        reverseY = False
    if p0.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
        Y = p0.geometry.vcoordinate.levels
    if Ycoordinate == None:
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
    fig, ax = plt.subplots(figsize=(6., 9.))
    if logscale:
        ax.set_yscale('log')
    if reverseY:
        plt.gca().invert_yaxis()
    i = 0
    for p in profiles:
        if len(p.validity) != 1:
            raise epygramError("plotprofiles can handle only profiles with one validity.")
        Y = numpy.array(p.geometry.vcoordinate.levels).flatten()
        if labels != None:
            label = labels[i]
        else:
            if fidkey != None:
                label = p.fid.get(fidkey, p.fid)
            else:
                label = str(p.fid)
        data = p.data
        if ema:
            mindata = numpy.inf
            maxdata = -numpy.inf
            alpha = 0.75
            templines = numpy.arange(round(min(data), -1) - 10,
                                     round(max(data), -1) + 10 + 50,
                                     10)
            for t in templines:
                plt.plot([t, t + (max(data) - min(data)) * alpha],
                         [max(Y), min(Y)],
                         color='grey', linestyle=':')
            ax.set_yticks(numpy.linspace(0, 1000, 11))
            data = data + abs(Y - Y[-1]) * (data.max() - data.min()) / \
                                           (Y.max() - Y.min()) * alpha
            unit = 'K'
            mindata = min(mindata, data.min())
            maxdata = max(maxdata, data.max())
        plt.plot(data.flatten(), Y.flatten(), label=label,
                 color=colors[i % len(colors)],
                 linestyle=linestyles[i // len(colors)])
        i += 1

    # Decoration
    if reverseY:
        ax.set_ylim(bottom=numpy.array(Y).max())
    else:
        ax.set_ylim(bottom=numpy.array(Y).min())
    if zoom != None:
        if 'ymin' in zoom.keys():
            ax.set_ylim(bottom=zoom['ymin'])
        if 'ymax' in zoom.keys():
            ax.set_ylim(top=zoom['ymax'])
    if title != None:
        ax.set_title(title)
    legend = ax.legend(loc='upper right', shadow=True)
    for label in legend.get_texts():
        label.set_fontsize('medium')
    ax.set_xlabel(r'$' + unit + '$')
    ax.set_ylabel(Ycoordinate)
    if ema:
        plt.grid(axis='y')
        plt.xlim(mindata - 10, maxdata + 10)
    else:
        plt.grid()

    return fig

def plotanimation(profile,
                 fidkey=None,
                 Ycoordinate=None,
                 unit='SI',
                 title='__auto__',
                 logscale=False,
                 ema=False,
                 zoom=None,
                 linecolor='black',
                 linestyle='-',
                 repeat=False,
                 interval=1000):
    """
    To plot a time-dependent profile as an animation. Returns a :mod:`matplotlib` *Figure*.

    Args: \n
    - *profile* being a :class:`epygram.fields.V1DField`.
    - *fidkey* = key of fid for labelling the curve with *fid[fidkey]*;
                  if *None*, labels with raw fid.
    - *Ycoordinate* = label for the Y coordinate.
    - *unit* = label for X coordinate.
    - *title* = title for the plot.
    - *logscale* = to set Y logarithmic scale
    - *ema* = to make emagram-like plots of Temperature
    - *zoom*: a dict containing optional limits to zoom on the plot. \n
      Syntax: e.g. {'ymax':500, ...}.
    - *linecolor*: color of the profile
    - *linestyle*: line style to use for plotting the profile
    - *repeat*: to repeat animation
    - *interval*: number of milliseconds between two validities
    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    plt.rc('font', family='serif')
    plt.rc('figure', autolayout=True)

    if len(profile.validity) == 1:
        raise epygramError("plotanimation can handle only profile with several validities.")

    if profile.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100):
        reverseY = True
    else:
        reverseY = False
    if profile.geometry.vcoordinate.typeoffirstfixedsurface in (118, 119):
        Y = profile.geometry.vcoordinate.levels
    if Ycoordinate == None:
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

    # Figure
    fig, ax = plt.subplots(figsize=(6., 9.))
    if logscale:
        ax.set_yscale('log')
    if reverseY:
        plt.gca().invert_yaxis()
    if profile.geometry.vcoordinate.typeoffirstfixedsurface == 100:
        Y = numpy.array(profile.geometry.vcoordinate.levels)  # / 100
    elif profile.geometry.vcoordinate.typeoffirstfixedsurface in (102, 103, 109):
        Y = numpy.array(profile.geometry.vcoordinate.levels)
    if fidkey != None:
        label = profile.fid.get(fidkey, profile.fid)
    else:
        label = str(profile.fid)
    data = profile.data
    if ema:
        raise epygramError("ema=True not tested for animations")
        mindata = numpy.inf
        maxdata = -numpy.inf
        alpha = 0.75
        templines = numpy.arange(round(min(data), -1) - 10,
                                 round(max(data), -1) + 10 + 50,
                                 10)
        for t in templines:
            plt.plot([t, t + (max(data) - min(data)) * alpha],
                     [max(Y), min(Y)],
                     color='grey', linestyle=':')
        ax.set_yticks(numpy.linspace(0, 1000, 11))
        data = data + abs(Y - Y[-1]) * (max(data) - min(data)) / \
                                       (max(Y) - min(Y)) * alpha
        unit = 'K'
        mindata = min(mindata, min(data))
        maxdata = max(maxdata, max(data))
    plot, = ax.plot(data[0], Y, label=label, color=linecolor, linestyle=linestyle)

    # Decoration
    if reverseY:
        ax.set_ylim(bottom=numpy.array(Y).max())
    else:
        ax.set_ylim(bottom=numpy.array(Y).min())
    if zoom != None:
        if 'ymin' in zoom.keys():
            ax.set_ylim(bottom=zoom['ymin'])
        if 'ymax' in zoom.keys():
            ax.set_ylim(top=zoom['ymax'])
    if title is not None:
        if title == '__auto__':
            _title = ''
        else:
            _title = title
        ax.set_title(_title + '\n' + profile.validity[0].get().isoformat())
    legend = ax.legend(loc='upper right', shadow=True)
    for label in legend.get_texts():
        label.set_fontsize('medium')
    ax.set_xlabel(r'$' + unit + '$')
    ax.set_ylabel(Ycoordinate)
    if ema:
        ax.grid(axis='y')
        ax.xlim(mindata - 10, maxdata + 10)
    else:
        ax.grid()

    def update(i):
        plot.set_xdata(data[i])
        if title is not None:
            ax.set_title(_title + '\n' + profile.validity[i].get().isoformat())
        return plot,

    ani = animation.FuncAnimation(fig, update, range(len(data)), interval=interval, repeat=repeat)

    return ani
