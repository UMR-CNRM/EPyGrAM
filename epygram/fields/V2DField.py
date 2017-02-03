#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class that handle a Vertical 2D field.
"""

import numpy

import footprints

from .D3Field import D3Field
from epygram import config, util, epygramError
from epygram.geometries import V2DGeometry

epylog = footprints.loggers.getLogger(__name__)



class V2DField(D3Field):
    """
    Vertical 2-Dimension (section) field class.
    A field is defined by its identifier 'fid',
    its data, its geometry, and its validity.

    At least for now, it is designed somehow like a collection of V1DFields.
    And so is V2DGeometry.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V2D'])),
            geometry=dict(
                type=V2DGeometry),
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
                  colorbar='vertical',
                  colorbar_over=None,
                  graphicmode='colorshades',
                  minmax=None,
                  levelsnumber=21,
                  center_cmap_on_0=False,
                  colormap='jet',
                  zoom=None,
                  title=None,
                  fidkey=None,
                  logscale=False,
                  minmax_in_title=True,
                  contourcolor='k',
                  contourwidth=1,
                  contourlabel=True,
                  x_is='distance',
                  mask_threshold=None):
        """
        Makes a simple (profile) plot of the field.

        Args: \n
        - *over* = any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several 
          subplots (axes); *ax* is the matplotlib axes on 
          which the drawing is done. When given (!= None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        - *colorbar_over*: an optional existing ax to plot the colorbar on.
        - *title* = title for the plot.
        - *fidkey* = type of fid for entitling the plot with *fid[fidkey]*,
                     if title is *None*;
                     if *None*, labels with raw fid.
        - *logscale* = to set Y logarithmic scale
        - *minmax*: defines the min and max values for the plot colorbar. \n
          Syntax: [min, max]. [0.0, max] also works. Default is min/max of the
          field.
        - *graphicmode*: among ('colorshades', 'contourlines').
        - *levelsnumber*: number of levels for contours and colorbar.
        - *colormap*: name of the **matplotlib** colormap to use.
        - *center_cmap_on_0*: aligns the colormap center on the value 0.
        - *colorbar*: if *False*, hide colorbar the plot; else, befines the
          colorbar orientation, among ('horizontal', 'vertical').
          Defaults to 'vertical'.
        - *zoom*: a dict containing optional limits to zoom on the plot. \n
          Syntax: e.g. {'ymax':500, ...}.
        - *minmax_in_title*: if True and minmax != None, adds min and max
          values in title
        - *contourcolor*: color or colormap to be used for 'contourlines'
          graphicmode. It can be either a legal html color name, or a colormap
          name.
        - *contourwidth*: width of contours for 'contourlines' graphicmode.
        - *contourlabel*: displays labels on contours.
        - *existingfigure*: to plot the field over an existing figure (e.g.
          contourlines over colorshades).
          Be aware that no check is done between the *existingfigure* geometry
          and either the field's one: there might be inconsistency.
        - *x_is*: abscissa to be among:
          - 'distance': distance from first point of transect
          - 'lon': longitude of points
          - 'lat': latitude of points
        - *mask_threshold*: dict with min and/or max value(s) to mask outside.

        Warning: requires **matplotlib**.
        """

        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)
        # User colormaps
        if colormap not in plt.colormaps():
            util.add_cmap(colormap)

        # Figure, ax
        fig, ax = util.set_figax(*over)
        # coords
        z = numpy.zeros((len(self.geometry.vcoordinate.levels),
                         self.geometry.dimensions['X']))
        for k in range(len(self.geometry.vcoordinate.levels)):
            z[k, :] = self.geometry.vcoordinate.levels[k]
        x = numpy.zeros((len(self.geometry.vcoordinate.levels),
                         self.geometry.dimensions['X']))
        transect = zip(*self.geometry.get_lonlat_grid())
        p0 = transect[0]
        plast = p0
        distance = 0
        for i in range(self.geometry.dimensions['X']):
            if x_is == 'distance':
                p = transect[i]
                distance += self.geometry.distance((plast[0], plast[1]),
                                                   (p[0], p[1]))
                x[:, i] = distance
                plast = p
            elif x_is == 'lon':
                x[:, i] = transect[i][0]
            elif x_is == 'lat':
                x[:, i] = transect[i][1]
        plast = transect[-1]
        mask_outside = {'min':-config.mask_outside,
                        'max':config.mask_outside}
        if mask_threshold is not None:
            mask_outside.update(mask_threshold)
        data = numpy.ma.masked_outside(self.getdata(),
                                       mask_outside['min'],
                                       mask_outside['max'])
        #FIXME: mask_threshold does not produce expected results
        if self.geometry.vcoordinate.typeoffirstfixedsurface in (119, 100, 160):
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
        if logscale:
            ax.set_yscale('log')
        ax.grid()
        if graphicmode == 'colorshades':
            pf = ax.contourf(x, z, data, levels, cmap=colormap,
                             vmin=vmin, vmax=vmax)
            if colorbar:
                position = 'right' if colorbar == 'vertical' else 'bottom'
                if colorbar_over is None:
                    cax = make_axes_locatable(ax).append_axes(position,
                                                              size="5%",
                                                              pad=0.2)
                else:
                    cax = colorbar_over
                cb = plt.colorbar(pf,
                                  orientation=colorbar,
                                  ticks=hlevels,
                                  cax=cax)
                if minmax_in_title:
                    cb.set_label(minmax_in_title)
        elif graphicmode == 'contourlines':
            pf = ax.contour(x, z, data, levels=levels, colors=contourcolor,
                            linewidths=contourwidth)
            if contourlabel:
                ax.clabel(pf, colors=contourcolor)
        if reverseY and not ax.yaxis_inverted():
            ax.invert_yaxis()
        # decoration
        surf = z[-1, :]
        bottom = max(surf) if reverseY else min(surf)
        ax.fill_between(x[-1, :], surf, numpy.ones(len(surf)) * bottom,
                        color='k')
        if self.geometry.vcoordinate.typeoffirstfixedsurface == 119:
            Ycoordinate = 'Level \nHybrid-Pressure \ncoordinate'
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 100:
            Ycoordinate = 'Pressure (hPa)'
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 102:
            Ycoordinate = 'Altitude (m)'
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 103:
            Ycoordinate = 'Height (m)'
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 118:
            Ycoordinate = 'Level \nHybrid-Height \ncoordinate'
        elif self.geometry.vcoordinate.typeoffirstfixedsurface == 109:
            Ycoordinate = 'Potential \nvortex \n(PVU)'
        else:
            Ycoordinate = 'unknown \ncoordinate'
        if x_is == 'distance':
            ax.set_xlabel('Distance from left-end point (m).')
        elif x_is == 'lon':
            ax.set_xlabel(u'Longitude (\u00B0).')
        elif x_is == 'lat':
            ax.set_xlabel(u'Latitude (\u00B0).')
        ax.set_ylabel(Ycoordinate)
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
                fid = self.fid
            else:
                fid = self.fid[fidkey]
            title = 'Section of ' + str(fid) + ' between \n' + \
                    '<- (' + str(p0[0]) + ', ' + \
                    str(p0[1]) + ')' + ' and ' + \
                    '(' + str(plast[0]) + ', ' + \
                    str(plast[1]) + ') -> \n' + \
                    str(self.validity.get())
        ax.set_title(title)

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
            title = title_prefix + '\n' + self.validity[0].get().isoformat(sep=' ')
        else:
            title_prefix = None
        field0 = self.deepcopy()
        field0.validity = self.validity[0]
        field0.setdata(self.getdata()[0, ...])
        mindata = self.getdata().min()
        maxdata = self.getdata().max()

        minmax = kwargs.get('minmax')
        if minmax is None:
            minmax = (mindata, maxdata)
        kwargs['minmax'] = minmax
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
                    title = title_prefix + '\n' + fieldi.validity.get().isoformat(sep=' ')
                fieldi.plotfield(title=title,
                                 **kwargs)

        anim = animation.FuncAnimation(fig, update,
                                       fargs=[ax, self, field0, title_prefix, kwargs],
                                       frames=range(len(self.validity) + 1),  # AM: don't really understand why but needed for the last frame to be shown
                                       interval=interval,
                                       repeat=repeat)

        return anim
