#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to for Point (0D == 1 value) fields.
"""

import numpy

from footprints import proxy as fpx
from bronx.graphics.axes import set_figax, set_nice_time_axis

from .D3Field import D3Field
from epygram import config, epygramError
from epygram.geometries import Geometry, VGeometry, UnstructuredGeometry
from epygram.util import ifNone_emptydict


class PointField(D3Field):
    """
    0-Dimension (point) field class.
    A field is defined by its identifier 'fid',
    its data, its geometry, and its validity.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['Point'])),
            geometry=dict(
                type=Geometry),
        )
    )

    def plotfield(self,
                  over=(None, None),
                  title='__timerange__',
                  unit=None,
                  datefmt=None,
                  showgrid=True,
                  figsize=None,
                  plot_kwargs=None,
                  legend_kwargs=None,
                  rcparams=None):
        """
        Plot a PointField if only it has a time dimension.

        :param over: any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several
          subplots (axes); *ax* is the matplotlib axes on
          which the drawing is done. When given (is not None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        :param title: title of the plot.
          No title if *None*, time range if '__timerange__'.
        :param unit: unit for Y axis
        :param datefmt: date format to use, e.g. "%Y-%m-%d %H:%M:%S %Z"
        :param showgrid: True/False to show grid or not
        :param figsize: figure sizes in inches, e.g. (5, 8.5).
                        Default figsize is config.plotsizes.
        :param plot_kwargs: arguments to be passed to matplotlib.axes.Axes.plot()
        :param legend_kwargs: arguments to be passed to matplotlib.axes.Axes.legend()
        :param rcparams: list of (*args, **kwargs) to be passed to pyplot.rc()
                         defaults to [(('font',), dict(family='serif')),]
        """
        import matplotlib.pyplot as plt
        from matplotlib import dates
        if rcparams is None:
            rcparams = [(('font',), {'family':'serif'}),]
        for args, kwargs in rcparams:
            plt.rc(*args, **kwargs)
        if figsize is None:
            figsize = config.plotsizes

        plot_kwargs = ifNone_emptydict(plot_kwargs)
        legend_kwargs = ifNone_emptydict(legend_kwargs)
        if 'label' not in plot_kwargs:
            plot_kwargs['label'] = str(self.fid)

        assert len(self.validity) > 1, 'only time-dimensioned PointField can be plot.'

        fig, ax = set_figax(*over, figsize=figsize)

        validities = [v.get() for v in self.validity]
        basis = [v.getbasis() for v in self.validity]
        if len(set(validities)) == 1:
            if len(set(basis)) == 1:
                raise epygramError('all validities in Point time dimension are identical')
            elif len(set(basis)) == len(self.validity):
                xaxis_label = 'Basis'
                validities = basis
            else:
                raise epygramError('inconsistent length of validities basis.')
        elif len(set(validities)) == len(self.validity):
            xaxis_label = 'Validity'
        else:
            raise epygramError('inconsistent length of validities.')

        values = self.getdata()

        ax.plot(validities, values, **plot_kwargs)
        ax.legend(**legend_kwargs)

        xmin = dates.num2date(ax.axis()[0]).replace(tzinfo=None)
        xmax = dates.num2date(ax.axis()[1]).replace(tzinfo=None)
        set_nice_time_axis(ax, 'x',
                           showgrid=showgrid, datefmt=datefmt)

        if title is not None:
            if title == '__timerange__':
                title = '[' + xmin.isoformat() + ' to ' + xmax.isoformat() + ']'
            ax.set_title(title)
        ax.set_xlabel(xaxis_label)
        if unit is not None:
            ax.set_ylabel(unit)

        return fig, ax


def gimme_one_point(longitude, latitude,
                    field_args=None,
                    geometry_args=None,
                    vertical_geometry_args=None,
                    geometry_class=UnstructuredGeometry):
    """
    Builds an empty PointField at given (longitude, latitude).

    :param field_args: to be passed to field constructor, e.g. a validity
    :param geometry_args: to be passed to geometry constructor
    :param vertical_geometry_args: to be passed to vertical_geometry constructor
    :param geometry_class: geometry class tu use (default: UnstructuredGeometry)
    """
    field_args = ifNone_emptydict(field_args)
    geometry_args = ifNone_emptydict(geometry_args)
    vertical_geometry_args = ifNone_emptydict(vertical_geometry_args)

    kwargs_vcoord = {'position_on_grid':'mass',
                     'levels':[0],
                     'grid':{'gridlevels':[]},
                     'typeoffirstfixedsurface':255}
    kwargs_vcoord.update(vertical_geometry_args)
    vcoordinate = VGeometry(**kwargs_vcoord)

    kwargs_geom = {'name':'unstructured',
                   'grid':{'longitudes':[longitude],
                           'latitudes':[latitude]},
                   'dimensions':{'X':1, 'Y':1},
                   'vcoordinate':vcoordinate,
                   }
    kwargs_geom.update(geometry_args)
    geom = geometry_class(**kwargs_geom)

    kwargs_field = {'structure':'Point',
                    'geometry':geom,
                    'fid':{}}
    kwargs_field.update(field_args)
    field = fpx.field(**kwargs_field)
    zeros = numpy.zeros(len(field.validity))
    field.setdata(zeros)

    return field
