#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class to for Point (0D == 1 value) fields.
"""

import datetime

from footprints import proxy as fpx

from .D3Field import D3Field
from epygram import config
from epygram.geometries import PointGeometry, VGeometry
from epygram.base import FieldValidity
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
                type=PointGeometry),
        )
    )

    def plotfield(self,
                  existingfigure=None,
                  title='__timerange__',
                  plot_kwargs=None,
                  legend_kwargs=None):
        """
        Plot a PointField if only it has a time dimension.
        
        Args:
        - *existingfigure*: to plot over an existing figure, must be composed
          of (fig, ax).
        - *title*: of the plot.
          No title if *None*, time range if '__timerange__'.
        - *plot_kwargs*: arguments to be passed to matplotlib.axes.Axes.plot() 
        - *legend_kwargs*: arguments to be passed to matplotlib.axes.Axes.legend()
        """

        import matplotlib.pyplot as plt
        from matplotlib import dates
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)

        plot_kwargs = ifNone_emptydict(plot_kwargs)
        legend_kwargs = ifNone_emptydict(legend_kwargs)
        if 'label' not in plot_kwargs.keys():
            plot_kwargs['label'] = str(self.fid)

        assert len(self.validity) > 1, 'only time-dimensioned PointField can be plot.'

        if existingfigure is None:
            (f, ax) = plt.subplots()
        else:
            (f, ax) = (existingfigure[0], existingfigure[1])

        validities = [v.get() for v in self.validity]
        values = self.getdata()

        ax.plot(validities, values, **plot_kwargs)
        ax.legend(**legend_kwargs)

        def set_locators(datetimerange):
            if datetimerange <= datetime.timedelta(2):
                major_locator = dates.HourLocator(interval=6)
                minor_locator = dates.HourLocator(interval=1)
            elif datetimerange <= datetime.timedelta(7):
                major_locator = dates.DayLocator(interval=1)
                minor_locator = dates.HourLocator(interval=6)
            else:
                major_locator = dates.AutoDateLocator()
                minor_locator = None
            ax.xaxis.set_major_locator(major_locator)
            ax.xaxis.set_major_formatter(dates.AutoDateFormatter(major_locator))
            ax.grid(True)
            if minor_locator is not None:
                ax.xaxis.set_minor_locator(minor_locator)
                ax.grid(True, which='minor', axis='x', color='grey')

        xmin = dates.num2date(ax.axis()[1]).replace(tzinfo=None)
        xmax = dates.num2date(ax.axis()[0]).replace(tzinfo=None)
        set_locators(xmax - xmin)

        if title is not None:
            if title == '__timerange__':
                title = '[' + xmax.isoformat() + ' to ' + xmin.isoformat() + ']'
            ax.set_title(title)

        return f, ax

def gimme_one_point(longitude, latitude,
                    field_args=None,
                    geometry_args=None,
                    vertical_geometry_args=None):
    """
    Builds an empty PointField at given (longitude, latitude).
    
    Arguments *field_args*, *geometry_args*, *vertical_geometry_args*
    can be provided as dictionaries containing kwargs for the corresponding
    objects.
    """

    field_args = ifNone_emptydict(field_args)
    geometry_args = ifNone_emptydict(geometry_args)
    vertical_geometry_args = ifNone_emptydict(vertical_geometry_args)

    kwargs_vcoord = {'position_on_grid':'mass',
                     'levels':[0],
                     'grid':{'gridlevels':[]},
                     'structure':'V',
                     'typeoffirstfixedsurface':255}
    kwargs_vcoord.update(vertical_geometry_args)
    vcoordinate = VGeometry(**kwargs_vcoord)

    kwargs_geom = {'structure':'Point',
                   'name':'unstructured',
                   'grid':{'longitudes':[longitude],
                           'latitudes':[latitude],
                           'LAMzone':None},
                   'dimensions':{'X':1, 'Y':1},
                   'vcoordinate':vcoordinate,
                   }
    kwargs_geom.update(geometry_args)
    geom = fpx.geometry(**kwargs_geom)

    kwargs_field = {'structure':'Point',
                    'geometry':geom,
                    'fid':{}}
    kwargs_field.update(field_args)
    field = fpx.field(**kwargs_field)

    return field
