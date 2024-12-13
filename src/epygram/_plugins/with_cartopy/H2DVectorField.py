#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Extend H2DVectorField with plotting methods using cartopy.
"""

import numpy
import copy
import cartopy.crs as ccrs

import footprints

from epygram import config, epygramError

epylog = footprints.loggers.getLogger(__name__)

def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([cartoplot, cartoimage],
                               plugin_name)
    from epygram.fields import H2DVectorField
    # defaults arguments for cartopy plots
    H2DVectorField.default_streamplot_kw = dict(color='k', linewidth=2)
    H2DVectorField._cartoplot_set_figure_and_module = _cartoplot_set_figure_and_module
    H2DVectorField._cartoplot_mask = _cartoplot_mask
    H2DVectorField._cartoplot_subsample = classmethod(_cartoplot_subsample)
    H2DVectorField._cartoplot_reproject_wind = _cartoplot_reproject_wind
    H2DVectorField._cartoplot_actual_plot = classmethod(_cartoplot_actual_plot)
    H2DVectorField._cartoimage_actual_plot = classmethod(_cartoimage_actual_plot)
    H2DVectorField.cartoplot = cartoplot
    H2DVectorField.cartoimage = cartoimage


def _cartoplot_set_figure_and_module(self,
                                     map_factor_correction,
                                     **module_cartoplot_kwargs):
    """Set figure, geometry, cartography, and plot module if requested."""
    module = self.to_module()
    if 'gauss' in self.geometry.name and self.geometry.grid['dilatation_coef'] != 1.:
        if map_factor_correction:
            module.operation_with_other('*', self.geometry.map_factor_field())
        else:
            epylog.warning('check carefully *map_factor_correction* w.r.t. dilatation_coef')
    fig, ax, projection = module.cartoplot_fig_init(module_cartoplot_kwargs.pop('fig', None),
                                                    module_cartoplot_kwargs.pop('ax', None),
                                                    module_cartoplot_kwargs.pop('projection', None),
                                                    module_cartoplot_kwargs.pop('figsize', config.plotsizes),
                                                    module_cartoplot_kwargs.pop('rcparams', config.default_rcparams),
                                                    # get: to be passed to cartoplot()
                                                    module_cartoplot_kwargs.get('set_global'))
    result = module.cartoplot(fig=fig,
                              ax=ax,
                              projection=projection,
                              takeover=True,
                              **{k:v for (k, v) in module_cartoplot_kwargs.items() if k != 'takeover'})
    result = {(k if k in ('fig', 'ax') else 'module_' + k):v for (k, v) in result.items()}
    return result


def _cartoplot_mask(self,
                    mask_threshold,
                    subzone):
    """Consistently mask data, according to thresholds."""
    mask_outside = {'min':-config.mask_outside,
                    'max':config.mask_outside}
    if mask_threshold is not None:
        mask_outside.update(mask_threshold)
    if not self.geometry.grid.get('LAMzone', False):
        subzone = None
    [u, v] = [numpy.ma.masked_outside(data,
                                      mask_outside['min'],
                                      mask_outside['max']) for data in
              self.getdata(subzone=subzone)]
    # share mask with lons, lats
    (lons, lats) = self.geometry.get_lonlat_grid(subzone=subzone)
    if any([isinstance(d, numpy.ma.masked_array) for d in (u,v)]):
        assert isinstance(u, numpy.ma.masked_array) == isinstance(u, numpy.ma.masked_array)
        common_mask = u.mask + v.mask
        u.mask = common_mask
        v.mask = common_mask
        lons = numpy.ma.masked_where(common_mask, lons)
        lats = numpy.ma.masked_where(common_mask, lats)
    return u, v, lons, lats


def _cartoplot_subsample(cls,
                         u, v,
                         lons, lats,
                         subsampling):
    """Subsample data for vectors plotting."""
    if u.ndim == 1:
        u = u[::subsampling]
        v = v[::subsampling]
        lons = lons[::subsampling]
        lats = lats[::subsampling]
    else:
        u = u[::subsampling, ::subsampling]
        v = v[::subsampling, ::subsampling]
        lons = lons[::subsampling, ::subsampling]
        lats = lats[::subsampling, ::subsampling]
    return u, v, lons, lats


def _cartoplot_reproject_wind(self,
                              u, v,
                              lons, lats,
                              components_are_projected_on,
                              map_factor_correction):
    """If necessary, reproject wind on lon/lat axes before plotting."""
    # calculate the orientation of the vectors
    assert components_are_projected_on in ('grid', 'lonlat')
    if (components_are_projected_on == 'lonlat' or
        (components_are_projected_on == 'grid' and
         self.geometry.name == 'regular_lonlat')):
        u_ll, v_ll = u, v
    else:
        # wind is projected on a grid that is not lonlat: rotate to lonlat
        if numpy.ma.is_masked(lons):
            # protection within mask
            u = copy.deepcopy(u)
            v = copy.deepcopy(v)
            lons = copy.deepcopy(lons)
            lats = copy.deepcopy(lats)
            for a in (lons, lats, u, v):
                a.data[a.mask] = -999.
        # ask geometry object to reproject
        (u_ll, v_ll) = self.geometry.reproject_wind_on_lonlat(
            u.flatten(), v.flatten(),
            lons.flatten(), lats.flatten(),
            map_factor_correction=map_factor_correction)
        u_ll = u_ll.reshape(u.shape)
        v_ll = v_ll.reshape(v.shape)
    return u_ll, v_ll


def _cartoplot_actual_plot(cls,
                           ax,
                           lons, lats,
                           u, v,
                           vector_plot_method,
                           vector_plot_kwargs):
    """Actual plot of vectors."""
    if vector_plot_method == 'quiver':
        if vector_plot_kwargs is None:
            vector_plot_kwargs = dict()
        elements = ax.quiver(lons, lats, u, v, transform=ccrs.PlateCarree(),
                             **vector_plot_kwargs)
    elif vector_plot_method == 'barbs':
        if vector_plot_kwargs is None:
            vector_plot_kwargs = dict()
        elements = ax.barbs(lons, lats, u, v, transform=ccrs.PlateCarree(),
                            **vector_plot_kwargs)
    elif vector_plot_method == 'streamplot':
        if vector_plot_kwargs is None:
            vector_plot_kwargs = cls.default_streamplot_kw
        if vector_plot_kwargs.get('linewidth') == '__module__':
            vector_plot_kwargs['linewidth'] = 4 * numpy.sqrt(u ** 2 + v ** 2) / min(u.max(), v.max())
        elements = ax.streamplot(lons, lats, u, v, transform=ccrs.PlateCarree(),
                                 **vector_plot_kwargs)
    else:
        raise NotImplementedError('vector_plot_method={}'.format(vector_plot_method))
    return elements


def _cartoimage_actual_plot(cls,
                            ax,
                            colorspace,
                            crs, extent,
                            data
                           ):
    if colorspace not in ('RGB', 'RGBA'):
        import PIL.Image
        data = PIL.Image.fromarray(data.astype(numpy.uint8)) # astype needed especially if interpolation occurred before
        data = data.convert('RGBA')
        data = numpy.array(data.getdata()).reshape(data.size[1], data.size[0], 4)
        
    elements = ax.imshow(data / 255., transform=crs, extent=extent, origin='lower') # division because imshow does not recognize numpy.int64 as int!
    return elements


def cartoplot(self,
              map_factor_correction=True,
              subsampling=1,
              components_are_projected_on=None,
              vector_plot_method='quiver',
              vector_plot_kwargs=None,
              quiverkey=None,
              # takeover
              takeover=False,
              **module_plot_kwargs):
    """
    Makes a simple plot of the vector field, with a number of options.

    Arguments to build the figure, geometry, cartography and to plot the
    module of vector shall be passed as additional arguments,
    cf. H2DField.cartoplot() arguments.

    To deactivate module plotting, add *plot_method=None*.

    :param map_factor_correction: if True, applies a correction of magnitude
        to vector due to map factor.
    :param subsampling: to subsample the number of gridpoints to plot.
        Ex: *subsampling* = 10 will only plot one gridpoint upon 10.
    :param components_are_projected_on: inform the plot on which axes the
        vector components are already projected on ('grid' or 'lonlat').
        If None, look for information in the field, and raise error if not found.
    :param vector_plot_method: the matplotlib Axes method to be used to
        plot, among ('quiver', 'barbs', 'streamplot').
    :param vector_plot_kwargs: arguments to be passed to the associated
        plot method.
    :param quiverkey: to activate quiverkey, in case
        vector_plot_method='quiver': may be a dict containing arguments
        to be passed to pyplot.quiverkey().
        
    Takeover:
    
    :param takeover: give the user more access to the objects used in the
        plot, by returning a dict containing them instead of only fig/ax
    """
    if self.spectral:
        raise epygramError("please convert to gridpoint with sp2gp()" +
                           " method before plotting.")
    if components_are_projected_on is None:
        components_are_projected_on = self.components_are_projected_on
        if components_are_projected_on is None:
            raise epygramError("Unable to plot if arg **components_are_projected_on/wind_projected_on** is undefined")
    # 1/ Ask module to prepare figure, and plot module if required
    result = self._cartoplot_set_figure_and_module(map_factor_correction,
                                                   **module_plot_kwargs)
    fig, ax = result['fig'], result['ax']
    # 2/ mask data
    u, v, lons, lats = self._cartoplot_mask(module_plot_kwargs.get('mask_threshold'),
                                            module_plot_kwargs.get('subzone'))
    # 3/ subsample
    u, v, lons, lats = self._cartoplot_subsample(u, v,
                                                 lons, lats,
                                                 subsampling)
    # 4/ reproject wind on lonlat
    u, v = self._cartoplot_reproject_wind(u, v,
                                          lons, lats,
                                          components_are_projected_on,
                                          map_factor_correction)
    # 5/ actual plot
    elements = self._cartoplot_actual_plot(ax,
                                           lons, lats,
                                           u, v,
                                           vector_plot_method=vector_plot_method,
                                           vector_plot_kwargs=vector_plot_kwargs)
    result['vector_plot_elements'] = elements
    if vector_plot_method == 'quiver' and quiverkey:
        if quiverkey is True:
            quiverkey = {'X':0.95, 'Y':1.01, 'U':5, 'label':'5m/s'}
        ax.quiverkey(elements, **quiverkey)
    if module_plot_kwargs.get('title') is None:
        ax.set_title(str(self.fid) + "\n" + str(self.validity.get()))
    else:
        ax.set_title(module_plot_kwargs.get('title'))
    return result if takeover else (fig, ax)

def cartoimage(self,
               # takeover
               takeover=False,
               **plot_kwargs):
    """
    Project an image, with a number of options.
        
    :param takeover: give the user more access to the objects used in the
        plot, by returning a dict containing them instead of only fig/ax

    Arguments to build the figure, geometry, cartography and to plot the
    image shall be passed as additional arguments,
    cf. H2DField.cartoplot() arguments.

    """
    if self.spectral:
        raise epygramError("please convert to gridpoint with sp2gp()" +
                           " method before plotting.")
    # 1/ Ask first component to prepare figure
    kwargs = plot_kwargs.copy()
    kwargs['plot_method'] = None
    kwargs['takeover'] = True
    result = self._cartoplot_set_figure_and_module(False,
                                                   **kwargs)
    fig, ax = result['fig'], result['ax']
    # 2/ get crs and extension
    crs = self.geometry.default_cartopy_CRS()
    extent = self.geometry.get_cartopy_extent(subzone=plot_kwargs.get('subzone'))
    
    # 3/ determine color space
    if not all(['color' in c.fid for c in self.components]):
        raise epygramError("cartoimage can be called only on special vector fields " + \
                            "containing the different channel of a color space")
    colorspace = tuple([c.fid['color'] for c in self.components])
    colorspace = {tuple(['R', 'G', 'B']):'RGB',
                  tuple(['C', 'M', 'Y', 'K']):'CMYK',
                  tuple(['Y', 'Cb', 'Cr']):'YCbCr',
                  tuple(['L*', 'a*', 'b*']):'CIELab'}[colorspace]
    data = numpy.moveaxis(numpy.ma.array(self.getdata()), 0, -1)
    
    # 4/ actual plot
    elements = self._cartoimage_actual_plot(ax,
                                            colorspace,
                                            crs, extent,
                                            data
                                           )
    result['image_plot_elements'] = elements

    if plot_kwargs.get('title') is None:
        ax.set_title(str(self.fid) + "\n" + str(self.validity.get()))
    else:
        ax.set_title(plot_kwargs.get('title'))
    return result if takeover else (fig, ax)
