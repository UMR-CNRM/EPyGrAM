#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for a Horizontal 2D field.

Plus a function to create a Vector field from 2 scalar fields.
"""

import numpy
import sys

from footprints import proxy as fpx, FPList

from epygram import config, epygramError, util, epylog
from epygram.base import Field, FieldValidityList
from . import H2DField



def make_vector_field(fX, fY):
    """
    Creates a new :class:`epygram.H2DVectorField` from two
    :class:`epygram.H2DField` *fX, fY* representing resp.
    the X and Y components of the vector in the field geometry.
    """

    if not isinstance(fX, H2DField) or not isinstance(fY, H2DField):
        raise epygramError("'fX', 'fY' must be H2DField.")
    if fX.geometry.dimensions != fY.geometry.dimensions:
        raise epygramError("'fX', 'fY' must be share their gridpoint" + \
                           " dimensions.")
    if fX.spectral_geometry != fY.spectral_geometry:
        raise epygramError("'fX', 'fY' must be share their spectral" + \
                           " geometry.")
    if fX.structure != fY.structure:
        raise epygramError("'fX', 'fY' must share their structure.")

    f = fpx.field(fid={'op':'make_vector()'},
                  structure=fX.structure,
                  validity=fX.validity,
                  processtype=fX.processtype,
                  vector=True,
                  components=[fX, fY])

    return f



class H2DVectorField(Field):
    """
    Horizontal 2-Dimensions Vector field class.

    This is a wrapper to a list of H2DField(s), representing the components
    of a vector projected on its geometry (the grid axes).

    The attribute H2DVectorField.data is a proxy to the data of the components.
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                info="Type of Field geometry.",
                values=set(['H2D'])),
            vector=dict(
                info="Intrinsic vectorial nature of the field.",
                type=bool,
                values=set([True])),
            validity=dict(
                info="Validity of the field.",
                type=FieldValidityList,
                optional=True,
                access='rwx',
                default=FieldValidityList()),
            components=dict(
                info="List of Fields that each compose a component of the \
                      vector.",
                type=FPList,
                optional=True,
                default=FPList([])),
            processtype=dict(
                optional=True,
                info="Generating process.")
        )
    )

##############
# ABOUT DATA #
##############

    @property
    def spectral_geometry(self):
        return self.components[0].spectral_geometry

    @property
    def spectral(self):
        """Returns True if the field is spectral."""
        return self.spectral_geometry != None

    @property
    def geometry(self):
        return self.components[0].geometry

    def attach_components(self, *components):
        """
        Attach components of the vector to the VectorField.
        *components* must be a series of H2DField.
        """

        for f in components:
            if not isinstance(f, H2DField):
                raise epygramError("*components* must be H2DField(s).")
        for f in components:
            self.components.append(f)

    def sp2gp(self):
        """
        Transforms the spectral field into gridpoint, according to its spectral
        geometry. Replaces data in place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *sp2gp()* method.
        """

        for f in self.components:
            f.sp2gp()

    def gp2sp(self, spectral_geometry=None):
        """
        Transforms the gridpoint field into spectral space, according to the
        spectral geometry mandatorily passed as argument. Replaces data in
        place.

        The spectral transform subroutine is actually included in the spectral
        geometry's *gp2sp()* method.
        """

        for f in self.components:
            f.gp2sp(spectral_geometry=spectral_geometry)

    def getdata(self, subzone=None, **kwargs):
        """
        Returns the field data, with 2D shape if the field is not spectral,
        1D if spectral, as a tuple with data for each component.

        - subzone: optional, among ('C', 'CI'), for LAM fields only, returns
          the data resp. on the C or C+I zone.
          Default is no subzone, i.e. the whole field.

        Shape of 2D data: (x (0) being the X component, y (1) the Y one) \n
        - Rectangular grids:\n
          grid[0,0,x] is SW, grid[-1,-1,x] is NE \n
          grid[0,-1,x] is SE, grid[-1,0,x] is NW
        - Gauss grids:\n
          grid[0,:Nj,x] is first (Northern) band of latitude, masked
          after Nj = number of longitudes for latitude j \n
          grid[-1,:Nj,x] is last (Southern) band of latitude (idem).
        """
        return [f.getdata(subzone=subzone, **kwargs) for f in self.components]

    @property
    def data(self):
        return self.getdata(d4=True)  #TOBECHECKED: overload of a footprints attribute seems not to work

    def setdata(self, data):
        """
        Sets data, checking it to be:

        - 2D if not spectral,
        - 1D if spectral.
        
        data = (data_i for i components)
        """

        if len(data) != len(self.components):
            raise epygramError("data must have as many components as VectorField.")
        for i in range(len(self.components)):
            self.components[i].setdata(data[i])

    def to_module(self):
        """
        Returns a :class:`epygram.H2DField` whose data is the module of the
        Vector field.
        """

        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata()
        else:
            datagp = self.getdata()
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
        else:
            loc_sqrt = numpy.sqrt
        module = loc_sqrt(datagp[0] ** 2 + datagp[1] ** 2)
        f = fpx.field(geometry=self.geometry,
                      structure=self.structure,
                      fid={'op':'H2DVectorField.to_module()'},
                      validity=self.validity,
                      processtype=self.processtype)
        f.setdata(module)
        if self.spectral:
            f.gp2sp(self.spectral_geometry)

        return f

    def compute_direction(self):
        """
        Returns a :class:`epygram.H2DField` whose data is the direction of the
        Vector field, in degrees. E.g. 45° is a South-Westerly vector.
        """
        #TOBECHECKED: not tested !
        epylog.warning("this method has not been tested ! It may not produce expected results.")
        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata()
        else:
            datagp = self.getdata()
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
            loc_arctan2 = numpy.ma.arctan2
        else:
            loc_sqrt = numpy.sqrt
            loc_arctan2 = numpy.arctan2
        module = loc_sqrt(datagp[0] ** 2 + datagp[1] ** 2)
        direction = -loc_arctan2(datagp[1] / module, datagp[0] / module) + numpy.pi / 2.
        direction *= 180. / numpy.pi + 360.
        direction = numpy.mod(direction, 360.)
        f = fpx.field(geometry=self.geometry,
                      structure=self.structure,
                      fid={'op':'H2DVectorField.compute_direction()'},
                      validity=self.validity,
                      processtype=self.processtype)
        f.setdata(direction)
        if self.spectral:
            f.gp2sp(self.spectral_geometry)

        return f

    def reproject_wind_on_lonlat(self,
                                 map_factor_correction=True,
                                 reverse=False):
        """
        Reprojects a wind vector (u, v) from the grid axes onto real
        sphere, i.e. with components on true zonal/meridian axes.
        
        If *map_factor_correction*, applies a
        correction of magnitude due to map factor.
        
        If *reverse*, apply the reverse reprojection.
        """

        #assert self.geometry.name == 'rotated_reduced_gauss', \
        #       "method only available for Gauss geometry."
        (lon, lat) = self.geometry.get_lonlat_grid()
        u = self.components[0].getdata()
        v = self.components[1].getdata()
        if self.geometry.name == 'rotated_reduced_gauss':
            (u, v) = self.geometry.reproject_wind_on_lonlat(u.compressed(),
                                                            v.compressed(),
                                                            lon.compressed(),
                                                            lat.compressed(),
                                                            map_factor_correction=map_factor_correction,
                                                            reverse=reverse)
        else:
            (u, v) = self.geometry.reproject_wind_on_lonlat(u, v, lon, lat,
                                                            map_factor_correction=map_factor_correction,
                                                            reverse=reverse)
        u = self.geometry.reshape_data(u, len(self.validity))
        v = self.geometry.reshape_data(v, len(self.validity))
        self.setdata([u, v])

    def map_factorize(self, reverse=False):
        """Multiply the field by its map factor. If *reverse, divide."""
        if self.spectral:
            spgeom = self.spectral_geometry
            self.sp2gp()
            was_spectral = True
        else:
            was_spectral = False
        m = self.geometry.map_factor_field()
        if reverse:
            op = '/'
        else:
            op = '*'
        self.components[0].operation_with_other(op, m)
        self.components[1].operation_with_other(op, m)
        if was_spectral:
            self.gp2sp(spgeom)

    def compute_vordiv(self, divide_by_m=False):
        """
        Compute vorticity and divergence fields from the vector field.
        
        If *divide_by_m* is True, apply f = f/m beforehand, where m is the map
        factor.
        """

        if divide_by_m:
            field = self.deepcopy()
            field.map_factorize(reverse=True)
        else:
            field = self
        (dudx, dudy) = field.components[0].compute_xy_spderivatives()
        (dvdx, dvdy) = field.components[1].compute_xy_spderivatives()
        vor = dvdx - dudy
        div = dudx + dvdy
        vor.fid = {'derivative':'vorticity'}
        div.fid = {'derivative':'divergence'}
        vor.validity = dudx.validity
        div.validity = dudx.validity

        return (vor, div)

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

    def plotfield(self,
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
                  map_factor_correction=True):
        """
        Makes a simple plot of the field, with a number of options.

        Options: \n
        - *over*: to plot the vectors over an existing figure
          (e.g. colorshades).
          Any existing figure and/or ax to be used for the
          plot, given as a tuple (fig, ax), with None for
          missing objects. *fig* is the frame of the
          matplotlib figure, containing eventually several 
          subplots (axes); *ax* is the matplotlib axes on 
          which the drawing is done. When given (!= None),
          these objects must be coherent, i.e. ax being one of
          the fig axes.
        - *subzone*: among ('C', 'CI'), for LAM fields only, plots the data
          resp. on the C or C+I zone. \n
          Default is no subzone, i.e. the whole field.
        - *gisquality*: among ('c', 'l', 'i', 'h', 'f') -- by increasing
          quality. Defines the quality for GIS elements (coastlines, countries
          boundaries...). Default is 'i'. Cf. 'basemap' doc for more details.
        - *specificproj*: enables to make basemap on the specified projection,
          among: 'kav7', 'cyl', 'ortho', ('nsper', {...}) (cf. Basemap doc). \n 
          In 'nsper' case, the {} may contain:\n
          - 'sat_height' = satellite height in km;
          - 'lon' = longitude of nadir in degrees;
          - 'lat' = latitude of nadir in degrees. \n
          Overwritten by *zoom*.
        - *zoom*: specifies the lon/lat borders of the map, implying hereby
          a 'cyl' projection.
          Must be a dict(lonmin=, lonmax=, latmin=, latmax=).\n
          Overwrites *specificproj*.
        - *use_basemap*: a basemap.Basemap object used to handle the
          projection of the map. If given, the map projection
          options (*specificproj*, *zoom*, *gisquality* ...)
          are ignored, keeping the properties of the
          *use_basemap* object. (because making Basemap is the most
          time-consuming step).
        - *drawrivers*: to add rivers on map.
        - *departments*: if True, adds the french departments on map (instead 
          of countries).
        - *boundariescolor*: color of lines for boundaries (countries,
          departments, coastlines)
        - *drawcoastlines*: to add coast lines on map.
        - *drawcountries*: to add countries on map.
        - *title*: title for the plot. Default is field identifier.
        - *meridians* and *parallels* enable to fine-tune the choice of lines to
          plot, with either:\n
          - 'auto': automatic scaling to the basemap extents
          - 'default': range(0,360,10) and range(-90,90,10)
          - a list of values
          - a grid step, e.g. 5 to plot each 5 degree.
          - None: no one is plot
          - *meridian* == 'greenwich' // 'datechange' // 'greenwich+datechange'
            *parallel* == 'equator' // 'polarcircles' // 'tropics' or any
            combination (+) will plot only these.
        - *subsampling*: to subsample the number of gridpoints to plot.
          Ex: *subsampling* = 10 will only plot one gridpoint upon 10.
        - *symbol*: among ('barbs', 'arrows', 'stream')
        - *symbol_options*: a dict of options to be passed to **barbs** or
          **quiver** method.
        - *plot_module*: to plot module as colorshades behind vectors.
        - *plot_module_options*: options (dict) to be passed to module.plotfield().
        - *bluemarble*: if > 0.0 (and <=1.0), displays NASA's "blue marble"
          as background. The numerical value sets its transparency.
        - *background*: if True, set a background color to 
          continents and oceans.
        - *quiverkey*: to activate quiverkey; must contain arguments to be
          passed to pyplot.quiverkey(), as a dict.
        - *components_are_projected_on*: inform the plot on which axes the
          vector components are projected on ('grid' or 'lonlat').
        - *map_factor_correction*: if True, applies a correction of magnitude
          to vector due to map factor.
        
        This method uses (hence requires) 'matplotlib' and 'basemap' libraries.
        """
        import matplotlib.pyplot as plt
        plt.rc('font', family='serif')
        plt.rc('figure', figsize=config.plotsizes)

        plot_module_options = util.ifNone_emptydict(plot_module_options)
        quiver_options = util.ifNone_emptydict(quiver_options)

        if self.spectral:
            raise epygramError("please convert to gridpoint with sp2gp()" + \
                               " method before plotting.")

        # 1. Figure, ax
        if not plot_module:
            fig, ax = util.set_figax(*over)

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
                    module.operation_with_other('*', self.geometry.map_factor_field())
                fig, ax = module.plotfield(use_basemap=bm,
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
        (lons, lats) = self.geometry.get_lonlat_grid(subzone=subzone)
        lons = lons[::subsampling, ::subsampling]
        lats = lats[::subsampling, ::subsampling]
        x, y = bm(lons, lats)
        data = [numpy.ma.masked_outside(data,
                                        - config.mask_outside,
                                        config.mask_outside) for data in
                self.getdata(subzone=subzone)]
        u = data[0][::subsampling, ::subsampling]
        v = data[1][::subsampling, ::subsampling]

        # Calculate the orientation of the vectors
        assert components_are_projected_on in ('grid', 'lonlat')
        if components_are_projected_on == 'grid' and 'gauss' not in self.geometry.name \
           and (specificproj == None and zoom == None):
            # map has same projection than components: no rotation necessary
            u_map = u
            v_map = v
        else:
            # (1or2) rotation(s) is(are) necessary
            if components_are_projected_on == 'lonlat' or self.geometry.name == 'regular_lonlat':
                (u_ll, v_ll) = (u, v)
            else:
                # wind is projected on a grid that is not lonlat: rotate to lonlat
                (u_ll, v_ll) = self.geometry.reproject_wind_on_lonlat(util.stretch_array(u), util.stretch_array(v),
                                                                      util.stretch_array(lons), util.stretch_array(lats),
                                                                      map_factor_correction=map_factor_correction)
            # rotate from lonlat to map projection
            (u_map, v_map) = bm.rotate_vector(u_ll, v_ll, util.stretch_array(lons), util.stretch_array(lats))
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
                xf = util.stretch_array(x)
                yf = util.stretch_array(y)
                u = util.stretch_array(u_map)
                v = util.stretch_array(v_map)
        if symbol == 'barbs':
            bm.barbs(xf, yf, u, v, ax=ax, **symbol_options)
        elif symbol == 'arrows':
            q = bm.quiver(xf, yf, u, v, ax=ax, **symbol_options)
            if quiverkey:
                ax.quiverkey(q, **quiverkey)
        elif symbol == 'stream':
            bm.streamplot(xf, yf, u, v, ax=ax, linewidth=speed_width, **symbol_options)
        if title == None:
            ax.set_title(str(self.fid) + "\n" + str(self.validity.get()))
        else:
            ax.set_title(title)

        return (fig, ax)

    def getvalue_ij(self, *args, **kwargs):
        """
        Returns the value of the different components of the field from indexes.
        """
        return [f.getvalue_ij(*args, **kwargs) for f in self.components]

    def getvalue_ll(self, *args, **kwargs):
        """
        Returns the value of the different components of the field from coordinates.
        """
        return [f.getvalue_ll(*args, **kwargs) for f in self.components]

    def min(self, subzone=None):
        """Returns the minimum value of data."""
        return [f.min(subzone=subzone) for f in self.components]

    def max(self, subzone=None):
        """Returns the maximum value of data."""
        return [f.max(subzone=subzone) for f in self.components]

    def mean(self, subzone=None):
        """Returns the mean value of data."""
        return [f.mean(subzone=subzone) for f in self.components]

    def std(self, subzone=None):
        """Returns the standard deviation of data."""
        return [f.std(subzone=subzone) for f in self.components]

    def quadmean(self, subzone=None):
        """Returns the quadratic mean of data."""
        return [f.quadmean(subzone=subzone) for f in self.components]

    def nonzero(self, subzone=None):
        """
        Returns the number of non-zero values (whose absolute
        value > config.epsilon).
        """
        return [f.nonzero(subzone=subzone) for f in self.components]

    def global_shift_center(self, longitude_shift):
        """
        For global RegLLGeometry grids only !
        Shifts the center of the geometry (and the data accordingly) by
        *longitude_shift* (in degrees). *longitude_shift* has to be a multiple
        of the grid's resolution in longitude.
        """

        if self.geometry.name != 'regular_lonlat':
            raise epygramError("only for regular lonlat geometries.")
        for f in self.components:
            f.global_shift_center(longitude_shift)

    def what(self, out=sys.stdout,
             vertical_geometry=True,
             cumulativeduration=True):
        """
        Writes in file a summary of the field.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *vertical_geometry*: if True, writes the vertical geometry of the
          field.
        """
        for f in self.components:
            f.what(out,
                   vertical_geometry=vertical_geometry,
                   cumulativeduration=cumulativeduration)

#############
# OPERATORS #
#############

    def _check_operands(self, other):
        """
        Internal method to check compatibility of terms in operations on fields.
        """

        if not 'vector' in other._attributes:
            raise epygramError("cannot operate a Vector field with a" + \
                               " non-Vector one.")
        else:
            if isinstance(other, self.__class__):
                if len(self.components) != len(other.components):
                    raise epygramError("vector fields must have the same" + \
                                       " number of components.")
            super(H2DVectorField, self)._check_operands(other)

    def __add__(self, other):
        """
        Definition of addition, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'+'} and null validity.
        """

        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] + other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] + other for i in range(len(self.components))]
        newfield = self._add(other,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __mul__(self, other):
        """
        Definition of multiplication, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'*'} and null validity.
        """

        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] * other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] * other for i in range(len(self.components))]
        newfield = self._mul(other,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __sub__(self, other):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """

        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] - other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] - other for i in range(len(self.components))]
        newfield = self._sub(other,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __div__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """

        if isinstance(other, self.__class__):
            newcomponents = [self.components[i] / other.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [self.components[i] / other for i in range(len(self.components))]
        newfield = self._div(other,
                             vector=True,
                             components=newcomponents)
        return newfield
