#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for a 3D field.

Plus a function to create a Vector field from 2 scalar fields.
"""

import numpy
import sys

from footprints import proxy as fpx, FPList

from epygram import epygramError
from epygram.base import Field, FieldValidityList
from . import D3Field, D3VirtualField

def make_vector_field(*components):
    """
    Creates a new :class:`epygram.D3VectorField` or subclass from
    several :class:`epygram.D3Field` or subclass representing
    the components of the vector in the field geometry.
    """
    if len(components) < 2:
        raise epygramError("One need at least two components to make a vector")
    if not all([isinstance(c, (D3Field, D3VirtualField)) for c in components]):
        raise epygramError("All components must be (subclass of) D3Field.")
    if any([components[0].geometry.dimensions != c.geometry.dimensions for c in components[1:]]):
        raise epygramError("All components must be share their gridpoint" +
                           " dimensions.")
    if any([components[0].spectral_geometry != c.spectral_geometry for c in components[1:]]):
        raise epygramError("All components must be share their spectral" +
                           " geometry.")
    if any([components[0].structure != c.structure for c in components[1:]]):
        raise epygramError("'fX', 'fY' must share their structure.")

    f = fpx.field(fid={'op':'make_vector()'},
                  structure=components[0].structure,
                  validity=components[0].validity.copy(),
                  processtype=components[0].processtype,
                  vector=True,
                  components=components)
    return f


def psikhi2uv(psi, khi):
    """
    Compute wind (on the grid) as a D3VectorField (or subclass)
    from streamfunction **psi** and velocity potential **khi**.
    """
    (dpsidx, dpsidy) = psi.compute_xy_spderivatives()
    (dkhidx, dkhidy) = khi.compute_xy_spderivatives()
    u = dkhidx - dpsidy
    v = dkhidy + dpsidx
    u.fid = {'derivative':'u-wind'}
    v.fid = {'derivative':'v-wind'}
    u.processtype = psi.processtype
    v.processtype = psi.processtype
    u.validity = psi.validity
    v.validity = psi.validity
    return make_vector_field(u, v)

def uv2psikhi(u, v):
    """
    Compute streamfunction **psi** and velocity potential **khi**
    as a D3VectorField (or subclass) from wind (on the grid).
    """
    # Check space
    if not (u.spectral and v.spectral):
        raise epygramError("Wind must be in spectral space to compute psi/khi.")

    # Compute voriticity / divergence
    uv = epygram.fields.make_vector_field(u, v)
    vor, div = uv.compute_vordiv()
    vor.gp2sp(u.spectral_geometry)
    div.gp2sp(u.spectral_geometry)

    if u.geometry.projected_geometry:
        # Prepare elliptic truncation
        M = u.spectral_geometry.truncation["in_X"]
        N = u.spectral_geometry.truncation["in_Y"]
        ellips = np.zeros((M+1),dtype=int)
        ellips[0] = N
        for jm in range(0,M-1):
            ellips[jm+1] = int(float(N)/float(M)*np.sqrt(float(M**2-(jm+1)**2))+1.0e-10)
        ellips[M] = 0
        nsefre = 0
        for jm in range(0,M+1):
            nsefre += 4*(ellips[jm]+1)
        if nsefre != np.size(u.data):
            raise epygramError("Inconsistent elliptic truncation.")

        # Prepare laplacian and inverse laplacian
        exwn = 2.0*np.pi/(u.geometry.grid["X_resolution"]*float(u.geometry.dimensions["X"]))
        eywn = 2.0*np.pi/(u.geometry.grid["Y_resolution"]*float(u.geometry.dimensions["Y"]))
        lapdir = np.zeros(nsefre)
        lapinv = np.zeros(nsefre)
        inm = 0
        for jm in range(0,M+1):
            for jn in range(0,ellips[jm]+1):
                for jq in range(0,4):
                    if jm == 0 and jn == 0:
                        lapdir[inm] = 0.0
                        lapinv[inm] = 0.0
                    elif jn == 0:
                        lapdir[inm] = -float(jm**2)*exwn**2
                        lapinv[inm] = 1.0/lapdir[inm]
                    else:
                        lapdir[inm] = -(float(jm**2)*exwn**2+float(jn**2)*eywn**2)
                        lapinv[inm] = 1.0/lapdir[inm]
                    inm += 1
    else:
        raise epygramError("uv2psikhi implemented for projected geometry only so far.")

    # Apply inverse laplacian to compute stream function and velocity potential
    psi = copy.deepcopy(u)
    khi = copy.deepcopy(v)
    psi.data = vor.data*lapinv
    khi.data = div.data*lapinv
    psi.sp2gp()
    khi.sp2gp()
    psi.fid["FA"] = u.fid["FA"].replace("WIND.U.PHYS","FONC.COURANT")
    khi.fid["FA"] = u.fid["FA"].replace("WIND.U.PHYS","POT.VITESSE")

    return psi, khi

class D3VectorField(Field):
    """
    3-Dimensions Vector field class.

    This is a wrapper to a list of D3Field(s), representing the components
    of a vector projected on its geometry (the grid axes).
    """

    _collector = ('field',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                info="Type of Field geometry.",
                values=set(['3D'])),
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
                info="List of Fields that each compose a component of the vector.",
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
        return all([c.spectral for c in self.components])

    @property
    def geometry(self):
        return self.components[0].geometry

    @property
    def components_are_projected_on(self):
        components = [f.misc_metadata.get('uvRelativeToGrid', None) for f in self.components]
        proj = components[0]
        if len(set(components)) > 1:
            proj = None
        elif proj == 0:
            proj = 'lonlat'
        elif proj == 1:
            proj = 'grid'
        return proj

    def attach_components(self, *components):
        """
        Attach components of the vector to the VectorField.
        *components* must be a series of D3Field.
        """
        for f in components:
            if not isinstance(f, D3Field):
                raise epygramError("*components* must inherit from D3Field(s).")
            if f.structure != self.structure:
                raise epygramError("*components* must share the same strucuture.")
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
        *spectral_geometry* mandatorily passed as argument. Replaces data in
        place.

        :param spectral_geometry: instance of SpectralGeometry, actually
                                  containing spectral transform subroutine (in
                                  in its own *gp2sp()* method).
        """
        for f in self.components:
            f.gp2sp(spectral_geometry=spectral_geometry)

    def getdata(self, subzone=None, **kwargs):
        """
        Returns the field data, with 2D shape if the field is not spectral,
        1D if spectral, as a tuple with data for each component.

        :param subzone: optional, among ('C', 'CI'), for LAM fields only, returns
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

    def setdata(self, data):
        """
        Sets data to its components.

        :param data: [data_i for i components]
        """
        if len(data) != len(self.components):
            raise epygramError("data must have as many components as VectorField.")
        for i in range(len(self.components)):
            self.components[i].setdata(data[i])

    def deldata(self):
        """Empties the data."""
        for i in range(len(self.components)):
            self.components[i].deldata()

    data = property(getdata, setdata, deldata, "Accessor to the field data.")

    def getlevel(self, level=None, k=None):
        """
        Returns a level of the field as a new field.

        :param level: the requested level expressed in coordinate value (Pa, m...)
        :param k: the index of the requested level
        """
        components = [comp.getlevel(level=level, k=k) for comp in self.components]
        return fpx.field(fid={'op':'getlevel()'},
                         structure=components[0].structure,
                         validity=components[0].validity.copy(),
                         processtype=components[0].processtype,
                         vector=True,
                         components=components)

    def to_module(self):
        """
        Returns a :class:`epygram.D3Field` (or subclass) whose data
        is the module of the Vector field.
        """
        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata(d4=True)
        else:
            datagp = self.getdata(d4=True)
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
        else:
            loc_sqrt = numpy.sqrt
        module = 0.
        for i in range(len(self.components)):
            module += datagp[i] ** 2
        module = loc_sqrt(module)
        f = fpx.field(geometry=self.geometry.copy(),
                      structure=self.structure,
                      fid={'op':'VectorField.to_module()'},
                      validity=self.validity.copy(),
                      processtype=self.processtype)
        f.setdata(module)
        if self.spectral:
            f.gp2sp(self.spectral_geometry)

        return f

    def compute_direction(self):
        """
        Returns a :class:`epygram.D3Field` or subclass whose data
        is the direction of the horizontal part of the Vector field
        (the two firsts components), in degrees.
        """
        if self.spectral:
            fieldcopy = self.deepcopy()
            fieldcopy.sp2gp()
            datagp = fieldcopy.getdata()
        else:
            datagp = self.getdata()
        if isinstance(datagp[0], numpy.ma.MaskedArray):
            loc_sqrt = numpy.ma.sqrt
            loc_arccos = numpy.ma.arccos
        else:
            loc_sqrt = numpy.sqrt
            loc_arccos = numpy.arccos
        module = loc_sqrt(datagp[0] ** 2 + datagp[1] ** 2)
        module_cal = numpy.where(module < 1.E-15, 1.E-15, module)
        u_norm = -datagp[0] / module_cal
        v_norm = -datagp[1] / module_cal
        numpy.clip(v_norm, -1, 1, out=v_norm)
        dd1 = loc_arccos(v_norm)
        dd2 = 2. * numpy.pi - dd1
        direction = numpy.degrees(numpy.where(u_norm >= 0., dd1, dd2))
        direction = numpy.where(module < 1.E-15, 0., direction)
        f = fpx.field(geometry=self.geometry.copy(),
                      structure=self.structure,
                      fid={'op':'VectorField.compute_direction()'},
                      validity=self.validity.copy(),
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
        Other components are kept untouched.

        :param map_factor_correction: if True, apply a correction of magnitude
                                      due to map factor.
        :param reverse: if True, apply the reverse reprojection.
        """
        (lon, lat) = self.geometry.get_lonlat_grid()
        assert not self.spectral
        u = self.components[0].getdata(d4=True)
        v = self.components[1].getdata(d4=True)
        if self.geometry.name == 'rotated_reduced_gauss':
            for t in range(u.shape[0]):
                for k in range(u.shape[1]):
                    (newu, newv) = self.geometry.reproject_wind_on_lonlat(u[t, k, ...].compressed(),
                                                                          v[t, k, ...].compressed(),
                                                                          lon.compressed(),
                                                                          lat.compressed(),
                                                                          map_factor_correction=map_factor_correction,
                                                                          reverse=reverse)
                    u[t, k, ...][~u[t, k, ...].mask] = newu
                    v[t, k, ...][~v[t, k, ...].mask] = newv
        else:
            for t in range(u.shape[0]):
                for k in range(u.shape[1]):
                    (u[t, k, ...], v[t, k, ...]) = self.geometry.reproject_wind_on_lonlat(u[t, k, ...], v[t, k, ...],
                                                                                          lon, lat,
                                                                                          map_factor_correction=map_factor_correction,
                                                                                          reverse=reverse)
        self.setdata([u, v] + [c.getdata(d4=True) for c in self.components[2:]])

    def map_factorize(self, reverse=False):
        """
        Multiply the field by its map factor.
        Only the first two components are affected.

        :param reverse: if True, divide.
        """
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

        :param divide_by_m: if True, apply f = f/m beforehand, where m is the
                            map factor.
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

    def remove_level(self, *args, **kwargs):
        """Cf. D3Field.remove_level()"""
        for component in self.components:
            component.remove_level(*args, **kwargs)

    def extract_subdomain(self, *args, **kwargs):
        """Cf. D3Field.extract_subdomain()"""
        result = make_vector_field(self.components[0].extract_subdomain(*args, **kwargs),
                                   self.components[1].extract_subdomain(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subdomain(*args, **kwargs))
        return result

    def extract_zoom(self, *args, **kwargs):
        """Cf. D3Field.extract_zoom()"""
        result = make_vector_field(self.components[0].extract_zoom(*args, **kwargs),
                                   self.components[1].extract_zoom(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_zoom(*args, **kwargs))
        return result

    def extract_subarray(self, *args, **kwargs):
        """Cf. D3Field.extract_subarray()"""
        result = make_vector_field(self.components[0].extract_subarray(*args, **kwargs),
                                   self.components[1].extract_subarray(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subarray(*args, **kwargs))
        return result

    def extract_subsample(self, *args, **kwargs):
        """Cf. D3Field.extract_subsample"""
        result = make_vector_field(self.components[0].extract_subsample(*args, **kwargs),
                                   self.components[1].extract_subsample(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.extract_subsample(*args, **kwargs))
        return result

    def resample(self, *args, **kwargs):
        """Cf. D3Field.resample()"""
        result = make_vector_field(self.components[0].resample(*args, **kwargs),
                                   self.components[1].resample(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.resample(*args, **kwargs))
        return result

    def resample_on_regularll(self, *args, **kwargs):
        """Cf. D3Field.resample_on_regularll()"""
        result = make_vector_field(self.components[0].resample_on_regularll(*args, **kwargs),
                                   self.components[1].resample_on_regularll(*args, **kwargs))
        for component in self.components[2:]:
            result.attach_components(component.resample_on_regularll(*args, **kwargs))
        return result

    def center(self, *args, **kwargs):
        """Cf. D3Field.center()"""
        for component in self.components:
            component.center(*args, **kwargs)

    def select_subzone(self, *args, **kwargs):
        """Cf. D3Field.select_subzone()"""
        for component in self.components:
            component.select_subzone(*args, **kwargs)

    def use_field_as_vcoord(self, *args, **kwargs):
        """Cf. D3Field.use_field_as_vcoord()"""
        for component in self.components:
            component.use_field_as_vcoord(*args, **kwargs)

###################
# PRE-APPLICATIVE #
###################
# (but useful and rather standard) !
# [so that, subject to continuation through updated versions,
#  including suggestions/developments by users...]

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
        Shifts the center of the geometry (and the data accordingly) by
        *longitude_shift* (in degrees). *longitude_shift* has to be a multiple
        of the grid's resolution in longitude.

        For global RegLLGeometry grids only.
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

        :param out: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        :param vertical_geometry: if True, writes the vertical geometry of the
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
        if 'vector' not in other._attributes:
            raise epygramError("cannot operate a Vector field with a" +
                               " non-Vector one.")
        else:
            if isinstance(other, self.__class__):
                if len(self.components) != len(other.components):
                    raise epygramError("vector fields must have the same" +
                                       " number of components.")
            super(D3VectorField, self)._check_operands(other)

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

        newid = {'op':'+'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
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
        newid = {'op':'*'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
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
        newid = {'op':'-'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
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
        newid = {'op':'/'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, other):
        """
        Definition of substraction, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'-'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [other.components[i] - self.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [other - self.components[i] for i in range(len(self.components))]
        newid = {'op':'-'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield

    def __rdiv__(self, other):
        """
        Definition of division, 'other' being:
        - a scalar (integer/float)
        - another Field of the same subclass.
        Returns a new Field whose data is the resulting operation,
        with 'fid' = {'op':'/'} and null validity.
        """
        if isinstance(other, self.__class__):
            newcomponents = [other.components[i] / self.components[i] for i in range(len(self.components))]
        else:
            newcomponents = [other / self.components[i] for i in range(len(self.components))]
        newid = {'op':'/'}
        newfield = fpx.field(fid=newid,
                             structure=self.structure,
                             validity=self.validity,
                             processtype=self.processtype,
                             vector=True,
                             components=newcomponents)
        return newfield
