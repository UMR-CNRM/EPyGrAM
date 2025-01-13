#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for Vertical geometry of fields.
"""

import numpy
import sys

from footprints import proxy as fpx

from epygram import profiles, epygramError, config
from epygram.util import (RecursiveObject, write_formatted, separation_line,
                          ifNone_emptydict)
from epygram.geometries import VGeometry


class VGeometry(RecursiveObject):
    """
    Handles the vertical geometry for fields.

    Here, the grid defines the vertical position of each level
    between a bottom and a top positions.

    The position of points w/r to the vertical grid (mass or flux points),
    is interpreted as:\n
    - mass: points are located on same levels as the grid points.
    - flux: points are located on half-levels, hence are N+1.

    levels is a list with one item for each level represented in data.
    Each item can be:\n
        - a scalar (constant value for all the data point),
        - an array with the horizontal geographic shape (level constant in time
          but varying on the horizontal),
        - an array with the first dimension corresponding to the validity
          lengthy and other dimensions to represent the horizotal.
    It is not allowed to have a level varying in time and constant on the geographic domain.
    """

    def __init__(self, typeoffirstfixedsurface, levels,
                 typeofsecondfixedsurface=None,
                 toplevel=None, bottomlevel=None,
                 grid=None, position_on_grid='__unknown__'):
        """
        :param str structure: must be 'V'
        :param int typeoffirstfixedsurface: Type of horizontal level, as of GRIB2 norm
                                            (inspired from GRIB_API).
        :param list levels: Effective levels on which data is available.
        :param dict grid: Handles description of the vertical grid.
        :param str position_on_grid: Position of points w/r to the vertical grid.
                                     (among ['mass', 'flux', '__unknown__'])
        """
        self.add_attr_int('typeoffirstfixedsurface')
        self.add_attr_int('typeofsecondfixedsurface')
        self.add_attr_list('levels')
        self.add_attr_dict('grid')
        self.add_attr_inlist('position_on_grid', ['mass', 'flux', '__unknown__'])

        self.typeoffirstfixedsurface = typeoffirstfixedsurface
        if typeofsecondfixedsurface is not None:
            self.typeofsecondfixedsurface = typeofsecondfixedsurface
        self.levels = levels
        self.grid = ifNone_emptydict(grid)
        self.position_on_grid = position_on_grid

    def what(self, out=sys.stdout, levels=True):
        """
        Writes in file a summary of the geometry.

        :param out: the output open file-like object
        :param levels: if True, writes the levels of the geometry
        """
        out.write("#########################\n")
        out.write("### VERTICAL GEOMETRY ###\n")
        out.write("#########################\n")
        write_formatted(out, "Number of vertical levels",
                        len(self.levels))
        write_formatted(out, "Type of first fixed surface",
                        self.typeoffirstfixedsurface)
        if levels:
            if len(self.levels) == 1:
                write_formatted(out, "Level", self.levels[0])
            else:
                write_formatted(out, "Levels", self.levels)
        if self.grid is not None and self.grid != {}:
            if self.typeoffirstfixedsurface in (118, 119):
                if self.typeoffirstfixedsurface == 119:
                    name = "Hybrid-pressure"
                elif self.typeoffirstfixedsurface == 118:
                    name = "Hybrid-height"
                out.write(name + " coord. A coefficients\n")
                Ai = [level[1]['Ai'] for level in self.grid['gridlevels']]
                for c in range(len(Ai)):
                    write_formatted(out, str(c + 1), Ai[c],
                                    align='>')
                out.write(separation_line)
                out.write(name + " coord. B coefficients\n")
                Bi = [level[1]['Bi'] for level in self.grid['gridlevels']]
                for c in range(len(Bi)):
                    write_formatted(out, str(c + 1), Bi[c],
                                    align='>')
            else:
                if len(self.grid['gridlevels']) == 1:
                    write_formatted(out, "Grid", self.grid['gridlevels'][0])
                else:
                    write_formatted(out, "Grid", self.grid['gridlevels'])
        out.write(separation_line)
        out.write("\n")

    def to_vgrid(self, **kwargs):
        """
        Convert the VGeometry object to a :mod:`vgrid` one.

        Other kwargs passed to HybridPressureVGrid() constructor.
        """
        if self.typeoffirstfixedsurface == 119:
            from pyvgrid import make_HybridPressureVGrid
            vg = make_HybridPressureVGrid(self, **kwargs)
        else:
            raise NotImplementedError('to_vgrid(): self.typeoffirstfixedsurface==' + str(self.typeoffirstfixedsurface))
        return vg


########################
# CONVERSION FUNCTIONS #
########################
def hybridP2pressure(hybridP_geometry, Psurf, vertical_mean,
                     gridposition=None):
    """
    Converts a 'hybrid_pressure' VGeometry to a 'pressure' (in hPa) VGeometry.

    :param VGeometry hybridP_geometry: the initial vertical coordinate
    :param float Psurf: the surface pressure in Pa, needed for integration of Ai and Bi.
    :param str vertical_mean: defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
    :param str gridposition: (= 'mass' or 'flux') is the target grid position. By
      default the data position in the origin geometry is taken.

    :rtype: VGeometry
    """
    assert isinstance(hybridP_geometry, VGeometry), "*hybridP_geometry* must be of type VGeometry."
    if hybridP_geometry.typeoffirstfixedsurface != 119:
        raise epygramError("*hybridP_geometry.typeoffirstfixedsurface* must be 119.")

    if hybridP_geometry.grid['ABgrid_position'] != 'flux':
        raise NotImplementedError("A and B must define flux levels.")
    if gridposition is None:
        gridposition = hybridP_geometry.position_on_grid
    # compute pressures
    A = [level[1]['Ai'] for level in hybridP_geometry.grid['gridlevels']][1:]
    B = [level[1]['Bi'] for level in hybridP_geometry.grid['gridlevels']][1:]
    if gridposition == 'mass':
        levels = profiles.hybridP2masspressure(A, B, Psurf, vertical_mean)
    elif gridposition == 'flux':
        levels = profiles.hybridP2fluxpressure(A, B, Psurf)
    else:
        raise epygramError("gridposition != 'mass' or 'flux'.")
    levels = [l / 100. for l in levels.squeeze()]
    kwargs_vcoord = {'typeoffirstfixedsurface': 100,
                     'position_on_grid': hybridP_geometry.position_on_grid,
                     # 'grid': {'gridlevels':levels},
                     'levels': levels
                     }
    return VGeometry(**kwargs_vcoord)


def hybridH2pressure(hybridH_geometry, P, position):
    """
    Converts a hybrid_height coordinate grid into pressure (in hPa).

    :param P: the vertical profile of pressure to use
    :param position: the position of P values on the grid ('mass' or 'flux')

    :rtype: VGeometry
    """
    assert isinstance(hybridH_geometry, VGeometry), "*hybridH_geometry* must be of type VGeometry."
    if hybridH_geometry.typeoffirstfixedsurface != 118:
        raise epygramError("*hybridH_geometry.typeoffirstfixedsurface* must be 118.")

    if position == hybridH_geometry.position_on_grid:
        levels = P
    elif position == 'flux' and hybridH_geometry.position_on_grid == 'mass':
        raise NotImplementedError('this set of positions is not implemented 1.')
    elif position == 'mass' and hybridH_geometry.position_on_grid == 'flux':
        Pnew = numpy.zeros(len(P))
        for k in range(len(P)):
            if k == 0:
                Pnew[k] = P[k]
            else:
                Pnew[k] = 0.5 * (P[k - 1] + P[k])
        levels = Pnew
    else:
        raise NotImplementedError("grid positions can only be 'mass' or 'flux'.")
    levels = levels[numpy.array(hybridH_geometry.levels) - 1]
    levels = [l / 100 for l in levels.squeeze()]
    kwargs_vcoord = {'typeoffirstfixedsurface': 100,
                     'position_on_grid': hybridH_geometry.position_on_grid,
                     # 'grid':{'gridlevels':levels},
                     'levels': levels
                     }
    return VGeometry(**kwargs_vcoord)


def hybridP2altitude(hybridP_geometry, R, T, Psurf, vertical_mean,
                     Pdep=None, Phi_surf=None):
    """
    Converts a hybrid_pressure coordinate grid into altitude of mass levels.

    :param VGeometry hybridP_geometry: the initial vertical coordinate
    :param list,numpy.ndarray R: the profile of specific gas constant (J/kg/K).
    :param list,numpy.ndarray T: the profile of temperature (K).
    :param float Psurf: the surface pressure, needed for integration of Ai and Bi.
    :param str vertical_mean: defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
    :param list,numpy.ndarray Pdep: the optional profile of NH pressure departures.
    :param float Phi_surf: the optional surface geopotential.
      If given, the final coordinate is altitude above sea level,
      else height above ground surface.

    :rtype: VGeometry
    """
    assert isinstance(hybridP_geometry, VGeometry), "*hybridP_geometry* must be of type VGeometry."
    if hybridP_geometry.typeoffirstfixedsurface != 119:
        raise epygramError("*hybridP_geometry.typeoffirstfixedsurface* must be 119.")
    A = [level[1]['Ai'] for level in hybridP_geometry.grid['gridlevels']][1:]
    B = [level[1]['Bi'] for level in hybridP_geometry.grid['gridlevels']][1:]
    if hybridP_geometry.grid['ABgrid_position'] == 'flux':
        # compute alt
        alt = profiles.hybridP2altitude(A, B, R, T, Psurf,
                                        vertical_mean,
                                        Pdep=Pdep,
                                        Phi_surf=Phi_surf,
                                        Ptop=config.default_Ptop)
        # and update info
        if Phi_surf is None:
            coordinate = 103
        else:
            coordinate = 102
        levels = alt
    elif hybridP_geometry.grid['ABgrid_position'] == 'mass':
        raise NotImplementedError("hybrid-pressure grid at mass-levels.")
    levels = levels[numpy.array(hybridP_geometry.levels) - 1]
    kwargs_vcoord = {'typeoffirstfixedsurface': coordinate,
                     'position_on_grid': hybridP_geometry.position_on_grid,
                     # 'grid':{'gridlevels':list(levels)},
                     'levels': list(levels)
                     }
    return VGeometry(**kwargs_vcoord)


def hybridH2altitude(hybridH_geometry, Zsurf,
                     gridposition=None, conv2height=False):
    """
    Converts a hybrid_height coordinate grid into altitude.

    :param Zsurf: the surface height, needed for integration of Ai and Bi.
    :param gridposition: if given ('mass' or 'flux'), the target grid is
      computed accordingly.
      By default the data position in the origin geometry is taken.
    :param conv2height: if True, conversion into height is performed instead of
      altitude.

    :rtype: VGeometry
    """
    assert isinstance(hybridH_geometry, VGeometry), "*hybridH_geometry* must be of type VGeometry."
    if hybridH_geometry.typeoffirstfixedsurface != 118:
        raise epygramError("*hybridH_geometry.coordinate* must be 118.")
    if hybridH_geometry.grid['ABgrid_position'] != 'flux':
        raise NotImplementedError('A and B must define flux levels')
    if gridposition is None:
        gridposition = hybridH_geometry.position_on_grid
    # compute altitudes
    A = [level[1]['Ai'] for level in hybridH_geometry.grid['gridlevels']]
    B = [level[1]['Bi'] for level in hybridH_geometry.grid['gridlevels']]
    # profiles.hybridH2*height wait for Ai and Bi describing the flux level of extra level under the ground
    A = [-A[1]] + A
    B = [2 - B[1]] + B
    if gridposition == 'mass':
        levels = profiles.hybridH2massheight(A, B, Zsurf,
                                             conv2height=conv2height)
    elif gridposition == 'flux':
        levels = profiles.hybridH2fluxheight(A, B, Zsurf,
                                             conv2height=conv2height)
    else:
        raise epygramError("gridposition != 'mass' or 'flux'.")
    levels = levels[numpy.array(hybridH_geometry.levels)]
    kwargs_vcoord = {'typeoffirstfixedsurface': 103 if conv2height else 102,
                     'position_on_grid': hybridH_geometry.position_on_grid,
                     # 'grid':{'gridlevels':list(levels)},
                     'levels': list(levels)
                     }
    return VGeometry(**kwargs_vcoord)


def height2altitude(height_geometry, Zsurf):
    """
    Converts a height coordinate grid into altitude.

    :param Zsurf: the surface height.

    :rtype: VGeometry
    """
    assert isinstance(height_geometry, VGeometry), "*height_geometry* must be of type VGeometry."
    if height_geometry.typeoffirstfixedsurface != 103:
        raise epygramError("*height_geometry.coordinate* must be 103.")
    levels = list(height_geometry.levels)
    for k in range(len(levels)):
        levels[k] = levels[k] + Zsurf
    kwargs_vcoord = {'typeoffirstfixedsurface': 102,
                     'position_on_grid': height_geometry.position_on_grid,
                     'levels': levels
                     }
    return VGeometry(**kwargs_vcoord)


def altitude2height(altitude_geometry, Zsurf):
    """
    Converts an altitude coordinate grid into height.

    :param Zsurf: the surface height.

    :rtype: VGeometry
    """
    assert isinstance(altitude_geometry, VGeometry), "*altitude_geometry* must be of type VGeometry."
    if altitude_geometry.typeoffirstfixedsurface != 102:
        raise epygramError("*altitude_geometry.coordinate* must be 102.")
    levels = list(altitude_geometry.levels)
    for k in range(len(levels)):
        levels[k] = levels[k] - Zsurf
    kwargs_vcoord = {'typeoffirstfixedsurface': 103,
                     'position_on_grid': altitude_geometry.position_on_grid,
                     'levels': levels
                     }
    return VGeometry(**kwargs_vcoord)


def pressure2altitude(pressure_geometry, R, T, vertical_mean,
                      Pdep=0., Phi_surf=0.):
    """
    Converts a pressure coordinate grid (on mass or flux levels) to
    altitude on mass levels).

    :param R: the profile of specific gas constant (J/kg/K).
    :param T: the profile of temperature (K).
    :param Pdep: the optional profile of NH pressure departures.
    :param Phi_surf: the optional surface geopotential.
      If given, the final coordinate is altitude above sea level,
      else height above ground surface.
    :param vertical_mean: defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
    """
    raise NotImplementedError("Must be modified")
    # La fonction n'a pas été revue depuis la modification sur les géométries verticales.
    # Je n'ai pas prévu le cas d'une grille définie sur des niveaux pression avec des niveaux
    # qui peuvent être différents (au moins à cause des mass/flux).
    # Pour intégrer la possibilité, il faudrait peut-être mettre une grille dans tous les cas
    # et que levels ne soit, à chaque fois, qu'une liste d'entiers indiquant les niveaux dans la grille.
    if pressure_geometry.grid['gridposition'] == 'flux':
        # compute alt
        levels = profiles.pressure2altitude(R, T, vertical_mean,
                                            pi_tilde=pressure_geometry.grid['levels'] * 100,
                                            Pdep=Pdep,
                                            Phi_surf=Phi_surf)
    elif pressure_geometry.grid['gridposition'] == 'mass':
        # compute alt
        levels = profiles.pressure2altitude(R, T, vertical_mean,
                                            pi=pressure_geometry.grid['levels'] * 100,
                                            Pdep=Pdep,
                                            Phi_surf=Phi_surf)
    # and update info
    if Phi_surf is None:
        coordinate = 103
    else:
        coordinate = 102
    geom = pressure_geometry.deppcopy()
    grid = pressure_geometry.grid.copy()
    grid.update({'gridposition':pressure_geometry.gridposition,
                 'levels':list(levels)})
    return geom


def hybridP_coord_and_surfpressure_to_3D_pressure_field(
        hybridP_geometry, Psurf, vertical_mean,
        gridposition=None):
    """
    From a hybridP VGeometry and a surface pressure (in Pa) H2D field,
    compute a 3D field containing the pressure (in hPa) at each hybridP level
    for each gridpoint.

    :param VGeometry hybridP_geometry: the hybridP VGeometry
    :param H2DField Psurf: the surface pressure H2DField in Pa, needed for
        integration of Ai and Bi.
    :param str vertical_mean: defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
    :param str gridposition: (= 'mass' or 'flux') is the target grid position. By
      default the data position in the origin geometry is taken.

    :rtype: D3Field
    """
    from epygram.fields import H2DField
    assert isinstance(Psurf, H2DField)
    if Psurf.spectral:
        Psurf.sp2gp()
    pressures = hybridP2pressure(hybridP_geometry, Psurf.data, vertical_mean,
                                 gridposition=gridposition).levels
    vgeom = hybridP_geometry.deepcopy()
    vgeom.levels = list(range(1, len(pressures) + 1))
    geom = Psurf.geometry.deepcopy()
    geom.vcoordinate = vgeom
    d3pressure = fpx.field(fid={'computed':'pressure'},
                           structure='3D',
                           geometry=geom,
                           units='hPa')
    if isinstance(pressures[0], numpy.ma.masked_array):
        pressures = numpy.ma.array(pressures)
    else:
        pressures = numpy.array(pressures)
    d3pressure.setdata(pressures)
    return d3pressure


def hybridP_coord_to_3D_altitude_field(
        hybridP_geometry, Psurf, vertical_mean,
        t3D, q3D,
        ql3D=None, qi3D=None, qr3D=None, qs3D=None, qg3D=None,
        Pdep3D=None,
        Phi_surf=None):
    """
    From a hybridP VGeometry, a surface pressure (in Pa) H2D field,
    and temperature and specific humidity 3D fields,
    compute a 3D field containing the altitude (in m) at each hybridP level
    for each gridpoint.

    Hydrometeors 3D fields can be provided for more accurate R computation.

    :param VGeometry hybridP_geometry: the hybridP VGeometry
    :param H2DField Psurf: the surface pressure H2DField in Pa, needed for
        integration of Ai and Bi.
    :param str vertical_mean: defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
    :param D3Field t3D: temperature D3Field
    :param D3Field q3D: specific humidity D3Field
    :param D3Field ql3D: liquid water content D3Field
    :param D3Field qi3D: ice water content D3Field
    :param D3Field qr3D: rain water content D3Field
    :param D3Field qs3D: snow water content D3Field
    :param D3Field qg3D: graupel water content D3Field
    :param D3Field Pdep3D: non-hydrostatic pressure departure D3Field
    :param H2DField Phi_surf: surface geopotential H2DField (for altitude vs. height)

    :rtype: D3Field
    """
    from epygram.fields import D3Field, H2DField
    from bronx.meteo.conversion import q2R
    assert isinstance(t3D, D3Field)
    assert isinstance(q3D, D3Field)
    if Phi_surf is not None:
        assert isinstance(Phi_surf, H2DField)
    # preparations
    if t3D.spectral:
        t3D.sp2gp()
    if Phi_surf is not None:
        if Phi_surf.spectral:
            Phi_surf.sp2gp()
        Phi_surf = Phi_surf.data
    hydrometeors = {}
    for qk in ('ql', 'qi', 'qr', 'qs', 'qg'):
        if locals()[qk + '3D'] is not None:
            hydrometeors[qk] = locals()[qk + '3D'].data
    if isinstance(Pdep3D, D3Field):
        Pdep3D = Pdep3D.data
    # computations
    p3D = hybridP_coord_and_surfpressure_to_3D_pressure_field(
        hybridP_geometry, Psurf, vertical_mean,
        gridposition='mass')
    p3D.operation('*', 100)
    R3D = q2R(q3D.data, **hydrometeors)
    alt3D = profiles.pressure2altitude(R3D, t3D.data,
                                       vertical_mean='geometric',
                                       pi=p3D.data,
                                       Phi_surf=Phi_surf,
                                       Pdep=Pdep3D)
    # arrange output field
    p3D.setdata(alt3D)
    alt3D = p3D
    if Phi_surf is None:
        alt3D.fid['computed'] = 'height'
    else:
        alt3D.fid['computed'] = 'altitude'
    alt3D.units = 'm'
    return alt3D
