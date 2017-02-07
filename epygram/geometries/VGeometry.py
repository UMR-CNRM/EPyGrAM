#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the class for Vertical geometry of fields.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import numpy
import sys

import footprints
from footprints import FootprintBase, FPDict, FPList, proxy as fpx

from epygram import profiles, epygramError, config
from epygram.util import RecursiveObject, write_formatted, separation_line


class VGeometry(RecursiveObject, FootprintBase):
    """
    Handles the vertical geometry for fields.

    Here, the grid defines the vertical position of each level
    between a bottom and a top positions.

    The position of points w/r to the vertical grid (mass or flux points),
    is interpreted as:\n
    - mass: points are located on same levels as the grid points.
    - flux: points are located on half-levels, hence are N+1.

    levels is a list with one item for each level represented in data.
    Each item can be a scalar (constant value for all the data point),
                     an array with the horizontal geographic shape (level constant in time but varying on the horizontal),
                     an array with the first dimension corresponding to the validity lengthy and other dimensions to represent the horizotal.
    It is not allowed to have a level varying in time and constant on the geographic domain.
    """

    _collector = ('geometry',)
    _footprint = dict(
        attr=dict(
            structure=dict(
                values=set(['V'])),
            typeoffirstfixedsurface=dict(
                info="Type of horizontal level, as of GRIB2 norm \
                      (inspired from GRIB_API).",
                type=int),
            levels=dict(
                type=FPList,
                info="Effective levels on which data is available."),
            grid=dict(
                type=FPDict,
                optional=True,
                default=FPDict({}),
                info="Handles description of the vertical grid."),
            position_on_grid=dict(
                optional=True,
                info="Position of points w/r to the vertical grid.",
                values=set(['mass', 'flux', '__unknown__']),
                default='__unknown__')
        )
    )

    def what(self, out=sys.stdout,
             levels=True):
        """
        Writes in file a summary of the geometry.

        Args: \n
        - *out*: the output open file-like object (duck-typing: *out*.write()
          only is needed).
        - *levels*: if True, writes the levels of the geometry
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


########################
# CONVERSION FUNCTIONS #
########################
def hybridP2pressure(hybridP_geometry, Psurf, vertical_mean,
                     gridposition=None):
    """
    Converts a 'hybrid_pressure' VGeometry to a 'pressure' VGeometry.

    *Psurf* is the surface pressure in Pa, needed for integration of Ai and Bi.

    *gridposition* (= 'mass' or 'flux') is the target grid position. By
    default the data position in the origin geometry is taken.

    *vertical_mean* defines the kind of averaging done on the vertical
    to compute half-levels from full-levels, or inverse: 'geometric' or
    'arithmetic'.
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
    kwargs_vcoord = {'structure':'V',
                     'typeoffirstfixedsurface': 100,
                     'position_on_grid': hybridP_geometry.position_on_grid,
                     'grid': {'gridlevels':levels},
                     'levels': levels
                     }

    return fpx.geometry(**kwargs_vcoord)


def hybridH2pressure(hybridH_geometry, P, position):
    """
    Converts a hybrid_height coordinate grid into pressure.
    - *P* is the vertical profile of pressure to use
    - *position* is the position of P values on the grid ('mass' or 'flux')
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
    kwargs_vcoord = {'structure':'V',
                     'typeoffirstfixedsurface': 100,
                     'position_on_grid': hybridH_geometry.position_on_grid,
                     'grid':{'gridlevels':levels},
                     'levels': levels
                     }

    return fpx.geometry(**kwargs_vcoord)


def hybridP2altitude(hybridP_geometry, R, T, Psurf, vertical_mean,
                     Pdep=0., Phi_surf=None):
    """
    Converts a hybrid_pressure coordinate grid into altitude of mass levels.

    - *R* is the profile of specific gas constant (J/kg/K).
    - *T* is the profile of temperature (K).
    - *Psurf* is the surface pressure, needed for integration of Ai and Bi.
    - *Pdep* is the optional profile of NH pressure departures.
    - *Phi_surf* is the optional surface geopotential.
      If given, the final coordinate is altitude above sea level,
      else height above ground surface.
    - *vertical_mean* defines the kind of averaging done on the vertical
      to compute half-levels from full-levels, or inverse: 'geometric' or
      'arithmetic'.
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
        if Phi_surf is None or \
           numpy.all(numpy.array(abs(Phi_surf)) < config.epsilon):
            coordinate = 103
        else:
            coordinate = 102
        levels = alt
    elif hybridP_geometry.grid['ABgrid_position'] == 'mass':
        raise NotImplementedError("hybrid-pressure grid at mass-levels.")
    levels = levels[numpy.array(hybridP_geometry.levels) - 1]
    kwargs_vcoord = {'structure':'V',
                     'typeoffirstfixedsurface': coordinate,
                     'position_on_grid': hybridP_geometry.position_on_grid,
                     'grid':{'gridlevels':list(levels)},
                     'levels': list(levels)
                     }

    return fpx.geometry(**kwargs_vcoord)


def hybridH2altitude(hybridH_geometry, Zsurf,
                     gridposition=None, conv2height=False):
    """
    Converts a hybrid_height coordinate grid into altitude.

    *Zsurf* is the surface pressure, needed for integration of Ai and Bi.

    If *gridposition* is given ('mass' or 'flux'), the target grid is
    computed accordingly.
    By default the data position in the origin geometry is taken.
    If conv2height is True, conversion into height is performed instead of
    altitude.
    """

    assert isinstance(hybridH_geometry, VGeometry), "*hybridH_geometry* must be of type VGeometry."
    if hybridH_geometry.typeoffirstfixedsurface != 118:
        raise epygramError("*hybridH_geometry.cooridnate* must be 118.")

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

    levels = levels[numpy.array(hybridH_geometry.levels) - 1]

    kwargs_vcoord = {'structure':'V',
                     'typeoffirstfixedsurface': 103 if conv2height else 102,
                     'position_on_grid': hybridH_geometry.position_on_grid,
                     'grid':{'gridlevels':list(levels)},
                     'levels': list(levels)
                     }

    return fpx.geometry(**kwargs_vcoord)


def pressure2altitude(pressure_geometry, R, T, vertical_mean,
                      Pdep=0., Phi_surf=0.):
    """
    Converts a pressure coordinate grid (on mass or flux levels) to
    altitude on mass levels).

    - *R* is the profile of specific gas constant (J/kg/K).
    - *T* is the profile of temperature (K).
    - *Pdep* is the optional profile of NH pressure departures.
    - *Phi_surf* is the optional surface geopotential.
      If given, the final coordinate is altitude above sea level,
      else height above ground surface.
    - *vertical_mean* defines the kind of averaging done on the vertical
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
    if abs(Phi_surf) < config.epsilon:
        coordinate = 103
    else:
        coordinate = 102
    grid = pressure_geometry.grid.copy()
    grid.update({'gridposition':pressure_geometry.gridposition,
                 'levels':list(levels)})

    return fpx.geometry(structure='V1D',
                        coordinate=coordinate,
                        grid=grid,
                        hlocation=pressure_geometry.hlocation,
                        position_on_grid=pressure_geometry.position_on_grid)


footprints.collectors.get(tag='geometrys').fasttrack = ('format',)
