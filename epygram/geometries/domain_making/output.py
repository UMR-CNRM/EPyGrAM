#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains functions for outing a LAM domain to namelists, plot, summary.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import six

from bronx.datagrip import namelist

from epygram.config import epsilon
from epygram.geometries.SpectralGeometry import truncation_from_gridpoint_dims

from .util import Ezone_minimum_width
from .build import compute_lonlat_included, build_CIE_field, build_lonlat_field


def lam_geom2namelists(geometry,
                       truncation='linear',
                       orography_subtruncation='quadratic',
                       Ezone_in_pgd=False,
                       Iwidth_in_pgd=False):
    """
    From the geometry, build the namelist blocks for the necessary namelists.

    :param truncation: the kind of truncation of spectral geometry
                       to generate, among ('linear', 'quadratic', 'cubic').
    :param orography_subtruncation: additional subtruncation for orography to
                                    be generated.
    """
    namelists = {}

    # compute additionnal parameters
    truncation = truncation_from_gridpoint_dims(geometry.dimensions,
                                                grid=truncation)
    subtruncation = truncation_from_gridpoint_dims(geometry.dimensions,
                                                   grid=orography_subtruncation)

    # PGD namelist
    nam = namelist.NamelistSet()
    namelists['namel_buildpgd'] = nam
    nam.add(namelist.NamelistBlock('NAM_CONF_PROJ'))
    nam.add(namelist.NamelistBlock('NAM_CONF_PROJ_GRID'))
    nam['NAM_CONF_PROJ']['XLON0'] = geometry.projection['reference_lon'].get('degrees')
    nam['NAM_CONF_PROJ']['XLAT0'] = geometry.projection['reference_lat'].get('degrees')
    nam['NAM_CONF_PROJ']['XRPK'] = geometry.projection['reference_lat'].get('cos_sin')[1]
    nam['NAM_CONF_PROJ']['XBETA'] = 0.  # geometry.projection['reference_lon'].get('degrees') - geometry.getcenter()[0].get('degrees')
    nam['NAM_CONF_PROJ_GRID']['XLONCEN'] = geometry.getcenter()[0].get('degrees')
    nam['NAM_CONF_PROJ_GRID']['XLATCEN'] = geometry.getcenter()[1].get('degrees')
    nam['NAM_CONF_PROJ_GRID']['NIMAX'] = geometry.dimensions['X_CIzone']
    nam['NAM_CONF_PROJ_GRID']['NJMAX'] = geometry.dimensions['Y_CIzone']
    nam['NAM_CONF_PROJ_GRID']['XDX'] = geometry.grid['X_resolution']
    nam['NAM_CONF_PROJ_GRID']['XDY'] = geometry.grid['Y_resolution']
    if Iwidth_in_pgd:
        nam['NAM_CONF_PROJ_GRID']['IWIDTH_I_X'] = geometry.dimensions['X_Iwidth']
        nam['NAM_CONF_PROJ_GRID']['IWIDTH_I_Y'] = geometry.dimensions['Y_Iwidth']
    if Ezone_in_pgd:
        nam['NAM_CONF_PROJ_GRID']['ILONE'] = geometry.dimensions['X'] - geometry.dimensions['X_CIzone']
        nam['NAM_CONF_PROJ_GRID']['ILATE'] = geometry.dimensions['Y'] - geometry.dimensions['Y_CIzone']

    # c923 namelist
    nam = namelist.NamelistSet()
    namelists['namel_c923'] = nam
    nam.add(namelist.NamelistBlock('NAMCT0'))
    nam.add(namelist.NamelistBlock('NAMDIM'))
    nam.add(namelist.NamelistBlock('NEMDIM'))
    nam.add(namelist.NamelistBlock('NEMGEO'))
    nam['NAMCT0']['LRPLANE'] = True
    nam['NAMDIM']['NDLON'] = geometry.dimensions['X']
    nam['NAMDIM']['NDLUXG'] = geometry.dimensions['X_CIzone']
    nam['NAMDIM']['NDGLG'] = geometry.dimensions['Y']
    nam['NAMDIM']['NDGUXG'] = geometry.dimensions['Y_CIzone']
    nam['NAMDIM']['NMSMAX'] = truncation['in_X']
    nam['NAMDIM']['NSMAX'] = truncation['in_Y']
    nam['NEMDIM']['NBZONL'] = geometry.dimensions['X_Iwidth']
    nam['NEMDIM']['NBZONG'] = geometry.dimensions['Y_Iwidth']
    nam['NEMGEO']['ELON0'] = geometry.projection['reference_lon'].get('degrees')
    nam['NEMGEO']['ELAT0'] = geometry.projection['reference_lat'].get('degrees')
    nam['NEMGEO']['ELONC'] = geometry.getcenter()[0].get('degrees')
    nam['NEMGEO']['ELATC'] = geometry.getcenter()[1].get('degrees')
    nam['NEMGEO']['EDELX'] = geometry.grid['X_resolution']
    nam['NEMGEO']['EDELY'] = geometry.grid['Y_resolution']

    # subtruncated grid namelist
    nam = namelist.NamelistSet()
    namelists['namel_c923_orography'] = nam
    nam.add(namelist.NamelistBlock('NAMDIM'))
    nam['NAMDIM']['NMSMAX'] = subtruncation['in_X']
    nam['NAMDIM']['NSMAX'] = subtruncation['in_Y']

    # couplingsurf namelist
    nam = namelist.NamelistSet()
    namelists['namel_e927'] = nam
    nam.add(namelist.NamelistBlock('NAMFPD'))
    nam.add(namelist.NamelistBlock('NAMFPG'))
    nam['NAMFPD']['NLON'] = geometry.dimensions['X']
    nam['NAMFPD']['NFPLUX'] = geometry.dimensions['X_CIzone']
    nam['NAMFPD']['NFPBZONL'] = geometry.dimensions['X_Iwidth']
    nam['NAMFPD']['NLAT'] = geometry.dimensions['Y']
    nam['NAMFPD']['NFPGUX'] = geometry.dimensions['Y_CIzone']
    nam['NAMFPD']['NFPBZONG'] = geometry.dimensions['Y_Iwidth']
    nam['NAMFPD']['RLONC'] = geometry.getcenter()[0].get('degrees')
    nam['NAMFPD']['RLATC'] = geometry.getcenter()[1].get('degrees')
    nam['NAMFPD']['RDELX'] = geometry.grid['X_resolution']
    nam['NAMFPD']['RDELY'] = geometry.grid['Y_resolution']
    nam['NAMFPG']['FPLON0'] = geometry.projection['reference_lon'].get('degrees')
    nam['NAMFPG']['FPLAT0'] = geometry.projection['reference_lat'].get('degrees')
    nam['NAMFPG']['NMFPMAX'] = truncation['in_X']
    nam['NAMFPG']['NFPMAX'] = truncation['in_Y']

    return namelists


def regll_geom2namelists(geometry, domain_name='__YOUR_DOM_NAME__'):
    """
    From the regular LonLat geometry, build the namelist blocks for the
    necessary namelists.
    """
    namelists = {}

    # PGD
    nam = namelist.NamelistSet()
    namelists['namel_buildpgd'] = nam
    nam.add(namelist.NamelistBlock('NAM_PGD_GRID'))
    nam.add(namelist.NamelistBlock('NAM_LONLAT_REG'))
    corners = geometry.gimme_corners_ll()
    nam['NAM_PGD_GRID']['CGRID'] = 'LONLAT REG'
    nam['NAM_LONLAT_REG']['XLONMIN'] = corners['ul'][0] - geometry.grid['X_resolution'].get('degrees') / 2.
    nam['NAM_LONLAT_REG']['XLONMAX'] = corners['lr'][0] + geometry.grid['X_resolution'].get('degrees') / 2.
    nam['NAM_LONLAT_REG']['XLATMIN'] = corners['lr'][1] - geometry.grid['Y_resolution'].get('degrees') / 2.
    nam['NAM_LONLAT_REG']['XLATMAX'] = corners['ul'][1] + geometry.grid['Y_resolution'].get('degrees') / 2.
    nam['NAM_LONLAT_REG']['NLON'] = geometry.dimensions['X']
    nam['NAM_LONLAT_REG']['NLAT'] = geometry.dimensions['Y']
    assert all([nam['NAM_LONLAT_REG']['XLATMIN'] >= -90.,
                nam['NAM_LONLAT_REG']['XLATMAX'] <= 90.]), \
        'latitude borders too close to poles (less than half a resolution)'

    # c923
    nam = namelist.NamelistSet()
    namelists['namel_c923'] = nam
    nam.add(namelist.NamelistBlock('NAMCT0'))
    nam.add(namelist.NamelistBlock('NAMDIM'))
    nam.add(namelist.NamelistBlock('NEMGEO'))
    nam['NAMCT0']['LRPLANE'] = False
    nam['NAMDIM']['NDGUXG'] = geometry.dimensions['Y']
    nam['NAMDIM']['NDGLG'] = geometry.dimensions['Y']
    nam['NAMDIM']['NDLUXG'] = geometry.dimensions['X']
    nam['NAMDIM']['NDLON'] = geometry.dimensions['X']
    nam['NEMGEO']['ELAT0'] = 0.
    nam['NEMGEO']['ELON0'] = 0.
    nam['NEMGEO']['ELONC'] = geometry.getcenter()[0].get('degrees')
    nam['NEMGEO']['ELATC'] = geometry.getcenter()[1].get('degrees')
    nam['NEMGEO']['EDELX'] = geometry.grid['X_resolution'].get('degrees')
    nam['NEMGEO']['EDELY'] = geometry.grid['Y_resolution'].get('degrees')

    # FullPos
    nam = namelist.NamelistSet()
    namelists['namel_fpos'] = nam
    nam.add(namelist.NamelistBlock('NAMFPC'))
    nam.add(namelist.NamelistBlock('NAMFPD'))
    nam.add(namelist.NamelistBlock('NAMFPF'))
    nam['NAMFPC']['CFPFMT'] = 'LALON'
    nam['NAMFPC']['CFPDOM'] = domain_name.upper()
    nam['NAMFPD']['NLON(1)'] = geometry.dimensions['X']
    nam['NAMFPD']['NLAT(1)'] = geometry.dimensions['Y']
    nam['NAMFPD']['RLONC(1)'] = geometry.getcenter()[0].get('degrees')
    nam['NAMFPD']['RLATC(1)'] = geometry.getcenter()[1].get('degrees')
    nam['NAMFPD']['RDELX(1)'] = geometry.grid['X_resolution'].get('degrees')
    nam['NAMFPD']['RDELY(1)'] = geometry.grid['Y_resolution'].get('degrees')
    # nam['NAMFPF']['NFPMAX(1)'] = int(max(geometry.dimensions['X'],
    #                                      geometry.dimensions['Y'])) / 2
    return namelists


def write_namelists(namelists, out=None, prefix='', suffix='geoblocks'):
    """
    Write out namelists blocks.

    :param namelists: dict of NamelistSet (from bronx.datagrip.namelist)
    :param out: if given as a str: write all in one file,
                else in separate files:
                     either as filename if out==dict(namelist:filename, ...)
                     or if None following syntax: "prefix.namelist.suffix".
    :param prefix: prefix for output names
    :param suffix: prefix for output names
    """

    if isinstance(out, six.string_types):
        out.write("# Namelists blocks #\n")
        out.write("  ================\n")
        for n in sorted(namelists.keys(), reverse=True):
            out.write(namelists[n].dumps())
    else:
        for n, nam in namelists.items():
            if isinstance(out, dict):
                nm = out[n]
            else:
                nm = [n, suffix]
                if prefix != '':
                    nm.insert(0, prefix)
                nm = '.'.join(nm)
            with open(nm, 'w') as out:
                out.write(nam.dumps())


def write_geometry_as_namelists(geometry, allinone=False):
    """
    Write out namelists blocks from a geometry.

    :param allinone: if True, write all in one file, else in separate files.
    """
    namelists_blocks = lam_geom2namelists(geometry)

    if allinone:
        outputfilename = "new_domain.namelists_blocks"
        with open(outputfilename, 'w') as out:
            out.write(summary(geometry) + '\n')
            write_namelists(namelists_blocks, out)
    else:
        outputfilename = "new_domain.summary"
        with open(outputfilename, 'w') as out:
            out.write(summary(geometry) + '\n')
        write_namelists(namelists_blocks, out=None)


def summary(geometry):
    """
    Returns a summary of geometry as a character string.

    :param geometry: a H2DGeometry instance
    """
    invprojections = {'lambert':'L', 'mercator':'M', 'polar_stereographic':'PS'}
    print_projections = {'L':'Lambert (conformal conic)', 'M':'Mercator', 'PS':'Polar Stereographic'}

    disp = ""
    disp += "# Geometry Summary #" + '\n'
    disp += "  ================  " + '\n'
    disp += "Center Longitude: " + str(geometry.getcenter()[0].get('degrees')) + '\n'
    disp += "Center Latitude:  " + str(geometry.getcenter()[1].get('degrees')) + '\n'
    disp += "Tilting:          " + \
        str(geometry.projection['reference_lon'].get('degrees') - geometry.getcenter()[0].get('degrees')) + '\n'
    disp += "  => Reference longitude: " + str(geometry.projection['reference_lon'].get('degrees')) + '\n'
    disp += "Projection: " + print_projections[invprojections[geometry.name]] + '\n'
    disp += "  Reference latitude: " + str(geometry.projection['reference_lat'].get('degrees')) + '\n'
    disp += "Resolution: " + str(geometry.grid['X_resolution']) + '\n'
    mapfactor = geometry.map_factor_field().getdata(subzone='CI')
    mapfactor_range = [mapfactor.min(), mapfactor.max()]
    disp += "Map factor range on C+I domain: [" + \
          '{:.{precision}{type}}'.format(mapfactor_range[0], type='E', precision=3) + " -> " + \
          '{:.{precision}{type}}'.format(mapfactor_range[1], type='E', precision=3) + "]" + '\n'
    disp += "---" + '\n'
    disp += "Dimensions    " + '{:^{width}}'.format("C+I", width=10) + \
          '{:^{width}}'.format("C+I+E", width=10) + \
          '{:^{width}}'.format("E-zone", width=10) + '\n'
    Ezone_Xwidth = geometry.dimensions['X'] - geometry.dimensions['X_CIzone']
    Ezone_Ywidth = geometry.dimensions['Y'] - geometry.dimensions['Y_CIzone']
    if Ezone_Xwidth == Ezone_minimum_width:
        disp += "X:            " + '{:^{width}}'.format(str(geometry.dimensions['X_CIzone']), width=10) + \
                '{:^{width}}'.format(str(geometry.dimensions['X']), width=10) + \
                '{:^{width}}'.format(str(Ezone_Xwidth), width=10) + \
                ' (optimal)' + '\n'
    else:
        disp += "X:            " + '{:^{width}}'.format(str(geometry.dimensions['X_CIzone']), width=10) + \
                '{:^{width}}'.format(str(geometry.dimensions['X']), width=10) + \
                '{:^{width}}'.format(str(Ezone_Xwidth), width=10) + \
                '/ ' + str(Ezone_minimum_width) + ' (optimal) => advised C+I zonal (X) dimension =' + \
                '{:^{width}}'.format(str(geometry.dimensions['X'] - Ezone_minimum_width), width=10) + '\n'
    if Ezone_Ywidth == Ezone_minimum_width:
        disp += "Y:            " + '{:^{width}}'.format(str(geometry.dimensions['Y_CIzone']), width=10) + \
                '{:^{width}}'.format(str(geometry.dimensions['Y']), width=10) + \
                '{:^{width}}'.format(str(Ezone_Ywidth), width=10) + \
                ' (optimal)' + '\n'
    else:
        disp += "Y:            " + '{:^{width}}'.format(str(geometry.dimensions['Y_CIzone']), width=10) + \
                '{:^{width}}'.format(str(geometry.dimensions['Y']), width=10) + \
                '{:^{width}}'.format(str(Ezone_Ywidth), width=10) + \
                '/ ' + str(Ezone_minimum_width) + ' (optimal) => advised C+I meridian (Y) dimension =' + \
                '{:^{width}}'.format(str(geometry.dimensions['Y'] - Ezone_minimum_width), width=10) + '\n'
    disp += "I zone width: " + str(geometry.dimensions['X_Iwidth']) + '\n'
    if (geometry.projection['reference_lon'].get('degrees') -
        geometry.getcenter()[0].get('degrees')) <= epsilon:
        disp += "---" + '\n'
        disp += "The domain contains (at least) the following lon/lat regular area:" + '\n'
        ll_included = compute_lonlat_included(geometry)
        disp += "Longitudes: " + '{:.{precision}{type}}'.format(ll_included['lonmin'], type='F', precision=4) + \
                " <--> " + '{:.{precision}{type}}'.format(ll_included['lonmax'], type='F', precision=4) + '\n'
        disp += "Latitudes:  " + '{:.{precision}{type}}'.format(ll_included['latmax'], type='F', precision=4) + '\n'
        disp += "            " + "  ^" + '\n'
        disp += "            " + "  |" + '\n'
        disp += "            " + "  v" + '\n'
        disp += "            " + '{:.{precision}{type}}'.format(ll_included['latmin'], type='F', precision=4) + '\n'
    disp += "--------------------------------------------------\n"

    return disp


def plot_geometry(geometry,
                  lonlat_included=None,
                  out=None,
                  background=True,
                  departments=False,
                  **_):
    """
    Plot the built geometry, along with lonlat included domain if given.

    :param lonlat_included: parameters of the lonlat domain to plot
    :param out: filename (.png) if not None (else interactive pyplot.show())
    :param background: if True, set a background color to continents and oceans.
    :param departments: if True, adds the french departments on map (instead
                         of countries).
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # plot
    CIEdomain = build_CIE_field(geometry)
    sat_height = geometry.distance(geometry.gimme_corners_ll()['ll'],
                                   geometry.gimme_corners_ll()['ur'])
    center = geometry.getcenter()
    projection = ccrs.NearsidePerspective(central_longitude=center[0].get('degrees'),
                                          central_latitude=center[1].get('degrees'),
                                          satellite_height=sat_height)
    if background:
        cartopy_features = [cfeature.LAND, cfeature.OCEAN,
                            cfeature.LAKES, cfeature.RIVERS]
    else:
        cartopy_features = []
    fig, ax = CIEdomain.cartoplot(projection=projection,
                                  colorsnumber=6,
                                  minmax=[-1.0, 3.0],
                                  colorbar=False,
                                  title='Domain: C+I+E',
                                  meridians=None,
                                  parallels=None,
                                  cartopy_features=cartopy_features,
                                  epygram_departments=departments,
                                  set_global=True)
    if lonlat_included is not None:
        ll_domain = build_lonlat_field(lonlat_included)
        ll_domain.cartoplot(fig=fig,
                            ax=ax,
                            projection=projection,
                            plot_method='contour',
                            title='Domain: C+I+E \n Blue contour: required lon/lat',
                            colorsnumber=2,
                            contourcolor='blue',
                            contour_kw=dict(contourwidth=4,),
                            contourlabel=False,
                            set_global=True,
                            meridians=None,
                            parallels=None)
    if out is not None:
        fig.savefig(out, bbox_inches='tight')
    else:
        import matplotlib.pyplot as plt
        plt.show()
