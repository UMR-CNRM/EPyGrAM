#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
The module is used throughout the package to share constant parameters.

The standard default parameters lie below. They can be overwritten by
the user in the User config file ``userconfig.py`` to be found under
:attr:`userlocaldir`.
"""

import os
import sys

from . import __version__


# INSTALL #
###########
#: Directory of epygram package install
installdir = __file__[:-(len(os.path.basename(__file__)) + 1)]
home = os.getenv('HOME')
#: User customization directory
userlocaldir = os.path.join(home, '.epygram')
#: epygram Colormaps
epygram_colormaps = {'aspect':os.path.join(installdir, 'data', 'colormaps', 'aspect.json'),
                     'gaspect':os.path.join(installdir, 'data', 'colormaps', 'gaspect.json'),
                     'radar':os.path.join(installdir, 'data', 'colormaps', 'radar.json'),
                     'rr1h':os.path.join(installdir, 'data', 'colormaps', 'rr1h.json'),
                     'rr6h':os.path.join(installdir, 'data', 'colormaps', 'rr6h.json'),
                     'rr24h':os.path.join(installdir, 'data', 'colormaps', 'rr24h.json'),
                     'ptype':os.path.join(installdir, 'data', 'colormaps', 'ptype.json'),
                     'ptype0':os.path.join(installdir, 'data', 'colormaps', 'ptype0.json'),
                     }

# PARAMETERS #
##############
#: Ptop: pressure @ summit of atmosphere. For vertical coordinates conversions.
default_Ptop = 0.
#: Epsilon
epsilon = sys.float_info.epsilon
#: Rounding decimal in FP operations
rounding_decimal = 12
#: Maximum number of truncations handled (setup spectral transforms)
KNUMMAXRESOL = 10
#: Plots sizes (in inches)
plotsizes = (16., 12.)
#: Default output for apptools
default_graphical_output = 'png'
#: Cache directory for internet requests
internet_cache_dir = None

# FORMATS #
###########
#: List of implemented/activated formats
#: for the formats_factory, the list is ordered by specificity for those of the
#: same family (e.g. FA before LFI, DDHLFA before LFA...)
#: Removing one of these (in userconfig file) may allow an incomplete install
#: of epygram, disabling one format.
#: In loading order to deal with sub-libraries compatibilities (e.g. h5py after
#: netCDF4, cf. https://github.com/Unidata/netcdf4-python/issues/1343)
implemented_formats = ['netCDFSAF', 'netCDFMNH', 'netCDF', 'GRIB', 'GeoPoints',
                       'TIFFMF', 'HDF5SAF', 'FA', 'LFI', 'DDHLFA', 'LFA']

#: FA default compression parameters
FA_default_compression = {'KNGRIB': 123, 'KDMOPL': 0, 'KPUILA': 0, 'KSTRON': 0,
                          'KNBPDG': 16, 'KNBCSP': 16}
#: FA default KNGRIB for spectral fields
FA_default_KNGRIB_spectral = 103
#: Default reference pressure coefficient for converting hybrid A coefficients
#: in FA
FA_default_reference_pressure = 101325.
#: geoid of FA files in pyproj syntax
FA_default_geoid = {'a':6371229., 'b':6371229.}
#: FA field dictionaries
FA_field_dictionaries_csv = {'default':os.path.join(installdir,
                                                    'data',
                                                    'Field_Dict_FA.csv'),
                             'user':os.path.join(userlocaldir,
                                                 'user_Field_Dict_FA.csv')}
#: FA (write) date & time precision: use FANDAR (minute) or FANDAX (second,
#: cy40t1 onwards)
FA_fandax = True
#: To avoid re-computing lons/lats of Gauss Grids from FA each time needed:
#: makes ARPEGE profiles and section extraction times acceptable (< 1min).
FA_buffered_gauss_grid = True
#: Allow MOCAGE fields to have multiple validities in file; in which case the
#: term is decoded from the fid[2:4]
FA_allow_MOCAGE_multivalidities = False
#: Maximum recommended encoding for FA (KNBPDG)
FA_max_encoding = 30
#: Mute messages from FAIPAR in FA4py
FA_mute_FA4py = False
#: FA limits
FA_limits = {
             #'JPXTRO':,  # Maximum truncation
             #'JPXLAT':,  # Maximum number of latitudes
             #'JPXNIV':,  # Maximum number of levels
             #'JPNXFA':,  # Maximum number of open files
             #'JPNXCA':,  # Maximum number of headers
             #'JPXRATIO':,  # Maximum X/Y ratio
             }

#: LFA maximum number of fields
LFA_max_num_fields = 1000
#: LFA maximum length of strings (as in most LFA useage)
LFA_maxstrlen = 200

#: GeoPoints write precision for lon/lat
GeoPoints_lonlat_precision = 4
#: GeoPoints write precision for other floats
GeoPoints_precision = 6
#: GeoPoints write width of columns
GeoPoints_col_width = 12

#: GRIB lowlevel API library to be used, among ('eccodes', 'grib_api', 'gribapi')
GRIB_lowlevel_api = 'eccodes'
#: GRIB default production centre -- write mode
GRIB_default_centre = 85  # Météo-France
#: GRIB samples from epygram (treated as clone from file)
GRIB_epygram_samples_path = installdir + '/data/grib_samples'
#: satellites local GRIB2 encoding
satellites_local_GRIB2 = {'METEOSAT7':192,
                          'METEOSAT8':193,
                          'METEOSAT9':194,
                          'GOES11':195,
                          'GOES12':196,
                          'MTSAT1':197}
#: sensors local GRIB2 encoding
sensors_local_GRIB2 = {'MVIRI':192,
                       'SEVIRI':193,
                       'IMAGER':194}
#: GRIB: errors while setting packing are fatal
GRIB_packing_fatal = True
#: GRIB: ignore errors while trying to decode validity
GRIB_ignore_validity_decoding_errors = False
#: Maximum recommended encoding for GRIB (bitsPerValue)
GRIB_max_bitspervalue = 30
#: Force bitspervalue to GRIB_max_bitspervalue if requested higher
GRIB_force_bitspervalue = False
#: Use temporary links to workaround a bug in grib indexes
#: set either a directory to use for that purpose (e.g. /tmp) or
#: False not to use this option
GRIB_safe_indexes = '/tmp'

#: LFI field dictionaries
LFI_field_dictionaries_csv = {'default':os.path.join(installdir,
                                                     'data',
                                                     'Field_Dict_LFI.csv'),
                              'user':os.path.join(userlocaldir,
                                                  'user_Field_Dict_LFI.csv')}
#: geoid of LFI files in pyproj syntax
LFI_default_geoid = {'a':6371229., 'b':6371229.}

#: netCDFMNH field dictionaries
netCDFMNH_field_dictionaries_csv = {'default':os.path.join(installdir,
                                                           'data',
                                                           'Field_Dict_netCDFMNH.csv'),
                                    'user':os.path.join(userlocaldir,
                                                        'user_Field_Dict_netCDFMNH.csv')}
#: geoid of netCDFMNH files in pyproj syntax
netCDFMNH_default_geoid = {'a':6371229., 'b':6371229.}

#: netCDF standard dimensions
netCDF_standard_dimensions = ['N_dimension',  # numerotation (obs, profile, ...)
                              'T_dimension',  # time
                              'X_dimension',  # X-axis (cartesian projection or longitude)
                              'Y_dimension',  # Y-axis (cartesian projection or latitude)
                              'Z_dimension']  # Z-axis (vertical)
#: netCDF usual names for standard dimensions
netCDF_usualnames_for_standard_dimensions = {'N_dimension':('N', 'n', 'transect', 'obs', 'profile',
                                                            'Number_of_points', 'gridpoints_number'),
                                             'T_dimension':('time', 'T', 't', 'validity', 'time_counter'),
                                             'X_dimension':('X', 'x', 'xx', 'nx',
                                                            'LON', 'lon', 'Nbre_lon', 'longitude',
                                                            'max_lon_number'),
                                             'Y_dimension':('Y', 'y', 'yy', 'ny',
                                                            'LAT', 'lat', 'Nbre_lat', 'latitude',
                                                            'lat_number'),
                                             'Z_dimension':('Z', 'z', 'zz', 'nz',
                                                            'level', 'nlev', 'nlevp1',
                                                            'pressure', 'hybrid-pressure',
                                                            'height', 'altitude', 'depth')}
#: netCDF usual names for storing lon, lat grids
netCDF_usualnames_for_lonlat_grids = {'X':['LON', 'lon', 'lons', 'longitude', 'longitudes'],
                                      'Y':['LAT', 'lat', 'lats', 'latitude', 'latitudes']}
#: netCDF default behaviour
#: the behaviour can be updated with specific netCDF_usualnames_for_standard_dimensions,
#: e.g. 'X_dimension':'nb_lon'
#: and eventual according grid, e.g. 'X_grid':'longitudes'
netCDF_default_behaviour = {'reverse_Yaxis':False,
                            # writing behaviours:
                            'flatten_horizontal_grids':False,
                            'write_lonlat_grid':True,
                            'H1D_is_H2D_unstructured':False,
                            }
#: netCDF default compression
netCDF_default_compression = 4
#: netCDF, replace dots in variable names by...
netCDF_replace_dot_in_variable_names = '.'
#: netCDF default standard global attributes
netCDF_default_global_attributes = {'made_with':'epygram-' + __version__}
#: netCDF variables data type
netCDF_default_variables_dtype = 'f8'
#: netCDF metavariables data type
netCDF_default_metavariables_dtype = 'f8'
#: netCDF variables fill value.
#: None will make netCDF ignore existence of a _FillValue
netCDF_default_variables_fill_value = -999999.9


# OPTIONS #
###########
#: A classical default geoid for *pyproj*
default_geoid = {'ellps':'WGS84'}
#: Protect unhappy writes: ask before opening a new file on an existing path
protect_unhappy_writes = False
#: Threshold on field absolute values to mask meaningless field values
mask_outside = 1e19
#: pyproj coordinates invalid values
pyproj_invalid_values = 1.e30
#: To hide footprints warnings...
hide_footprints_warnings = True
#: To raise an error if the memory needed for Legendre transforms exceeds
#: this percentage of the available memory.
prevent_swapping_legendre = 0.75
#: Use footprints.proxy builder to generate a field.
#: True: more flexible, False: faster
footprints_proxy_as_builder = False  # CLEANME: remove this, less useful since fasttrack ?
#: Vector graphical symbol
vector_symbol = 'barbs'
#: Default quality for figures
default_figures_dpi = 150
#: ordering of spectral coefficients, with regards to arpifs spaces:
#: 'model' or 'FA'.
#: => 'model': read/write by wfacilo/wfaieno,
#:    'FA': read/write by wfacile/wfaienc
spectral_coeff_order = 'model'
#: To call epygram.init_env() automatically at import
#: ! Should not be True if using Vortex !
init_at_import = False
#: hide messages when guessing format
#: True is dangerous, causes troubles in logging
silent_guess_format = False
#: Number or margin within C-zone to generate a lonlat-included domain
margin_points_within_Czone = 3
#: Defaults for matplotlib rcparams
default_rcparams = [(('font',), dict(family='serif')), ]
#: Plugins to be activated by default
activate_plugins = ['with_vtk', 'with_cartopy']


# USER CUSTOMIZATION #
######################
# usercolormaps should also remain empty here
#: In userconfig, this should be a dict whose keys are the colormap name and
#: values the source absolute path of the colormap definition;
#: e.g. {'aspect', '/home/mary/.epygram/aspect.json'}.
usercolormaps = {}


# OVERWRITE WITH USER CONFIG #
##############################
if os.path.exists(os.path.join(userlocaldir, 'userconfig.py')):
    sys.path.insert(0, userlocaldir)
    from userconfig import *
    sys.path.remove(userlocaldir)
#: colormaps gathers epygram and user colormaps
colormaps = {}
colormaps.update(epygram_colormaps)
colormaps.update(usercolormaps)
