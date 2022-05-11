:mod:`epygram.config` --- Configuration: parameters, options
============================================================

.. automodule:: epygram.config
   :synopsis: Configuration: parameters, options used throughout the package.

-----------------------------------------------------------

Installation
-------------

.. autodata:: installdir
.. autodata:: userlocaldir
.. autodata:: epygram_colormaps

-----------------------------------------------------------

Parameters
----------

.. autodata:: default_Ptop
.. autodata:: epsilon
.. autodata:: KNUMMAXRESOL
.. autodata:: plotsizes
.. autodata:: default_graphical_output

-----------------------------------------------------------

Formats Specifications
----------------------

.. autodata:: implemented_formats

FA
..

.. autodata:: FA_default_compression
.. autodata:: FA_default_reference_pressure
.. autodata:: FA_default_geoid
.. autodata:: FA_field_dictionaries_csv
.. autodata:: FA_fandax
.. autodata:: FA_buffered_gauss_grid
.. autodata:: FA_allow_MOCAGE_multivalidities
.. autodata:: FA_max_encoding
.. autodata:: FA_mute_FA4py

LFA
...

.. autodata:: LFA_max_num_fields
.. autodata:: LFA_maxstrlen

GeoPoints
.........

.. autodata:: GeoPoints_lonlat_precision
.. autodata:: GeoPoints_precision
.. autodata:: GeoPoints_col_width

GRIB
....

.. autodata:: GRIB_lowlevel_api
.. autodata:: GRIB_default_centre
.. autodata:: GRIB_epygram_samples_path
.. autodata:: satellites_local_GRIB2
.. autodata:: sensors_local_GRIB2
.. autodata:: GRIB_packing_fatal
.. autodata:: GRIB_ignore_validity_decoding_errors
.. autodata:: GRIB_max_bitspervalue
.. autodata:: GRIB_force_bitspervalue
.. autodata:: GRIB_safe_indexes

LFI
...

.. autodata:: LFI_field_dictionaries_csv
.. autodata:: LFI_default_geoid

netCDF
......

.. autodata:: netCDF_standard_dimensions
.. autodata:: netCDF_usualnames_for_standard_dimensions
.. autodata:: netCDF_usualnames_for_lonlat_grids
.. autodata:: netCDF_default_behaviour
.. autodata:: netCDF_default_compression
.. autodata:: netCDF_replace_dot_in_variable_names
.. autodata:: netCDF_default_global_attributes
.. autodata:: netCDF_default_variables_dtype
.. autodata:: netCDF_default_metavariables_dtype
.. autodata:: netCDF_default_variables_fill_value

Options
-------

.. autodata:: default_geoid
.. autodata:: protect_unhappy_writes
.. autodata:: mask_outside
.. autodata:: pyproj_invalid_values
.. autodata:: hide_footprints_warnings
.. autodata:: prevent_swapping_legendre
.. autodata:: footprints_proxy_as_builder
.. autodata:: vector_symbol
.. autodata:: default_figures_dpi
.. autodata:: spectral_coeff_order
.. autodata:: init_at_import
.. autodata:: silent_guess_format
.. autodata:: margin_points_within_Czone
.. autodata:: default_rcparams

-----------------------------------------------------------

User modules
------------

.. autodata:: usermodules

-----------------------------------------------------------

Colormaps
---------

See :ref:`Colormaps <epy-colormaps>` for illustration of the EPyGrAM colormaps.

.. autodata:: colormaps
.. autodata:: usercolormaps
