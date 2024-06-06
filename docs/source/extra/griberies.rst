:mod:`griberies` --- utilities around GRIB
==========================================

.. automodule:: epygram.extra.griberies
   :synopsis: utilities around GRIB.

-----------------------------------------------------------

.. autoclass:: GribDef
   :show-inheritance:
   :members:
   :member-order: alphabetical

-----------------------------------------------------------

functions
---------

.. autofunction:: complete_grib_paths
.. autofunction:: complete_grib_samples_paths
.. autofunction:: complete_grib_definition_paths
.. autofunction:: set_definition_path
.. autofunction:: get_samples_paths
.. autofunction:: get_definition_paths
.. autofunction:: parse_GRIBstr_todict
.. autofunction:: read_gribdef

-----------------------------------------------------------

tables
------

.. autodata:: epygram.extra.griberies.tables.productionStatusOfProcessedData_dict
.. autodata:: epygram.extra.griberies.tables.typeOfGeneratingProcess_dict
.. autodata:: epygram.extra.griberies.tables.pyproj_geoid_shapes
.. autodata:: epygram.extra.griberies.tables.statistical_processes
.. autodata:: epygram.extra.griberies.tables.typeoffixedsurface2sample

-----------------------------------------------------------

defaults
--------

.. autodata:: epygram.extra.griberies.defaults.GRIB2_keyvalue
.. autodata:: epygram.extra.griberies.defaults.GRIB2_metadata_to_embark

.. autodata:: epygram.extra.griberies.defaults.GRIB1_sample
.. autodata:: epygram.extra.griberies.defaults.GRIB1_packing
.. autodata:: epygram.extra.griberies.defaults.GRIB1_ordering
.. autodata:: epygram.extra.griberies.defaults.GRIB1_keyvalue
