:mod:`griberies` --- utilities around GRIB
==========================================

.. automodule:: griberies
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

.. autodata:: griberies.tables.productionStatusOfProcessedData_dict
.. autodata:: griberies.tables.typeOfGeneratingProcess_dict
.. autodata:: griberies.tables.pyproj_geoid_shapes
.. autodata:: griberies.tables.statistical_processes
.. autodata:: griberies.tables.typeoffixedsurface2sample

-----------------------------------------------------------

defaults
--------

.. autodata:: griberies.defaults.GRIB2_keyvalue
.. autodata:: griberies.defaults.GRIB2_metadata_to_embark

.. autodata:: griberies.defaults.GRIB1_sample
.. autodata:: griberies.defaults.GRIB1_packing
.. autodata:: griberies.defaults.GRIB1_ordering
.. autodata:: griberies.defaults.GRIB1_keyvalue
