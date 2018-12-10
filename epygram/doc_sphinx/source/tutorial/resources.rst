Resources
=========

.. highlight:: python

What are the different formats implemented in my configuration ?
(some may be deactivated in config for simpler installation)

>>> import epygram
>>> epygram.init_env()
>>> epygram.config.implemented_formats
['netCDF', 'GRIB', 'GeoPoints', 'TIFFMF', 'FA', 'LFI', 'DDHLFA', 'LFA']

There is a *proxy* to open an existing resource without making an explicit reference to its format, e.g.:

>>> a_resource = epygram.formats.resource('ICMSHAROM+0042', 'r')

(where 'r' stands for 'read' opening mode and 'ICMSHAROM+0042' is the filename
of the resource), that uses the underneath guessing function,
that you may find useful (?)

>>> epygram.formats.guess('ICMSHAROM+0042')
'FA'

For opening a new resource in writing mode, specifying the format is though necessary:

>>> a_resource = epygram.formats.resource('ICMSHAROM+0042_new', 'w', fmt='FA')

All resources of different formats share a set of common methods, such as
:meth:`listfields`, :meth:`readfield()` and :meth:`writefield`, that behave
similarly: :meth:`readfield()` always return the same kind of objects:
:class:`Fields`.

-----------------------------------------------------------

Explore resources
-----------------
Let's open an historic FA file from AROME and explore it.

>>> r = epygram.formats.resource('ICMSHAROM+0042', 'r')
>>> r.format
'FA'
>>> type(r)
<class 'epygram.formats.FA.FA'>
>>> r.isopen
True
>>> r.openmode
'r'

But I may want later to add a field in this resource ? Let's reopen it with 'append' mode.

>>> r.close()
>>> r.isopen
False
>>> r.open(openmode='a')
>>> r.openmode
'a'

Where is my resource stored ?

>>> r.container.absdir
'/home/mary/worktmp/'
>>> r.container.basename
'ICMSHAROM+0042'
>>> print r.container.abspath
'/home/mary/worktmp/ICMSHAROM+0042'

Now let's explore what's inside this FA file.

>>> r.empty
False
# there are fields in there
>>> r.listfields()
['SURFTENS.TURB.ZO', 'SURFTENS.TURB.ME',
...
'S090WIND.V.PHYS']
>>> len(r.listfields())
1237

Some FA additional properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The case of FA is a bit specific in that the temporal and geometric metadata is
common to the whole resource and shared by all fields.
Therefore, the resource "has a geometry":

>>> print r.geometry
D3ProjectedGeometry containing:
    _center_lon: Angle containing:
        _degrees: 2.0
        _radians: 0.0349065850399
        _origin_value: 0.0349065850399
        _origin_unit: radians
    projection: 
        reference_lat: Angle containing:
            _cos_sin: (0.6858183529273761, 0.7277727576572106)
            _degrees: 46.7
            _radians: 0.815068760681
            _origin_value: 0.815068760681
            _origin_unit: radians
        reference_lon: Angle containing:
...

and optionally a spectral geometry as well:

>>> print r.spectral_geometry
SpectralGeometry containing:
    truncation: 
        in_X: 719
        in_Y: 767
        shape: elliptic
    space: bi-fourier

and also a validity (embedded in a list of 1 element, because time can be
a dimension of fields):

>>> print r.validity[0]
FieldValidity containing:
    _basis: 2014-12-01 00:00:00
    _date_time: 2014-12-02 18:00:00
    _statistical_process_on_duration: None
    _cumulativeduration: 3:00:00

Here, we can see the validity is 2014-12-02 18:00:00,
starting from *basis* 2014-12-01 00:00:00,
which mean the term of the resource is 18h.

Also, has been included a function to look for fields with a generic *seed*,
e.g.:

>>> r.find_fields_in_resource('*RAY*')
['SOMMFLU.RAY.SOLA', 'SURFFLU.RAY.SOLA', 'SOMMFLU.RAY.THER', 'SURFFLU.RAY.THER',
'S001RAYT SOL CL', 'S090RAYT SOL CL', 'S001RAYT THER CL', 'S090RAYT THER CL',
'SURFRAYT DIR SUR', 'TOPRAYT DIR SOM', 'SURFRAYT SOLA DE', 'SURFRAYT THER DE', 
'SOMMRAYT.SOLAIRE', 'SURFRAYT.SOLAIRE', 'SOMMRAYT.TERREST', 'SURFRAYT.TERREST']
>>> r.find_fields_in_resource('S06[1-3]WIND.?.PHYS')
['S061WIND.U.PHYS', 'S061WIND.V.PHYS', 'S062WIND.U.PHYS', 'S062WIND.V.PHYS',
'S063WIND.U.PHYS', 'S063WIND.V.PHYS']
>>> r.find_fields_in_resource(['S090TEMP*', 'SURF*'])
...

The encoding of fields is also available:

- on request:

  >>> r.fieldencoding('SURFTEMPERATURE')
  {'spectral': False, 'KSTRON': 0, 'KPUILA': 0, 'KNGRIB': 2, 'KNBITS': 16}
  >>> r.fieldencoding('SPECSURFGEOPOTEN')
  {'spectral': True, 'KSTRON': 0, 'KPUILA': 0, 'KNGRIB': 0, 'KNBITS': 0}

- and stored by time of reading:
  
  >>> r.readfield('S001TEMPERATURE')
  >>> r.fieldscompression
  {'S001TEMPERATURE': {'KNBPDG': 18, 'KSTRON': 106, 'KPUILA': 1, 'KNGRIB': 2, 'KNBCSP': 18},
  ...
  }

-----------------------------------------------------------

Field identifier (**fid**)
--------------------------

FA fields are identified by a character string name. Other formats may identify
fields differently, for instance GRIB with a set of **key:value** pairs.

As an example for GRIB, the :mod:`epygram.formats.GRIB.GRIB.listfields` method
returns a list of dicts:

>>> g = epygram.formats.resource('GRIDHSTFRANGP0025+0003', 'r')
>>> g.format
'GRIB'
>>> g.listfields()
[{'typeOfLevel': 'surface', 'indicatorOfTypeOfLevel': 1, 'name': 'Temperature',
'level': 0, 'table2Version': 1, 'editionNumber': 1, 'shortName': 't',
'paramId': 130, 'indicatorOfParameter': 11},
...
]

Field identifiers as an attribute of :doc:`../library/fields` objects will be
detailed in section :ref:`Field identifier <tuto-fid>` of the tutorial.

-----------------------------------------------------------

Juggling with resources
-----------------------

Transferring a field from one resource to another is almost as simple as
telling it:

>>> source_r = epygram.formats.resource('ICMSHAROM+0042', 'r')
>>> dest_r = epygram.formats.resource('ICMSHAROM+0042_bis', 'a')
>>> f = source_r.readfield('SURFTEMPERATURE')
>>> type(f)
<class 'epygram.fields.H2DField.H2DField'>
>>> dest_r.writefield(f)
# [2016/05/04-15:20:53][epygram.formats.FA][writefield:0980][INFO]: there
already is a field with the same name in this FA: overwrite.


-----------------------------------------------------------

Resource modifiers
---------------------

Sometimes, fields in a resource does not take the appropriate form.
For example, 3D plots require 3D fields; animations require fields with
time evolution; complex treatments can require to work on a subdomain;
and some applications can need variables that are not in the resource
but can be computed from it.
To deal with this problems, one can implement the transforms or use one of
the resource modifiers provided by epygram:

- CombineLevelsResource takes one resource and tries to expose 3D fields built
  from the H2D fields actually present in the resource
- MultiValiditiesResource takes several resources. Each resource must contain
  the same fields but for different validities; the new resource join the different
  fields to return a field with a time dimension
- SubdomainResource takes one resource and return the fields on
  a sub-domain defined by a geometry or by indexes bounds
- DiagnosticsResource takes one resource and tries to compute new fields from
  the fields contained in the resource. For example, one can request the
  temperature field and the resource returns it if it is already
  in the resource or computes it from potential temperature.

All this modifiers should work in a pipeline if needed (one can compute
diagnostics on a multivalidities resource for example).

Here are some examples to build a resource modifier:

>>> from footprints import proxy as fpx
>>> r = epygram.formats.resource(filename, 'r')
>>> rDiag = fpx.resource_modificator(name='Diagnostics', resource=r, openmode='r', ...)
>>> r3d = fpx.resource_modificator(name='CombineLevels', resource=r)
>>> r_time = fpx.resource_modificator(name='MultiValidities', resources=[r1, r2, ...])
>>> r_subdo = fpx.resource_modificator(name='Subdomain', resource=r, geometry=geom)

