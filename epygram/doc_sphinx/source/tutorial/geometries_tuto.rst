Geometries
==========

.. highlight:: python

Again, fields objects (and their :attr:`geometry` attribute) are rather
independent from the resource they originate. In the following, we handle
geometries and fields coming from FA files, for convenience only.

>>> import epygram
>>> epygram.init_env()
>>> arp = epygram.formats.resource('icmsharpe+0000', 'r')  # an ARPEGE file
>>> aro = epygram.formats.resource('ICMSHAROM+0042', 'r')  # an AROME file

Horizontal geometries
---------------------

The most usual handled fields are horizontal, 2D fields.
``epygram`` derive geometries in different classes depending on the kind of grid
(Gauss grids, regular Lon/Lat, regular grids on plane projections, or 
unstructured grids).

>>> t_arp = arp.readfield('CLSTEMPERATURE')
>>> type(t_arp.geometry)
<class 'epygram.geometries.H2DGeometry.H2DGaussGeometry'>
>>> t_aro = arp.readfield('CLSTEMPERATURE')
>>> type(t_aro.geometry)
<class 'epygram.geometries.H2DGeometry.H2DProjectedGeometry'>

>>> g_arp = t_arp.geometry; g_aro = t_aro.geometry  # shortcuts
>>> g_arp.rectangular_grid  # reduced Gauss grids have a variable number of longitudes depending on the latitude
False
>>> g_aro.rectangular_grid
True
>>> g_arp.name
'rotated_reduced_gauss'
>>> g_aro.name
'lambert'


Noticeable attributes
^^^^^^^^^^^^^^^^^^^^^

- :attr:`geometry.dimensions` : grid dimensions (including C+I+E zones for LAM)
- :attr:`geometry.grid` : grid other parameters (positionning, definition, resolution)
- :attr:`geometry.projection` : parameters of projection, if appropriate
- :attr:`geometry.vcoordinate` : the vertical characteristics of a horizontal field

Note that lon/lat angles as parameters of a geometry are stored in
:class:`epygram.util.Angle` objects, in order to:

- ensure to keep bit-reproducibility of original values (as read in a given format)
- enable easy unit conversion (degrees, radians, (cos, sin) tuples...)

>>> g_aro.projection['reference_lat']._origin_unit
'radians'
>>> g_aro.projection['reference_lat'].get('degrees')
46.70000000000001
>>> g_aro.projection['reference_lat']._origin_value
0.815068760681
>>> numpy.radians(numpy.degrees(g_aro.projection['reference_lat']._origin_value))
0.81506876068135203

Most useful methods
^^^^^^^^^^^^^^^^^^^

- get grid, eventually directly a subzone of (LAM)

  >>> (lons, lats) = g_aro.get_lonlat_grid()
  >>> corners = g_aro.gimme_corners_ll(subzone='C')  # get the corners lon/lat of the C zone

- lon/lat <=> grid indexes

  >>> g_aro.ll2ij(1.25, 45.666)  # lon/lat to grid indexes
  (669.16059153010667, 673.76289192747686)
  >>> g_aro.ij2ll(669, 674)  # grid indexes to lon/lat
  (1.2472761341493812, 45.668753485393296)

- create a Basemap object (:mod:`mpl_toolkits.basemap`), as much appropriate
  to the geometry or with requested specificities

  >>> g_aro.make_basemap()
  >>> # or
  >>> g_arp.make_basemap(specificproj='kav7')

...

Vertical geometries
-------------------

Vertical kind of geometries are usually built from either resources or
3D fields:

>>> p = aro.extractprofile('S*TEMPERATURE', 1.26, 45.3)
>>> type(p)
<class 'epygram.fields.V1DField.V1DField'>
>>> p.geometry.vcoordinate.typeoffirstfixedsurface
119
>>> # levels originate from FA historic file : 119 = GRIB2 code for hybrid-pressure vertical coordinate
>>> from epygram.geometries.VGeometry import hybridP2pressure  # vertical coordinate conversion functions
>>> vert_coord_as_pressure = hybridP2pressure(p.geometry.vcoordinate, Psurf=102500, vertical_mean='geometric')
>>> # and eventually
>>> p.geometry.vcoordinate = vert_coord_as_pressure  # but be careful: no consistency check is done here
>>> # actually, this kind of transforms is integrated:
>>> p = aro.extractprofile('S*TEMPERATURE', 1.26, 45.3, vertical_coordinate=100)
>>> p.geometry.vcoordinate.typeoffirstfixedsurface
100

Cf. http://apps.ecmwf.int/codes/grib/format/grib2/ctables/4/5 for vertical coordinate types.
