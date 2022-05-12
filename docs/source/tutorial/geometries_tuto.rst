Geometries
==========

.. _tuto-geometry:

.. highlight:: python

Again, fields objects (and their :attr:`geometry` attribute) are rather
independent from the resource they originate. In the following, we handle
geometries and fields coming from FA files, for convenience only.

>>> import epygram
>>> epygram.init_env()
>>> arp = epygram.formats.resource('icmsharpe+0000', 'r')  # an ARPEGE file
>>> aro = epygram.formats.resource('ICMSHAROM+0042', 'r')  # an AROME file

-----------------------------------------------------------

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

- ensure to keep bit-reproducibility of original values (as read in a given
  format, for instance)
- enable easy unit conversion (degrees, radians, (cos, sin)...)

>>> g_aro.projection['reference_lat']._origin_unit
'radians'
>>> g_aro.projection['reference_lat'].get('degrees')
46.70000000000001
>>> g_aro.projection['reference_lat']._origin_value
0.815068760681
>>> g_aro.projection['reference_lat'].get('radians')
0.815068760681


Most useful methods
^^^^^^^^^^^^^^^^^^^

.. note::
    *Autocompletion, in interactive (i)Python session or smart editors,
    may be an even better (than doc/tuto) way to explore the available methods
    of objects.*

- get grid, eventually directly a subzone of (LAM)

  >>> (lons, lats) = g_aro.get_lonlat_grid()
  >>> corners = g_aro.gimme_corners_ll(subzone='C')  # get the corners lon/lat of the C zone

- lon/lat <=> grid indexes

  >>> g_aro.ll2ij(1.25, 45.666)  # lon/lat to grid indexes
  (669.16059153010667, 673.76289192747686)
  >>> g_aro.ij2ll(669, 674)  # grid indexes to lon/lat
  (1.2472761341493812, 45.668753485393296)

...

-----------------------------------------------------------

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

These V1D fields are also *plotable*, by the way:

>>> p.plotfield(title='A simple profile') 

-----------------------------------------------------------

3D plot
-------
We can plot, whis a 3D rendering, a lat/lon image on a surface.
The surface is described by the vertical coordinate of the geometry.
Two methods are available to plot the NOAA's bluemarble image or
tiles from a map tiles server.

>>> import vtk #We need to import vtk before epygram even if do not use it directly in the script
>>> import epygram
Vortex 1.4.0 loaded ( Friday 19. October 2018, at 14:40:31 )
>>> epygram.init_env() #initialisation of environment, for FA/LFI and spectrals transforms sub-libraries
>>> r = epygram.formats.resource(filename, 'r')
>>> 
>>> #We need a geometry containing the altitude of the ground
... zs = r.readfield('SPECSURFGEOPOTEN') #surface geopotential
>>> zs.sp2gp() #convert spectral data into grid points
>>> zs.setdata(zs.getdata() / 9.8) #Convert geopotential height into height
>>> zs.use_field_as_vcoord(zs, 103) #We replace the vertical coordinate of the field by the height values
>>> 
>>> #Set-up of the view
... offset = zs.geometry.gimme_corners_ll()['ll'] #We translate the domain
>>> hCoord = 'll' #We use lat/lon on the horizontal
>>> z_factor = 0.0005 #0.0005 horizontal degree of lat/lon is represented by the same length as one meter on the vertical
>>> ren = epygram.util.vtk_set_window((0.5, 0.5, 0.5), (800, 800))
>>> 
>>> zs.geometry.plot3DMaptiles(ren,                                                   #window to plot on
...                            "https://a.tile.openstreetmap.org/${z}/${x}/${y}.png", #url of the map tiles server
...                            2,                                                     #ratio between field and tiles resolutions
...                            interpolation='linear',                                #interpolation method
...                            hCoord=hCoord, z_factor=z_factor, offset=offset)       #helpful only on the first 3D plot
>>> 
>>> ren['interactor'].Start()



