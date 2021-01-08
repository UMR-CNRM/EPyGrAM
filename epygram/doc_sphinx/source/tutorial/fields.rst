Fields
======

.. highlight:: python



>>> import epygram
>>> epygram.init_env()
>>> r = epygram.formats.resource('ICMSHAROM+0042', 'r')
>>> field = r.readfield('S058WIND.U.PHYS')

**Components of a meteorological field**

A meteorological field has (mainly):

- a :ref:`fid <tuto-fid>`, that identifies the field nature
- a :ref:`geometry <tuto-geometry>`, that defines the position of the of the
  field's data; detailed in :doc:`../library/geometries_lib`
- a :ref:`validity <tuto-validity>`, that defines the temporal validity of the
  field's data
- a :ref:`spectral geometry <tuto-spectral>`, if the field is spectral.

>>> field.structure
'H2D'
>>> type(field.geometry)
<class 'epygram.geometries.H2DGeometry.H2DProjectedGeometry'>
>>> type(field.validity)
<class 'epygram.base.FieldValidityList'>

In the following, we focus on H2DField, but all kind of fields behave similarly. 

-----------------------------------------------------------

.. _tuto-fid:

Field identifier (fid)
----------------------

The fid of a field handle its *nature*, i.e. the physical parameter and the 
material (surface, soil, atmosphere) it represents (often with a detailed 
vertical location).
It is handled as a dict() which keys are formats names, because fields in files
are stored under different fid depending on the format, e.g.

>>> field.fid
{'FA':'S058WIND.U.PHYS',
 'GRIB1':{'indicatorOfParameter':33,
          'indicatorOfTypeOfLevel':109,
          'level':58,
          'table2Version':1,}
 'GRIB2':{'discipline':0,
          'parameterCategory':3,
          'parameterNumber':2,
          'typeofFirstFixedSurface':119,
          'level':58,}
}


Take care that for a field to be written in a resource, the only-but-mandatory
required field characteristics is that its fid in the resource format must be
present!

In other words, a field whose fid is ``{'FA':'SURFTEMPERATURE', 'GRIB':{...}}``
will be writeable either in a GRIB or a FA file, but not in a LFI or any other
format...

-----------------------------------------------------------

.. _tuto-validity:

Validity
--------
The ``validity`` attribute of a field is a list of
:class:`epygram.base.FieldValidity`, because a field can have a temporal dimension. 

>>> field.validity[0].getbasis()
datetime.datetime(2014, 12, 1, 0, 0)
>>> field.validity[0].getbasis(fmt='IntStr')
20141201000000
>>> field.validity[0].term()
datetime.timedelta(1, 64800)

Some fields may have a cumulative duration, e.g. precipitation fields are accumulated:

>>> field = r.readfield('SURFACCGRAUPEL')
>>> field.validity[0].statistical_process_on_duration()
'accumulation'
>>> field.validity[0].cumulativeduration()
datetime.timedelta(0, 10800)

The ``statistical_process_on_duration``, among which one can find min, max,
average and so on, is coded as GRIB2 norm
(http://apps.ecmwf.int/codes/grib/format/grib2/ctables/4/10), if available.

-----------------------------------------------------------

Fields useful methods
---------------------

.. note::
    *Autocompletion, in interactive (i)Python session or smart editors,
    may be an even better (than doc/tuto) way to explore the available methods
    of objects.*

.. _tuto-spectral:

Spectralness
^^^^^^^^^^^^

Spectral transforms are done ("in place") through the two methods
:meth:`Field.sp2gp() <epygram.fields.D3Field.D3Field.sp2gp()>` and
:meth:`Field.gp2sp(a_sp_geom) <epygram.fields.D3Field.D3Field.gp2sp()>`:

>>> field = r.readfield('S090TEMPERATURE')
>>> field.spectral
True

The :class:`epygram.geometries.SpectralGeometry` contains the kind of spectral
space (bi-Fourier//LAM or Legendre//global), the truncation(s), and the actual
spectral transforms routines. 

>>> spgeom = field.spectral_geometry
>>> field.sp2gp()
>>> field.spectral
False
>>> field.spectral_geometry
None
>>> field.gp2sp(spgeom)
>>> field.spectral
True


Data
^^^^

Some methods have been implemented to ease a comprehensive access to the data:

- basic statistics:

  >>> field.sp2gp()
  >>> field.stats()
  {'std': 6.1556479204416092, 'nonzero': 2211840, 'quadmean': 280.99889135008499, 'min': 259.33480158698774, 'max': 293.21439026360537, 'mean': 280.93145950331785}

- field value at some lon/lat point: 

  >>> field.getvalue_ll(1.5, 45.6) # default interpolation = 'nearest'
  274.278588481978
  >>> field.getvalue_ll(1.5, 45.6, neighborinfo=True) # get info about the nearest neighbor gridpoint used
  (274.278588481978, (1.4988145157970028, 45.60001936720281))
  >>> field.getvalue_ll(1.5, 45.6, interpolation='linear')
  274.25775674717687

Also, although the field's data is accessible through its attribute ``data``,
it is strongly advised to access the data through the method :meth:`Field.getdata`,
because the internal storage of the data may differ from expected by the user.

- modifying the field data should resemble:

  >>> data = field.getdata()
  >>> type(data)
  <type 'numpy.ndarray'>
  >>> data.shape
  (1536, 1440)
  >>> data[100:800,500:600] += 10*numpy.random.rand(700,100)
  >>> field.setdata(data)

  after what the field can of course be re-written in a resource.

- some patterned operations on fields are facilitated through
  the :meth:`Field.operation() <epygram.base.Field.operation>` method: any of the four basic operations (+,*,-,/)
  with scalars or any :mod:`numpy` function (exp, sin, log...):

  >>> field.operation('-', 273.15)  # e.g. go from K to °C
  >>> field.operation('sin')  # does field.data = numpy.sin(field.data) 

- of spectral fields can also be computed horizontal derivatives:

  >>> t = r.readfield('S045TEMPERATURE')
  >>> t.spectral
  True
  >>> (dx, dy) = t.compute_xy_spderivatives()
  >>> type(dx)
  <class 'epygram.fields.H2DField.H2DField'>
  >>> dx.spectral
  False
  >>> dx.max()
  0.0051387105385038408

- of 2D fields can be computed spectra (:class:`epygram.spectra.Spectrum`):

  >>> t.sp2gp()
  >>> s = t.dctspectrum()
  >>> type(s)
  <class 'epygram.spectra.Spectrum'>

Operations between fields
^^^^^^^^^^^^^^^^^^^^^^^^^

Operations between fields can be done in two ways:

- standard Python syntax; in case a new Field object is created,
  with uninitialized ``validity`` (what is the validity of an operation between
  two fields of potential different validity ?) and ``fid``:
  
  >>> field90 = r.readfield('S090TEMPERATURE')
  >>> field89 = r.readfield('S089TEMPERATURE')
  >>> field_diff = field90 - field89
  
- the :meth:`Field.operation() <epygram.base.Field.operation>` method;
  in case the field values are modified "in place":
  
  >>> field90 = r.readfield('S090WIND.U.PHYS')
  >>> field89 = r.readfield('S089WIND.U.PHYS')
  >>> field90.operation('+', field89)
  
In any case, a simple consistency check is done on the fields' geometry,
basically on their dimensions.

-----------------------------------------------------------

Building Vector Fields
----------------------

Wind fields (for instance) can be re-assembled from their U/V components
into :doc:`H2DVectorField <../library/H2DVectorField>` or
:doc:`D3VectorField <../library/D3VectorField>` for more integrated functionalities
(re-projection, computation of derivatives or direction/module, plotting and
so on...).

>>> u = r.readfield('S090WIND.U.PHYS')
>>> v = r.readfield('S090WIND.V.PHYS')
>>> wind = epygram.fields.make_vector_field(u,v)
>>> wind.sp2gp()

- reprojection: FA wind fields are projected on the grid axes (here, a Lambert
  projection); let's get the wind components on true zonal/meridian axes:
  
  >>> wind.getvalue_ij(0,0)
  [0.5525041298918116, -2.8212975453933336]
  >>> wind.reproject_wind_on_lonlat()
  >>> wind.getvalue_ij(0,0)
  [0.9307448483516318, -2.759376908801778]
  
- derivatives: just as the :meth:`Field.compute_xy_spderivatives` method
  enable to compute derivatives of spectral fields, the
  :meth:`H2DVectorField.compute_vordiv() <epygram.fields.H2DVectorField.compute_vordiv>` 
  method enable to compute vorticity and divergence of a spectral wind field:
  
  >>> wind.gp2sp(r.spectral_geometry)
  >>> (vor, div) = wind.compute_vordiv()
  >>> type(vor)
  <class 'epygram.fields.H2DField.H2DField'>
  
- direction/module: to compute a wind direction or wind module field from
  vectors: 

  >>> wind.sp2gp()
  >>> ff = wind.to_module()
  >>> type(ff)
  <class 'epygram.fields.H2DField.H2DField'>

-----------------------------------------------------------

Plots
-----

Cf. *Plots* sections of the `Gallery <../gallery/index.rst>`_.

-----------------------------------------------------------

3D Plots
--------
Fields can be plotted in a 3D view in three ways:

- contour with :meth:`Field.plot3DContour` method (which becomes a
  classical contour plot if field is 2D, vertical or horizontal)
- volume with :meth:`Field.plot3DVolume`
- color with :meth:`Field.plot3DColor` (corresponding to the matplotlib contourf method)

A vector field can be plotted using :meth:`Field.plot3DVector` (to plot arrows)
or :meth:`Field.plot3DStream` to plot (stream lines or tubes).

>>> import vtk #We need to import vtk before epygram even if do not use it directly in the script
>>> import epygram
>>> epygram.init_env() #initialisation of environment, for FA/LFI and spectrals transforms sub-libraries
>>> r = epygram.formats.resource(filename, 'r', true3d=True)
>>> 
>>> CF = r.readfield('S---CLOUD_FRACTI')
>>> 
>>> #Set-up of the view
... offset = CF.geometry.gimme_corners_ll()['ll'] #We translate the domain
>>> hCoord = 'll' #We use lat/lon on the horizontal
>>> z_factor = 0.1 #0.1 horizontal degree of lat/lon is represented by the same length as one model level on the vertical
>>> ren = epygram.util.vtk_set_window((0.5, 0.5, 0.5), (800, 800))
>>> 
>>> CF.plot3DContour(ren, [1.], color='White', hCoord=hCoord, offset=offset, z_factor=z_factor)
((vtkRenderingOpenGL2Python.vtkOpenGLActor)0x7f899c230390, (vtkRenderingOpenGL2Python.vtkOpenGLPolyDataMapper)0x7f89b8750a78)
>>> 
>>> ren['interactor'].Start()
