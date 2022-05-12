Some examples
=============

*These examples are meant to illustrate the relative simplicity of
handling and use of the library, as well as a quick glance at its
possibilities.*

For a more comprehensive introduction to the package, please refer to the :ref:`Tutorial <tuto-index>`.

For more advanced examples, please refer to the :ref:`Gallery <tuto-index>`.

About Resources
---------------

- Open a resource and list its contents:

  .. code-block:: python

    >>> import epygram
    >>> epygram.init_env()
    >>> fcst5 = epygram.formats.resource(filename='../ICMSHAROM+0005', openmode='r')

  And, get the list of its fields:

  .. code-block:: python

    >>> fcst5.listfields()
    ['SPECSURFGEOPOTEN', 'SURFPRESSION', 'S001WIND.U.PHYS', 'S002WIND.U.PHYS', ...]

- Read a field in a resource:
  
  .. code-block:: python

    >>> surftemp = fcst5.readfield('SURFTEMPERATURE')
    >>> print surftemp
    H2DField containing:
        comment: __unknown__
        fid: 
            FA: SURFTEMPERATURE
        geometry: H2DGeometry containing:
            projection:
               ...
        processtype: forecast
        spectral_geometry: __unknown__
        data: [[ 279.93986033   ...,  296.24581933]]
        validity: FieldValidity containing:
            ...

And many other features depending on the format...

--------------------------------------------------------------------------------

FA: Some specific features
..........................

- Open a new FA sharing its geometry and validity with an existing one:

  .. code-block:: python

    >>> import epygram
    >>> epygram.init_env()
    >>> fcst5 = epygram.formats.resource(filename='../ICMSHAROM+0005', openmode='r', fmt='FA')
    >>> new_fcst5 = epygram.formats.resource(filename='../ICMSHAROM+0005_bis', openmode='w',
                                             fmt='FA',
                                             headername=fcst5.headername, validity=fcst5.validity)

- If you want it to be even more similar, also add arguments:

  .. code-block:: python

    >>> new_fcst5 = epygram.formats.resource(filename='../ICMSHAROM+0005_bis', openmode='w',
                                             fmt='FA',
                                             headername=fcst5.headername, validity=fcst5.validity,
                                             cdiden=fcst5.cdiden,
                                             default_compression=fcst5.default_compression,
                                             processtype=fcst5.processtype)

--------------------------------------------------------------------------------

About Fields
------------

- Read a spectral field, convert it to gridpoint, modify values, and get back
  to spectral space:

  .. code-block:: python

    >>> s30temp = fcst5.readfield('S030TEMPERATURE')
    >>> s30temp.spectral
    True
    >>> sp_geom = s30temp.spectral_geometry # save info about spectral geometry (lost after conversion)
    >>> s30temp.sp2gp() # conversion (in place) to gridpoint
    >>> s30temp.spectral
    False
    >>> s30temp.mean(subzone='CI')
    268.2468305845095
    >>> s30temp.spectral_geometry
    None
    >>> data = s30temp.data
    >>> data[10:-10,20:-20] = data[10:-10,20:-20] + 2 # heat the "center" by 2K
    >>> s30temp.setdata(data)
    >>> s30temp.mean(subzone='CI')
    270.14124437511464
    >>> s30temp.gp2sp(sp_geom)
    >>> s30temp.spectral
    True

- Compute the magnitude of wind rotation between two levels (it's just an 
  example...):

  .. code-block:: python

    >>> wfields = fcst5.readfields('S03[0-1]WIND.?.PHYS')
    >>> type(wfields)
    <class 'epygram.base.FieldSet'>
    >>> wfields.listfields('FA')
    ['S030WIND.U.PHYS', 'S031WIND.U.PHYS', 'S030WIND.V.PHYS', 'S031WIND.V.PHYS']
    >>> for f in wfields:
    ...     f.sp2gp()
    ... 
    >>> du = wfields[0]-wfields[1]
    >>> dv = wfields[2]-wfields[3]
    >>> rot = du*du + dv*dv
    >>> rot.mean()
    1.3670194410302572
    >>> rot.plotfield() # rot is still a Field object

- Compute the windspeed:

  .. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> u = fcst5.readfield('S030WIND.U.PHYS')
    >>> v = fcst5.readfield('S030WIND.V.PHYS')
    >>> vectwind = epygram.fields.make_vector_field(u, v)
    >>> vectwind.sp2gp()
    >>> FF = vectwind.to_module()

--------------------------------------------------------------------------------

About Geometry
--------------

- Get the lon/lat coordinates of the whole grid of a field, the corresponding
  max map factor, ask for the lower-right corner and checks whether a point is inside grid:

  .. code-block:: python

    >>> (lons, lats) = s30temp.geometry.get_lonlat_grid()
    (array([[ -8.3539362, ..., 12.6589942],
                          ...
            [-11.8154972, ..., 16.2188262]]),
     array([[ 37.3330505, ..., 37.3009714],
                          ...
            [ 53.2651904, ..., 53.2224197]]))
    >>> s30temp.geometry.map_factor_field().max()
    1.0108440620737038
    >>> s30temp.geometry.gimme_corners_ll()['lr']
    (12.658994623143194, 37.300971400173346)
    >>> s30temp.geometry.point_is_inside_domain(15.0, 37.0)
    False

