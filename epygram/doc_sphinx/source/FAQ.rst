FAQ
===
.. highlight:: python

The Library
-----------

I have rather esoterical errors when I try to read a FA file...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   => the :func:`epygram.init_env` may be initialized at the beginning of your
   script(s), to be sure the :mod:`arpifs4py` library do not use MPI, or uses
   OpenMP with an initialized environment.

There is a field in my FA file, unknown by ``epygram`` (or it seems to be recognized as a MiscField); how can I manage to read it ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   => edit ``$HOME/.epygram/user_Field_Dict_FA.csv``, and add your field in the
   manner of the examples given therein.
   And warn the ``epygram`` team about it, so that it will enter next version
   default ``Field_Dict_FA.csv``.

I have built a new format class ``FMFILE`` (or whatever else) for ``epygram``, and I want it to be "fully" integrated in the package locally on my platform (so that the :func:`epygram.formats.resource` can return it). How can I do ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   => edit ``$HOME/.epygram/userconfig.py``, and:
     - add
     
       >>> usermodules = [{'name':'FMFILE', 'abspath':'/path/to/FMFILE.py'}]
     - copy the variable :obj:`epygram.config.implemented_formats` in it, adding ``FMFILE``

How do ``epygram`` and ``vortex`` interconnect ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   => Well ! Here is an example (having ``vortex`` installed),
   that fetches a resource on the archive and returns an ``epygram`` resource:
   
   >>> import epygram
   >>> epygram.init_env()
   >>> import usevortex
   >>> r = usevortex.get_resource(getmode='epygram',   # for the function to return the resource as an epygram object
                                  experiment='864G',                # XPID
                                  block='forecast',                 # Olive 'forecast' block (directory in archive)
                                  kind='gridpoint',                 # post-processed fields
                                  nativefmt='grib',                 # GRIBbed files
                                  date='2015041500',                # initial date and time
                                  term=3,                           # forecast term
                                  geometry='frangp0025',            # BDAP domain
                                  local='fcst_[term].[nativefmt]')  # local filename, once fetched
      
   Other resource descriptors are available, cf. :func:`usevortex.get_resource` documentation.
      
   ``vortex`` also is able to use ``epygram`` in order to handle a file's content, cf. Vortex doc.

I want to add new features or methods to a class, :class:`epygram.fields.H2DField` for instance, and be sure that my modifications will not be overwritten at the next upgrade of ``epygram``...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  

   => build your class ``myH2DField`` in ``/path/to/myH2DField.py``, making it inherit from :class:`epygram.fields.H2DField`, as follows:
  
     .. code-block:: python
     
       #!/usr/bin/env python
       # -*- coding: utf-8 -*-
       import copy
       import footprints
       footprints.priorities.set_after('default','user')
       from epygram.fields import H2DField
       
       class myH2DField(H2DField):
           _footprint = dict(
               priority = dict(
                   level = footprints.priorities.top.level('user')
               )
           )

     For this class to be used by ``epygram``, you simply have to add it in ``$HOME/.epygram/userconfig.py``:
     
     >>> usermodules = [{'name':'myH2DField', 'abspath':'/path/to/myH2DField.py'}]
    
    Anyway, if your modifications may be useful to others, propose to the ``epygram`` team its integration in the next version !
    
I want to add a personal colormap to be used by ``epygram``.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   => write your colormap as RGB as below
    
    .. code-block:: python
    
        r1,g1,b1;
        r2,g2,b2;
        ...
        rn,gn,bn
    
    into file ``mycolormap.cmap``. You may need the help of http://colormap.org
    
    Then in ``$HOME/.epygram/userconfig.py`` add:
    
    >>> usercolormaps = {'mycolormap', '/path/to/mycolormap.cmap'}
    
    and the colormap is now accesible to ``epygram``. 


**(to be continued...)**



Applicative tools
-----------------

**Option -h is your best friend !**

