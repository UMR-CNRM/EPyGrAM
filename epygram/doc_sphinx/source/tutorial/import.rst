Import
======

.. highlight:: python

First of all, import ``epygram`` ! and call the initialization function
(to position environment variables for some compiled sub-modules). 

>>> import epygram
>>> epygram.init_env()

Then you may check your version and config (:mod:`epygram.config`) ?

>>> epygram.__version__
0.6.9
>>> epygram.showconfig()
######################
### epygram.config ###
######################
- Cpd = 1004.709
- FA_buffered_gauss_grid = True
- FA_default_compression = {'KNGRIB': 2, 'KDMOPL': 0, 'KPUILA': 0, 'KSTRON': 0, 'KNBPDG': 16, 'KNBCSP': 16}
... #(and so on) ...

Create objects
^^^^^^^^^^^^^^

In a general way in ``epygram``, users should not have to generate an object
using its class name and constructor (except in the library development, inside
internal methods). First, because it may be painful to do so (objects can be
rather complex...). Second, because there surely is a method to do it in an
integrated, comprehensive manner (except if you are doing very exotic things) !

Explore objects
^^^^^^^^^^^^^^

Generally, when wondering about the internal components of an object, just type

>>> print <my_object>

Another, better-looking way of inspecting objects (fields, geometries,
resources and validities) is the :meth:`what` method, that dumps info into
a file-like object:

>>> my_object.what(open('dump_of_my_object.txt', 'w'))
>>> # or
>>> import sys
>>> my_object.what(sys.stdout)

Customize installation
^^^^^^^^^^^^^^^^^^^^^^

The ``epygram`` library is designed to be customized easily by a given user:
modifying defaults, bringing modifications to existing classes (Fields, Formats,
Geometries) or adding new (fields, formats) classes to be used "on the fly" by
the library. 

Yet, for a cleaner partitioning of installations and avoid unwanted overwriting
of files, the library is able to look for source files in *user* directory(ies).

Declaration:

- *configuration*: values of :mod:`epygram.config` can be modified in
  ``$HOME/.epygram/userconfig.py``
- *classes* (modified or new): fill the **usermodules** variable in
  ``$HOME/.epygram/userconfig.py``, with the below formalism:
>>> usermodules = [{'name':'module1', 'abspath':'/path/to/module1.py'},
>>>                {'name':'module2', 'abspath':'/path/to/module2.py'}
>>>                ]
