Introduction
============

.. highlight:: python

First of all, import ``epygram`` ! and call the initialization function
(to position environment variables for some compiled sub-modules). 

>>> import epygram
>>> epygram.init_env()

Then you may check your version and config (:mod:`epygram.config`) ?

>>> epygram.__version__
1.1.6
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
rather complex...). Second, because there surely is a method to do it with an
integrated, comprehensive manner (except if you are doing very exotic things) !

Explore objects
^^^^^^^^^^^^^^^

Generally, when wondering about the internal components of an object, just type

>>> print my_object

the representation method of these objects dig into its attributes
(which can be verbose...).

Another, better-looking way of inspecting objects (fields, geometries,
resources and validities) is the :meth:`what` method, that dumps organized,
formatted info to the stdout (or to a file-like object):

>>> my_object.what() # to stdout

