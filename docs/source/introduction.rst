Introduction
============

The ``epygram`` library package is a set of Python classes and functions
designed to handle meteorological fields in Python, as well as interfacing
their storage in various usual (or not) data formats.



For who, for what ?
-------------------

The purpose of the library is to provide user-friendly interfaces to formats
of data resources as FA, LFI, LFA (Météo France historical formats), GRIB/1-2
(...) and handy classes for manipulating meteorological fields, whatever their
geometry is.

Therefore, it has been conceived basically for people who need to process data
from Numerical Weather Prediction (NWP) or Climate Modelisation, whatever they
need to do with, as soon as they build their application in Python.

One strong will during package development was to give the most simple access
to data and meta-data (its description in time and space, basically). The
object-oriented design seemed to be the most appropriate manner to reach this
goal. Thus, a set of features has been developped within the ``epygram`` classes,
so that the user should: 

- 1/ not really need to dig into the objects, only their methods;

- 2/ wonder (and hence search the doc) if there is not already a
  functionality dealing with the data/meta-data processing he faces;
 
- 3/ print the object to have a recursively indented overview of its components.



.. _general-design:

General design
--------------

There are 3 basic concepts used in ``epygram``: **fields, geometries, resources**.

- a **field** is the union of a data and its meta-data: its identification, its 
  temporal validity and geographical description, and eventually further
  documentation. There is no *a priori* about the geometry of the field,
  it can be either a horizontal surface, a vertical profile, a transect, a
  vertical section, a single point, or a 3D cloud of points.

  The field can even be represented in spectral space.

  In any case, a set of basic features has been intimately attached to the 
  field, providing it handy manipulation.

- the **geometry** of a field is a set of parameters and methods that enables,
  basically, to know precisely what is the 3D-earth-round localization of any 
  point of the field it describes. As listed above, a geometry can be of 
  several natures and dimensions.

- a **resource** is an aggregation of fields, stored in a given data and meta-data
  format. The resource is dissociated with its container, *i.e.* the 
  system and physical support is it written on (*e.g.* a file on disk, a memory
  address, remote database...).

A 4th important element is the **fid** (field identifier) of a field,
which identifies its nature and can index it inside resources.

The ``epygram`` package hence provides as much as possible easy read/write of
fields from/to resources, as well as basic features on fields and geometries,
described further in the documentation.



Genesis
-------

The ``epygram`` conception comes from the observation that more and more people
in Météo France / CNRM-GAME were coming to Python in order to process NWP/CM 
data. Rather than letting everyone face the same obstacles with geometries and
historical data formats, raised the good idea to develop a Python package once 
for all, releasing time for scientific activities rather than technical common
problems...


