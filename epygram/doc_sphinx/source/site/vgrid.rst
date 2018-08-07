:mod:`vgrid` --- Contains handling of Hybrid-Pressure vertical grid generation and plotting.
============================================================================================

.. automodule:: vgrid
   :synopsis: Contains handling of Hybrid-Pressure vertical grid generation and plotting.

-----------------------------------------------------------

Classes
-------

.. autoclass:: HybridPressureVGrid
   :members:
   :member-order: alphabetical

-----------------------------------------------------------

Standard atmosphere
-------------------

.. automodule:: vgrid.standard_atmosphere

.. autofunction:: pressure_at
.. autofunction:: altitude_at

It uses underlying Fortran routines that may need to be recompiled (to .so)
platform-wise. Sources and Makefile are provided in the sub-directory
``vertical_discretization``.

-----------------------------------------------------------

Tools
-----

``mkvgrid.py`` is a command-line tool used to generate a HybridPressure Vertical
grid according to a set of requirements in terms of:
  
  - number of levels
  - height of first level, pressure of last level
  - number of levels by slices of atmosphere
  - thickness of these levels
  - dynamics options on how to compute levels

It uses an underlying Fortran program that may need to be recompiled
platform-wise. Sources and Makefile are provided in the sub-directory
``vertical_discretization``.

The tool is autodocumented with option ``-h``.

These functionalities are also accessible from python:

.. code-block:: python
	
	from vgrid.mkvgrid import generate_vertical_grid
	generate_vertical_grid(...)

or

.. code-block:: python
		
	from vgrid.mkvgrid import main as mkvgrid
	mkvgrid(...)

.. automodule:: vgrid.mkvgrid

.. autofunction:: generate_vertical_grid
.. autofunction:: main
 