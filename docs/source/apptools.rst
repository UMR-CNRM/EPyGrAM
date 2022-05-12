Applicative tools
=================

.. highlight:: python

A series of applicative tools, to be used in shell command-line, is also
developped and distributed along epygram.
They are available in directory ``apptools``.

Each tool is self-documented: ``epy_tool.py -h`` 

Currently available:

- ``epy_what.py``: explore the contents of a file
- ``epy_conv.py``: format converter of fields
- ``epy_movefield.py``: move fields from one file to another, with optional
  *on-the-fly* basic operations
- ``epy_delfield.py``: delete fields from a file
- ``epy_profile.py``: extract and plot a vertical profile from 3D fields
- ``epy_section.py``: extract and plot a vertical section from 3D fields
- ``epy_stats.py``: compute basic statistics on fields
- ``epy_plot.py``: plot horizontal 2D fields using basemap package (DEPRECATED)
- ``epy_cartoplot.py``: plot horizontal 2D fields using cartopy package
- ``epy_hist.py``: compute and plot a histogram of fields values
- ``epy_point.py``: extract the value of fields on a given location
- ``epy_spectrum.py``: compute and plot the DCT spectrum of a field
- ``fa_sp2gp.py``: convert to gridpoint the spectral fields of a FA file
- ``domain_maker.py``: interactive construction and display of a new LAM domain
- ``ddhlfa_plot.py``: basic plots of LFA-formatted outputs from DDH
- ``run_epyweb.py``: local fields-plotting html server (requires ``vortex``)
