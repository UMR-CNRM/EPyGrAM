EPyGrAM
=======

__*Enhanced Python for Graphics and Analysis of Meteorological fields*__

---

The epygram library package is a set of Python classes and functions designed to handle meteorological fields in Python, as well as interfacing their storage in various usual (or not) data formats.

Dependencies
------------

EPyGrAM dependencies are available from Pypi (pip install ...), including the packages [`footprints`](https://pypi.org/project/footprints/) and [`bronx`](https://pypi.org/project/bronx/) from the Vortex project.

Installation
------------

Installation can be done manually by setting up path to `epygram` package and `site` subdirectories in your PYTHONPATH.
Installation can be done in local user site-packages with the help of the `_install/setup_epygram.py` script.
This script usage is autodocumented (option `-h`) about its options.

Then follow instructions printed by the script to complete installation.

NOTE: from version 1.5.0 (not released yet) it is planned to set EPyGrAM on PyPI so as to install it with `pip`.

Tests
-----

To run tests, cf. [`tests/README.md`](tests/README.md).

Documentation
-------------

To generate Sphinx doc: `make doc`. It will be generated in `docs/build/html`.
Online doc of the latest release on `master` branch is available at https://umr-cnrm.github.io/EPyGrAM-doc

License
-------

This software is governed by the open-source [CeCILL-C](http://www.cecill.info) license under French law, cf. LICENSE.txt.
Downloading and using this code means that you have had knowledge of the CeCILL-C license and that you accept its terms.

