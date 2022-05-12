EPyGrAM
=======

__*Enhanced Python for Graphics and Analysis of Meteorological fields*__

---

The epygram library package is a set of Python classes and functions designed to handle meteorological fields in Python, as well as interfacing their storage in various usual (or not) data formats.

Dependencies
------------

EPyGrAM dependencies are mostly available from Pypi (pip install ...), except a few site packages from the [Vortex](https://opensource.umr-cnrm.fr/projects/vortex) project, mainly the `footprints` and `bronx` packages.

Installation
------------

Installation can be done in local user site-packages with the help of the `_install/setup_epygram.py` script.
This script usage is autodocumented (option `-h`) about its options.

Then follow instructions printed by the script to complete installation.

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

