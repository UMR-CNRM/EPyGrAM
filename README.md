EPyGrAM
=======

__*Enhanced Python for Graphics and Analysis of Meteorological fields*__

---

The epygram library package is a set of Python classes and functions designed to handle meteorological fields in Python, as well as interfacing their storage in various usual (or not) data formats.

Dependencies
------------

EPyGrAM dependencies are available from Pypi (pip install ...), and listed in `pyproject.toml`.
Some packages are mandatory, others are optional, only necessary for the use of specific functionalities or formats.
Formats for which the import of the according underlying package fails are deactivated at runtime.

Installation
------------

`pip install epygram`

or

`pip3 install epygram`

To use specific functionalities which dependencies are not covered by default,
you may need to manually pip install the according package(s).

You can also install all optional dependencies using:

`pip install epygram[all]`

or

more specifically one of the extra group of dependencies:

`pip install epygram[<opt>]`

with `<opt>` among (`graphics`, `docs`, `features`, `extra_formats`), cf. `pyproject.toml`.

Tests
-----

To run tests, cf. [`tests/README.md`](tests/README.md).

Documentation
-------------

To generate Sphinx doc: `make doc`. It will be generated in `docs/build/html`.
Online doc of the latest release on `master` branch is available at https://umr-cnrm.github.io/EPyGrAM-doc

Applicative tools
-----------------

Some applicative tools in command line are provided and installed by pip.

These tools are available through a single command line `epygram` with sub-commands:

- `epygram -h` to list available sub-commands
- `epygram <sub-command> -h` for auto-documentation of each tool/sub-command.

or as `epy_<sub-command>` (pip should have placed them in your `$PATH`).

Example, to plot a field:

- `epygram cartoplot <file> -f <field>`

or

- `epy_cartoplot <file> -f <field>`

are equivalent to

- `epy_cartoplot.py <file> -f <field>` in versions prior to 1.6.0

License
-------

This software is governed by the open-source [CeCILL-C](http://www.cecill.info) license under French law, cf. LICENSE.txt.
Downloading and using this code means that you have had knowledge of the CeCILL-C license and that you accept its terms.

