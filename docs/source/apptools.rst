Command-Line applicative tools
==============================

.. highlight:: python

Some applicative tools in command line are provided and installed by pip.

These tools are available through a single command line ``epygram`` with sub-commands:

- ``epygram -h`` to list available sub-commands
- ``epygram <sub-command> -h`` for auto-documentation of each tool/sub-command.

or as ``epy_<sub-command>`` (pip should have placed them in your ``$PATH``).

Example, to plot a field:

- ``epygram cartoplot <file> -f <field>``

or

- ``epy_cartoplot <file> -f <field>``

are equivalent to 

- ``epy_cartoplot.py <file> -f <field>`` in versions prior to 1.6.0

See :ref:`Commented list of commands <subcommands>`
