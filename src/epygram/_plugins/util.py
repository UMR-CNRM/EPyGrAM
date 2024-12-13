#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

def notify_doc_requires_plugin(functions, plugin):
    """Adds a note in the doc of functions/methods, specifying the need for the plugin."""
    plugin = plugin.split('.')[-1]
    for f in functions:
        f.__doc__ = ".. note:: Requires plugin: **{}** (config.activate_plugins)\n".format(plugin) + f.__doc__
