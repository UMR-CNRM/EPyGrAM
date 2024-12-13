#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Plugins: add functionalities to standard epygram classes.

Plugins must be packages which name starts with "with_*" and contain an
activate() function in addition to the definition of extensions.
"""

import os

# actual import
fatal = False
implemented = [m for m in os.listdir(os.path.dirname(os.path.abspath(__file__)))
               if m.startswith('with_')]
_successful_import = {}
_failed_import = {}
for p in implemented:
    try:
        import importlib
        pkg = importlib.import_module('.' + p, __name__)
    except Exception as e:  # Exception: we need to catch any exception raised when an import fails
        if fatal:
            raise e
        else:
            _failed_import[p] = e
    else:
        _successful_import[p] = pkg


def why_plugins_failed():
    """Tells the reason why eventual plugins import failed."""
    for p, e in _failed_import.items():
        print("Plugin:{} -> Error:{}".format(p, str(e)))


def available():
    """Return the list of available plugins."""
    return tuple(_successful_import.keys())


def activate(plugin):
    """Activate the required plugin."""
    if plugin not in implemented:
        raise NotImplementedError("plugin '{}'".format(plugin))
    elif plugin in _failed_import:
        raise ImportError(("An error ({}) occurred trying to import the '{}' plugin " +
                           "(probably because of a missing dependency).").format(str(_failed_import[plugin]), plugin))
    else:
        _successful_import[plugin].activate()
