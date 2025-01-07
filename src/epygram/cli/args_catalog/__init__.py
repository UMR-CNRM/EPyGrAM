#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains a catalog of reusable command-line arguments for argparse,
with a function *add_arg_to_parser* to import them to an *argparse* parser.

The arguments are classified by categories, each category being a dict.

Each argument is a list composed as follows:
['name', 'alternate_name1', 'alternate_name2', ..., dict(kw1, kw2, kw3, ...)]
where the keywords 'kwi' are argparse.ArgumentParser.add_argument() optional
arguments.
"""

from epygram import config


def add_arg_to_parser(parser, arg, **flychanges):
    """
    Wrapper to add one item *arg* of the following dictionaries to a *parser*.

    *flychanges* enable to change argument options on the fly.
    """

    arg[-1].update(flychanges.items())
    parser.add_argument(*arg[:-1], **arg[-1])

_defaults = {}
_defaults.update(config.__dict__)

from .domain_maker import d as domain_maker_args
from .extraction import d as extraction_args
from .fields import d as fields_args
from .files import d as files_args
from .graphical import d as graphical_args
from .misc import d as misc_args
from .operational import d as operational_args
from .output import d as output_args
from .runtime import d as runtime_args

