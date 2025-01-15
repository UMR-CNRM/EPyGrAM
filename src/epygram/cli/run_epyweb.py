#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""A web interface for plotting with epygram and vortex."""

import argparse

import epygram
from epygram.extra import epyweb
from . import epilog
from .args_catalog import add_arg_to_parser, runtime_args

_description = __doc__


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        epylog.setLevel('INFO')
    else:
        epylog.setLevel('WARNING')
    epyweb.main(open_browser=args.open_browser,
                port=args.port,
                verbose=args.verbose)


def get_args():
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    parser.add_argument('-o',
                        action='store_true',
                        dest='open_browser',
                        help='opens a web browser tab with epyweb interface.',
                        default=False)
    parser.add_argument('-p',
                        dest='port',
                        type=int,
                        help='port to be used, default to 8080.',
                        default=8080)
    add_arg_to_parser(parser, runtime_args['verbose'])
    return parser.parse_args()

