#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""EPyGrAM Command-Line Interface"""
import os
import argparse
import sys
import importlib
import epygram

commands = sorted([f[:-3] for f in os.listdir(os.path.dirname(__file__))
                   if f.endswith('.py') and f not in ('__init__.py',)])
epilog = 'End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')'
_description = __doc__


def main():
    cmd = sys.argv[1]
    if cmd == "-h":
        args = get_args()
    # this is a bit sketchy:
    sys.argv.pop(1)
    sys.argv[0] = sys.argv[0] + ' ' + cmd
    # import module and call main():
    module = importlib.import_module('.' + cmd, __name__)
    module.main()


def get_args():
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    parser.add_argument('command',
        help=' '.join(['Command to be executed: will call the "main" function from epygram.cli.<command> module.',
                       'Each command is then auto-documented: `epygram <command> -h`.']),
        choices=commands)
    return parser.parse_args()

