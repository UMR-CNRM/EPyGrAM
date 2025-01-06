#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import os
import argparse
import sys
import importlib
import epygram


commands = [f.strip('.py') for f in os.listdir(os.path.dirname(__file__)) if f.endswith('.py') and f not in ('__init__.py')]


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
    parser = argparse.ArgumentParser(description='EPyGrAM main command-line interface.',
                                     epilog='End of help for: %(prog)s (EPyGrAM-' + epygram.__version__ + ')')
    parser.add_argument('command',
            help=' '.join(['Command to be executed: will call the "main" function from epygram.cli.<command> module.',
                           'Available commands are: {}.'.format(commands),
                           'Each command is then auto-documented with option -h.']))
    return parser.parse_args()

