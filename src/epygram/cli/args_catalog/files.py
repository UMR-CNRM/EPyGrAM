#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Arguments dealing with files
"""

d = {
    'principal_file':[
        'filename',
        dict(type=str,
             help='name of the file to be processed.')],
    'several_files':[
        'filenames',
        dict(type=str,
             nargs='+',
             help='names of the files to be processed.')],
    'file_to_refer_in_diff':[
        '-d', '--diff',
        dict(type=str,
             dest='refname',
             help='name of the 2nd (reference) file to be processed, to which\
                   comparison is done.')],
    'file_to_refer_in_diffonly':[
        '-D', '--diffonly',
        dict(type=str,
             dest='Drefname',
             help='same as -d/--diff, but only fields difference is processed.')],
    'source_file':[
        '-s', '--source',
        dict(type=str,
             dest='refname',
             help='name of the 2nd file from which to extract the fields.',
             required=True)],
    'replace_by_diff':[
        '-d', '--diff',
        dict(action='store_const',
             const='diff',
             dest='replace_op',
             help='replaces the fields by the difference fields\
                   (filename - source).',
             default=False)],
    'replace_by_reversediff':[
        '-r', '--reversediff',
        dict(action='store_const',
             const='reversediff',
             dest='replace_op',
             help='replaces the fields by the reverse difference fields\
                   (source - filename).',
             default=False)],
    'replace_by_addition':[
        '-a', '--add',
        dict(action='store_const',
             const='add',
             dest='replace_op',
             help='replaces the fields by the addition fields\
                   (filename + source).',
             default=False)],
    'replace_by_product':[
        '-m', '--multiply',
        dict(action='store_const',
             const='multiply',
             dest='replace_op',
             help='replaces the fields by the product fields\
                   (filename * source).',
             default=False)],
    'in_place':[
        '-i', '--in_place',
        dict(action='store_true',
             help='the operation on fields is done "in place" on the file,\
                   not on a new file.',
             default=False)],
                    }

