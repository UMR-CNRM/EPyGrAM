#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
List the external (potential) dependancies of EPyGrAM,
and append the source/library/epygram.rst with list.

Some of these dependancies may be dynamic, only needed for actual call to some
methods or functionalities.
"""
from __future__ import print_function, absolute_import, unicode_literals, division
import glob
import epygram
import os
import io

inner_packages = [
    'arpifs4py',
    'ctypesForFortran',
    'epygram',
    'epyweb',
    'grib_utilities',
    'pyexttiff',
    'userconfig',
    'usevortex',
    'vgrid',]

vortex_site_packages = [
    'bronx',
    'footprints',
    'taylorism',]

stdlib_packages = [
    '__future__',
    'argparse',
    'collections',
    'contextlib',
    'copy',
    'ctypes',
    '_ctypes',
    'datetime',
    'datetime',
    'ftplib',
    'glob',
    'hashlib',
    'imp',
    'importlib',
    'inspect',
    'io',
    'itertools',
    'json',
    'math',
    'mmap',
    'multiprocessing',
    'os',
    'os',
    'platform',
    're',
    'resource',
    'shutil',
    'socket',
    'string',
    'subprocess',
    'sys',
    'tempfile',
    'threading',
    'time',
    'uuid',
    'webbrowser',]


def look_for_packages(lines):
    """Filter package name from import lines."""
    # filter out comments
    pckgs = [l.strip() for l in lines if not l.strip().startswith('#')]
    # filter out sub-modules
    pckgs = [l for l in pckgs if not l.startswith('from .')]
    # keep only real imports
    pckgs = [l for l in pckgs if l.startswith('from ') or l.startswith('import ')]
    # keep only package
    # "from package.module import sthg" // "import package.shtg"
    pckgs = [l.split(' ')[1].split('.')[0] for l in pckgs if ' ' in l]
    # filter stdlib and inner packages/modules
    pckgs = [p for p in pckgs if p not in stdlib_packages]
    pckgs = [p for p in pckgs if p not in inner_packages]
    pckgs = [p for p in pckgs if p not in vortex_site_packages]
    # unique and sorted
    pckgs = sorted(set(pckgs))
    return pckgs


def look_for_imports(filepattern, exclude=[]):
    """
    Find import lines in a series of files to be glob-expanded from filepattern.
    """
    found = []
    for f in glob.glob(filepattern):
        if len([e for e in exclude if e in f]) == 0:
            with io.open(f, 'r') as f:
                for l in f:
                    if 'import ' in l:
                        found.append(l.rstrip('\n'))
    return found


def main(printout=True, write_to_file=True):
    all_dependancies = {}
    install_dir = os.path.dirname(os.path.dirname(os.path.abspath(epygram.__file__)))
    # epygram
    subdir = 'epygram'
    found = []
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*.py')))
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*/*.py')))
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*/*/*.py')))
    packages = look_for_packages(found)
    all_dependancies[subdir] = packages
    # site
    subdir = 'site'
    found = []
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*.py'),
                                  exclude=['usevortex.py',]))
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*/*.py')))
    packages = look_for_packages(found)
    all_dependancies[subdir] = packages
    # apptools
    subdir = 'apptools'
    found = []
    found.extend(look_for_imports(os.path.join(install_dir, subdir, '*.py')))
    packages = look_for_packages(found)
    all_dependancies[subdir] = packages

    if printout:
        for k, packages in all_dependancies.items():
            print('### {} ###'.format(k))
            for p in packages:
                print(p)
    if write_to_file:
        with open(os.path.join(install_dir, 'dependancies.rst'), 'w') as f:
            f.write('Dependancies\n')
            f.write('============\n\n')
            f.write('.. _dependancies:\n\n')
            f.write('Some of which may be dynamically needed only, ' +
                    'i.e. needed at actual use of some methods.\n\n')
            alldep = []
            for packages in all_dependancies.values():
                alldep.extend(packages)
            for p in sorted(set(alldep)):
                f.write('- :mod:`' + p + '`\n')
            f.write('\n')
            f.write('VORTEX site-packages (distributed with epygram)\n')
            f.write('-----------------------------------------------\n\n')
            for p in sorted(vortex_site_packages):
                f.write('- :mod:`' + p + '`\n')


if __name__ == '__main__':
    main()