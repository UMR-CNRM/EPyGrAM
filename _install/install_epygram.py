#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
import os
import shutil
import io
import argparse
from collections import defaultdict

epygram_repositories = {
    'cnrm':'/home/common/epygram',
    'bullx':'/home/gmap/mrpe/mary/public',
    'dsidev':'/soprano/home/marp999/epygram'}
userconfigs = defaultdict(lambda:'userconfig_no_arpifs4py.py',  # default
                          cnrm='userconfig_empty.py',
                          bullx='userconfig_empty.py')
linkname = 'src'
epygram_home = os.path.join(os.environ['HOME'], '.epygram')
profile = os.path.join(epygram_home, 'profile')

if any([h in os.environ.get('HOSTNAME', '') for h in
        ['beaufix', 'prolix']]):
    localhost = 'bullx'
elif any([h in os.environ.get('HOSTNAME', '') for h in
          ['alose', 'pagre', 'orphie', 'rason', 'guppy']]):
    localhost = 'dsidev'
else:
    localhost = 'cnrm'
epygram_repo = epygram_repositories.get(localhost,
                                        epygram_repositories['cnrm'])
userconfig = userconfigs[localhost]
_install_profile = localhost + '_profile'


def main(version='',
         fromdir=epygram_repo,
         update_epygram_profile=False,
         update_bash_profile=False):
    """
    Link to **version** from **fromdir**, copy adequate profile and
    make .bash_profile source it.
    """
    if version != '':
        if version.startswith('EPyGrAM'):
            version = version[7:]
        elif not version.startswith('-'):
            version = '-' + version
    if not os.path.exists(epygram_home):
        os.mkdir(epygram_home)
    os.chdir(epygram_home)
    if os.path.islink(linkname):
        os.remove(linkname)
    os.symlink(os.path.join(fromdir, 'EPyGrAM' + version),
               linkname)
    if update_epygram_profile or not os.path.exists(profile):
        shutil.copy(os.path.join(linkname, '_install', _install_profile),
                    profile)
    if not os.path.exists('userconfig.py'):
        shutil.copy(os.path.join(linkname, '_install', userconfig),
                    'userconfig.py')
    ufdf = 'user_Field_Dict_FA.csv'
    if not os.path.exists(ufdf):
        shutil.copy(os.path.join(linkname, '_install', ufdf), ufdf)
    if update_bash_profile:
        with io.open(os.path.join(os.environ['HOME'], '.bash_profile'), 'a') as pf:
            pf.write('\n#\n')
            pf.write('# epygram & vortex environment\n')
            pf.write('if [ -f {} ]; then\n'.format(profile))
            pf.write('  . {}\n'.format(profile))
            pf.write('fi\n')
    print("Local installation complete in: {}".format(epygram_home))
    print("To use it, restart session or source {}".format(profile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Helper to install or update EPyGrAM @ CNRM')
    parser.add_argument('-v', '--version_to_be_linked',
                        help=' '.join(['version to be linked, within available'
                                       'versions in --from directory']),
                        required=False,
                        default='')
    parser.add_argument('-f', '--from',
                        help=' '.join(['absolute path to directory in which to',
                                       'find the required version, defaults to',
                                       '{}']).format(epygram_repo),
                        default=epygram_repo,
                        dest='fromdir')
    parser.add_argument('-e', '--epygram_profile',
                        help='update epygram profile {}'.format(profile),
                        action='store_true',
                        default=False)
    parser.add_argument('-b', '--bash_profile',
                        help=' '.join(['update bash_profile, making it source'
                                       '{}']).format(profile),
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    main(args.version_to_be_linked,
         fromdir=args.fromdir,
         update_epygram_profile=args.epygram_profile,
         update_bash_profile=args.bash_profile)
