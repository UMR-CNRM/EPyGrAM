#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
import os
import shutil
import io
import re
import argparse
from collections import defaultdict
import sys
import site

EPyGrAM = 'EPyGrAM'

epygram_repositories = {
    'cnrm':'/home/common/epygram/public/EPyGrAM',
    'bullx':'/home/gmap/mrpe/mary/public/EPyGrAM',
    'bull_sequana':'/home/gmap/mrpe/mary/public/EPyGrAM',
    'dsidev':'/soprano/home/marp999/epygram/EPyGrAM',
    'ecmwf_cc':'/home/ms/fr/rm9/public',  # TODO: move /EPyGrAM
    'ecgate':'/home/ms/fr/rm9/public',  # TODO: move /EPyGrAM
    }
vortex_repositories = {
    'cnrm':'/home/common/sync/vortex',
    'bullx':'/home/mf/dp/marp/verolive/vortex',
    'bull_sequana':'/home/mf/dp/marp/verolive/vortex',
    'dsidev':'/soprano/home/marp999/vortex',
    'ecmwf_cc':'/home/ms/fr/sos/vortex',
    'ecgate':'/home/ms/fr/sos/vortex'}
py2_eccodes_installdir = {
    'cnrm':'/home/common/epygram/ext/eccodes/lib64/python2.7/site-packages',
    'bullx':'/opt/softs/libraries/ICC16.1.150/eccodes-2.7.0-b80884e7ca77a8f8ead5b4b1a2bd9011448b961e/lib/python2.7/site-packages',
    'dsidev':'/usr/local/sopra/eccodes/lib64/python2.6/site-packages'}
userconfigs = defaultdict(lambda:'userconfig_no_arpifs4py.py',  # default
                          cnrm='userconfig_empty.py',
                          bullx='userconfig_empty.py',
                          bull_sequana='userconfig_empty.py',
                          ecmwf_cc='userconfig_empty.py')
profiles = defaultdict(lambda:'.bash_profile',
                       ecmwf_cc='.user_profile',
                       ecgate='.user_profile')
VORTEX_LINKNAME = 'vortex'
epygram_home = os.path.join(os.environ['HOME'], '.epygram')
usersite = site.getusersitepackages()
if isinstance(usersite, list):
    usersite = usersite[0]
profile = os.path.join(epygram_home, 'profile')

hostname = os.environ.get('HOSTNAME', '')
if any([hostname.startswith(h) for h in
        ['beaufix', 'prolix']]):
    localhost = 'bullx'
elif any([hostname.startswith(h) for h in
          ['epona', 'belenos', 'taranis']]):
    localhost = 'bull_sequana'
elif any([h in hostname for h in
          ['alose', 'pagre', 'orphie', 'rason', 'guppy'] +
          ['sotrtm{}-sidev'.format(n) for n in range(31, 35)]]):
    localhost = 'dsidev'
elif any([hostname.startswith(h) for h in
          ['cca', 'ccb']]):
    localhost = 'ecmwf_cc'
elif any([hostname.startswith(h) for h in
          ['ecgb',]]):
    localhost = 'ecgate'
else:
    localhost = 'cnrm'

epygram_repo = epygram_repositories.get(localhost,
                                        epygram_repositories['cnrm'])
vortex_repo = vortex_repositories.get(localhost,
                                      vortex_repositories['cnrm'])
userconfig = userconfigs[localhost]
_install_profile = localhost + '_profile'
_vortex_install_profile = 'vortex_profile'


def _version_():
    abspath = os.path.abspath(__file__)
    version = abspath.split(os.path.sep)[-3]
    return version


def main(version='stable',
         fromdir=epygram_repo,
         update_epygram_profile=False,
         update_bash_profile=False,
         install_vortex=True,
         vortex_version='olive',
         vortex_from=vortex_repo,
         link_eccodes=False):
    """
    Link to **version** from **fromdir**, copy adequate profile and
    make .bash_profile source it.
    
    :param version: must be a label (e.g. 'dev') or an actual version number
        (e.g. '1.4.5')
    """
    # EPyGrAM version
    # link sources in site-packages
    if not os.path.exists(usersite):
        os.makedirs(usersite)
    os.chdir(usersite)
    if os.path.islink(EPyGrAM):
        os.remove(EPyGrAM)
    os.symlink(os.path.join(fromdir, version),
               EPyGrAM)
    epygram_in_site = os.path.join(usersite, EPyGrAM)
    shutil.copy(os.path.join(EPyGrAM, '_install', EPyGrAM + '.pth'),
                EPyGrAM + '.pth')
    print("EPyGrAM has been successfully installed under: {}".format(usersite))
    # epygram_home local directory
    if not os.path.exists(epygram_home):
        os.mkdir(epygram_home)
    os.chdir(epygram_home)
    # vortex (if needed)
    if install_vortex:
        if vortex_version != '':
            if not vortex_version.startswith('-'):
                vortex_version = '-' + vortex_version
        if os.path.islink(VORTEX_LINKNAME):
            os.remove(VORTEX_LINKNAME)
        os.symlink(os.path.join(vortex_from, 'vortex' + vortex_version),
                   VORTEX_LINKNAME)
    # profile
    if update_epygram_profile or not os.path.exists(profile):
        with open(os.path.join(epygram_in_site, '_install', _install_profile), 'r') as p:
            lines = p.readlines()
        lines.insert(0, "export EPYGRAM_DIR={}\n\n".format(epygram_in_site))
        lines.append("\n# EPyGrAM apptools\n")
        lines.append('export PATH=$PATH:$EPYGRAM_DIR/{}\n'.format('apptools'))
        if install_vortex:
            with open(os.path.join(epygram_in_site, '_install', _vortex_install_profile), 'r') as p:
                lines.extend(p.readlines())
        with open(profile, 'w') as p:
            for l in lines:
                p.write(l)
    # user customization
    if not os.path.exists('userconfig.py'):
        shutil.copy(os.path.join(epygram_in_site, '_install', userconfig),
                    'userconfig.py')
    for example in ('sfxflddesc_mod.F90', 'gribapi.def.0'):
        if not os.path.exists(example):
            source = os.path.join(epygram_in_site, '_install', example)
            if os.path.isdir(source):
                shutil.copytree(source, example)
            else:
                shutil.copy(source, example)
    print("EPyGrAM user customization has been setup under: {}".format(epygram_home))
    # bash_profile
    if update_bash_profile:
        with io.open(os.path.join(os.environ['HOME'], profiles[localhost]), 'a') as pf:
            pf.write('\n#\n')
            pf.write('# epygram & vortex environment\n')
            pf.write('if [ -f {} ]; then\n'.format(profile))
            pf.write('  . {}\n'.format(profile))
            pf.write('fi\n')
    # eccodes
    if link_eccodes and sys.version_info.major == 2:
        targetdir = py2_eccodes_installdir.get(localhost, None)
        if targetdir is None:
            raise NotImplementedError("eccodes linking on this kind of platform: {}".format(localhost))
        for lib in ('eccodes', 'gribapi'):
            link = os.path.join(usersite, lib)
            if not os.path.exists(link):
                os.symlink(os.path.join(targetdir, lib), link)
            else:
                print("!!! Link eccodes failed: link already exists: {}".format(link))
    print("To use EPyGrAM, restart session (if option -b) or source {}".format(profile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Helper to install or update EPyGrAM @ CNRM')
    parser.add_argument('-v', '--version_to_be_linked',
                        help=' '.join(['version to be linked, within available'
                                       'versions in --from directory']),
                        required=False,
                        default=_version_())
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
                        help=' '.join(['update bash_profile, making it source',
                                       '{}']).format(profile),
                        action='store_true',
                        default=False)
    parser.add_argument('--nv', '--no_vortex',
                        help='do not install vortex (supposedly you already have it your own way)',
                        action='store_false',
                        dest='install_vortex',
                        default=True)
    parser.add_argument('--vv', '--vortex_version',
                        help=' '.join(['Vortex version to be linked']),
                        dest='vortex_version',
                        required=False,
                        default='olive')
    parser.add_argument('--link_eccodes',
                        help='link eccodes python2 interface installation in .local (CNRM workstations and beaufix/prolix only)',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    main(args.version_to_be_linked,
         fromdir=args.fromdir,
         update_epygram_profile=args.epygram_profile,
         update_bash_profile=args.bash_profile,
         install_vortex=args.install_vortex,
         vortex_version=args.vortex_version,
         link_eccodes=args.link_eccodes)
