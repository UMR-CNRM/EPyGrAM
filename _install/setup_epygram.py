#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals

import os
import shutil
import io
import argparse
from collections import defaultdict
import sys
import site

EPyGrAM = 'EPyGrAM'


class _LocalHost(object):
    """Object to handle locations specific to host platform."""

    @property
    def user_site(self):
        """Python user site, where packages are looked for"""
        user_site = site.getusersitepackages()
        if isinstance(user_site, list):
            user_site = user_site[0]
        if not os.path.exists(user_site):
            os.makedirs(user_site)
        return user_site

    @property
    def kind(self):
        """Kind == category of host."""
        return self._kind()

    @staticmethod
    def _kind():
        """Kind == category of host."""
        hostname = os.environ.get('HOSTNAME', '')
        if any([hostname.startswith(h) for h in
                ['beaufix', 'prolix']]):
            host_kind = 'bullx'
        elif any([hostname.startswith(h) for h in
                  ['belenos', 'taranis']]):
            host_kind = 'atos_sequana'
        elif any([h in hostname for h in
                  ['alose', 'pagre', 'orphie', 'rason', 'guppy'] +
                  ['sotrtm{}-sidev'.format(n) for n in range(31, 35)]]):
            host_kind = 'dsidev'
        elif any([hostname.startswith(h) for h in
                  ['cca', 'ccb']]):
            host_kind = 'ecmwf_cc'
        elif any([hostname.startswith(h) for h in
                  ['ecgb', ]]):
            host_kind = 'ecgate'
        else:
            host_kind = 'default'
        return host_kind

    @property
    def eccodes_python2_installdir(self):
        """Location of installation of eccodes for python2 on host platform."""
        install_dirs = {
            'default':'/home/common/epygram/ext/eccodes/lib64/python2.7/site-packages',
            'atos_sequana':'/opt/softs/libraries/ICC_2018.5.274/eccodes-2.17.0/lib/python2.7/site-packages',
            'bullx':os.path.join('/opt/softs/libraries/ICC16.1.150/',
                                 'eccodes-2.7.0-b80884e7ca77a8f8ead5b4b1a2bd9011448b961e',
                                 '/lib/python2.7/site-packages'),
            'dsidev':'/usr/local/sopra/eccodes/lib64/python2.6/site-packages'}
        if self.kind in install_dirs:
            return install_dirs[self.kind]
        else:
            raise NotImplementedError("eccodes/python2 install dir is unknown for this kind of host: {}".format(
                self.kind))

    @property
    def vortex_installdir(self):
        """Location of installation of Vortex on host platform."""
        install_dirs = {
            'default': '/home/common/sync/vortex/vortex',
            'bullx': '/home/mf/dp/marp/verolive/vortex/vortex-olive',
            'atos_sequana': '/home/mf/dp/marp/verolive/vortex/vortex-olive',
            'dsidev': '/soprano/home/marp999/vortex/vortex',
            'ecmwf_cc': '/home/ms/fr/sos/vortex/vortex-olive',
            'ecgate': '/home/ms/fr/sos/vortex/vortex-olive'}
        if self.kind in install_dirs:
            return install_dirs[self.kind]
        else:
            raise NotImplementedError("Vortex install dir is unknown for this kind of host: {}".format(self.kind))

    @property
    def general_profile(self):
        """Path to .bash_profile or equivalent, depending on host platform."""
        profiles = defaultdict(lambda: '.bash_profile',
                               ecmwf_cc='.user_profile',
                               ecgate='.user_profile')
        return os.path.join(os.environ['HOME'], profiles[self.kind])

    def include_in_general_profile(self, header, profile_path):
        """Utility to include sourcing of a **profile_path**, signaled with a **header**, into the general_profile."""
        with io.open(self.general_profile, 'a') as pf:
            pf.write('\n#\n')
            pf.write('# ' + header + '\n')
            pf.write('if [ -f {} ]; then\n'.format(profile_path))
            pf.write('  . {}\n'.format(profile_path))
            pf.write('fi\n')
        print("* {} has been included in {}".format(header, self.general_profile))


class EpygramInstaller(object):
    """Helper to install EPyGrAM."""

    #: epygram user customization directory
    epygram_home = os.path.join(os.environ['HOME'], '.epygram')

    def __init__(self):
        self.host = _LocalHost()
        # this file is /path/to/EPyGrAM/_install/setup_epygram.py
        self.this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # user custom dir: ~/.epygram
        if not os.path.exists(self.epygram_home):
            os.mkdir(self.epygram_home)
        self.in_site = os.path.join(self.host.user_site, 'EPyGrAM')
        # profile
        self.profile = os.path.join(self.epygram_home, 'profile')

    def link_sources_in_site_packages(self, uninstall_before=False):
        """
        Link sources in site-packages python directory.
        If **uninstall_before**, uninstall any preexisting EPyGrAM therein.
        """
        os.chdir(self.host.user_site)
        # uninstall if present
        if os.path.exists(EPyGrAM):
            if uninstall_before:
                self.uninstall_from_site_packages()
            else:
                raise RuntimeError("EPyGrAM is already installed in {} : please uninstall (remove).".format(
                    self.host.user_site))
        # link
        os.symlink(self.this_dir, EPyGrAM)
        pth = EPyGrAM + '.pth'
        shutil.copy(os.path.join(EPyGrAM, '_install', pth),
                    pth)
        print("* EPyGrAM has been successfully linked in: {}".format(self.host.user_site))

    def uninstall_from_site_packages(self):
        """Uninstall an EPyGrAM from site-packages."""
        os.chdir(self.host.user_site)
        if os.path.islink(EPyGrAM):
            os.remove(EPyGrAM)
        elif os.path.isdir(EPyGrAM):
            shutil.rmtree(EPyGrAM)
        pth = EPyGrAM + '.pth'
        if os.path.exists(pth):
            os.remove(pth)

    def init_user_customization(self):
        """Initialize epygram user customization directory."""
        os.chdir(self.epygram_home)
        # customizable samples
        for customizable in ('sfxflddesc_mod.F90', 'gribapi.def.0', 'userconfig.py'):
            if not os.path.exists(customizable):
                source = os.path.join(self.this_dir, '_install', customizable)
                if os.path.isdir(source):
                    shutil.copytree(source, customizable)
                else:
                    shutil.copy(source, customizable)
        print("* EPyGrAM user customization has been setup under: {}".format(self.epygram_home))

    def set_profile(self):
        """Set epygram profile file, to be sourced."""
        os.chdir(self.epygram_home)
        install_profile = self.host.kind + '_profile'
        with open(os.path.join(self.this_dir, '_install', install_profile), 'r') as p:
            lines = p.readlines()
        lines.insert(0, "export EPYGRAM_DIR={}\n\n".format(self.in_site))
        lines.append("\n# PATH for apptools\n")
        lines.append('export PATH=$PATH:$EPYGRAM_DIR/{}\n'.format('apptools'))
        with open(self.profile, 'w') as p:
            p.writelines(lines)
        print("* EPyGrAM profile has been written in {} : source it or include in ~/.bash_profile".format(self.profile))

    def include_profile_in_general_profile(self):
        """Include sourcing of epygram profile in the general profile."""
        self.host.include_in_general_profile('EPyGrAM profile', self.profile)


class VortexInstaller(object):

    vortex_home = os.path.join(os.environ['HOME'], '.vortexrc')

    def __init__(self):
        self.host = _LocalHost()
        if not os.path.exists(self.vortex_home):
            os.mkdir(self.vortex_home)
        self.profile = os.path.join(self.vortex_home, 'profile')

    def set_profile(self, vortex_install_dir=None):
        """Set Vortex profile file, to be sourced."""
        os.chdir(self.vortex_home)
        if vortex_install_dir is None:
            vortex_install_dir = self.host.vortex_installdir
        lines = ["VORTEX_INSTALL_DIR={}".format(vortex_install_dir),
                 "export PYTHONPATH=$VORTEX_INSTALL_DIR:$VORTEX_INSTALL_DIR/src:$VORTEX_INSTALL_DIR/site:$PYTHONPATH",
                 "export PATH=$VORTEX_INSTALL_DIR/bin:$VORTEX_INSTALL_DIR/site/arpifs_listings/bin:$PATH"]
        with open(self.profile, 'w') as p:
            p.writelines([line + '\n' for line in lines])
        print("* Vortex profile has been written in {} : source it or include in ~/.bash_profile".format(self.profile))

    def include_profile_in_general_profile(self):
        """Include sourcing of Vortex profile in the general profile."""
        self.host.include_in_general_profile('Vortex (& site packages) profile', self.profile)


class EccodesInstaller(object):
    """Helper to install eccodes for python2."""

    def __init__(self):
        self.host = _LocalHost()

    def setup(self, uninstall_before=False):
        """Setup eccodes for python2."""
        if sys.version_info.major == 2:
            self._link_eccodes_for_py2_in_site_packages(uninstall_before=uninstall_before)
        else:
            print("For Python3, please install eccodes using pip: pip install --user eccodes")

    def _link_eccodes_for_py2_in_site_packages(self, uninstall_before=False):
        """
        Link eccodes/python2 in site-packages python directory.
        If **uninstall_before**, uninstall any preexisting eccodes therein.
        """
        os.chdir(self.host.user_site)
        for lib in ('gribapi', 'eccodes'):
            if os.path.exists(lib):
                if uninstall_before:
                    self.uninstall_from_site_packages(lib)
                else:
                    raise Warning("(eccodes already installed in : {}".format(self.host.user_site))
            if not os.path.exists(lib):
                os.symlink(os.path.join(self.host.eccodes_python2_installdir, lib), lib)
                if lib == 'eccodes':
                    print("* eccodes/python2 has been successfully linked in {}".format(self.host.user_site))

    def uninstall_from_site_packages(self, lib):
        """Uninstall an eccodes from site-packages."""
        os.chdir(self.host.user_site)
        if os.path.islink(lib):
            os.remove(lib)
        elif os.path.isdir(lib):
            shutil.rmtree(lib)


def main(include_profile_in_general_profile=False,
         link_eccodes_for_py2_in_site_packages=False,
         set_vortex_profile=False,
         vortex_install_dir=None,
         uninstall_before=False):
    """
    Setup EPyGrAM and optionally Vortex and eccodes/python2.

    :param include_profile_in_general_profile: update bash_profile to make it source EPyGrAM profile
        (& Vortex if requested)
    :param link_eccodes_for_py2_in_site_packages: link eccodes/python2 interface
    :param set_vortex_profile: create or overwrite a Vortex profile
    :param vortex_install_dir: Vortex install directory to be used
    :param uninstall_before: uninstall packages if already present, before reinstalling
    """
    installer = EpygramInstaller()
    installer.link_sources_in_site_packages(uninstall_before=uninstall_before)
    installer.init_user_customization()
    installer.set_profile()
    if include_profile_in_general_profile:
        installer.include_profile_in_general_profile()
    # eccodes
    if link_eccodes_for_py2_in_site_packages:
        eccodes_installer = EccodesInstaller()
        eccodes_installer.setup(uninstall_before=uninstall_before)
    # Vortex
    if set_vortex_profile:
        vortex_installer = VortexInstaller()
        vortex_installer.set_profile(vortex_install_dir)
        if include_profile_in_general_profile:
            vortex_installer.include_profile_in_general_profile()
    else:
        print("Vortex site-packages 'footprints' & 'bronx' are required for EPyGrAM:",
              "assuming pre-installed (cf. options).")


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Helper to install or update EPyGrAM')
    parser.add_argument('-b', '--bash_profile',
                        help='update bash_profile to make it source EPyGrAM profile (& Vortex if requested)',
                        action='store_true',
                        default=False)
    parser.add_argument('-v', '--vortex_profile',
                        help='create or overwrite a Vortex profile in {}'.format(VortexInstaller.vortex_home),
                        action='store_true',
                        default=False)
    parser.add_argument('--vortex_install_dir',
                        help='Vortex install directory to be used',
                        default=None)
    parser.add_argument('-e', '--link_eccodes',
                        help='link eccodes/python2 interface (for python3, use: pip3 install --user eccodes)',
                        action='store_true',
                        default=False)
    parser.add_argument('-u', '--uninstall',
                        help='uninstall packages if already present, before reinstalling',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    main(include_profile_in_general_profile=args.bash_profile,
         link_eccodes_for_py2_in_site_packages=args.link_eccodes,
         set_vortex_profile=args.vortex_profile,
         vortex_install_dir=args.vortex_install_dir,
         uninstall_before=args.uninstall)
