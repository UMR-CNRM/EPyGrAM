#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2016-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import os
import getpass
import sys
import socket
import web

import matplotlib
matplotlib.use("Agg")

#from . import web

import epygram
from epygram import epylog

from . import util

__all__ = ['main']
__authors__ = ['Ghislain Faure', ]
__contributors__ = ['Alexandre Mary', ]


###########
# Workdir #####################################################################
###########
# Vortex cache
location_of_vortex_cache = 'MTOOLDIR'
vortex_cache_dir = os.getenv(location_of_vortex_cache)
if not vortex_cache_dir:
    location_of_vortex_cache = 'FTDIR'
    vortex_cache_dir = os.getenv(location_of_vortex_cache)
    if not vortex_cache_dir:
        location_of_vortex_cache = 'WORKDIR'
        vortex_cache_dir = os.getenv(location_of_vortex_cache)
        if not vortex_cache_dir:
            location_of_vortex_cache = 'TMPDIR'
            vortex_cache_dir = os.getenv(location_of_vortex_cache)
if not vortex_cache_dir:
    raise ValueError('no rootdir has been defined for the Vortex cache: please export $MTOOLDIR.')
if location_of_vortex_cache == 'MTOOLDIR':
    vortex_cache = os.path.join(vortex_cache_dir, 'cache')
elif location_of_vortex_cache in ('FTDIR', 'WORKDIR', 'TMPDIR'):
    vortex_cache = os.path.join(vortex_cache_dir, 'mtool', 'cache')
if not os.path.exists(vortex_cache):
    os.makedirs(vortex_cache)

# Epyweb workdir for tmp files (basemap pickle, resources hardlinks and figures)
epyweb_workdir = os.path.join(vortex_cache, 'epyweb', getpass.getuser())
# epyweb_workdir = os.path.join(vortex_cache, 'epyweb', os.getlogin())
basemap_pickle_path = os.path.join(epyweb_workdir, 'basemap.pickle')

# Debug mode
all_fatal_exceptions = False

#############
# Web stuff ###################################################################
#############
shared_urls = ('/epyweb', 'Epyweb',
               '/getfieldsasjson', 'GetFieldsAsJSON',
               '/getminmax', 'GetMinMax',
               '/getdomain', 'GetDomain',
               '/getFile', 'GetFile',
               '/getCacheSize', 'GetCacheSize',
               '/getGeometries', 'GetGeometries',
               '/getPNG/(.+)', 'GetPNG',
               '/myplot', 'MyPlot',
               '/myplot_overlay', 'MyPlotOverlay',
               '/myplot_diff', 'MyPlotDiff',
               )

render = web.template.render('templates', base='base')

from .actions import *


########
# MAIN ########################################################################
########
def main(open_browser=False,
         port=8080,
         verbose=True):
    """
    Run the 'epyweb' local server.

    :param open_browser: if True, open a web browser tab with 'epyweb' interface
    :param port: the port to be used for the server can be specified.
    :param verbose: verbosity
    """
    epygram.init_env()
    epylog.setLevel('WARNING')
    if verbose:
        epylog.setLevel('INFO')

    if len(sys.argv) > 1:
        sys.argv = sys.argv[:1]  # workaround to a bug in web.py: command-line arguments can be other than port
    epyweb_url = 'http://' + socket.gethostname() + ':' + str(port) + '/epyweb'

    util.init_workdir(epyweb_workdir)
    # to avoid any path issues, "cd" to the web root. #FIXME: ? needed for the templates
    web_root = os.path.abspath(os.path.dirname(__file__))
    if os.getcwd() != web_root:
        os.chdir(web_root)

    print("=====================")
    print("EPYWEB is running within epygram version:", epygram.__version__)
    print("EPYWEB Interface =>", epyweb_url)
    print("EPYWEB Workdir   =>", epyweb_workdir)
    print("VORTEX Cache     =>", vortex_cache)
    print("(based on $" + location_of_vortex_cache + "=" + vortex_cache_dir + ")")
    if location_of_vortex_cache != 'MTOOLDIR':
        print("(the *vortex* cache location is accessible by priority order through:")
        print("$MTOOLDIR, $FTDIR, $WORKDIR, $TMPDIR")
    if location_of_vortex_cache == 'TMPDIR':
        epylog.warning(' '.join(['the use of $TMPDIR as rootdir for the Vortex',
                                 'cache is hazardous. You should define a',
                                 'better rootdir using $MTOOLDIR.']))
    print("To close the server: Ctrl-C")
    print("=====================")

    if open_browser:
        import threading
        t = threading.Thread(target=util.func_open_browser,
                             args=[epyweb_url],
                             kwargs={'delay':1.})
        t.start()  # FIXME: no join ?!!
    app = util.PortApplication(shared_urls, globals())
    try:
        app.run(port)
    except socket.error:
        print("!!! ERROR : EPYWEB is already running !!!")
        raise
    finally:
        util.clean_workdir(epyweb_workdir)
