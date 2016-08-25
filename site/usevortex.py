#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains proxies to reach resources using Vortex.

Of course, this module need to have a proper Vortex installation !
"""

import uuid
import tempfile
import ftplib
import netrc
import os

from vortex import toolbox
import common
import olive

def set_defaults(**defaults):
    """
    Set defaults key/value pairs for get_resource().
    """
    if len(defaults) == 0:
        defaults = dict(model='arome',
                        cutoff='prod',
                        vapp='arome',
                        vconf='france',
                        geometry='franmgsp',
                        origin='hst')
    toolbox.defaults(**defaults)

def get_resources(getmode='epygram', uselocalcache=False, **description):
    """
    Get resources, given their description.
    
    *getmode*: 'epygram' return the epygram resources
               'locate' return the physical resolved location of the resource
               'exist' return the physical resolved location of the resource and its existence
               'fetch' fetches the resource in local, as filename *local*
               'vortex' return the vortex data handler of the resource
               'prestaging' puts a prestaging request for the resolved resources,
                            on the archive system
    *uselocalcache*: if True, store resources in a local cache (not
                     automatically cleaned, take care) defined either (and by
                     priority order) in $MTOOL_STEP_CACHE, $MTOOLDIR, $FTDIR,
                     $WORKDIR, $TMPDIR.

    Examples:
    
    - for the analysis of AROME-france from experiment 864G on 2015/08/15/00,
      description will look like:
      
      + experiment='864G',  # the experiment id
      + model='arome',
      + block='analysis',  # the OLIVE block
      + kind='analysis',  # the kind of resource
      + date='2015081500',  # the initial date and time
      + geometry='franmgsp',  # the name of the model domain
      + local='analysis_[experiment]',  # the local filename of the resource, once fetched.
        
    - for the model state at term 18h of ALADIN-reunion oper, production
      cutoff, on 2015/08/15/00,
      description will look like:
      
      + suite='oper',  # the suite // 'dble' = e-suite
      + kind='historic',  # model state
      + date='2015081500',  # the initial date and time
      + term=18,  # the forecast term
      + geometry='reunionsp',  # the name of the model domain
      + local='ICMSHALAD_[term]',  # the local filename of the resource, once fetched.
      + cutoff='prod',  # type of cutoff
      + vapp='aladin',  # type of application in operations namespace
      + vconf='reunion',  # name of config in operation namespace
      + model='aladin',  # name of the model
        
    - for the GRIB post-processed output on MASCA025 BDAP domain at terms
      22->24h of ALADIN-reunion oper, production cutoff, on 2015/08/15/00,
      description will look like:
      
      + suite='oper',  # the suite // 'dble' = e-suite
      + kind='gridpoint',  # model state
      + date='2015081500',  # the initial date and time
      + term=[22,23,24],  # the forecast terms
      + geometry='masca025',  # the name of the post-processing domain
      + local='[geometry]_[term].grb',  # the local filename of the resource, once fetched.
      + nativefmt='grib'  # to fetch the gribs and not the FA post-processed files
      + cutoff='prod',  # type of cutoff
      + vapp='aladin',  # type of application in operations namespace
      + vconf='reunion',  # name of config in operation namespace
      + model='aladin',  # name of the model
    """
    import common.util.usepygram

    # alias xp = experiment
    xp = description.pop('xp', None)
    if xp:
        description['experiment'] = xp

    # set namespace according to status xp/oper...
    if uselocalcache:
        _domain = 'multi'
    else:
        _domain = 'archive'
    if description.get('experiment', None) and 'namespace' not in description.keys():
        description['namespace'] = '.'.join(['olive', _domain, 'fr'])
    elif description.get('suite', None) in ('oper', 'dble') and 'namespace' not in description.keys():
        description['namespace'] = '.'.join(['oper', _domain, 'fr'])
    if not description.get('local', None):
        description['local'] = str(uuid.uuid4()) + "_[term]"
    if description.get('kind', None) in ('analysis', 'historic'):
        description['nativefmt'] = 'fa'
    elif description.get('nativefmt', None) == 'grib':
        description['kind'] = 'gridpoint'
    if description.get('model', None) and not description.get('vapp', None):
        description['vapp'] = description['model']

    # resolve resource description
    if getmode == 'check':
        check_desc = {}
        for k in description.keys():
            if isinstance(description[k], list):
                check_desc[k] = description[k][0]
            else:
                check_desc[k] = description[k]
            resolved = toolbox.rload(**check_desc)
    else:
        if getmode == 'prestaging':
            # force namespace archive for prestaging
            description['namespace'] = '.'.join([description['namespace'].split('.')[0],
                                                 'archive',
                                                 'fr'])
        resolved = toolbox.rload(**description)

    # and complete reaching resource according to getmode
    if getmode == 'locate':
        resources = [r.locate() for r in resolved]
    elif getmode == 'exist':
        resources = []
        for r in resolved:
            try:
                resources.append((r.locate(), bool(r.check())))
            except Exception:
                resources.append((r.locate(), False))
    elif getmode in ('epygram', 'fetch'):
        for r in resolved:
            ok = r.get()
            if not ok:
                raise IOError("fetch failed: " + r.locate())
        if getmode == 'epygram':
            resources = [r.contents.data for r in resolved]
        elif getmode == 'fetch':
            resources = [r.container.abspath for r in resolved]
    elif getmode == 'vortex':
        resources = resolved
    elif getmode == 'check':
        resources = resolved[0].complete
    elif getmode == 'prestaging':
        # connect to archive
        archive_name = 'hendrix'
        try:
            (_login, _, _passwd) = netrc.netrc().authenticators(archive_name)
        except TypeError:
            if netrc.netrc().authenticators(archive_name) is None:
                raise IOError("host " + archive_name + " is unknown in .netrc")
        ftp = ftplib.FTP(archive_name)
        ftp.login(_login, _passwd)
        # build request
        tmpdir = '/dev/shm'
        stagedir = '/DemandeMig/ChargeEnEspaceRapide'
        if 'mail' in description.keys():
            request = ["#MAIL=" + description['mail'] + '\n', ]
        else:
            request = []
        request += [r.locate().split('hendrix.meteo.fr:')[1] + '\n' for r in resolved]
        with open(tempfile.mkstemp(prefix='.'.join([_login, 'staging_request', '']),
                                   dir=tmpdir,
                                   suffix='.MIG')[1], 'w') as f:
            f.writelines(request)
            staging_request = os.path.basename(f.name)
            staging_request_fullname = f.name
        # transfer request to archive
        with open(staging_request_fullname, 'rb') as f:
            ftp.cwd(stagedir)
            ftp.storbinary('STOR ' + staging_request, f)
        os.remove(staging_request_fullname)
        ftp.quit()
        # send back request identifier
        resources = ['/'.join([stagedir, staging_request])]
    else:
        raise ValueError('*getmode* unknown: ' + getmode)

    return resources

#################
### SHORTCUTS ###
#################
# AROME

def get_gribfc_arome_oper(date, term, geometry='FRANGP0025', **others):
    """
    Proxy for AROME oper GRIBs.
    """

    _others = dict(suite='oper',
                   model='arome',
                   cutoff='prod',
                   vapp='arome',
                   vconf='france',
                   origin='hst',
                   kind='gridpoint',
                   nativefmt='grib')
    _others.update(others)
    return get_resources(date=date, term=term, geometry=geometry,
                         **_others)

def get_histfc_arome_oper(date, term, **others):
    """
    Proxy for AROME oper historic FAs.
    """

    _others = dict(suite='oper',
                   model='arome',
                   cutoff='prod',
                   vapp='arome',
                   vconf='france',
                   geometry='franmgsp',
                   kind='historic',
                   nativefmt='fa')
    _others.update(others)
    return get_resources(date=date, term=term, **_others)

def get_gribfc_arome_xp(xp, date, term, geometry='FRANGP0025',
                        **others):
    """
    Proxy for AROME experiment GRIBs.
    """

    _others = dict(model='arome',
                   cutoff='prod',
                   block='forecast',
                   vapp='arome',
                   vconf='france',
                   origin='hst',
                   kind='gridpoint',
                   nativefmt='grib')
    _others.update(others)
    return get_resources(experiment=xp, date=date, term=term, geometry=geometry,
                         **_others)

def get_histfc_arome_xp(xp, date, term, **others):
    """
    Proxy for AROME experiment historic FAs.
    """

    _others = dict(model='arome',
                   cutoff='prod',
                   block='forecast',
                   vapp='arome',
                   vconf='france',
                   geometry='franmgsp',
                   kind='historic',
                   nativefmt='fa')
    _others.update(others)
    return get_resources(experiment=xp, date=date, term=term,
                         **_others)

# ARPEGE ######################################################################

def get_gribfc_arpege_oper(date, term, geometry='FRANX01', **others):
    """
    Proxy for ARPEGE oper GRIBs.
    """

    _others = dict(suite='oper',
                   model='arpege',
                   cutoff='prod',
                   vapp='arpege',
                   vconf='france',
                   origin='hst',
                   kind='gridpoint',
                   nativefmt='grib')
    _others.update(others)
    return get_resources(date=date, term=term, geometry=geometry,
                         **_others)

def get_histfc_arpege_oper(date, term, **others):
    """
    Proxy for ARPEGE oper historic FAs.
    """

    _others = dict(suite='oper',
                   model='arpege',
                   cutoff='prod',
                   vapp='arpege',
                   vconf='france',
                   geometry='global',
                   kind='historic',
                   nativefmt='fa')
    _others.update(others)
    return get_resources(date=date, term=term,
                         **_others)

def get_gribfc_arpege_xp(xp, date, term, geometry='FRANX01', **others):
    """
    Proxy for ARPEGE experiment GRIBs.
    """

    _others = dict(model='arpege',
                   cutoff='prod',
                   block='forecast',
                   vapp='arpege',
                   vconf='france',
                   origin='hst',
                   kind='gridpoint',
                   nativefmt='grib')
    _others.update(others)
    return get_resources(experiment=xp, date=date, term=term, geometry=geometry,
                         **_others)

def get_histfc_arpege_xp(xp, date, term, **others):
    """
    Proxy for ARPEGE experiment historic FAs.
    """

    _others = dict(model='arpege',
                   cutoff='prod',
                   block='forecast',
                   vapp='arpege',
                   vconf='france',
                   geometry='global',
                   kind='historic',
                   nativefmt='fa')
    _others.update(others)
    return get_resources(experiment=xp, date=date, term=term, **_others)
