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
import os
import datetime
import copy
from contextlib import contextmanager

import footprints
from bronx.system.mf import prestage
from bronx.stdtypes.date import daterange
import taylorism

import vortex
from vortex import toolbox
import common  # @UnusedImport : for footprints to load classes
import olive  # @UnusedImport : for footprints to load classes


def list_vortex_geometries():
    """Proxy for epyweb not to be vortex directly dependant."""
    return list(vortex.data.geometries.keys())


def set_defaults(**defaults):
    """Set defaults key/value pairs for get_resource()."""
    if len(defaults) == 0:
        defaults = dict(model='arome',
                        cutoff='prod',
                        vapp='arome',
                        vconf='france',
                        geometry='franmgsp',
                        origin='hst')
    toolbox.defaults(**defaults)


@contextmanager
def quiet_get(loggers=['vortex.data.stores',
                       'vortex.data.handlers',
                       'vortex.tools.net']):
    """
    Shut off some loggers (set level to ERROR) while executing action,
    then restore initial level.

    Example of use::

        with quiet_get():
            get_resources(...)
    """
    old_levels = {l:footprints.loggers.getLogger(l).getEffectiveLevel()
                  for l in loggers}
    for l in loggers:
        footprints.loggers.getLogger(l).setLevel('ERROR')
    yield
    for k,v in old_levels.items():
        footprints.loggers.getLogger(k).setLevel(v)


def get_resources(getmode='epygram',
                  uselocalcache=False,
                  **description):
    """
    Get resources, given their description.

    :param getmode: controls the action and returned object:\n
                    - 'check' return True only if the description is complete
                    - 'epygram' return the epygram resources
                    - 'locate' return the physical resolved location of the
                      resource
                    - 'exist' return the physical resolved location of the
                      resource and its existence
                    - 'fetch' fetches the resource in local, as filename *local*
                    - 'vortex' return the vortex data handler of the resource
                    - 'prestaging' puts a prestaging request for the resolved
                      resources, on the archive system
    :param uselocalcache: if True, store resources in a local cache (not
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
      + origin='hst',  # origin of post-processed files: historic
      + vapp='aladin',  # type of application in operations namespace
      + vconf='reunion',  # name of config in operation namespace
      + model='aladin',  # name of the model
    """
    import common.util.usepygram  # @UnusedImport : to load Epygram FormatAdapters
    t = vortex.ticket()

    if uselocalcache:
        _domain = 'multi'
    else:
        _domain = 'archive'

    # completion of description...
    xp = description.pop('xp', None)
    if xp:
        description['experiment'] = xp
    if description.get('experiment', None) in ('oper', 'dble') and 'suite' not in description.keys():
        description['suite'] = description['experiment']  # for oper/dble, alias experiment to suite
    if description.get('experiment', None) and 'namespace' not in description.keys():
        description['namespace'] = '.'.join(['olive', _domain, 'fr'])  # set namespace if not given
    elif description.get('suite', None) in ('oper', 'dble') and 'namespace' not in description.keys():
        description['namespace'] = '.'.join(['oper', _domain, 'fr'])  # set namespace if not given
    if not description.get('local', None):
        description['shouldfly'] = True
    if description.get('kind', None) in ('analysis', 'historic'):
        description['nativefmt'] = 'fa'  # set nativefmt if recognized
    elif description.get('nativefmt', None) == 'grib':
        description['kind'] = 'gridpoint'  # set kind if recognized
    if description.get('vapp', None) and not description.get('model', None):
        description['model'] = description['vapp']  # set model (useless?) if not given
    if getmode == 'prestaging':
        # force namespace archive for prestaging
        description['namespace'] = '.'.join([description['namespace'].split('.')[0],
                                             'archive',
                                             'fr'])

    if getmode == 'check':
        # just check completion of the resource
        try:
            t.context.record_off()  # avoid overconsumption of memory
            toolbox.rload(**description)
            t.context.record_on()
            resources = True
        except (toolbox.VortexToolboxDescError, footprints.FootprintException):
            resources = False
    else:
        # resolve resource description: raise an error if the description is not complete
        t.context.record_off()  # avoid overconsumption of memory
        resolved = toolbox.rload(**description)
        t.context.record_on()
        # and complete reaching resource according to getmode
        if getmode == 'vortex':
            resources = resolved
        elif getmode == 'locate':
            resources = [r.locate() for r in resolved]
        elif getmode == 'exist':
            with t.sh.ftppool():  # unique ftp connection for the whole request
                resources = [(r.locate(), bool(r.check())) for r in resolved]
        elif getmode == 'prestaging':
            with t.sh.ftppool():  # unique ftp connection for the whole request
                resources = [prestage([r.locate().split('hendrix.meteo.fr:')[1]
                                       for r in resolved],
                                      description.get('mail', None))]
        elif getmode in ('fetch', 'epygram'):
            ok = []
            with t.sh.ftppool():  # unique ftp connection for the whole request
                ok = [r.get() for r in resolved]
            if all(ok):
                if getmode == 'fetch':
                    resources = [r.container.abspath for r in resolved]
                elif getmode == 'epygram':
                    resources = [r.contents.data for r in resolved]
            else:
                raise IOError("fetch failed, at least: " + resolved[ok.index(False)].locate())
        else:
            raise ValueError('*getmode* unknown: ' + getmode)

    return resources


default_nc_global_attributes = dict(
    institution='Météo France',
    history='Extracted on {} using epygram.'.format(datetime.datetime.now().date().isoformat()))


def extractor(vortex_description,
              coords,
              start_cutoff,
              end_cutoff,
              start_term,
              end_term,
              points_fields=None,
              profiles_FAfields=None,
              everycutoff_in_hours=24,
              everyterm_in_hours=1,
              # options
              vertical_coordinate=None,
              extraction_options=None,
              use_local_cache=False,
              prestaging=False,
              outdir=os.getcwd(),
              progressmode=None,
              nc_global_attributes=default_nc_global_attributes):
    """
    Extract temporal series of points (and profiles for FA historical files) to
    netCDF files.

    Series are either:
      - one series == one file per cutoff if **start_term* != **end_term**
      - one series == one file containing all cutoffs if **start_term* == **end_term**

    Shall be extended to series of profiles for other formats later on.

    :param vortex_description: vortex description of the resources,
                               incl. providers, excepted dates/terms
                               Example:
                               {'model':'arome',
                                'vapp':'arome',
                                'vconf':'3dvarfr',
                                'kind':'historic',
                                'suite':'oper',
                                'geometry':'franmgsp',
                                'cutoff':'prod',
                                #'namespace':'oper.archive.fr',
                                #'experiment':'86U8',
                                #'block':'forecast',
                                #'nativefmt':'fa',
                                #'member':8,}
    :param coords: either a list of lon/lat coordinates or a dict of,
                   in which case key would be a "name of the point",
                   used in output filename
    :param start_cutoff: first cutoff of the series (datetime.datetime instance)
    :param end_cutoff: last cutoff of the series
    :param start_term: first term of the series (datetime.datetime instance)
    :param end_term: last term of the series
    :param points_fields: dictionnary of fields to extract with options, e.g.:\n
                          {fid1:{'nc_name':nc_fid1,
                                 'nc_attributes':{'unit':...,
                                                  'comment':...,
                                                  },
                                 'decumulate':True,
                                 'operation':('-',273.15) or ('exp',),
                                 'reproject_wind_on_lonlat':False,
                                 },
                           fid2:...
                           }\n
                          - 'decumulate' can be absent==False, True,
                            or 'center' for centered decumulation
                            (cf. D3Field.decumulate for more details);
                          - 'operation' is done after decumulation, if both
                            present
                          - 'reproject_wind_on_lonlat: FA only: reproject winds
                            on lonlat axes (True by default)
    :param profiles_FAfields: like points_fields but for profiles,
                              but only for FA historic resources !
                              fids supposed to be 'S*PARAMETER'
    :param everycutoff_in_hours: step between 2 cutoffs
    :param everyterm_in_hours: step between 2 terms
    :param vertical_coordinate: vertical coordinate requested for profiles,
                                among ('pressure', 'height', 'altitude', None)
    :param extraction_options: options to FA.extractprofile() and
                               H2DField.extract_point()
                               concerning interpolation, as dict...
    :param use_local_cache: use a local Vortex cache, in which files are stored
    :param prestaging: triggers prestaging for the requested resources (Hendrix)
    :param outdir: directory in which to store output files
    :param progressmode: sets the verbosity of the task among
                         (None, 'percentage', 'verbose')
    :param nc_global_attributes: given as a dict. Defaults are above in
                                 module variable default_nc_global_attributes
    """
    import epygram
    from epygram.base import FieldValidity
    from epygram.fields import gimme_one_point
    vd = {'pressure':epygram.geometries.Pressure,
          'altitude':epygram.geometries.Altitude,
          'height':epygram.geometries.Height}
    points_fields = epygram.util.ifNone_emptydict(points_fields)
    profiles_FAfields = epygram.util.ifNone_emptydict(profiles_FAfields)
    extraction_options = epygram.util.ifNone_emptydict(extraction_options)

    # 1. Prepare stuff
    dt_cutoffs = daterange(start_cutoff, end_cutoff,
                           'PT{}H'.format(everycutoff_in_hours))
    cutoffs = [dt.isoformat(sep=str(' ')).replace('-', '').replace(' ', '')[:10]
               for dt in dt_cutoffs]
    terms = range(start_term, end_term + 1, everyterm_in_hours)
    if isinstance(coords, list):
        coords = {'{}E,{}N'.format(*c):c for c in coords}
    metadata = {pf.get('nc_name', str(f)):pf.get('nc_attributes', {})
                for f, pf in points_fields.items()}
    metadata.update({pf.get('nc_name', str(f)):pf.get('nc_attributes', {})
                     for f, pf in profiles_FAfields.items()})
    decumulate = {pf.get('nc_name', str(f)):pf.get('decumulate', False)
                  for f, pf in points_fields.items()}
    decumulate.update({pf.get('nc_name', str(f)):pf.get('decumulate', False)
                       for f, pf in profiles_FAfields.items()})
    operations = {pf.get('nc_name', str(f)):pf.get('operation', False)
                  for f, pf in points_fields.items()}
    operations.update({pf.get('nc_name', str(f)):pf.get('operation', False)
                       for f, pf in profiles_FAfields.items()})
    vortex_description.update({'getmode':'epygram',
                               'uselocalcache':use_local_cache})
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # 2. definition of functions
    # extraction function from resource to a dict of extracted variables
    def get_from_to(res, variables, term, dt_cutoff):
        for f in sorted(points_fields.keys()):
            if progressmode == 'verbose':
                print(f)
            nc_name = points_fields[f].get('nc_name', str(f))
            if f in res.listfields():
                fld = res.readfield(f)
                if fld.spectral:
                    fld.sp2gp()
                if points_fields[f].get('reproject_wind_on_lonlat', True):
                    try:
                        d, other = epygram.formats.fafields.find_wind_pair(f)
                    except epygram.epygramError:
                        pass
                    else:
                        other = res.readfield(other)
                        if d == 'y':  # other is v
                            wind = (fld, other)
                        elif d == 'x':  # other is u
                            wind = (other, fld)
                        wind = epygram.fields.make_vector_field(*wind)
                        if 'WIND.U.PHYS' in fld.fid['FA']:
                            map_factor_correction = False
                        else:
                            map_factor_correction = True
                        wind.reproject_wind_on_lonlat(map_factor_correction=map_factor_correction)
                        if d == 'y':
                            fld = wind.components[0]
                        else:
                            fld = wind.components[1]
                for loc, cll in coords.items():
                    point = fld.extract_point(*cll)
                    point.fid['netCDF'] = nc_name
                    if point.fid['netCDF'] not in variables[loc].keys():
                        variables[loc][point.fid['netCDF']] = point
                    else:
                        variables[loc][point.fid['netCDF']].extend(point)
            elif term == 0:  # exception for term 0, cumulated fields are not therein
                for loc, cll in coords.items():
                    val = FieldValidity(basis=dt_cutoff,
                                        date_time=dt_cutoff)
                    point = gimme_one_point(*cll,
                                            field_args={'fid':{'netCDF':nc_name},
                                                        'validity':val})
                    if point.fid['netCDF'] not in variables[loc].keys():
                        variables[loc][point.fid['netCDF']] = point
                    else:
                        variables[loc][point.fid['netCDF']].extend(point)
            else:
                raise RuntimeError('field {} not found in resource at {}, term={}'.format(
                    f, cutoffs[i], term))
        for f in profiles_FAfields.keys():
            if progressmode == 'verbose':
                print(f)
            nc_name = profiles_FAfields[f].get('nc_name', str(f))
            for loc, cll in coords.items():
                profile = res.extractprofile(f, *cll,
                                             vertical_coordinate=vd.get(vertical_coordinate, None),
                                             **extraction_options)
                profile.fid['netCDF'] = nc_name
                if profile.fid['netCDF'] not in variables[loc].keys():
                    variables[loc][profile.fid['netCDF']] = profile
                else:
                    variables[loc][profile.fid['netCDF']].extend(profile)

    # decumulation/operation
    def finalize_field(fld):
        if decumulate.get(fld.fid['netCDF']):
            center = decumulate.get(fld.fid['netCDF']) == 'center'
            fld.decumulate(center=center)
        if operations.get(fld.fid['netCDF']):
            fld.operation(*operations.get(fld.fid['netCDF']))

    # write to file
    def write_to(variables, nc, loc):
        for var in sorted(variables[loc].keys()):
            fld = variables[loc][var]
            finalize_field(fld)
            nc.writefield(fld, metadata=metadata.get(fld.fid['netCDF']))
        nc.set_global_attributes(**nc_global_attributes)
        nc.close()

    # 3. process data
    if prestaging:
        vortex_description['getmode'] = 'prestaging'
        vortex_description['date'] = cutoffs
        vortex_description['term'] = terms
        to_return = get_resources(**vortex_description)
    else:
        to_return = []
        if len(terms) > 1:
            for i in range(len(cutoffs)):
                vortex_description['date'] = cutoffs[i]
                variables = {loc:{} for loc in coords.keys()}
                for term in terms:
                    tmpfile = str(uuid.uuid4())
                    vortex_description['term'] = term
                    vortex_description['local'] = tmpfile
                    r = get_resources(**vortex_description)[0]
                    get_from_to(r, variables, term, dt_cutoffs[i])
                    r.close()
                    os.remove(r.container.abspath)
                for loc in coords.keys():
                    nc = epygram.formats.resource(
                        os.path.join(
                            outdir,
                            '_'.join([loc, 'T'.join([cutoffs[i][:8],
                                                     cutoffs[i][8:]]) + '.nc'])),
                        'w',
                        fmt='netCDF')
                    write_to(variables, nc, loc)
                    to_return.append(nc.container.abspath)
        else:  # 1 term only: series of dates at frozen term
            variables = {loc:{} for loc in coords.keys()}
            for i in range(len(cutoffs)):
                vortex_description['date'] = cutoffs[i]
                tmpfile = str(uuid.uuid4())
                vortex_description['term'] = terms[0]
                vortex_description['local'] = tmpfile
                r = get_resources(**vortex_description)[0]
                get_from_to(r, variables, terms[0], dt_cutoffs[i])
                r.close()
                os.remove(r.container.abspath)
            for loc in coords.keys():
                out = '-'.join(['T'.join([cutoffs[0][:8], cutoffs[0][8:]]),
                                'T'.join([cutoffs[-1][:8], cutoffs[-1][8:]])])
                nc = epygram.formats.resource(
                    os.path.join(
                        outdir,
                        '_'.join([loc, out + '.nc'])),
                    'w',
                    fmt='netCDF')
                write_to(variables, nc, loc)
                to_return.append(nc.container.abspath)
    return to_return
# End of extractor()


class Extractor(taylorism.Worker):
    """Independant extractor."""

    _footprint = dict(
        info="Run extractor().",
        attr=dict(
            directives=dict(
                info="Contains all arguments to extractor(...).",
                type=footprints.FPDict)
        )
    )

    def _task(self):
        directives = copy.copy(self.directives)
        return extractor(directives.pop('vortex_description'),
                         directives.pop('coords'),
                         directives.pop('start_cutoff'),
                         directives.pop('end_cutoff'),
                         directives.pop('start_term'),
                         directives.pop('end_term'),
                         **directives)


#############
# SHORTCUTS #
#############
# AROME
def get_gribfc_arome_oper(date, term, geometry='FRANGP0025', **others):
    """Proxy for AROME oper GRIBs."""

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
    """Proxy for AROME oper historic FAs."""

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


def get_gribfc_arome_xp(xp, date, term, geometry='FRANGP0025', **others):
    """Proxy for AROME experiment GRIBs."""

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
    """Proxy for AROME experiment historic FAs."""

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
    """Proxy for ARPEGE oper GRIBs."""

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
    """Proxy for ARPEGE oper historic FAs."""

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
    """Proxy for ARPEGE experiment GRIBs."""

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
    """Proxy for ARPEGE experiment historic FAs."""

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
