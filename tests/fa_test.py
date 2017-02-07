#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import epygram
import matplotlib.pyplot as plt
import os
epygram.init_env()

def check_geom(filename):
    print '-->', filename
    thefa = epygram.formats.resource(filename, openmode='r', fmt='FA')
    fieldname = 'SURFIND.TERREMER'
    field = thefa.readfield(fieldname)
    field.plotfield(levelsnumber=3)
    plt.show()
    if mapfactor:
        mf = thefa.geometry.map_factor_field()
        mf.plotfield()
        plt.show()

def check_fields(filename):
    print '-->', filename
    infa = epygram.formats.resource(filename, openmode='r', fmt='FA')
    outfile = epygram.containers.File(filename=filename + '.copy')
    if os.path.exists(outfile.abspath):
        os.remove(outfile.abspath)
    outfa = epygram.formats.resource(filename + '.copy', openmode='w', fmt='FA',
                                     cdiden=infa.cdiden,
                                     headername=infa.headername,
                                     validity=infa.validity)
    for f in infa.listfields():
        field = infa.readfield(f)
        if epygram.formats.FA.inquire_field_dict(f)['type'] == 'H2D':
            if not field.spectral and not outfa.isopen:
                outfa.open(geometry=infa.geometry,
                           spectral_geometry=infa.spectral_geometry,
                           validity=infa.validity)
            fieldcompression = infa.fieldscompression[f]
            outfa.writefield(field, fieldcompression)
        else:
            if not outfa.isopen:
                outfa.open(geometry=infa.geometry,
                           spectral_geometry=infa.spectral_geometry,
                           validity=infa.validity)
            outfa.writefield(field)
    outfa.close()
    print "--- R/W ended; begin R/diff..."
    outfa = epygram.formats.resource(filename + '.copy', openmode='r', fmt='FA')
    significant_error = False
    for f in infa.listfields():
        infield = infa.readfield(f)
        outfield = outfa.readfield(f)
        if epygram.formats.FA.inquire_field_dict(f)['type'] == 'meteo':
            diff = infield - outfield
            stats = diff.stats()
            if max(abs(stats['min']), abs(stats['max'])) > epygram.config.epsilon:
                print "Difference (max) on", f, ":", max(abs(stats['min']), abs(stats['max']))
                significant_error = True
    if significant_error: print "-- End of R/W/R: numerical differences."
    else: print "-- End of R/W/R: reproducibility @", epygram.config.epsilon

def check_spectral(filename):
    print '-->', filename
    thefa = epygram.formats.resource(filename, openmode='r', fmt='FA')
    spfieldname = 'SPECSURFGEOPOTEN'
    spfield = thefa.readfield(spfieldname)
    spfield2 = thefa.readfield(spfieldname)
    spfield2.sp2gp()
    gpfield = thefa.readfield('SURFGEOPOTENTIEL')
    diff = gpfield - spfield2
    print "Spectral geopotential vs. Gridpoint geopotential"
    print diff.stats()
    spfield2.gp2sp(spfield.spectral_geometry)
    diff = spfield2 - spfield
    print "Difference after spectral round-trip (" + spfieldname + "):"
    print diff.stats()



###########
# Monitor #
###########
geom = True
fields = False
spectral = False
mapfactor = False
#footprints.logger.setLevel('DEBUG')

# Geometries
if geom:
    geomroot = '/home/mary/EPyGrAM/bench_files/FA/geometries/clim_'
    geomfileset = ['SP_HN_t', 'SP_HS_t', 'L_HN_t', 'L_HS_t', 'M_HN', 'M_HS', 'RG', 'RRG', 'LL_large', 'LL_small']
    for g in geomfileset[:]:
        check_geom(geomroot + g)

# Fields RWR (different processes, different fields)
if fields:
    fieldsroot = '/home/mary/EPyGrAM/bench_files/FA/fields/'
    fieldsfileset = ['arpege_previ3', 'arpege_analyse', 'arpege_canari', 'arpege_clim',
                     'arome_previ3', 'arome_analyse', 'arome_cpl', 'arome_clim',
                     'aromesurfex_previ3', 'aromesurfex_canari', 'aromesurfex_pgd',
                     'fullpos_lonlat_glob', 'fullpos_lonlat_lam',
                     'clim_arome_2.5.fa']
    for f in fieldsfileset[:1]:
        check_fields(fieldsroot + f)

# Spectral round-trip
if spectral:
    spectralroot = '/home/mary/EPyGrAM/bench_files/FA/spectral/clim_'
    spectralfileset = ['LAM', 'arpege']
    for s in spectralfileset[1:2]:
        check_spectral(spectralroot + s)

