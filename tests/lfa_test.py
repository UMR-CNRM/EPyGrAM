#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import sys
import matplotlib.pyplot as plt
import epygram
import footprints
#footprints.logger.setLevel('INFO')
epygram.init_env()
filename = '/home/mary/EPyGrAM/bench_files/LFA/DHFDLARPE+0024'
simple_lfa_file = footprints.proxy.file(filename=filename)

"""
print simple_lfa

liste = simple_lfa.listfields()
print liste

f = simple_lfa.readfield('DOCFICHIER')
print f
f = simple_lfa.readfield('DOCD002')
print f
f = simple_lfa.readfield('DATE')
print f
f = simple_lfa.readfield('ECHEANCE')
print f
f = simple_lfa.readfield('VPP1')
print len(f.data)/22

simple_lfa.close()"""

myddh = epygram.formats.resource(filename, openmode='r', fmt='DDHLFA')
#myddh.what(sys.stdout)
print myddh.listfields()
"""for f in myddh.listfields():
    try:
        s = myddh.readfield(f)#[0].geometry.structure
        print f, s
    except: pass"""
print myddh.readfield('G01')[0]
f = myddh.readfield('VUU1')[0]
print f.geometry
f.plotfield()
plt.show()
#print myddh.readfield('FCTTUR')[0]
#print myddh.readfield('VUU1')[0]
