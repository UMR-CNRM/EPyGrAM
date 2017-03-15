#!/usr/bin/env python

import epygram
from footprints import proxy as fpx
import numpy
import matplotlib.pyplot as plt

epygram.init_env()

#Test configuration
test_plot = False #To plot some figures during the test
doTest = range(18) #list of tests to execute
terms = range(4) #terms of the run to use
baseDir = "/cnrm/mesonh/data1/riette/Public/epygram/smallRunFa" #Directory in which ICMSHAROM+* files are stored

#Results
allResults = {}

#Common resources used by different tests
resources = [epygram.formats.resource(baseDir + '/ICMSHAROM+' + str(term).zfill(4), 'r') for term in terms]
resources3D = [fpx.resource_modificator(name='CombineLevels', resource=r, openmode='r') for r in resources]
mult = fpx.resource_modificator(name='MultiValidities', resources=resources, openmode='r')
one = resources[0]
one3D = resources3D[0]
mult3D_1 = fpx.resource_modificator(name='CombineLevels', resource=mult, openmode='r') #combinelevels(multivalidities)
mult3D_2 = fpx.resource_modificator(name='MultiValidities', resources=resources3D, openmode='r') #multivalidities(combinelevels)
profile1 = one.extractprofile('S*TEMPERATURE', 8.77, 41.9)
section1 = one.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4))

if test_plot:
    #H2D anim
    field = mult.readfield('S060TEMPERATURE')
    field.sp2gp()
    anim=field.plotanimation(repeat=True, interval=100)
    plt.show()

    #Test H1D extraction and plot
    section1.getlevel(k=5).plotfield()
    sections = mult3D_2.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4))
    anim=sections.getlevel(k=1).plotanimation(repeat=True, interval=500)
    plt.show()

if 1 in doTest:
    #Build a virtualField from H1D - one validity
    print("Test 1")
    fields = epygram.base.FieldSet()
    for k in range(len(section1.geometry.vcoordinate.levels)):
        field = section1.getlevel(k=k)
        fields.append(field)
    section1bis = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='V2D', fieldset=fields)
    allResults[(1, 1)] = numpy.all(section1.getdata()==section1bis.getdata())
    allResults[(1, 2)] = section1.geometry==section1bis.geometry

if 2 in doTest:
    #Build a virtualField from H1D - several validities
    print("Test 2")
    sections = mult3D_2.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4))
    fields = epygram.base.FieldSet()
    for k in range(len(sections.geometry.vcoordinate.levels)):
        field = sections.getlevel(k=k)
        fields.append(field)
    sectionsbis = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='V2D', fieldset=fields)
    allResults[(2, 1)] = numpy.all(sections.getdata()==sectionsbis.getdata())
    allResults[(2, 2)] = sections.geometry==sectionsbis.geometry
    
if 3 in doTest:
    #Build a virtualField from H1D - one validity - sigma levels converted in P
    print("Test 3")
    newsection1 = one.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4), vertical_coordinate=100)
    fields = epygram.base.FieldSet()
    for k in range(len(newsection1.geometry.vcoordinate.levels)):
        field = newsection1.getlevel(k=k)
        fields.append(field)
    section1bis = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='V2D', fieldset=fields)
    if test_plot:
        plot = section1bis.plotfield()
        plt.show()
    allResults[(3, 1)] = numpy.all(newsection1.getdata()==section1bis.getdata())
    allResults[(3, 2)] =  newsection1.geometry==section1bis.geometry
    
sections = None
if 4 in doTest:
    #Build a virtualField from H1D - several validities - sigma levels converted in P
    print("Test 4")
    sections = mult3D_2.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4), vertical_coordinate=100, interpolation='nearest')
    if test_plot:
        anim = sections.plotanimation(repeat=True, interval=500, zoom={'ymax':900, 'ymin':1020})
        plt.show() #we can see mountains movement
    profiles = sections.as_profiles()[0]
    lon, lat = sections.geometry.get_lonlat_grid()
    profilesBis = mult.extractprofile('S*TEMPERATURE', lon[0], lat[0], vertical_coordinate=100, interpolation='nearest')
    allResults[(4, 1)] = numpy.all(profiles.getdata()==profilesBis.getdata())
    allResults[(4, 2)] = profiles.geometry==profilesBis.geometry
    if test_plot:
        anim = profiles.plotanimation(repeat=True, interval=500, zoom={'ymax':900, 'ymin':1020})
        plt.show()
    
if 5 in doTest:
    print("Test 5")
    if sections is None:
        sections = mult3D_2.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4), vertical_coordinate=100, interpolation='nearest')
    for i in range(len(resources)):
        sect1 = resources[i].extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4), vertical_coordinate=100, interpolation='nearest')
        sect2 = sections.getvalidity(i)
        allResults[(5, i, 1)] = numpy.all(sect1.getdata()==sect2.getdata())
        allResults[(5, i, 2)] = sect1.geometry == sect2.geometry
    
if 6 in doTest:
    print("Test 6")
    if sections is None:
        sections = mult3D_2.extractsection('S*TEMPERATURE', (8.77, 41.9), (8.9, 41.4), vertical_coordinate=100, interpolation='nearest')
    fields = epygram.base.FieldSet()
    for k in range(len(sections.geometry.vcoordinate.levels)):
        field = sections.getlevel(k=k)
        fields.append(field)
    sectionsbis = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='V2D', fieldset=fields)
    allResults[(6, 1)] = numpy.all(sections.getdata()==sectionsbis.getdata())
    allResults[(6, 2)] = sections.geometry==sectionsbis.geometry
    
# XXXXXXXX Test Subdomain, CombineLevels and D3VirtualField
print "Test Subdomain, CombineLevels and D3VirtualField"
one3D_profile = fpx.resource_modificator(name="Subdomain", resource=one3D, openmode='r', geometry=profile1.geometry)
profile2 = one3D_profile.readfield({'discipline': 0, 'parameterCategory': 0, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 0})
    
if 7 in doTest:
    print "Test 7 comp 1 and 2"
    allResults[(7, 1)] = numpy.all(profile1.data == profile2.data)
    allResults[(7, 2)] =  profile1.geometry == profile2.geometry
    
if 8 in doTest:
    print "Test 8 comp 1 and 3"
    temp_fids = [fid for fid in resources[0].listfields() if (fid[0:2] == 'S0' and fid[4:] == 'TEMPERATURE')]
    one3D_virtual = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='3D', resource=resources[0], resource_fids=temp_fids)
    one3D_virtual.sp2gp()
    profile3 = one3D_virtual.extract_subdomain(profile1.geometry)
    allResults[(8, 1)] = numpy.all(profile1.data == profile3.data)
    allResults[(8, 2)] = profile1.geometry == profile3.geometry
    
if 9 in doTest:
    print "Test 9 comp 1 and 3bis"
    temp_fids = [fid for fid in resources[0].listfields() if (fid[0:2] == 'S0' and fid[4:] == 'TEMPERATURE')]
    fields = epygram.base.FieldSet()
    for fid in temp_fids:
        field = resources[0].readfield(fid)
        fields.append(field)
    one3D_virtual = fpx.field(fid={'FA':'S*TEMPERATURE'}, structure='3D', fieldset=fields)
    one3D_virtual.sp2gp()
    profile3bis = one3D_virtual.extract_subdomain(profile1.geometry)
    allResults[(9, 1)] = numpy.all(profile1.data == profile3bis.data)
    allResults[(9, 2)] = profile1.geometry == profile3bis.geometry
    
if 10 in doTest:
    print "Test 10 comp 1 and 4"
    one3D_virtual_real = one3D_virtual.as_real_field()
    profile4 = one3D_virtual_real.extract_subdomain(profile1.geometry)
    allResults[(10, 1)] = numpy.all(profile1.data == profile4.data)
    allResults[(10, 2)] = profile1.geometry == profile4.geometry
    
if 11 in doTest:
    print "Test 11 comp 1 and 5"
    one3DVirtual = fpx.resource_modificator(name='CombineLevels', resource=one, openmode='r', virtual=True)
    profile5 = one3DVirtual.readfield({'discipline': 0, 'parameterCategory': 0, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 0})
    profile5.sp2gp()
    profile5 = profile5.extract_subdomain(profile1.geometry)
    allResults[(11, 1)] = numpy.all(profile1.data == profile5.data)
    allResults[(11, 2)] = profile1.geometry == profile5.geometry
    if test_plot:
        profile1.plotfield()
        plt.show()
    
if 12 in doTest:
    print("Test 12")
    fid = {'discipline': 0, 'parameterCategory': 1, 'level': 2, 'typeOfFirstFixedSurface': 103, 'parameterNumber': 1}
    f1 = one3DVirtual.readfield(fid)
    if f1.spectral: f1.sp2gp()
    f2 = one3D.readfield(fid)
    if f2.spectral: f2.sp2gp()
    f3 = one.readfield('CLSHUMI.RELATIVE')
    if f3.spectral: f3.sp2gp()
    allResults[(12, 1)] = numpy.all(f1.data==f2.data)
    allResults[(12, 2)] = numpy.all(f1.data==f3.data)
    if test_plot:
        f1.plotfield()
        f2.plotfield()
        f3.plotfield()
        plt.show()

# XXXXXXXX Test Subdomain, CombineLevels and MultiValidities
print "Test Subdomain, CombineLevels and MultiValidities"
fid = {'discipline': 0, 'parameterCategory': 0, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 0}
for iTest, geometry in [(13, profile1.geometry.deepcopy()), (14, section1.geometry.deepcopy())]:
    if iTest in doTest:
        fields = []
        print("Test " + str(iTest))
        for _ in range(len(geometry.vcoordinate.levels)) : geometry.vcoordinate.levels.pop()
        
        #subdo(comb(mult([file])))
        chain1 = fpx.resource_modificator(name="Subdomain", resource=mult3D_1, openmode='r', geometry=geometry.deepcopy(), interpolation='nearest')
        fields.append(chain1.readfield(fid))
        
        if test_plot:
            field = fields[0]
            anim=field.plotanimation(repeat=True, interval=100)
            plt.show()
        
        #subdo(mult([comb(file)]))
        chain2 = fpx.resource_modificator(name="Subdomain", resource=mult3D_2, openmode='r', geometry=geometry.deepcopy(), interpolation='nearest')
        fields.append(chain2.readfield(fid))
        
        #comb(subdo(mult([file])))
        subdo_mult = fpx.resource_modificator(name="Subdomain", resource=mult, openmode='r', geometry=geometry.deepcopy(), interpolation='nearest')
        chain3 = fpx.resource_modificator(name='CombineLevels', resource=subdo_mult, openmode='r')
        fields.append(chain3.readfield(fid))
        
        #mult([subdo(comb(file))])
        subdo_comb = [fpx.resource_modificator(name="Subdomain", resource=r3, openmode='r', geometry=geometry.deepcopy(), interpolation='nearest') for r3 in resources3D]
        chain4 = fpx.resource_modificator(name='MultiValidities', resources=subdo_comb, openmode='r')
        fields.append(chain4.readfield(fid))
        
        #comb(mult([subdo(file)]))
        subdo = [fpx.resource_modificator(name="Subdomain", resource=r, openmode='r', geometry=geometry.deepcopy(), interpolation='nearest') for r in resources]
        mult_subdo = fpx.resource_modificator(name='MultiValidities', resources=subdo, openmode='r')
        chain5 = fpx.resource_modificator(name='CombineLevels', resource=mult_subdo, openmode='r')
        fields.append(chain5.readfield(fid))
        
        #mult([comb(subdo(file))])
        comb_subdo = [fpx.resource_modificator(name='CombineLevels', resource=r, openmode='r') for r in subdo]
        chain6 = fpx.resource_modificator(name='MultiValidities', resources=comb_subdo, openmode='r')
        fields.append(chain6.readfield(fid))
        
        #comb(mult([file])).extractprofile
        field = mult3D_1.readfield(fid)
        if field.spectral:
            field.sp2gp()
        fields.append(field.extract_subdomain(geometry.deepcopy(), interpolation='nearest'))
        
        #mult(comb([file])).extractprofile
        field = mult3D_2.readfield(fid)
        if field.spectral:
            field.sp2gp()
        fields.append(field.extract_subdomain(geometry.deepcopy(), interpolation='nearest'))
        
        test = [numpy.all(field.data == fields[0].data) for field in fields]
        allResults[(iTest, 1)] = all(test)#, test
        test = [field.geometry == fields[0].geometry for field in fields]
        allResults[(iTest, 2)] = all(test)#, test

#XXXXXXX Test profile animation
if test_plot:
    print "Test profile animation"
    fid = {'discipline': 0, 'parameterCategory': 0, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 0}
    geometry = profile1.geometry.deepcopy()
    for _ in range(len(geometry.vcoordinate.levels)) : geometry.vcoordinate.levels.pop()
    profile = fpx.resource_modificator(name="Subdomain", resource=mult3D_1, openmode='r', geometry=geometry).readfield(fid)
    ani = profile.plotfield(force_mode='animation')
    # Set up formatting for the movie files
    #import matplotlib.animation as animation
    #writer = animation.writers['ffmpeg'](fps=15, metadata=dict(artist='Me'), bitrate=1800)
    #ani.save('profile_animation.mp4', fps=1, writer='imagemagick')
    plt.show()

# XXXXXXX Test CombineLevels and MultiValidities
if 15 in doTest:
    print "Test 15 CombineLevels and MultiValidities"
    field3D=one3D.readfield({'discipline': 0, 'parameterCategory': 2, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 11})
    field3D.sp2gp()
    res = []
    for i in range(60):
        field = one.readfield('S'+str(i+1).zfill(3)+'VERTIC.DIVER')
        field.sp2gp()
        res.append(numpy.all(field.getdata() == field3D.getdata()[i]))
    allResults[(15,)] = all(res)
    
if 16 in doTest:
    print("Test 16")
    field_1 = mult3D_1.readfield({'discipline': 0, 'parameterCategory': 2, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 11})
    field_2 = mult3D_2.readfield({'discipline': 0, 'parameterCategory': 2, 'typeOfFirstFixedSurface': 119, 'parameterNumber': 11})
    field_2.fid.pop('MultiValidities')
    allResults[(16,)] = field_1 == field_2
    
if 17 in doTest:
    #Test multivialidities
    print "Test 17 Multivalidities"
    fieldname = 'SURFTEMPERATURE'
    fieldname = 'SPECSURFGEOPOTEN'
    
    print "for one resource"
    one_field = one.readfield(fieldname)
    if one_field.spectral: one_field.sp2gp()
    one_data = one_field.getdata(subzone='CI', d4=True)
    
    print "for all resources"
    mult_field = mult.readfield(fieldname)
    if mult_field.spectral: mult_field.sp2gp()
    mult_data = mult_field.getdata(subzone='CI')
    
    allResults[(17,)] = numpy.all(one_data[0]==mult_data[0])

print("Results")
for k, v in allResults.items():
    print "  ", k, v
print
print(str(sum(allResults.values())) + " tests over " + str(len(allResults)) + " are successful.")
