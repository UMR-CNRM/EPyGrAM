#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import argparse
import time
import os

import footprints
import taylorism

import epygram
from epygram import epylog
from epygram.util import printstatus
from epygram.args_catalog import add_arg_to_parser, \
                                 files_management, fields_management, \
                                 runtime_options, misc_options, \
                                 output_options



FA2GRIB1_dict = {'TEMPERATURE':{'table2Version':128, 'indicatorOfParameter':130, },
                 #'TEMPERATURE':{'table2Version':1, 'indicatorOfParameter':11, },
                 'HUMI.SPECIFI':{'table2Version':128, 'indicatorOfParameter':133, },
                 #'HUMI.SPECIFI':{'table2Version':1, 'indicatorOfParameter':51, },
                 'VOR':{'table2Version':128, 'indicatorOfParameter':138, },
                 #'VOR':{'table2Version':1, 'indicatorOfParameter':43, },
                 'DIV':{'table2Version':128, 'indicatorOfParameter':155, },
                 #'DIV':{'table2Version':1, 'indicatorOfParameter':44, },
                 'SURFPRESSION':{'table2Version':128, 'indicatorOfParameter':152, },
                 #'SURFPRESSION':{'table2Version':128, 'indicatorOfParameter':152, },
                 }

default_fieldseed = ['S[0-9][0-9][0-9]TEMPERATURE',
                     'S[0-9][0-9][0-9]HUMI.SPECIFI',
                     'SURFPRESSION',
                     'S[0-9][0-9][0-9]WIND.*.PHYS']
default_GRIB_options = {'centre':85}



class Femarser(taylorism.Worker):

    _footprint = dict(
        info="Compute difference between 2 FA files and store it as a spectral GRIB.",
        attr=dict(
            fileA=dict(
                info="First file (difference = FileA-FileB)."),
            fileB=dict(
                info="Second file (difference = FileA-FileB)."),
            fileout=dict(
                info="Name of output file."),
            fieldseed=dict(
                optional=True,
                type=footprints.FPList,
                default=footprints.FPList(default_fieldseed)),
            other_grib_options=dict(
                info="to specify other GRIB options (key/value pairs).",
                type=footprints.FPDict,
                optional=True,
                default=footprints.FPDict({}),
                access='rwx'),
            progressmode=dict(
                info="to print either progress 'percentage' of job, \
                      or field fids ('verbose').",
                optional=True,
                default=None)
        )
    )

    def _task(self):

        done = femars_2files(self.fileA, self.fileB,
                             fieldseed=self.fieldseed,
                             fileout=self.fileout,
                             other_GRIB_options=self.other_grib_options,
                             progressmode=self.progressmode)

        return done

def femars_2files(fileA, fileB,
                  fieldseed=default_fieldseed,
                  fileout=None,
                  other_GRIB_options=default_GRIB_options,
                  progressmode=None):
    """
    Computes fieldA-fieldB for the fields listed in *fieldseed*
    from *fileA*, *fileB*, transforms to spectral space if necessary,
    and writes the result into a GRIB file.
    """
    t0 = time.time()

    # set defaults
    packing = {'sphericalHarmonics':1,
               'bitsPerValue':24,
               'packingType':'spectral_simple'}
    sample = 'sh_ml_grib1'
    write_kwargs = {'grib_edition':1,
                    'packing':packing,
                    'sample':sample,
                    'other_GRIB_options':other_GRIB_options}
    if fileout is None:
        fileout = fileA + '-' + os.path.basename(fileB) + '.femars.out'

    # resources
    resourceA = epygram.formats.resource(fileA, openmode='r', fmt='FA')
    resourceB = epygram.formats.resource(fileB, openmode='r', fmt='FA')
    fileout = epygram.formats.resource(fileout, 'w', fmt='GRIB')
    assert resourceA.spectral_geometry == resourceB.spectral_geometry

    # set fieldslist
    fieldslistA = resourceA.find_fields_in_resource(fieldseed)
    fieldslistB = resourceB.find_fields_in_resource(fieldseed)
    assert fieldslistA == fieldslistB
    fieldslist = fieldslistA
    numfields = len(fieldslist)

    # separate wind if present, because of specific treatment: compute vor/div
    windfields = epygram.util.find_re_in_list('*WIND.?.PHYS', fieldslist)
    for f in windfields:
        fieldslist.pop(fieldslist.index(f))
    windfields = set([f.replace('U', '?').replace('V', '?') for f in windfields])
    windfields = sorted(list(windfields))

    n = 1
    # loop over fields
    for f in fieldslist:
        if progressmode == 'verbose':
            epylog.info(str(f))
        elif progressmode == 'percentage':
            printstatus(n, numfields)
            n += 1
        field = resourceA.readfield(f)
        if 'gauss' in field.geometry.name:
            raise NotImplementedError('Gauss geometries in femars: not yet !')
        field.operation_with_other('-', resourceB.readfield(f))
        if not field.spectral:
            field.gp2sp(resourceA.spectral_geometry)
        field.fid['GRIB1'] = FA2GRIB1_dict.get(field.fid['FA'], FA2GRIB1_dict.get(field.fid['FA'][4:]))
        field.validity.set(basis=field.validity.get())
        fileout.writefield(field, **write_kwargs)

    # U/V => vor/div particular case
    for f in windfields:
        if progressmode == 'verbose':
            epylog.info(str(f))
        elif progressmode == 'percentage':
            printstatus(n, numfields)
            n += 1
        (Ufid, Vfid) = resourceA.split_UV(f)
        # A
        UA = resourceA.readfield(Ufid[0])
        VA = resourceA.readfield(Vfid[0])
        vectA = epygram.fields.make_vector_field(UA, VA)
        vectA.map_factorize(reverse='True')
        if not vectA.spectral:
            vectA.gp2sp(resourceA.spectral_geometry)
        (vor, div) = vectA.compute_vordiv()
        # B
        UB = resourceB.readfield(Ufid[0])
        VB = resourceB.readfield(Vfid[0])
        vectB = epygram.fields.make_vector_field(UB, VB)
        vectB.map_factorize(reverse='True')
        if not vectB.spectral:
            vectB.gp2sp(resourceB.spectral_geometry)
        (vorB, divB) = vectB.compute_vordiv()
        # diff and back to spectral space
        vor.operation_with_other('-', vorB)
        div.operation_with_other('-', divB)
        vor.gp2sp(resourceA.spectral_geometry)
        div.gp2sp(resourceA.spectral_geometry)
        vor.fid['GRIB1'] = FA2GRIB1_dict['VOR']
        div.fid['GRIB1'] = FA2GRIB1_dict['DIV']
        fileout.writefield(vor, **write_kwargs)
        fileout.writefield(div, **write_kwargs)

    return ' '.join(['Difference', fileA, '-', fileB,
                     'computed and written to spectral GRIB in',
                     str(time.time() - t0), 's.'])

def femars_ensemble(filenames,
                    fieldseed=default_fieldseed,
                    other_GRIB_options=default_GRIB_options,
                    threads_number=2,
                    progressmode=None):
    """
    Computes field differences for a series of files, file to file,
    for the fields listed in *fieldseed*, transforms to spectral space if
    necessary, and writes the result into a GRIB file.
    
    Mandatory arguments:
    - *filenames*: list of names of the files to be processed 
    
    Optional arguments:
    - *fieldseed*: seed used to generate list of fields to be processed.
    - *other_GRIB_options*: to specify other GRIB options (key/value pairs).
    
    Technical named (optional) arguments:
    - *threads_number*: parallelisation of files processing
    - *progressmode*: among ('verbose', 'percentage', None)
    """

    filenames_ = filenames[1:] + [filenames[0]]
    if len(filenames) == 2:  # only one pair to be computed
        filenames = filenames[:-1]
        filenames_ = filenames_[:-1]

    common_instructions = {'other_grib_options':other_GRIB_options}
    if threads_number == 1:
        common_instructions['progressmode'] = progressmode
    if isinstance(fieldseed, str):
        common_instructions['fieldseed'] = footprints.FPList([fieldseed])
    individual_instructions = {'fileA':filenames_,
                               'fileB':filenames}
    #TODO: clean this
    temp_hack = True
    if temp_hack:
        fileouts = ['P' + filenames[i][7:10] + 'P' + filenames_[i][7:10] + '_gribdiff_1' for i in range(len(filenames))]
    else:
        fileouts = ['-'.join([filenames_[i], os.path.basename(filenames[i])])for i in range(len(filenames))]
    individual_instructions['fileout'] = fileouts

    # run !
    taylorism.batch_main(common_instructions, individual_instructions,
                         scheduler=taylorism.schedulers.MaxThreadsScheduler(threads_number),
                         verbose=(progressmode == 'verbose'))

def femars_avg(filenames,
               fieldseed=default_fieldseed,
               other_GRIB_options=default_GRIB_options,
               threads_number=2,
               progressmode=None):
    """
    Computes field differences to the average field for a series of files,
    for the fields listed in *fieldseed*, transforms to spectral space if
    necessary, and writes the result into a GRIB file.
    
    Mandatory arguments:
    - *filenames*: list of names of the files to be processed 
    
    Optional arguments:
    - *fieldseed*: seed used to generate list of fields to be processed.
    - *other_GRIB_options*: to specify other GRIB options (key/value pairs).
    
    Technical named (optional) arguments:
    - *threads_number*: parallelisation of files processing
    - *progressmode*: among ('verbose', 'percentage', None)
    """
    import tempfile

    # compute mean and save it (non parallel zone)
    tmpdir = tempfile.mkdtemp(dir=os.getenv('TMPDIR'))
    avg_tmp = '/'.join([tmpdir, 'avg_of_members.fa'])
    with epygram.formats.resource(filenames[0], 'r') as r:
        avg_resource = epygram.formats.resource(avg_tmp, 'w', fmt='FA',
                                                headername=r.headername,
                                                validity=r.validity,
                                                cdiden=r.cdiden,
                                                default_compression=r.default_compression,
                                                processtype=r.processtype)
    avg_resource.close()
    avg_resource.open(openmode='a')
    for i in range(len(filenames)):
        with epygram.formats.resource(filenames[i], 'r') as r:
            fieldslist = r.find_fields_in_resource(fieldseed)
            for f in fieldslist:
                field = r.readfield(f)
                if i == 0:  # only copy
                    avg_resource.writefield(field)
                else:  # add
                    avg = avg_resource.readfield(f)
                    avg.operation('+', field)
                    if i == len(filenames) - 1:  # last: divide by n !
                        avg.operation('/', float(len(filenames)))
                    avg_resource.writefield(avg)
    del avg_resource

    # classical femars against average (parallel zone)
    common_instructions = {'other_grib_options':other_GRIB_options,
                           'fileB':avg_tmp}
    if threads_number == 1:
        common_instructions['progressmode'] = progressmode
    if isinstance(fieldseed, str):
        common_instructions['fieldseed'] = footprints.FPList([fieldseed])
    individual_instructions = {'fileA':filenames}

    #TODO: clean this
    temp_hack = True
    if temp_hack:
        fileouts = ['P' + 'avg' + 'P' + filenames[i][7:10] + '_gribstd' for i in range(len(filenames))]
    else:
        fileouts = ['.'.join([f, 'femars_avg.out'])for f in filenames]
    individual_instructions['fileout'] = fileouts

    # run !
    taylorism.batch_main(common_instructions, individual_instructions,
                         scheduler=taylorism.schedulers.MaxThreadsScheduler(threads_number),
                         verbose=(progressmode == 'verbose'))

    # cleanings
    for f in os.listdir(tmpdir):
        os.remove('/'.join([tmpdir, f]))
    os.rmdir(tmpdir)

# end of main() ###############################################################



if __name__ == '__main__':

    ### 1. Parse arguments
    ######################
    parser = argparse.ArgumentParser(description='FEMARS based on epygram.',
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['several_files'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['FA_multiple_fields'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, misc_options['femars.diff_to_avg'])
    add_arg_to_parser(parser, output_options['outputfilename'], default=None)
    add_arg_to_parser(parser, output_options['GRIB_other_options'])
    add_arg_to_parser(parser, runtime_options['threads_number'])
    status = parser.add_mutually_exclusive_group()
    add_arg_to_parser(status, runtime_options['verbose'])
    #add_arg_to_parser(status, runtime_options['percentage'])

    args = parser.parse_args()

    ### 2. Initializations
    ######################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.verbose:
        progressmode = 'verbose'
    else:
        progressmode = None
    if args.GRIB_other_options is not None:
        other_GRIB_options = epygram.util.parse_str2dict(args.GRIB_other_options, try_convert=int)
    else:
        other_GRIB_options = default_GRIB_options

    # 2.2 list of fields to be processed
    if args.field != None:
        fieldseed = args.field
    elif args.listoffields != None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            fieldseed = l.readlines()
        for n in range(len(fieldseed)):
            fieldseed[n] = fieldseed[n].replace('\n', '').strip()
    else:
        fieldseed = default_fieldseed

    ### 3. Main
    ###########
    if args.diff_to_avg:
        femars_avg(args.filenames,
                   fieldseed=fieldseed,
                   other_GRIB_options=other_GRIB_options,
                   threads_number=args.threads_number,
                   progressmode=progressmode)
    else:
        femars_ensemble(args.filenames,
                        fieldseed=fieldseed,
                        other_GRIB_options=other_GRIB_options,
                        threads_number=args.threads_number,
                        progressmode=progressmode)

###########
### END ###
###########
