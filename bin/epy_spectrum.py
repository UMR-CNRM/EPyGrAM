#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from __future__ import print_function, absolute_import, unicode_literals, division
import six

import os
import sys
import numpy
import argparse
import copy

from footprints import FPDict
from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

# Automatically set the python path
package_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_path, 'src'))
import epygram
from epygram import epylog, epygramError
import epygram.spectra as esp
from epygram.args_catalog import (add_arg_to_parser,
                                  files_management, fields_management,
                                  misc_options, output_options,
                                  runtime_options, graphical_options)
import matplotlib.pyplot as plt


def main(filename,
         fieldseed,
         subzone=None,
         refname=None,
         diffonly=False,
         computewind=False,
         verbose=False,
         noplot=False,
         legend=None,
         superposeplots=False,
         slopes=[{'exp':-3, 'offset':1, 'label':'-3'}, {'exp':-5. / 3., 'offset':1, 'label':'-5/3'}],
         zoom=None,
         unit='SI',
         output=False,
         outputfilename=None,
         figures_dpi=epygram.config.default_figures_dpi):
    """
    Computes and plots spectra.

    :param filename: name of the file to be processed.
    :param fieldseed: either a fid or a list of fid, used as a seed for
                   generating the list of fields to be processed.
    :param subzone: LAM zone among ('C', 'CI', 'CIE').
    :param refname: name of the reference file to be compared to.
    :param diffonly: if True, only plots the difference field.
    :param computewind: from fieldseed, gets U and V components of wind, and
                     computes the module and the module's spectrum.
    :param verbose: if True, verbose mode on.
    :param noplot: if True, disable the plot of profiles.
    :param legend: legend of plot.
    :param superposeplots: if True, superpose fields spectra on one plot.
    :param slopes: list of dict(exp=x where x is exposant of a A*k**-x slope
                             offset=A where A is logscale offset in a A*k**-x
                                      slope; an offset=1 is fitted to intercept
                                      the first spectra at wavenumber = 2.
                             label=(optional label) appearing 'k = label' in
                                                    legend)
    :param zoom: a dict(xmin=,xmax=,ymin=,ymax=) to restrain the plot.
    :param interpolation: kind of interpolation from grid to coordinates, among
                       ('nearest', 'linear', 'cubic').
    :param unit: scientifical unit of plot.
    :param output: output format for graphics, among ('png', 'pdf', False).
    :param outputfilename: specify an output filename for the profile.
                        The graphical output, if requested, will be completed
                        by the requested output format.
    :param figures_dpi: quality of saved figures.
    """
    resource = epygram.formats.resource(filename, openmode='r')
    diffmode = refname is not None
    if diffmode:
        reference = epygram.formats.resource(refname, openmode='r')

    if not diffmode and not computewind:
        spectra = []
        spectraplots = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(f) for f in fidlist]
        for f in fidlist:
            epylog.info(str(f))
            field = resource.readfield(f)
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            if field.spectral:
                field.sp2gp()
            if not field.geometry.projected_geometry:
                raise NotImplementedError("cannot compute spectra on regular_lonlat or Gauss grids.")
            variances = esp.dctspectrum(field.getdata(subzone=subzone),
                                        log=epylog,
                                        verbose=verbose)
            spectra.append(esp.Spectrum(variances[1:],
                                        name=str(f),
                                        resolution=field.geometry.grid['X_resolution'] / 1000.,
                                        mean2=variances[0]))
            if not noplot:
                # plot
                if legend is not None:
                    title = legend
                else:
                    title = resource.container.basename
                if not superposeplots:
                    spectraplots[f] = esp.plotspectra(spectra[-1],
                                                      slopes=slopes,
                                                      zoom=zoom,
                                                      unit=unit,
                                                      title=title)
                if not output:
                    plt.show()
        if superposeplots:
            sortedspectra = esp.sort(spectra)
            spectraplots['superposed'] = esp.plotspectra(sortedspectra,
                                                         slopes=slopes,
                                                         zoom=zoom,
                                                         unit=unit,
                                                         title=title)
            if not output:
                plt.show()
    elif not diffmode and computewind:
        spectra = []
        spectraplots = {}
        (Ufid, Vfid) = (sorted(resource.find_fields_in_resource(seed=fieldseed[0], fieldtype='H2D')),
                        sorted(resource.find_fields_in_resource(seed=fieldseed[1], fieldtype='H2D')))
        if resource.format == 'GRIB':
            Ufid = [FPDict(f) for f in Ufid]
            Vfid = [FPDict(f) for f in Vfid]
        assert len(Ufid) == len(Vfid), "number of fields for U mismatch those for V."
        for l in range(0, len(Ufid)):
            epylog.info(str(Ufid[l]) + str(Vfid[l]))
            U = resource.readfield(Ufid[l])
            if U.spectral:
                U.sp2gp()
            V = resource.readfield(Vfid[l])
            if V.spectral:
                V.sp2gp()
            field = U * U + V * V
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            field.setdata(numpy.sqrt(field.getdata()))
            variances = esp.dctspectrum(field.getdata(subzone=subzone),
                                        log=epylog,
                                        verbose=verbose)
            if resource.format == 'GRIB':
                name = {k:v for k, v in Ufid[l].items()
                        if Ufid[l][k] == Vfid[l][k]}
                name.update({k:str(Ufid[l][k]) + str(Vfid[l][k])
                             for k in Ufid[l].keys()
                             if Ufid[l][k] != Vfid[l][k]})
                name = str(name)
            else:
                _u = [c for c in Ufid[l]]
                _v = [c for c in Vfid[l]]
                name = ''
                for i in range(len(_u)):
                    if _u[i] == _v[i]:
                        name += _u[i]
                    else:
                        name += '*'
            spectra.append(esp.Spectrum(variances[1:],
                                        name=name,
                                        resolution=field.geometry.grid['X_resolution'] / 1000.,
                                        mean2=variances[0]))
            if not noplot:
                # plot
                if legend is not None:
                    title = legend
                else:
                    title = resource.container.basename
                if not superposeplots:
                    spectraplots[name] = esp.plotspectra(spectra[-1],
                                                         slopes=slopes,
                                                         zoom=zoom,
                                                         unit=unit,
                                                         title=title)
                if not output:
                    plt.show()
        if superposeplots:
            sortedspectra = esp.sort(spectra)
            spectraplots['superposed'] = esp.plotspectra(spectra,
                                                         slopes=slopes,
                                                         zoom=zoom,
                                                         unit=unit,
                                                         title=title)
            if not output:
                plt.show()
    else:
        spectra = []
        refspectra = []
        diffspectra = []
        spectraplots = {}
        fidlist = resource.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        reffidlist = reference.find_fields_in_resource(seed=fieldseed, fieldtype='H2D')
        if resource.format == 'GRIB':
            fidlist = [FPDict(f) for f in fidlist]
            reffidlist = [FPDict(f) for f in reffidlist]
        unionfidlist = list(set(fidlist).union(set(reffidlist)))
        intersectionfidlist = list(set(fidlist).intersection(set(reffidlist)))
        unionfidlist.sort()
        for f in unionfidlist:
            epylog.info(str(f))
            if f in fidlist:
                epylog.info("- in " + resource.container.basename)
                field = resource.readfield(f)
                if not field.geometry.grid.get('LAMzone', False):
                    subzone = None
                if field.spectral:
                    field.sp2gp()
                if not diffonly:
                    variances = esp.dctspectrum(field.getdata(subzone=subzone),
                                                log=epylog,
                                                verbose=verbose)
                    spectrum = esp.Spectrum(variances[1:],
                                            name=str(f),
                                            resolution=field.geometry.grid['X_resolution'] / 1000.,
                                            mean2=variances[0])
                    spectra.append(spectrum)
            if f in reffidlist:
                epylog.info("- in " + reference.container.basename)
                reffield = reference.readfield(f)
                if not field.geometry.grid.get('LAMzone', False):
                    subzone = None
                if reffield.spectral:
                    reffield.sp2gp()
                if not diffonly:
                    variances = esp.dctspectrum(reffield.getdata(subzone=subzone),
                                                log=epylog,
                                                verbose=verbose)
                    refspectrum = esp.Spectrum(variances[1:],
                                               name=str(f),
                                               resolution=reffield.geometry.grid['X_resolution'] / 1000.,
                                               mean2=variances[0])
                    refspectra.append(refspectrum)
            if f in intersectionfidlist:
                epylog.info("- on difference")
                diff = field - reffield
                variances = esp.dctspectrum(diff.getdata(subzone=subzone),
                                            log=epylog,
                                            verbose=verbose)
                diffspectrum = esp.Spectrum(variances[1:],
                                            name=str(f),
                                            resolution=diff.geometry.grid['X_resolution'] / 1000.,
                                            mean2=variances[0])
                diffspectra.append(diffspectrum)
            # PLOTS
            if not noplot:
                spectratoplot = []
                if f in fidlist and not diffonly:
                    spectratoplot.append(copy.copy(spectra[-1]))
                    spectratoplot[-1].name = resource.container.basename
                if f in reffidlist and not diffonly:
                    spectratoplot.append(copy.copy(refspectra[-1]))
                    spectratoplot[-1].name = reference.container.basename
                if f in intersectionfidlist:
                    spectratoplot.append(copy.copy(diffspectra[-1]))
                    spectratoplot[-1].name = 'diff'
                if legend is not None:
                    title = legend
                else:
                    title = str(f)
                sortedspectra = esp.sort(spectratoplot)
                spectraplots[f] = esp.plotspectra(sortedspectra,
                                                  slopes=slopes,
                                                  zoom=zoom,
                                                  unit=unit,
                                                  title=title)
                if not output:
                    plt.show()

    # Output
    epylog.info("save output...")
    suffix = "spectrum.out"
    # spectra
    if not diffmode:
        spectratosave = spectra
        outputrootname = resource.container.abspath
    else:
        spectratosave = diffspectra
        outputrootname = resource.container.absdir + \
                         '.'.join(['diff',
                                   '-'.join([resource.container.basename,
                                             reference.container.basename])])
    for s in spectratosave:
        if outputfilename:
            filename = outputfilename
        else:
            filename = '.'.join([outputrootname,
                                 smooth_string(s.name),
                                 suffix])
        s.write(open(filename, 'w'))
    if output:
        for p in spectraplots.keys():
            if outputfilename:
                filename = outputfilename
            else:
                filename = '.'.join([outputrootname,
                                     smooth_string(p),
                                     suffix,
                                     output])
            spectraplots[p][0].savefig(filename, dpi=figures_dpi)


if __name__ == '__main__':

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description="An EPyGrAM tool for computing DCT spectrum of \
                                                  meteorological fields (or difference of fields) from a resource.",
                                     epilog='End of help for: %(prog)s (EPyGrAM v' + epygram.__version__ + ')')

    add_arg_to_parser(parser, files_management['principal_file'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_management['field'])
    add_arg_to_parser(flds, fields_management['list_of_fields'])
    add_arg_to_parser(parser, fields_management['windfieldU'])
    add_arg_to_parser(parser, fields_management['windfieldV'])
    multi = parser.add_mutually_exclusive_group()
    add_arg_to_parser(multi, files_management['file_to_refer_in_diff'])
    add_arg_to_parser(multi, files_management['file_to_refer_in_diffonly'])
    add_arg_to_parser(multi, graphical_options['superpose_spectra_plots'])
    add_arg_to_parser(parser, output_options['output'])
    add_arg_to_parser(parser, output_options['outputfilename'])
    add_arg_to_parser(parser, output_options['noplot'])
    add_arg_to_parser(parser, misc_options['LAMzone'])
    add_arg_to_parser(parser, graphical_options['legend'])
    add_arg_to_parser(parser, graphical_options['scientifical_unit'])
    add_arg_to_parser(parser, graphical_options['spectra_slopes'])
    add_arg_to_parser(parser, graphical_options['spectra_zoom'])
    add_arg_to_parser(parser, graphical_options['figures_dpi'])
    add_arg_to_parser(parser, runtime_options['verbose'])

    args = parser.parse_args()

    # 2. Initializations
    ####################
    epygram.init_env()
    # 2.0 logs
    epylog.setLevel('WARNING')
    if args.verbose:
        epylog.setLevel('INFO')

    # 2.1 options
    if args.Drefname is not None:
        refname = args.Drefname
        diffonly = True
    else:
        refname = args.refname
        diffonly = False
    if args.zone in ('C', 'CI'):
        subzone = args.zone
    else:
        subzone = None
    slopes = []
    for s in args.kindofslopes.split(','):
        s = s.split()
        if len(s) > 0:
            # exp
            exp = s[0].strip()
            try:
                exp = float(exp)
            except Exception:  # try conversion of fraction string to float
                try:
                    exp = float(exp.split('/')[0]) / float(exp.split('/')[1])
                except Exception:
                    try:
                        exp = float(exp.split('e')[0]) * 10 ** float(exp.split('e')[1])
                    except Exception:
                        raise ValueError("slope exposant has to be numeric, fraction, or in scientific notation (1.1e-2), received " + exp)
            # offset
            try:
                offset = s[1].strip()
            except Exception:
                offset = 1.0
            else:
                try:
                    offset = float(offset)
                except Exception:
                    try:
                        exp = float(exp.split('e')[0]) * 10 ** float(exp.split('e')[1])
                    except Exception:
                        raise ValueError("slope offset has to be numeric or in scientific notation (1.1e-2), received " + offset)
            # label
            try:
                label = s[2].strip()
            except Exception:
                label = str(exp)
            slopes.append({'exp':exp, 'offset':offset, 'label':label})
    if args.zoom is not None:
        zoom = str2dict(args.zoom, float)
    else:
        zoom = None

    # 2.2 list of fields to be processed
    computewind = False
    if args.field is not None:
        fieldseed = args.field
    elif args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        fieldseed = (args.Ucomponentofwind, args.Vcomponentofwind)
        if None in fieldseed:
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        computewind = True
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            fieldseed = l.readlines()
        for n in range(len(fieldseed)):
            fieldseed[n] = fieldseed[n].replace('\n', '').strip()
    else:
        fieldseed = None

    # 3. Main
    #########
    main(six.u(args.filename),
         fieldseed,
         subzone=subzone,
         refname=refname,
         diffonly=diffonly,
         computewind=computewind,
         verbose=args.verbose,
         noplot=args.noplot,
         legend=args.legend,
         superposeplots=args.superposeplots,
         slopes=slopes,
         zoom=zoom,
         unit=args.unit,
         output=args.output if args.output != 'X' else None,
         outputfilename=args.outputfilename,
         figures_dpi=args.figures_dpi)
