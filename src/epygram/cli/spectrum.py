#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""An EPyGrAM tool for computing DCT spectrum of meteorological fields (or difference of fields) from a resource."""

import numpy
import argparse
import copy

from footprints import FPDict
from bronx.syntax.parsing import str2dict
from bronx.syntax.pretty import smooth_string

import epygram
from epygram import epylog as logger
from epygram import epygramError
import epygram.spectra as esp
from . import epilog
from .args_catalog import (add_arg_to_parser,
                           files_args,
                           fields_args,
                           misc_args,
                           output_args,
                           runtime_args,
                           graphical_args)
import matplotlib.pyplot as plt
from epygram.geometries import GaussGeometry, SpectralGeometry

_description = __doc__


def get_spectral_geometry(field, resource, verbose=False):
    """
    Returns the SpectralGeometry object of the field or resource.

    If the field has no spectral geometry, return the spectral geometry of the resource.
    If the resource has no spectral geometry and the grid is a Gaussian grid, 
    return a spectralGeometry object assuming linear and triangular truncation.
    """
    spectral_geometry = None
    if field.spectral_geometry is not None:
        spectral_geometry = field.spectral_geometry
    elif hasattr(resource, "spectral_geometry"):
        spectral_geometry = resource.spectral_geometry
    elif isinstance(field.geometry, GaussGeometry):
        if verbose:
            print(
                "Build spectral geometry assuming linear and triangular truncation"
            )
        truncation = dict(
            max=field.geometry.dimensions["lat_number"] - 1,
            shape="triangular",
        )
        spectral_geometry = SpectralGeometry("legendre", truncation)
    if verbose:
        print("Spectral geometry is", spectral_geometry)
    return spectral_geometry


def main():
    epygram.init_env()
    args = get_args()
    if args.verbose:
        logger.setLevel('INFO')
    else:
        logger.setLevel('WARNING')
    spectrum(
        args.filename,
        args.fieldseed,
        subzone=args.subzone,
        refname=args.refname,
        diffonly=args.diffonly,
        computewind=args.computewind,
        verbose=args.verbose,
        noplot=args.noplot,
        legend=args.legend,
        superposeplots=args.superposeplots,
        slopes=args.slopes,
        zoom=args.zoom,
        unit=args.unit,
        output=args.output if args.output != 'X' else None,
        outputfilename=args.outputfilename,
        figures_dpi=args.figures_dpi)


def spectrum(
    filename,
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
            logger.info(str(f))
            field = resource.readfield(f)
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            spectral_geometry = get_spectral_geometry(field, resource, verbose=verbose)
            spectra.append(field.spectrum(f,
                                          spectral_geometry=spectral_geometry,
                                          subzone=subzone,
                                          verbose=verbose))
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
            logger.info(str(Ufid[l]) + str(Vfid[l]))
            U = resource.readfield(Ufid[l])
            if U.spectral:
                U.sp2gp()
            V = resource.readfield(Vfid[l])
            if V.spectral:
                V.sp2gp()
            field = U * U + V * V
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None
            field.setdata(numpy.ma.sqrt(field.getdata()))
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
            spectral_geometry = get_spectral_geometry(field, resource, verbose=verbose)
            spectra.append(field.spectrum(name,
                                          spectral_geometry=spectral_geometry,
                                          subzone=subzone,
                                          verbose=verbose))
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
            logger.info(str(f))
            if f in fidlist:
                logger.info("- in " + resource.container.basename)
                field = resource.readfield(f)
                if not field.geometry.grid.get('LAMzone', False):
                    subzone = None
                if not diffonly:
                    spectral_geometry = get_spectral_geometry(field, resource, verbose=verbose)
                    spectra.append(field.spectrum(f,
                                                  spectral_geometry=spectral_geometry,
                                                  subzone=subzone,
                                                  verbose=verbose))
            if f in reffidlist:
                logger.info("- in " + reference.container.basename)
                reffield = reference.readfield(f)
                if not field.geometry.grid.get('LAMzone', False):
                    subzone = None
                if not diffonly:
                    spectral_geometry = get_spectral_geometry(reffield, reference, verbose=verbose)
                    refspectra.append(reffield.spectrum(f,
                                                        spectral_geometry=spectral_geometry,
                                                        subzone=subzone,
                                                        verbose=verbose))
            if f in intersectionfidlist:
                logger.info("- on difference")
                diff = field - reffield
                spectral_geometry = get_spectral_geometry(diff, resource, verbose=verbose)
                diffspectra.append(diff.spectrum(f,
                                                 spectral_geometry=spectral_geometry,
                                                 subzone=subzone,
                                                 verbose=verbose))
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
    logger.info("save output...")
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


def get_args():

    # 1. Parse arguments
    ####################
    parser = argparse.ArgumentParser(description=_description, epilog=epilog)
    add_arg_to_parser(parser, files_args['principal_file'])
    flds = parser.add_mutually_exclusive_group()
    add_arg_to_parser(flds, fields_args['field'])
    add_arg_to_parser(flds, fields_args['list_of_fields'])
    add_arg_to_parser(parser, fields_args['windfieldU'])
    add_arg_to_parser(parser, fields_args['windfieldV'])
    multi = parser.add_mutually_exclusive_group()
    add_arg_to_parser(multi, files_args['file_to_refer_in_diff'])
    add_arg_to_parser(multi, files_args['file_to_refer_in_diffonly'])
    add_arg_to_parser(multi, graphical_args['superpose_spectra_plots'])
    add_arg_to_parser(parser, output_args['output'])
    add_arg_to_parser(parser, output_args['outputfilename'])
    add_arg_to_parser(parser, output_args['noplot'])
    add_arg_to_parser(parser, misc_args['LAMzone'])
    add_arg_to_parser(parser, graphical_args['legend'])
    add_arg_to_parser(parser, graphical_args['scientifical_unit'])
    add_arg_to_parser(parser, graphical_args['spectra_slopes'])
    add_arg_to_parser(parser, graphical_args['spectra_zoom'])
    add_arg_to_parser(parser, graphical_args['figures_dpi'])
    add_arg_to_parser(parser, runtime_args['verbose'])
    args = parser.parse_args()

    # 2. Initializations
    ####################
    # 2.1 options
    if args.Drefname is not None:
        args.refname = args.Drefname
        args.diffonly = True
    else:
        args.refname = args.refname
        args.diffonly = False
    if args.zone in ('C', 'CI'):
        args.subzone = args.zone
    else:
        args.subzone = None
    args.slopes = []
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
            args.slopes.append({'exp':exp, 'offset':offset, 'label':label})
    if args.zoom is not None:
        args.zoom = str2dict(args.zoom, float)
    # 2.2 list of fields to be processed
    args.computewind = False
    if args.field is not None:
        args.fieldseed = args.field
    elif args.Ucomponentofwind is not None or args.Vcomponentofwind is not None:
        args.fieldseed = (args.Ucomponentofwind, args.Vcomponentofwind)
        if None in args.fieldseed:
            raise epygramError("wind mode: both U & V components of wind must be supplied")
        args.computewind = True
    elif args.listoffields is not None:
        listfile = epygram.containers.File(filename=args.listoffields)
        with open(listfile.abspath, 'r') as l:
            args.fieldseed = [line.replace('\n', '').strip() for line in l.readlines()]
    else:
        args.fieldseed = None

    return args

