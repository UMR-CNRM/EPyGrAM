#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
This module contains:

- a class to handle variance spectrum;
- a function to compute DCT spectrum from a 2D field;
- a function to sort spectra with regards to their name;
- a function to plot a series of spectra;
- a function to read a dumped spectrum;
- a function to create a spectral geometry, e.g., for GRIBs.
"""

from bronx.graphics.axes import set_figax

import numpy
import copy

from epygram import epygramError, config
from epygram.geometries import GaussGeometry, SpectralGeometry
from epygram.util import RecursiveObject, write_formatted_table


_file_id = 'epygram.spectra.Spectrum'
_file_columns = ['#', 'lambda', 'variance']

def nlonlat_to_nsmax(nlon, nlat, stretching, trunctype):
    """
    Relationship between grid-point space and spectral space.
    nlat is the number of latitudes, nlon is the maximum number of longitudes.
    trunc_type should be one of "linear", "quadratic" or "cubic"
    Returns the maximum total dimensionless wavenumber.
    """
    trunctype2ratio = {"linear": 2, "quadratic": 3, "cubic": 4}
    if trunctype not in trunctype2ratio.keys():
        raise ValueError("trunc_type should be one of 'linear', 'quadratic' or 'cubic'")
    ratio = trunctype2ratio[trunctype]
    if stretching == 1.0:
        return int(numpy.floor((nlon - 1) / ratio))
    else:
        return int(numpy.floor(min(nlon - 1, 2 * nlat - 3) / ratio))


def make_spectral_geometry(geometry, trunctype="linear", verbose=False):
    """
    Returns a SpectralGeometry object consistent with the input grid-point
    geometry and with the required truncation type.
    This is only implemented for Gaussiand grids, with a triangular truncation.
    """
    if trunctype not in ("linear", "quadratic", "cubic"):
        raise ValueError("trunctype should be either 'linear', 'quadratic' or 'cubic'")
    if not isinstance(geometry, GaussGeometry):
        raise NotImplementedError(
            "No meaningful spectral transform implemented for " + geometry.name
        )
    if verbose:
        print(
            f"Build spectral geometry assuming {trunctype} and triangular truncation."
        )

    stretching = geometry.grid["dilatation_coef"]
    nlat = geometry.dimensions["lat_number"]
    nlon = geometry.dimensions["max_lon_number"]
    truncation = dict(
        max=nlonlat_to_nsmax(nlon, nlat, stretching, trunctype),
        shape="triangular",
    )
    spectral_geometry = SpectralGeometry("legendre", truncation)

    if verbose:
        print("Built spectral geometry:", spectral_geometry)

    return spectral_geometry


def get_spectral_geometry(field, resource, verbose=False):
    """
    Returns the SpectralGeometry object of the field or resource.
    If the field has no spectral geometry, return the spectral geometry of the resource.
    If the resource has no spectral geometry, returns None.
    """
    spectral_geometry = None
    if field.spectral_geometry is not None:
        spectral_geometry = field.spectral_geometry
    elif hasattr(resource, "spectral_geometry"):
        spectral_geometry = resource.spectral_geometry
    if verbose:
        print("Spectral geometry is", spectral_geometry)

    return spectral_geometry


def read_Spectrum(filename):
    """Read a Spectrum written in file and return it."""
    with open(filename, 'r') as _file:
        init_kwargs = {}
        # Spectrum file id
        assert _file.readline()[:-1] == _file_id, \
               ' '.join(["file:", filename, "does not contain a Spectrum."])
        # header: other stuff
        line = _file.readline()[:-1]
        while line[0] != '#':
            init_kwargs[line.split('=')[0].strip()] = line.split('=')[1].strip()
            line = _file.readline()[:-1]
        # columns description
        assert line.split() == _file_columns
        # values
        table = [line.split() for line in _file.readlines()]
        if int(table[0][0]) == 0:
            init_kwargs['mean2'] = float(table.pop(0)[2])
        elif not int(table[0][0]) == 1:
            raise epygramError("first wavenumber must be 0 or 1.")
        if 'resolution' in init_kwargs:
            init_kwargs['resolution'] = float(init_kwargs['resolution'])
        else:
            k = int(table[-1][0])
            init_kwargs['resolution'] = float(table[-1][1]) * k / (2 * (k + 1))
        variances = [float(line[2]) for line in table]

        return Spectrum(variances, **init_kwargs)


class Spectrum(RecursiveObject):
    """
    A spectrum can be seen as a quantification of a signal's variance with
    regards to scale.
    If the signal is defined in physical space on N points, its spectral
    representation will be a squared mean value (wavenumber 0) and variances for
    N-1 wavenumbers.
    For details and documentation, see
        Denis et al. (2002) : 'Spectral Decomposition of Two-Dimensional
        Atmospheric Fields on Limited-Area Domains
        Using the Discrete Cosine Transform (DCT)'
    """

    def __init__(self, variances,
                 name=None,
                 resolution=None,
                 mean2=None,
                 **kwargs):
        """
        :param variances: variances of the spectrum, from wavenumber 1 to N-1.
        :param name: an optional name for the spectrum.
        :param resolution: an optional resolution for the field represented by
                           the spectrum. It is used to compute the according
                           wavelengths.
                           Resolution unit is arbitrary, to the will of the
                           user.
        :param mean2: the optional mean^2 of the field, i.e. variance of
                      wavenumber 0 of the spectrum.
        """
        self.variances = numpy.array(variances)
        self.name = name
        self.resolution = resolution
        self.mean2 = mean2
        for k, v in kwargs.items():
            setattr(self, k, v)

    @property
    def wavenumbers(self):
        """Gets the wavenumbers of the spectrum."""
        return numpy.arange(1, len(self.variances) + 1)

    @property
    def wavelengths(self):
        """Gets the wavelengths of the spectrum."""
        K = len(self.variances) + 1
        return numpy.array([2. * self.resolution * K / k
                            for k in self.wavenumbers])

    def write(self, out):
        """
        Writes the spectrum with formatted output in *out*.

        :param out: must be an output open file-like object.
        """
        out.write(_file_id + '\n')
        if self.name is not None:
            out.write('name = ' + str(self.name) + '\n')
        if self.resolution is not None:
            out.write('resolution = ' + str(self.resolution) + '\n')
        table = [_file_columns,
                 [0, '-', self.mean2]]
        wn = self.wavenumbers
        wl = self.wavelengths
        var = self.variances
        for k in range(len(var)):
            table.append([wn[k], wl[k], var[k]])
        write_formatted_table(out, table)

    def dump(self, filename):
        """
        Writes the spectrum with formatted output in *filename*.
        """
        with open(filename, 'w') as _file:
            self.write(_file)

    def plotspectrum(self,
                     together_with=[],
                     over=(None, None),
                     slopes=[{'exp':-3, 'offset':1, 'label':'-3'},
                             {'exp':-5. / 3., 'offset':1, 'label':'-5/3'}],
                     zoom=None,
                     unit='SI',
                     title=None,
                     figsize=None,
                     takeover=False):
        """
        Plot the spectrum.

        :param together_with: another spectrum or list of spectra to plot on the
                              same ax.
        Cf. function plotspectra() of this module for other arguments.
        """
        if isinstance(together_with, Spectrum):
            together_with = [together_with]
        for s in together_with:
            assert isinstance(s, Spectrum)
        return plotspectra([self] + together_with,
                           over=over,
                           slopes=slopes,
                           zoom=zoom,
                           unit=unit,
                           title=title,
                           figsize=figsize)

##########
# internal
    def _check_operands(self, other):
        """Check compatibility of both spectra."""
        if isinstance(other, Spectrum):
            assert all((len(self.variances) == len(other.variances),
                        self.resolution == other.resolution or
                        None in (self.resolution, other.resolution))), \
                   "operations between spectra require that they share dimension and resolution."
        else:
            try:
                _ = float(other)
            except (ValueError, TypeError) as e:
                raise type(e)('*other* must be a Spectrum or a float-convertible.')
        if isinstance(other, Spectrum):
            othermean2 = other.mean2
            othername = other.name
            otherval = other.variances
        else:
            othermean2 = other
            othername = str(other)
            otherval = other

        return (othermean2, othername, otherval)

    def __add__(self, other):
        (othermean2, othername, otherval) = self._check_operands(other)
        mean2 = None if None in (self.mean2, othermean2) else self.mean2 + othermean2
        name = None if (self.name is othername is None) else str(self.name) + '+' + str(othername)
        return Spectrum(self.variances + otherval,
                        name=name,
                        resolution=self.resolution,
                        mean2=mean2)

    def __sub__(self, other):
        (othermean2, othername, otherval) = self._check_operands(other)
        mean2 = None if None in (self.mean2, othermean2) else self.mean2 - othermean2
        name = None if (self.name is othername is None) else str(self.name) + '+' + str(othername)
        return Spectrum(self.variances - otherval,
                        name=name,
                        resolution=self.resolution,
                        mean2=mean2)

    def __mul__(self, other):
        (othermean2, othername, otherval) = self._check_operands(other)
        mean2 = None if None in (self.mean2, othermean2) else self.mean2 * othermean2
        name = None if (self.name is othername is None) else str(self.name) + '+' + str(othername)
        return Spectrum(self.variances * otherval,
                        name=name,
                        resolution=self.resolution,
                        mean2=mean2)

    def __div__(self, other):
        (othermean2, othername, otherval) = self._check_operands(other)
        mean2 = None if None in (self.mean2, othermean2) else self.mean2 / othermean2
        name = None if (self.name is othername is None) else str(self.name) + '+' + str(othername)
        return Spectrum(self.variances / otherval,
                        name=name,
                        resolution=self.resolution,
                        mean2=mean2)


# FUNCTIONS FOR SPECTRA #
#########################
def sort(spectra):
    """ Sort a list of spectra with regards to their name. """
    untied_spectra = copy.copy(spectra)
    sortedspectra = []
    for f in sorted([s.name for s in untied_spectra], reverse=True):
        for s in range(len(untied_spectra)):
            if untied_spectra[s].name == f:
                sortedspectra.append(untied_spectra.pop(s))
                break
    return sortedspectra


def dctspectrum(x, verbose=False, log=None):
    """
    Function *dctspectrum* takes a 2D-array as argument and returns its 1D
    DCT ellipse spectrum.

    For details and documentation, see
        Denis et al. (2002) : 'Spectral Decomposition of Two-Dimensional
        Atmospheric Fields on Limited-Area Domains Using
        the Discrete Cosine Transform (DCT).'

    :param verbose: verbose mode
    :param log: an optional logging.Logger instance to which print info
                in *verbose* case.
    """
    import scipy.fftpack as tfm

    # compute transform
    if log is not None and verbose:
        log.info("dctspectrum: compute DCT transform...")
    norm = 'ortho'  # None
    y = tfm.dct(tfm.dct(x, norm=norm, axis=0), norm=norm, axis=1)

    # compute spectrum
    if log is not None and verbose:
        log.info("dctspectrum: compute variance spectrum...")
    N, M = y.shape
    N2 = N ** 2
    M2 = M ** 2
    MN = M * N
    K = min(M, N)
    variance = numpy.zeros(K)
    variance[0] = y[0, 0] ** 2 / MN
    for j in range(0, N):
        j2 = float(j) ** 2
        for i in range(0, M):
            var = y[j, i] ** 2 / MN
            k = numpy.sqrt(float(i) ** 2 / M2 + j2 / N2) * K
            k_inf = int(numpy.floor(k))
            k_sup = k_inf + 1
            weightsup = k - k_inf
            weightinf = 1.0 - weightsup
            if 0 <= k < 1:
                variance[1] += weightsup * var
            if 1 <= k < K - 1:
                variance[k_inf] += weightinf * var
                variance[k_sup] += weightsup * var
            if K - 1 <= k < K:
                variance[k_inf] += weightinf * var

    return variance


def global_spectrum(field):
    """
    Return variance spectrum of a global spectral field on a GaussGeometry,
    from its spectral coefficients.
    """
    assert field.geometry.isglobal
    assert isinstance(field.geometry, GaussGeometry)
    assert field.spectral_geometry.truncation["shape"] == "triangular"
    assert field.spectral
    data = field.getdata()
    nsmax = field.spectral_geometry.truncation["max"]
    assert data.size == (nsmax+1) * (nsmax+2)
    variances = numpy.zeros(nsmax + 1)
    jdata = 0
    # Loop on zonal wavenumber m
    for m in range(nsmax + 1):
        # Loop on total wavenumber n
        for n in range(m, nsmax + 1):
            squaredmodule = data[jdata] ** 2 + data[jdata+1] ** 2
            if m == 0:
                variances[n] += squaredmodule
                # Check that coefficient (m, n) == (0, n) is its own complex conjugate,
                # i.e. it has zero imaginary part
                assert data[jdata+1] == 0
            else:
                # Coefficient (-m, n) is the complex conjugate of (m, n).
                # It is not stored so (m, n) should be counted twice.
                variances[n] += 2 * squaredmodule
            jdata += 2
    return variances


def plotspectra(spectra,
                over=(None, None),
                slopes=[{'exp':-3, 'offset':1, 'label':'-3'},
                        {'exp':-5. / 3., 'offset':1, 'label':'-5/3'}],
                zoom=None,
                unit='SI',
                title=None,
                figsize=None,
                takeover=False):
    """
    To plot a series of spectra.

    :param over: any existing figure and/or ax to be used for the
                 plot, given as a tuple (fig, ax), with None for
                 missing objects. *fig* is the frame of the
                 matplotlib figure, containing eventually several
                 subplots (axes); *ax* is the matplotlib axes on
                 which the drawing is done. When given (is not None),
                 these objects must be coherent, i.e. ax being one of
                 the fig axes.
    :param spectra: a Spectrum instance or a list of.
    :param unit: string accepting LaTeX-mathematical syntaxes
    :param slopes: list of dict(
                   - exp=x where x is exposant of a A*k**-x slope
                   - offset=A where A is logscale offset in a A*k**-x slope;
                     a offset=1 is fitted to intercept the first spectra at wavenumber = 2
                   - label=(optional label) appearing 'k = label' in legend)
    :param zoom: dict(xmin=,xmax=,ymin=,ymax=)
    :param title: title for the plot
    :param figsize: figure sizes in inches, e.g. (5, 8.5).
                    Default figsize is config.plotsizes.
    :param takeover: give the user more access to the objects used in the
                     plot, by returning a dict containing them instead of only fig/ax
    """
    import matplotlib.pyplot as plt
    plt.rc('font', family='serif')
    if figsize is None:
        figsize = config.plotsizes

    fig, ax = set_figax(*over, figsize=figsize)
    result = {'fig':fig, 'ax':ax}

    if isinstance(spectra, Spectrum):
        spectra = [spectra]
    # prepare dimensions
    window = dict()
    window['ymin'] = min([min(s.variances) for s in spectra]) / 10
    window['ymax'] = max([max(s.variances) for s in spectra]) * 10
    window['xmax'] = max([max(s.wavelengths) for s in spectra]) * 1.5
    window['xmin'] = min([min(s.wavelengths) for s in spectra]) * 0.8
    if zoom is not None:
        for k, v in zoom.items():
            window[k] = v
    x1 = window['xmax']
    x2 = window['xmin']

    # colors and linestyles
    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'darkolivegreen',
              'yellow', 'salmon', 'black']
    linestyles = ['-', '--', '-.', ':']

    # axes
    if title is not None :
        ax.set_title(title)
    ax.set_yscale('log')
    ax.set_ylim(window['ymin'], window['ymax'])
    ax.set_xscale('log')
    ax.set_xlim(window['xmax'], window['xmin'])
    ax.grid()
    ax.set_xlabel('wavelength ($km$)')
    ax.set_ylabel(r'variance spectrum ($' + unit + '$)')

    # plot slopes
    # we take the second wavenumber (of first spectrum) as intercept, because
    # it is often better fitted with spectrum than the first one
    x_intercept = spectra[0].wavelengths[1]
    y_intercept = spectra[0].variances[1]
    i = 0
    result['slopes'] = []
    for slope in slopes:
        # a slope is defined by y = A * k**-s and we plot it with
        # two points y1, y2
        try:
            label = slope['label']
        except KeyError:
            # label = str(Fraction(slope['exp']).limit_denominator(10))
            label = str(slope['exp'])
        # because we plot w/r to wavelength, as opposed to wavenumber
        s = -slope['exp']
        A = y_intercept * x_intercept ** (-s) * slope['offset']
        y1 = A * x1 ** s
        y2 = A * x2 ** s
        line = ax.plot([x1, x2], [y1, y2], color='0.7',
                       linestyle=linestyles[i % len(linestyles)],
                       label=r'$k^{' + label + '}$')
        result['slopes'].append(line)
        i += 1

    # plot spectra
    i = 0
    result['spectra'] = []
    for s in spectra:
        line = ax.plot(s.wavelengths, s.variances, color=colors[i % len(colors)],
                       linestyle=linestyles[i // len(colors)], label=s.name)
        result['spectra'].append(line)
        i += 1

    # legend
    legend = ax.legend(loc='lower left', shadow=True)
    result['legend'] = legend
    for label in legend.get_texts():
        label.set_fontsize('medium')

    return result if takeover else (fig, ax)
