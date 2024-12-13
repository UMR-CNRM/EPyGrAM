#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Custom colormaps reading and tuning.
"""

import numpy
import json

from bronx.fancies import loggers

from . import config

#: No automatic export
__all__ = []

logger = loggers.getLogger(__name__)

_loaded_colormaps = {}
loaded_ColormapHelpers = {}


def register_colormap_from_json(filename):
    """
    Load colormap and metadata from file (json), register to matplotlib and
    return metadata.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    with open(filename, 'r') as f:
        asdict = json.load(f)
    # colormap itself
    colormap = asdict['name']
    colors = numpy.array(asdict['colors_RGB'], dtype=numpy.float64)
    if colors.max() > 1.:
        colors /= 255.
    asdict['colors_RGB'] = colors
    cmap = matplotlib.colors.ListedColormap(colors, name=colormap)
    asdict['cmap'] = cmap
    if colormap not in plt.colormaps():
        matplotlib.colormaps.register(cmap=cmap)
    else:
        raise ValueError('this colormap is already registered: {}'.format(colormap))
    _loaded_colormaps[filename] = asdict
    return asdict


def load_colormap(colormap):
    """Load colormap from epygram colormaps if needed."""
    import matplotlib.pyplot as plt
    if colormap not in plt.colormaps():
        if colormap in config.colormaps:
            cmapfile = config.colormaps[colormap]
            cmap = register_colormap_from_json(cmapfile)['cmap']
        else:
            raise ValueError("unknown colormap: {}".format(colormap))
    else:
        cmap = plt.get_cmap(colormap)
    return cmap


def get_ColormapHelper(colormap):
    if colormap in config.colormaps:
        colormap_file = config.colormaps[colormap]
    else:
        raise ValueError("unknown colormap: {}".format(colormap))
    return get_ColormapHelper_fromfile(colormap_file)


def get_ColormapHelper_fromfile(filename):
    """Get colormap from file (json) and build ad hoc ColormapHelper."""
    if filename in _loaded_colormaps:
        asdict = _loaded_colormaps[filename]
    else:
        asdict = register_colormap_from_json(filename)
    colormap = asdict['name']
    # colormap helper
    if asdict.get('colorcenters', False):
        ch = CenteredColormapHelper(colormap,
                                    explicit_colorcenters=asdict['colorcenters'],
                                    normalize=asdict.get('normalize', False),)
    elif asdict.get('colorbounds', False):
        ticks = asdict.get('ticks', None)
        if ticks == 'colorbounds':
            ticks = asdict['colorbounds']
        ch = ColormapHelper(colormap,
                            explicit_colorbounds=asdict['colorbounds'],
                            normalize=asdict.get('normalize', False),
                            explicit_ticks=ticks)
    else:
        ch = ColormapHelper(colormap, normalize=False)
    loaded_ColormapHelpers[colormap] = ch
    return ch


class ColormapHelper(object):
    """
    An integrated object helping for colormapping.
    """

    max_ticks = 15

    def __init__(self, colormap,
                 explicit_colorbounds=None,
                 normalize=False,
                 explicit_ticks=None):
        """
        A ColormapHelper is meant to help dealing with colormapping,
        especially mapping values to color changes, and prepares colormapping
        arguments for plotting functions.

        :param colormap: name of the colormap to be used
        :param explicit_colorbounds: to specify explicitly the colorbounds,
            i.e. values where colors need to change.
            Includes min value as first item, and max value as last item.
        :param normalize: if colors need to be normalized, i.e. that each color
            interval need to occupy the same space on the colorbar.
        :param explicit_ticks: to specify the ticks values to be shown
        """
        self.cmap_object = load_colormap(colormap)
        self.explicit_colorbounds = explicit_colorbounds
        self.normalize = normalize
        self.explicit_ticks = explicit_ticks
    
    @property
    def colormap(self):
        return self.cmap_object.name
    
    def colorbounds(self, minmax=None, number=None, step=None):
        """
        Get color bounds, i.e. values where colors change.
        Arguments are needed for implicit colorbounds only.

        :param minmax: (min, max) values of the values to be plotted
        :param number: number of different colors
        :param step: step in values from min to max, where to change colors

        Arguments number and step are mutually exclusive.
        """
        if self.explicit_colorbounds is not None:
            return self.explicit_colorbounds
        elif minmax is not None:
            assert None in (number, step), "Cannot provide both number and step"
            if step is None:
                if number is None:
                    number = 50
                return numpy.linspace(minmax[0], minmax[1], number)
            else:
                return numpy.arange(minmax[0], minmax[1] + step, step)
        elif minmax is None:
            raise ValueError("Must provide minmax if colorbounds is not known a priori.")

    def norm(self):
        """
        Normalize colormap, for each color occupy the same space on colorbar.
        Return the Norm object to be used by matplotlib.
        """
        import matplotlib
        assert self.explicit_colorbounds is not None, "Cannot compute norm if explicit_colorbounds are not known"
        colors = matplotlib.colors
        return colors.BoundaryNorm(boundaries=self.explicit_colorbounds,
                                   ncolors=self.cmap_object.N)

    def ticks_label(self, *args, **kwargs):
        """
        Return the labels of the ticks to be shown on colorbar.
        Arguments are passed to ticks_position(), in case of implicit ticks.
        """
        return self.ticks_position(*args, **kwargs)

    def ticks_position(self, *args, **kwargs):
        """
        Return the position of the ticks to be shown on colorbar.
        Arguments are passed to colorbounds(), in case of implicit ticks.
        """
        if self.explicit_ticks is not None:
            return self.explicit_ticks
        else:
            cbounds = self.colorbounds(*args, **kwargs)
            if len(cbounds) <= self.max_ticks:
                ticks = cbounds
            else:
                L = int((len(cbounds) - 1) // self.max_ticks) + 1
                ticks = [cbounds[i]
                         for i in range(len(cbounds) - (L // 3 + 1))
                         if i % L == 0] + [cbounds[-1]]
            return ticks

    def kwargs_for_plot(self, plot_method,
                        minmax=None,
                        center_cmap_on_0=False,
                        **colorbounds_kw):
        """
        Get kwargs for plot.

        :param plot_method: because arguments depend on plotting method.
        :param minmax: min and max values to be plot.
        :param center_cmap_on_0: if the colormap is to be centered on 0
            (diff plots).

        Other arguments are passed to colorbounds().
        """
        kwargs = dict(cmap=self.colormap)
        if plot_method not in ('scatter', 'pcolormesh'):
            kwargs['levels'] = self.colorbounds(minmax=minmax,
                                                **colorbounds_kw)
        if self.explicit_colorbounds is not None:
            if self.normalize:
                kwargs['norm'] = self.norm()
            kwargs['vmin'] = None
            kwargs['vmax'] = None
        else:
            if center_cmap_on_0:
                vmax = max(abs(minmax[0]), minmax[1])
                kwargs['vmin'] = -vmax
                kwargs['vmax'] = vmax
            else:
                kwargs['vmin'] = minmax[0]
                kwargs['vmax'] = minmax[1]
        return kwargs


class CenteredColormapHelper(ColormapHelper):

    def __init__(self, colormap, explicit_colorcenters, normalize=True):
        """
        A specific ColormapHelper, where the values of ticks are explicitly
        defined, and colors must surround each tick.

        :param colormap: colormap: name of the colormap to be used
        :param explicit_colorcenters: explicit values of the ticks and center
            values of each color.
        :param normalize: if colors need to be normalized, i.e. that each color
            interval need to occupy the same space on the colorbar.
        """
        self.cmap_object = load_colormap(colormap)
        colorbounds = [float(explicit_colorcenters[0]) -
                       0.5 * abs(explicit_colorcenters[0])]
        colorbounds += [float(explicit_colorcenters[i + 1] +
                              explicit_colorcenters[i]) / 2.
                        for i in range(len(explicit_colorcenters) - 1)]
        colorbounds += [float(explicit_colorcenters[-1]) +
                        0.5 * abs(explicit_colorcenters[-1])]
        self.explicit_colorbounds = colorbounds
        self.normalize = normalize
        self.explicit_ticks = explicit_colorcenters

    def ticks_label(self, *_, **__):
        """Return the labels of the ticks to be shown on colorbar."""
        return self.explicit_ticks

    def ticks_position(self, *_, **__):
        """Return the position of the ticks to be shown on colorbar."""
        return [(self.explicit_colorbounds[i] +
                 self.explicit_colorbounds[i + 1]) / 2.
                for i in range(len(self.explicit_colorbounds) - 1)]
