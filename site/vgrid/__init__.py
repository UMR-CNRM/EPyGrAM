#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package **vgrid**:

Contains handling of Hybrid-Pressure vertical grid generation and plotting.
"""
from __future__ import print_function, absolute_import, unicode_literals, division

import io
import os
import re
import six
from collections import OrderedDict
import itertools

try:
    import bokeh
    from bokeh import models as bkm
    from bokeh import io as bkio
    from bokeh import plotting as bkp
    from bokeh import core as bkc
    from bokeh.palettes import Category10
except ImportError:
    raise ImportError('Seems like bokeh is not installed: install locally with "pip install --user bokeh" !')

from bronx.datagrip import namelist
from epygram.geometries.VGeometry import VGeometry
from epygram.profiles import hybridP2fluxpressure, flux2masspressures

MKVGRID_BINARY = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                              'vertical_discretization',
                              'mkvgrid.x')


class HybridPressureVGrid(object):
    """Handle and represent information about a Hybrid-Pressure Vertical Grid."""

    _re_grid_labels = '\s+'.join(
        ['(ILa)',
         '(A)', '(B)',
         '(Sigma)', '(1 - Sigma)', '(Rap Si-Hy)',
         '(Prehyd\(lbar\))', '(Prehyd\(l\))', '(\[D Prehyd\]\(l\))',
         '(Alti\(l\))', '(Alti\(lbar\))', '(\[D Alti\]\(l\))'])
    colors = itertools.cycle(Category10[10])
    complete_tooltips = [('level', "@level"),
                         ('p', "@p{0.2f} hPa"),
                         ('z', "@z{0.2f} m"),
                         ('top', "@psup{0.2f} hPa | @zsup{0.2f} m"),
                         ('bottom', "@pinf{0.2f} hPa | @zinf{0.2f} m"),
                         ('thickness', "@pthickness{0.2f} hPa | @zthickness{0.2f} m")]

    @classmethod
    def next_color(cls):
        return bkc.properties.value(next(cls.colors))

    def __init__(self, source,
                 vertical_mean=None,
                 reference_pressure=101325.,
                 ptop=0.,
                 vgrid_name=None):
        """
        Build the grid object from a **source**, which may be

            - an infosup file from ``mkvgrid.x``
            - an epygram Hybrid-Pressure VGeometry; in which case the below
              arguments must/may be provided:

        Case of an epygram Hybrid-Pressure VGeometry **source**:

        :param vertical_mean: mandatory, among ('geometric', 'arithmetic', 'LAPRXPK=False')
        :param reference_pressure: at the surface, for building a standard atmosphere
        :param ptop: pressure at the top of the model (upper boundary of the upper layer)
        :param vgrid_name: name of the grid, for plot/saving purpose
        """
        if isinstance(source, six.string_types):
            self.name = os.path.basename(source).replace('.infosup', '')
            self._parse_infosup(source)
        elif isinstance(source, VGeometry):
            assert vertical_mean is not None, \
                "must provide a **vertical_mean** among ('geometric', 'arithmetic', 'LAPRXPK=False')"
            self._from_epygram(source,
                               vertical_mean=vertical_mean,
                               reference_pressure=reference_pressure,
                               ptop=ptop,
                               vgrid_name=vgrid_name)
        else:
            raise NotImplementedError('construction from else that source:' + str(type(source)))

    def _parse_infosup(self, infosup_file):
        """Parse an infosup file from ``mkvgrid.x`` program."""
        # read
        with io.open(infosup_file, 'r') as f:
            lines = [line.strip('\n').strip() for line in f.readlines()]
        lines = [line for line in lines if line != '']
        # dimensions
        self.dimensions = dict()
        for line in lines[0:5]:
            k,v = line.split('=')
            k = k.strip('*').strip()
            self.dimensions[k] = int(v)
        # domains
        self.domains = list()
        cursor = 5
        for line in lines[cursor:cursor + self.ndoms]:
            k,v = line.split('=')
            k = k.strip('*').strip()
            k = k.split()[-1]
            self.domains.append((k, int(v)))
        self.domains = self.domains[::-1]
        cursor += self.ndoms
        # domains altitudes and heights
        self.domains_altitudes = dict()
        self.domains_widths = dict()
        for line in lines[cursor:cursor + self.ndoms]:
            alt, height = line.split(';')
            k,v = alt.split('=')
            k = k.strip('*').strip()
            self.domains_altitudes[k] = float(v)
            k,v = height.split('=')
            k = k.strip()
            self.domains_widths[k] = float(v)
        cursor += self.ndoms
        alt, height = lines[cursor].split(';')
        k,v = alt.split(':')
        k = k.strip('*').strip()
        self.domains_altitudes[k] = float(v)
        k,v = height.split('=')
        k = k.strip()
        self.domains_widths[k] = float(v)
        cursor += 1
        # parameters
        self.params = dict()
        for line in lines[cursor:cursor + 5]:
            k,v = line.split('=')
            k = k.strip('*').strip()
            try:
                self.params[k] = float(v)
            except ValueError:
                self.params[k] = True if v.strip() == 'T' else False
        cursor += 5
        # levels
        self.levels = dict()
        match = re.match(self._re_grid_labels, lines[cursor])
        if match:
            labels = match.groups()
        else:
            raise SyntaxError('Labels of levels table !')
        cursor += 2
        for line in lines[cursor:cursor + self.nlevels]:
            line = line.split()
            i, vals = int(line[0]), [float(v) for v in line[1:]]
            self.levels[i] = dict(zip(labels[1:], vals))

    def _from_epygram(self, vgeometry,
                      vertical_mean,
                      reference_pressure=101325.,
                      ptop=0.,
                      vgrid_name=None):
        """
        Convert an epygram Hybrid-Pressure VGeometry object.

        :param vertical_mean: among ('geometric', 'arithmetic', 'LAPRXPK=False')
        :param reference_pressure: at the surface, for building a standard atmosphere
        :param ptop: pressure at the top of the model (upper boundary of the upper layer)
        :param vgrid_name: name of the grid, for plot/saving purpose
        """
        from . import standard_atmosphere
        assert vgeometry.typeoffirstfixedsurface == 119
        A = [level[1]['Ai'] for level in vgeometry.grid['gridlevels'][1:]]
        B = [level[1]['Bi'] for level in vgeometry.grid['gridlevels'][1:]]
        # compute standard atmosphere at levels
        p_lbar = hybridP2fluxpressure(A, B, reference_pressure).squeeze()
        p_l = flux2masspressures(p_lbar, vertical_mean, Ptop=ptop)
        a_l = standard_atmosphere.altitude_at(p_l)
        a_lbar = standard_atmosphere.altitude_at(p_lbar)
        # dimensions
        self.dimensions = {'Number of levels':len(p_l)}
        # parameters
        self.params = {
            'LLAPRXPK':vertical_mean != 'LAPRXPK=False',
            'Reference pressure at 0m: ZVP00PR':reference_pressure}
        if vgrid_name is not None:
            self.name = vgrid_name
        else:
            self.name = 'from_field_L' + str(self.nlevels)
        # levels
        self.levels = dict()
        for i in range(1, self.nlevels + 1):
            if i == 1:
                D_a = 2 * (a_lbar[0] - a_lbar[1])
            else:
                D_a = (a_l[i - 2] - a_lbar[i - 2]) + (a_l[i - 1] - a_lbar[i - 1])
            self.levels[i] = {
                'A':A[i - 1],
                'B':B[i - 1],
                'Alti(l)':a_l[i - 1],
                'Alti(lbar)':a_lbar[i - 1],
                '[D Alti](l)':D_a,
                'Prehyd(l)':p_l[i - 1],
                'Prehyd(lbar)':p_lbar[i - 1]}

    def _as_bokeh_ColumnDataSource(self):
        """Convert to a bokeh ColumnDataSource object."""
        return bkm.ColumnDataSource(self.grid)

    @property
    def nlevels(self):
        return self.dimensions['Number of levels']

    @property
    def ndoms(self):
        return self.dimensions['Number of vertical domains']

    def unit(self, yaxis):
        return {'p':'hPa', 'z':'m', 'l':'#'}[yaxis]

    @property
    def grid(self):
        """
        Get grid height and pressure of full(mass) levels and frontiers.
        Pressures are given in hPa.
        """
        grid = {'z':[self.levels[i]['Alti(l)'] for i in range(1, self.nlevels + 1)],
                'zinf':[self.levels[i]['Alti(lbar)'] for i in range(1, self.nlevels + 1)],
                'zsup':([self.levels[1]['Alti(l)'] + self.levels[1]['[D Alti](l)'] / 2.] +  # TOBECHECKED:
                        [self.levels[i - 1]['Alti(lbar)'] for i in range(2, self.nlevels + 1)]),
                'p':[self.levels[i]['Prehyd(l)'] / 100 for i in range(1, self.nlevels + 1)],
                'pinf':[self.levels[i]['Prehyd(lbar)'] / 100 for i in range(1, self.nlevels + 1)],
                'psup':[0.] + [self.levels[i - 1]['Prehyd(lbar)'] / 100 for i in range(2, self.nlevels + 1)],
                'level':range(1, self.nlevels + 1),
                }
        for y in ('p', 'z'):
            grid[y + 'thickness'] = [abs(grid[y + 'sup'][i] - grid[y + 'inf'][i])
                                     for i in range(self.nlevels)]
        return grid

    def domains_lower_limit(self):
        """
        Get the lower limit of each domain in pressure, height and level number.
        Pressures are given in hPa.
        Lower level is "last included level + 0.5".
        """
        domains_lower_limit = OrderedDict()
        i = 0
        for d, d_nlev in self.domains:  # these are ordered from top to bottom
            i += d_nlev
            d_lower_pressure = self.levels[i]['Prehyd(lbar)'] / 100
            d_lower_height = self.levels[i]['Alti(lbar)']
            domains_lower_limit[d] = {'z':d_lower_height,
                                      'p':d_lower_pressure,
                                      'level':float(i) + .5}
        return domains_lower_limit

    def domains_upper_limit(self):
        """
        Get the upper limit of each domain in pressure, height and level number.
        Pressures are given in hPa.
        Lower level is "last included level - 0.5".
        """
        domains_upper_limit = OrderedDict()
        domains_upper_limit[self.domains[0][0]] = {'z':self.grid['z'][0] + self.grid['zthickness'][0] / 2.,
                                                   'p':self.grid['p'][0] + self.grid['pthickness'][0] / 2.,
                                                   'level':0.5}
        domains = list(self.domains_lower_limit().keys())
        limits = list(self.domains_lower_limit().values())
        for i in range(1, len(domains)):
            # the upper limit of the below domain is the lower limit of the above domain
            domains_upper_limit[domains[i]] = limits[i - 1]
        return domains_upper_limit

    def get_AB(self):
        """Get the lists of A & B coefficients of the N+1 flux levels."""
        return {'A':[0.] + [self.levels[i]['A'] for i in range(1, self.nlevels + 1)],
                'B':[0.] + [self.levels[i]['B'] for i in range(1, self.nlevels + 1)]}

    def set_domains(self, domain_names, domain_levels_numbers):
        """
        Set a series of domains.

        :param domain_names: list of the name of each domain,
            ordered from top to bottom
        :param domain_levels_numbers: list of the number of levels in each domain,
            ordered from top to bottom
        """
        assert len(domain_names) == len(domain_levels_numbers)
        assert sum(domain_levels_numbers) == self.nlevels
        self.dimensions['Number of vertical domains'] = len(domain_names)
        self.domains = list()
        for i in range(len(domain_names)):
            self.domains.append((domain_names[i], domain_levels_numbers[i]))

# applicative methods
    def write_AB_to_namelist(self):
        """Write A and B coefficients to namelist(s)."""
        AB = self.get_AB()
        A = AB['A']
        B = AB['B']
        NAM = namelist.NamelistBlock('NAMFPG')
        NAM['FPVALH(0)'] = A
        NAM['FPVBH(0)'] = B
        with io.open('.'.join([self.name, NAM.name]), 'w') as o:
            o.write(NAM.dumps(sorting=namelist.FIRST_ORDER_SORTING))
        NAM = namelist.NamelistBlock('NAMVV1')
        NAM['DVALH(0)'] = A
        NAM['DVBH(0)'] = B
        with io.open('.'.join([self.name, NAM.name]), 'w') as o:
            o.write(NAM.dumps(sorting=namelist.FIRST_ORDER_SORTING))

    def back_to_mkvgrid_namelist(self, namel_filename=None):
        """Generate back a mkvgrid.x namelist."""
        if namel_filename is None:
            namel_filename = '.'.join([self.name, 'mkvgrid', 'nam'])
        # NAM_DIM
        assert hasattr(self, 'domains')
        NAM = namelist.NamelistBlock('NAM_DIM')
        NAM['JPN'] = self.nlevels
        NAM['JPNPRES'] = len([b for b in self.get_AB()['B'][1:] if b == 0.])
        NAM['JPNSIGM'] = len([b for b in self.get_AB()['B'][1:] if b == 1.])
        NAM['JPDOM'] = self.ndoms
        with io.open(namel_filename, 'w') as o:
            o.write(NAM.dumps())
        # NAM_REF
        NAM = namelist.NamelistBlock('NAM_REF')
        NAM['ZP1'] = self.grid['p'][0] * 100
        NAM['ZVP00'] = self.grid['pinf'][-1] * 100
        NAM['ZVP00PR'] = self.params['Reference pressure at 0m: ZVP00PR']
        NAM['ZALT_BOT'] = self.grid['z'][-1]
        with io.open(namel_filename, 'a') as o:
            o.write(NAM.dumps())
        # NAM_DOM
        NAM = namelist.NamelistBlock('NAM_DOM')
        for i in range(1, self.ndoms):
            dom, nlev = self.domains[self.ndoms - i]
            NAM['ZALT_DOM({})'.format(i)] = self.domains_upper_limit()[dom]['z']
            upper_level = int(self.domains_upper_limit()[dom]['level'] + 0.5)
            NAM['ZDALT_DOM({})'.format(i)] = self.grid['zthickness'][upper_level - 1]
            NAM['IT_DOM({})'.format(i)] = nlev
            NAM['CLNAM_DOM({})'.format(i)] = dom
        NAM['ZDALT_DOM({})'.format(self.ndoms)] = self.grid['zinf'][0] - self.grid['zinf'][1]
        dom, nlev = self.domains[0]
        NAM['IT_DOM({})'.format(self.ndoms)] = nlev
        NAM['CLNAM_DOM({})'.format(self.ndoms)] = dom
        with io.open(namel_filename, 'a') as o:
            o.write(NAM.dumps(sorting=namelist.FIRST_ORDER_SORTING))
        # NAM_PARAM
        NAM = namelist.NamelistBlock('NAM_PARAM')
        if 'LLAPRXPK' in self.params:
            NAM['LLAPRXPK'] = self.params['LLAPRXPK']
        with io.open(namel_filename, 'a') as o:
            o.write(NAM.dumps())
        return namel_filename

    def bokeh_plot_y_vs_x(self, yaxis, xaxis,
                          over=None,
                          hover_attachment='right',
                          bokeh_kwargs={},
                          fig_kwargs={},
                          legend_location='top_right',
                          plot_domains=None,
                          add_domains_options={}):
        """
        Plot grid and return the figure.

        :param yaxis: 'level', 'p' (pressure), 'z' (height) or 'pthickness'//'zthickness'
        :param xaxis: 'level', 'p' (pressure), 'z' (height) or 'pthickness'//'zthickness'
        :param over: an eventual existing bokeh figure on which to plot
        :param hover_attachment: position of the hover
        :param bokeh_kwargs: passed to bokeh's inner plot functions
        :param fig_kwargs: to be passed to bokeh's figure()
        :param plot_domains: plot a horizontal bar separating each vertical domain (True/False/None)
        :param add_domains_options: options to be passed to bokeh_add_domains() method
        """
        import uuid
        thisid = str(uuid.uuid4())
        cds = self._as_bokeh_ColumnDataSource()
        # options
        axes_domain = dict(xaxis=xaxis, yaxis=yaxis)
        kwargs = dict(source=cds,
                      name=thisid)
        if yaxis in ('z', 'p') and xaxis in ('z', 'p'):
            # p/z or z/p
            method = 'vbar'
            kwargs.update(x=xaxis,
                          top=yaxis + 'sup',
                          bottom=yaxis + 'inf',
                          width=xaxis + 'thickness')
        elif yaxis in ('z', 'p'):
            # z or p as y: thickness as bar
            kwargs.update(x=xaxis)
            if xaxis == 'level':
                method = 'vbar'
                kwargs.update(top=yaxis + 'sup',
                              bottom=yaxis + 'inf',
                              width=1.)
            else:
                method = 'circle'
                kwargs.update(y=yaxis,
                              line_color='black',
                              size=9)
        elif xaxis in ('z', 'p'):
            # z or p as x: thickness as bar
            kwargs.update(y=yaxis)
            if yaxis == 'level':
                method = 'hbar'
                kwargs.update(left=xaxis + ('inf' if xaxis == 'z' else 'sup'),
                              right=xaxis + ('sup' if xaxis == 'z' else 'inf'),
                              height=1.)
            else:
                method = 'circle'
                kwargs.update(x=xaxis,
                              line_color='black',
                              size=9)
        else:
            method = 'circle'
            kwargs.update(x=xaxis,
                          y=yaxis,
                          line_color='black',
                          size=9)
        color = bokeh_kwargs.get('color', self.next_color())
        legend = bokeh_kwargs.get('legend_label', self.name)
        kwargs.update(color=color,
                      legend_label=legend)
        kwargs.update(bokeh_kwargs)
        # figure
        if over is None:
            def axes(axis):
                return ({'z':'Altitude',
                         'p':'Pressure',
                         'level':'Level',
                         'zthickness':'Thickness',
                         'pthickness':'Thickness'}[axis] +
                        ' ({})'
                        ).format(self.unit(axis[0]))
            fig = bkp.figure(sizing_mode='stretch_both',
                             title="Vertical Coordinate Grid",
                             **fig_kwargs)
            fig.xaxis.axis_label = axes(xaxis)
            fig.yaxis.axis_label = axes(yaxis)
            # reverse pressure axis
            if yaxis == 'p':
                extra = max(cds.data[yaxis + 'inf']) * 0.05
                fig.y_range = bkm.Range1d(max(cds.data[yaxis + 'inf']) + extra,
                                          0. - extra)
            elif yaxis == 'level':
                extra = float(self.nlevels) * 0.05
                fig.y_range = bkm.Range1d(self.nlevels + extra,
                                          1. - extra)
        else:
            fig = over
        if not any([isinstance(t, bkm.tools.CrosshairTool) for t in fig.tools]):
            fig.add_tools(bkm.tools.CrosshairTool())
        # actual plot command
        getattr(fig, method)(**kwargs)
        # hover
        tooltips = self.complete_tooltips
        if 'legend_label' in kwargs:
            tooltips = [('GRID', kwargs['legend_label'])] + tooltips
        fig.add_tools(bkm.tools.HoverTool(tooltips=tooltips,
                                          mode={'circle':'mouse',
                                                'vbar':'hline',
                                                'hbar':'vline'}[method],
                                          attachment=hover_attachment,
                                          names=[thisid]))
        if plot_domains is None:
            if hasattr(self, 'domains'):
                plot_domains = True
            else:
                plot_domains = False
        if plot_domains:
            if 'thickness' in yaxis:
                axes_domain['yaxis'] = None
            if 'thickness' in xaxis:
                axes_domain['xaxis'] = None
            if hasattr(self, 'domains'):
                kwargs = dict(line_color=color, line_width=2)
                kwargs.update(axes_domain)
                kwargs.update(add_domains_options)
                self.bokeh_add_domains(fig, **kwargs)
            else:
                print("WARNING: **plot_domains** is True while grid has no domains information; **plot_domains** ignored.")
        fig.legend.location = legend_location
        return fig

    def bokeh_add_domains(self, fig, yaxis=None, xaxis=None,
                          **span_options):
        """
        Add domains limits to figure.

        :param fig: the figure on which to plot.
        :param yaxis: 'level', 'p' (pressure), 'z' (height) or 'pthickness'//'zthickness'
        :param xaxis: 'level', 'p' (pressure), 'z' (height) or 'pthickness'//'zthickness'

        Other options are passed to bokeh's Span to plot domains limits.
        """
        dll = self.domains_lower_limit()
        doms = []
        border_line_color = span_options.get('line_color', 'black')
        border_line_width = span_options.get('line_width', 2)
        for d, lim in dll.items():
            if yaxis is not None:
                # horizontal lines
                if yaxis in lim:
                    doms.append(bkm.Span(location=lim[yaxis], dimension='width',
                                         **span_options))
                    fig.add_layout(bkm.Label(x=1, y=lim[yaxis], text=d,
                                             text_align='left',x_units='screen',
                                             border_line_color=border_line_color,
                                             border_line_width=border_line_width))
            if xaxis is not None:
                # vertical lines
                if xaxis in lim:
                    text_align = {'level':'right', 'z':'left', 'p':'right'}[xaxis]
                    doms.append(bkm.Span(location=lim[xaxis], dimension='height',
                                         **span_options))
                    fig.add_layout(bkm.Label(y=0, x=lim[xaxis], text=d,
                                             text_align=text_align, y_units='screen',
                                             border_line_color=border_line_color,
                                             border_line_width=border_line_width))
        fig.renderers.extend(doms)
        return fig

    def bokeh_fig_to_html(self, fig, htmlname=None, show_figure=True):
        """
        Outputs figure to html file, and optionally display it.

        :param htmlname: name of html file; if None, name of file
            is auto-generated from self.name
        :param show_figure: to show the plot in a WebBrowser tab
        """
        if htmlname is None:
            htmlname = self.name + '.html'
        bkio.output_file(htmlname, title='VGrid')
        if show_figure:
            bkio.show(fig)
        else:
            bkio.save(fig)
