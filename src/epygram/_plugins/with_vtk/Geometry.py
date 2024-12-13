#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains the classes that handle the vtk plotting on geometries.
"""

import numpy
import os
import tempfile
import math

from footprints import proxy as fpx

from epygram.util import Angle, get_file, epygramError
from epygram.geometries import AcademicGeometry, VGeometry, RegLLGeometry, ProjectedGeometry


def activate():
    """Activate extension."""
    from . import __name__ as plugin_name
    from epygram._plugins.util import notify_doc_requires_plugin
    notify_doc_requires_plugin([make_vtkGrid, plot3DProjImage, plot3DBluemarble,
                                plot3DMaptiles],
                               plugin_name)
    from epygram.geometries import Geometry, _need_pyproj_geod
    Geometry.make_vtkGrid = _need_pyproj_geod(make_vtkGrid)
    Geometry.plot3DProjImage = plot3DProjImage
    Geometry.plot3DBluemarble = plot3DBluemarble
    Geometry.plot3DMaptiles = plot3DMaptiles


def make_vtkGrid(self, rendering, subzone=None):
    """
    Makes an empty grid to use with vtk

    :param rendering: a usevtk.Usevtk instance
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    """
    import vtk  # @UnresolvedImport
    from vtk.numpy_interface import dataset_adapter as dsa  # @UnresolvedImport
    from vtk.numpy_interface import algorithms as algs  # @UnresolvedImport

    if rendering.hCoord == 'geoid' and not self.vcoordinate.typeoffirstfixedsurface in (102, 103):
        raise ValueError("The 'geoid' option is only meaningful with " + \
                         "vertical coordinate expressed in meters!")
        
    # Compute x, y, z
    lons, lats = self.get_lonlat_grid(d4=True, nb_validities=1, subzone=subzone)
    z = self.get_levels(d4=True, nb_validities=1, subzone=subzone)[0, ...]
    shape = z.shape
    total_mask = None
    if any([isinstance(arr, numpy.ma.masked_array) for arr in (lons, lats, z)]):
        total_mask = numpy.ma.getmaskarray(lons.flatten())
        for arr in (lats, z):
            total_mask = numpy.logical_or(total_mask, numpy.ma.getmaskarray(arr.flatten()))
    lons, lats, z = (arr.filled() if isinstance(arr, numpy.ma.masked_array) else arr for arr in (lons, lats, z))
    x, y, z = rendering.proj3d(lons, lats, z)

    # sgrid_point
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(*shape[::-1])
    points = vtk.vtkPoints()
    x, y, z = [a.astype(numpy.float32) for a in (x, y, z)]
    coordinates = algs.make_vector(x, y, z.flatten())
    points.SetData(dsa.numpyTovtkDataArray(coordinates, None, array_type=vtk.VTK_FLOAT))
    grid.SetPoints(points)
    if total_mask is not None:
        for index in numpy.nonzero(total_mask)[0]:
                grid.BlankPoint(index)
    return grid


def plot3DProjImage(self, rendering,
                    filename, geometry, subzone=None,
                    interpolation='nearest', color_transform=None):
    """
    This method adds an image to the vtk rendering system. The image is projected
    on the vertical coordinate of this geometry.

    :param rendering: a usevtk.Usevtk instance
    :param filename: path to an image file
    :param geometry: projection associated to the geometry
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :param interpolation: interpolation method to use to compute levels at each pixel
    :param color_transform: if not None, must be a user-defined function with parameters
                            (rgb, alpha, (lons, lats)) with rgb an array representing
                            the rgb values on each of the image points, alpha the alpha
                            values on the same points and (lons, lats) a turple of the
                            points coordinate. The function must return (rgb, alpha). This
                            enables the color modification of the image.
                            For example, the following function transforms the color
                            into grey levels:
        def grey(rgb, alpha, ll):
            rgb[:, :, :] = rgb.mean(axis=2)[:, :, numpy.newaxis]
            return rgb, alpha
    first_corner and last_corner describe the coverage of the image from the left edge
    of the leftmost pixel to the right edge of the rightmost pixel.
    """
    from PIL import Image
    import vtk  # @UnresolvedImport
    from vtk.numpy_interface import dataset_adapter as dsa  # @UnresolvedImport

    # Image reading
    with Image.open(filename) as image:
        rgb = numpy.array(image.convert('RGB'))[::-1, ...]
    alpha = numpy.zeros_like(rgb[:, :, 0])  # not visible by default

    lons, lats = geometry.get_lonlat_grid()

    # We hide non-visible parts
    m = 1 if interpolation == 'linear' else 0
    inside = numpy.array(self.point_is_inside_domain_ll(lons.flatten(),
                                                        lats.flatten(),
                                                        margin=m,
                                                        subzone=subzone)).reshape(lons.shape)
    if numpy.count_nonzero(inside) == 0:
        return  (None, None)  # image is completely out of the domain
    alpha[inside] = 255  # Points inside domain are set visible

    # color transformation
    if color_transform is not None:
        rgb, alpha = color_transform(rgb, alpha, (lons, lats))

    # Projection of levels on the image geometry
    level_field = self.vcoord_as_field(2)
    z_inside = level_field.getvalue_ll(lons[inside], lats[inside],
                                       k=0, interpolation=interpolation)
    if isinstance(z_inside, numpy.ma.masked_array):
        z = numpy.ma.zeros(lons.shape, fill_value=z_inside.fill_value)
    else:
        z = numpy.zeros(lons.shape)
    z[inside] = z_inside
    del z_inside
    if isinstance(z, numpy.ma.masked_array):
        alpha[numpy.ma.getmaskarray(z)] = 0
        inside[numpy.ma.getmaskarray(z)] = False

    # Extrapolation to suppress vtk artifacts due to interpolation between
    # invisible and visible parts of the domain
    inside_extrapolation = inside.copy()
    for _ in range(1):  # currently only one point but code is ready to extrapolate on more points
        mask = numpy.logical_and(inside_extrapolation[1:, :], numpy.logical_not(inside_extrapolation[:-1, :]))
        j, i = numpy.where(mask)
        z[j, i] = z[j + 1, i]
        inside_extrapolation[j, i] = inside_extrapolation[j + 1, i]
        mask = numpy.logical_and(inside_extrapolation[:-1, :], numpy.logical_not(inside_extrapolation[1:, :]))
        j, i = numpy.where(mask)
        z[j + 1, i] = z[j, i]
        inside_extrapolation[j + 1, i] = inside_extrapolation[j, i]
        mask = numpy.logical_and(inside_extrapolation[:, 1:], numpy.logical_not(inside_extrapolation[:, :-1]))
        j, i = numpy.where(mask)
        z[j, i] = z[j, i + 1]
        inside_extrapolation[j, i] = inside_extrapolation[j, i + 1]
        mask = numpy.logical_and(inside_extrapolation[:, :-1], numpy.logical_not(inside_extrapolation[:, 1:]))
        j, i = numpy.where(mask)
        z[j, i + 1] = z[j, i]
        inside_extrapolation[j, i + 1] = inside_extrapolation[j, i]
    del inside_extrapolation

    # We cut the useless part to reduce the grid size to render
    j, i = numpy.where(inside)
    imin, imax, jmin, jmax = i.min(), i.max(), j.min(), j.max()
    rgb = rgb[jmin:jmax + 1, imin:imax + 1, :]
    alpha = alpha[jmin:jmax + 1, imin:imax + 1]
    z = z[jmin:jmax + 1, imin:imax + 1]

    # VTK grid
    new_geometry = geometry.make_subarray_geometry(imin, imax, jmin, jmax)
    new_geometry.vcoordinate.levels[0] = z
    
    grid = new_geometry.make_vtkGrid(rendering)

    # Filling grid
    if numpy.all(alpha == 255):
        data = rgb
        num_components = 3
    else:
        # use alpha-channel only when really needed
        data = numpy.append(rgb, alpha[..., numpy.newaxis], axis=2)
        num_components = 4
    arr = dsa.numpyTovtkDataArray(data.flatten(), "CellID")
    arr.SetNumberOfComponents(num_components)
    grid.GetPointData().SetScalars(arr)

    # Render image
    sgridGeom = vtk.vtkStructuredGridGeometryFilter()
    sgridGeom.SetInputData(grid)
    sgridGeomMap = vtk.vtkPolyDataMapper()
    sgridGeomMap.SetInputConnection(sgridGeom.GetOutputPort())
    sgridGeomMapActor = vtk.vtkActor()
    sgridGeomMapActor.SetMapper(sgridGeomMap)
    rendering.renderer.AddActor(sgridGeomMapActor)
    # sgridGeomMapActor.GetProperty().SetAmbient(1);
    # sgridGeomMapActor.GetProperty().SetDiffuse(0);
    # sgridGeomMapActor.GetProperty().SetSpecular(0)
    return (sgridGeomMapActor, sgridGeomMap)


def plot3DBluemarble(self, rendering,
                     subzone=None,
                     interpolation='nearest', color_transform=None):
    """
    This method adds the NOAA's bluemarble image to the vtk rendering system. The image
    is projected on the vertical coordinate of this geometry.

    :param rendering: a usevtk.Usevtk instance
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :param interpolation: interpolation method to use to compute levels at each pixel
    :param color_transform: if not None, must be a user-defined function with parameters
                            (rgb, alpha, (lons, lats)) with rgb an array representing
                            the rgb values on each of the image points, alpha the alpha
                            values on the same points and (lons, lats) a turple of the
                            points coordinate. The function must return (rgb, alpha). This
                            enables the color modification of the image.
                            For example, the following function transforms the color
                            into grey levels:
        def grey(rgb, alpha, ll):
            rgb[:, :, :] = rgb.mean(axis=2)[:, :, numpy.newaxis]
            return rgb, alpha

    Image can be found here:
        https://sos.noaa.gov/datasets/blue-marble-without-clouds/
        ftp://public.sos.noaa.gov/land/blue_marble/earth_vegetation/4096.jpg
    """
    from PIL import Image

    if isinstance(self, AcademicGeometry):
        raise epygramError("We cannot plot lon-lat images in an academic projection")

    path = tempfile.mkstemp()[1]
    get_file("ftp://public.sos.noaa.gov/land/blue_marble/earth_vegetation/4096.jpg",
             path, authorize_cache=True, subst=None)
    # Compute the geometry
    with Image.open(path) as image:
        size = image.size
    resol = 360. / size[0], 180. / size[1]
    grid = {'input_lon':Angle(-180. + resol[0] / 2., 'degrees'),
            'input_lat':Angle(-90. + resol[1] / 2., 'degrees'),
            'input_position':(0, 0),
            'X_resolution':Angle(resol[0], 'degrees'),
            'Y_resolution':Angle(resol[1], 'degrees')}
    kwargs_vcoord = {'typeoffirstfixedsurface': self.vcoordinate.typeoffirstfixedsurface,
                     'position_on_grid': self.vcoordinate.position_on_grid,
                     'grid': self.vcoordinate.grid,
                     'levels': [0]
                     }
    kwargs_geom = dict(name='regular_lonlat',
                       grid=grid,
                       dimensions=dict(X=size[0], Y=size[1]),
                       vcoordinate=VGeometry(**kwargs_vcoord),
                       position_on_horizontal_grid='center',
                       geoid=self.geoid)
    geometry = RegLLGeometry(**kwargs_geom)

    result = self.plot3DProjImage(rendering, path, geometry, subzone,
                                  interpolation, color_transform=color_transform)
    os.remove(path)
    return result


def plot3DMaptiles(self, rendering, url, resol_factor,
                   minzoom=1, maxzoom=18,
                   subzone=None,
                   interpolation='nearest',
                   color_transform=None):
    """
    This method adds tiles image to the vtk rendering system. The image
    is projected on the vertical coordinate of this geometry.

    :param rendering: a usevtk.Usevtk instance
    :param url: url to get the map tiles to use
                ex: https://a.tile.openstreetmap.org/${z}/${x}/${y}.png
                where ${z}, ${x} and ${y} are place holders for
                zoom, x and y position of the tile
    :param resol_factor: how many times the tile resolution must be approximately
                         better than the geometry resolution. A factor of 2 (resp. 1/2.)
                         induce using the next (resp. preceding) zoom level.
    :param minzoom/maxzoom: minimum and maximum zoom levels to use
    :param subzone: optional, among ('C', 'CI'), for LAM grids only, returns
                    the grid resp. for the C or C+I zone off the C+I+E zone. \n
                    Default is no subzone, i.e. the whole field.
    :param interpolation: interpolation method to use to compute levels at each pixel
    :param color_transform: if not None, must be a user-defined function with parameters
                            (rgb, alpha, (lons, lats)) with rgb an array representing
                            the rgb values on each of the image points, alpha the alpha
                            values on the same points and (lons, lats) a turple of the
                            points coordinate. The function must return (rgb, alpha). This
                            enables the color modification of the image.
                            For example, the following function transforms the color
                            into grey levels:
        def grey(rgb, alpha, ll):
            rgb[:, :, :] = rgb.mean(axis=2)[:, :, numpy.newaxis]
            return rgb, alpha

    Zoom computation is done supposing that each tile is made of 256 * 256 pixels.

    URL examples (these URLs have been found on internet but some are protected with
    copyrights and must not be used with this tool or kept in cache; please check
    before use):
        - osm: a (full?) list can be found on https://wiki.openstreetmap.org/wiki/Tile_servers
               or https://wiki.openstreetmap.org/wiki/Tiles
               The most common one is certainly: https://a.tile.openstreetmap.org/${z}/${x}/${y}.png
        - google maps: it seems to be forbidden to use this tiles outside of google products.
    """
    from PIL import Image
    
    if isinstance(self, AcademicGeometry):
        raise epygramError("We cannot plot lon-lat images in an academic projection")

    def deg2num(lon, lat, zoom):
        """
        Returns x and y of the tile containing the point (lon, lat) in degrees
        """
        lat_rad = math.radians(lat)
        n = 2.0 ** zoom
        xtile = int((lon + 180.0) / 360.0 * n)
        ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
        return (xtile, ytile)

    def num2deg(xtile, ytile, zoom):
        """
        Returns the upper left corner coordinates (in degrees)
        corresponding to the (x, y) tile
        """
        n = 2.0 ** zoom
        lon_deg = xtile / n * 360.0 - 180.0
        lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
        lat_deg = math.degrees(lat_rad)
        return (lon_deg, lat_deg)

    radius = 6378137.   # found on the web but with (normally) no impact
                        # as long as we use the same for the geoid
    def zoom2resol(zoom):
        """
        Returns width of a tile at equator (from zoom)
        To be divided by the number of pixel to obtain the resolution
        of the pixel
        """
        return 2 * math.pi * radius / 2**zoom
    
    # Looking for bounds
    lons, lats = self.get_lonlat_grid(subzone=subzone)
    lonmin, lonmax, latmin, latmax = lons.min(), lons.max(), lats.min(), lats.max()
    surface = (latmax - latmin) * (lonmax - lonmin) / lons.size  # mean pixel surface in geometry
    surface = surface / resol_factor ** 2  # mean pixel surface wanted for the tile
    surface = surface * 256 * 256  # mean tile surface wanted
    zoom = max(minzoom, min(maxzoom, int(math.log(360 * 180 / surface) / (2 * math.log(2)))))  # zoom level
    del lons, lats
    xmin, ymin = deg2num(lonmin, latmax, zoom)
    xmax, ymax = deg2num(lonmax, latmin, zoom)

    # Building complete image from tiles
    path = tempfile.mkstemp()[1]
    size = None
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            get_file(url, path, authorize_cache=True,
                     subst={'${z}': zoom, '${x}': x, '${y}': y})
            small_img = Image.open(path)
            if size is None:
                size = small_img.size
                new_im = Image.new('RGB', ((xmax + 1 - xmin) * size[0], (ymax + 1 - ymin) * size[1]))
            else:
                assert size == small_img.size, "All tiles must share the same size"
            new_im.paste(small_img, ((x - xmin) * size[0], (y - ymin) * size[1]))
    os.remove(path)
    path += '.png'
    new_im.save(path)
    new_im.close()
    
    # Building geometry corresponding to the image
    inputpos = num2deg(xmin, ymin, zoom) #position of the North-West image corner
    resol = zoom2resol(zoom) / size[0], zoom2resol(zoom) / size[1]
    
    # an offset is needed to go from the image corner to the pixel center
    grid = {'input_lon':Angle(inputpos[0], 'degrees'),
            'input_lat':Angle(inputpos[1], 'degrees'),
            'input_position':(0 - .5, (ymax + 1 - ymin) * size[1] + .5),
            'X_resolution':resol[0],
            'Y_resolution':resol[1],
            'LAMzone':None}
    kwargs_vcoord = {'typeoffirstfixedsurface': self.vcoordinate.typeoffirstfixedsurface,
                     'position_on_grid': self.vcoordinate.position_on_grid,
                     'grid': self.vcoordinate.grid,
                     'levels': [0]
                     }
    projection = {'reference_lon':Angle(0, 'radians'),
                  'reference_lat':Angle(0, 'radians'),
                  'rotation':Angle(0., 'radians')}
    kwargs_geom = dict(name='mercator',
                       grid=grid,
                       dimensions=dict(X=(xmax + 1 - xmin) * size[0], Y=(ymax + 1 - ymin) * size[1]),
                       projection=projection,
                       vcoordinate=VGeometry(**kwargs_vcoord),
                       position_on_horizontal_grid='center',
                       geoid={'a':radius, 'b':radius})
    geometry = ProjectedGeometry(**kwargs_geom)

    # Rendering
    result = self.plot3DProjImage(rendering, path, geometry, subzone,
                                interpolation, color_transform=color_transform)
    os.remove(path)
    
    return result

