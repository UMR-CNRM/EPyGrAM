#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2016-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import json
import datetime
import copy
import web

import matplotlib
matplotlib.use("Agg")

#from . import web

import epygram


####################
# Epyweb functions ############################################################
####################
def whichFields(fichier):
    """List and returns fields of selected file"""

    try:
        resource = epygram.formats.resource(fichier, 'r')
        listoffields = resource.listfields()
        #Cas de Surfex : on ne lit que les champs de type H2D (plus long que listfields)
        if resource.format == "FA" and listoffields[0].startswith("SFX."):
            listoffields = resource.find_fields_in_resource(fieldtype="H2D") #resource.listfields()
        resource.close()
        return listoffields
    except ValueError:
        raise Exception('whichField error')


def getAjaxArg(sArg, sDefault=''):
    """Picks out and returns a single value, regardless of GET or POST."""

    try:
        data = web.data()
        dic = None
        if data:
            dic = json.loads(data)
        else:
            dic = dict(web.input())
        if dic:
            if sArg in dic:
                if dic[sArg]:
                    return dic[sArg]
                else:
                    return sDefault
            else:
                return sDefault
        else:
            return sDefault
    except ValueError:
        raise Exception('getAjaxArg - no JSON arguments to decode. This method required a POST with JSON arguments.')


def getAjaxArgSmart(sArg, sDefault=''):
    """
    Picks out and returns a single value, regardless of GET or POST.
    Convert "" to None and true / false strings to True False, and string with starting { to dict
    USELESS when correct use of Javascript is done !!!!!!!!!!!!!!!!!!!!!!
    """

    try:
        data = web.data()
        dic = None
        if data:
            dic = json.loads(data)
        else:
            dic = dict(web.input())
        if dic:
            if sArg in dic:
                if dic[sArg] == '':
                    return None
                elif dic[sArg] == 'false':
                    return False
                elif dic[sArg] == 'true':
                    return True
                elif isinstance(dic[sArg], str) and dic[sArg].startswith('{'):
                    return json.loads(dic[sArg])
                else:
                    return dic[sArg]
            else:
                return sDefault
        else:
            return sDefault
    except ValueError:
        raise Exception('getAjaxArg - no JSON arguments to decode. This method required a POST with JSON arguments.')


def print_code(myplot_args, existingbasemap, existingfigure):
    try:
        zecode = ', '.join(["field.plotfield(subzone=''",
                            "gisquality='" + unicode(myplot_args["gisquality"]) + "'",
                            "specificproj='" + unicode(myplot_args["specificproj"]) + "'",
                            "minmax='" + unicode(myplot_args["minmax"]) + "'",
                            "graphicmode='" + unicode(myplot_args["graphicmode"]) + "'",
                            "levelsnumber=" + unicode(myplot_args["levelsnumber"]),
                            "colormap='" + unicode(myplot_args["colormap"]) + "'",
                            "center_cmap_on_0=" + unicode(myplot_args["center_cmap_on_0"]),
                            "drawrivers=" + unicode(myplot_args["drawrivers"]),
                            "meridians=" + unicode(myplot_args["meridians"]),
                            "parallels=" + unicode(myplot_args["parallels"]),
                            "bluemarble=" + unicode(myplot_args["bluemarble"]),
                            "departments=" + unicode(myplot_args["departments"]),
                            "existingfigure=" + unicode(existingfigure),
                            "pointsize=" + unicode(myplot_args["pointsize"])
                            ])
        return zecode
    except Exception as ex:
        print("Error print_code !")
        print(str(ex))


def get_common_args(fileid):
    myplot_args = {}
    mini = getAjaxArgSmart('min')
    maxi = getAjaxArgSmart('max')
    mypointsize = getAjaxArgSmart('pointsize')[fileid]

    myplot_args["minmax"] = (mini[fileid], maxi[fileid])
    #myplot_args["levelsnumber"] = getAjaxArgSmart('levelsnumber')[fileid]  # nombre/slidebar ?
    myplot_args["colorsnumber"] = getAjaxArgSmart('levelsnumber')[fileid]  # nombre/slidebar ?    
    myplot_args["colormap"] = getAjaxArg('colormap')[fileid]  # un menu déroulant avec des mini-images de chaque colormap ?
    #myplot_args["graphicmode"] = getAjaxArgSmart('graphicmode')[fileid]  # cases radiobutton [colorshades,contourlines,points]

    myplot_args["plot_method"] = getAjaxArgSmart('graphicmode')[fileid]  # cases radiobutton [colorshades,contourlines,points]
    if myplot_args["plot_method"] == "colorshades":
        myplot_args["plot_method"] = "contourf"
    if myplot_args["plot_method"] == "contourlines":
        myplot_args["plot_method"] = "contour"
    if myplot_args["plot_method"] == "points":
        myplot_args["plot_method"] = "scatter" #DOES NOT WORK!
    
    
    
    if (mypointsize != ""):
        myplot_args["pointsize"] = mypointsize  # nombre/slidebar ?
    myplot_args["contourcolor"] = getAjaxArgSmart('contourcolor')[fileid]
    #myplot_args["vectorcolor"] = getAjaxArgSmart('vectorcolor')[fileid]
    myplot_args["subzone"] = getAjaxArgSmart('subzone')[fileid]  # cases radiobutton [C,CI,CIE]
    myplot_args["gisquality"] = getAjaxArgSmart('gisquality')  # pour régler la finesse des traits de côte : cases radiobutton [c,l,i,h,f]
    myplot_args["specificproj"] = getAjaxArgSmart('specificproj')  # pour utiliser une projection de la carte particulière [kav7,cyl,ortho,nsperXXXX]
    myplot_args["center_cmap_on_0"] = getAjaxArgSmart('center_cmap_on_0')  # checkbox
    myplot_args["drawrivers"] = getAjaxArgSmart('drawrivers')  # checkbox
    myplot_args["meridians"] = getAjaxArgSmart('meridians')  # nombre/slidebar
    myplot_args["parallels"] = getAjaxArgSmart('parallels')
    myplot_args["bluemarble"] = getAjaxArgSmart('bluemarble')
    myplot_args["epygram_departments"] = getAjaxArgSmart('departments')  # checkbox

    return myplot_args


def datex(start, end=None, step=None):
    """
    Extended date expansion : YYYYMMDDHH-YYYYMMDDHH-HH
    """
    rangevalues = list()
    arguments = start.decode().split('-')
    start_arg = arguments[0]

    if len(arguments) == 1:
        end_arg = start_arg
        delta_arg = 24
    elif len(arguments) == 2:
        end_arg = arguments[1]
        delta_arg = 24
    elif len(arguments) == 3:
        end_arg = arguments[1]
        delta_arg = arguments[2]
    else:
        print("Uncorrect date range")

    start_date = arg2date(start_arg)
    end_date = arg2date(end_arg)
    delta = datetime.timedelta(days=0, seconds=int(delta_arg) * 3600)

    d = start_date
    while d <= end_date:
        rangevalues.append(d.strftime("%Y%m%d%H"))
        d += delta

    return (rangevalues)


def arg2date(myarg):
    out = datetime.datetime(int(myarg[0:4]),
                            int(myarg[4:6]),
                            int(myarg[6:8]),
                            int(myarg[8:10]))
    return out

def make_my_plot2(resource, field, cle, champ, champ_v, FF, vecteurs,
                 over, vectors_subsampling, myplot_args, monzoom = None):
    """
    Generic method for plotting (ok for plot, plot_both, overlay, but not for diff).
    """

    '''
    try:
        vectorcolor = myplot_args["vectorcolor"]
        myplot_args.pop("vectorcolor", None)
    except:
        vectorcolor = "black"
    '''
    print("MakeMyPlot starts ")

    # On enleve des options en cas de tracé de vecteur
    myplot_args_vect = copy.copy(myplot_args)
    myplot_args_vect.pop("pointsize", None)
    myplot_args_vect.pop("colormap", None)
    #myplot_args_vect.pop("levelsnumber", None)
    #myplot_args_vect.pop("graphicmode", None)
    myplot_args_vect.pop("minmax", None)
    myplot_args_vect.pop("center_cmap_on_0", None)
    #myplot_args_vect.pop("departments", None)
    myplot_args_vect.pop("contourcolor", None)
    myplot_args_vect.pop("colorbar", None)

    #TMP GF DEBUG
    myplot_args.pop("pointsize", None)
    myplot_args.pop("gisquality", None)
    myplot_args.pop("specificproj", None)
    myplot_args.pop("drawrivers", None)
    myplot_args.pop("bluemarble", None)

    if FF[cle] or vecteurs[cle]:
        field_v = resource.readfield(champ_v[cle])
        if field_v.spectral:
            field_v.sp2gp()
        if monzoom is not None:
            #pass
            #field = field.extract_zoom(monzoom)
            field_v = field_v.extract_zoom(monzoom)
        vectwind = epygram.fields.make_vector_field(field,
                                                    field_v)
        if FF[cle]:
            FF_field = vectwind.to_module()
            myplot = FF_field.cartoplot(title=str(champ[cle]) + '\n' + str(field.validity.get()),**myplot_args)
            if vecteurs[cle]:
                myplot = vectwind.plotfield(over=myplot,
                                            title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                            subsampling=vectors_subsampling[cle],
                                            plot_module=False,
                                            symbol_options={'color': vectorcolor},
                                            **myplot_args_vect)
            del field_v
            del vectwind
        elif vecteurs[cle]:
            myplot = vectwind.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                        over=over,
                                        subsampling=vectors_subsampling[cle],
                                        plot_module=False,
                                        symbol_options={'color': vectorcolor},
                                        **myplot_args_vect)
            del field_v
            del vectwind
    else:

        myplot = field.cartoplot(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                 **myplot_args)

    return myplot


'''
def make_my_plot(resource, field, cle, champ, champ_v, FF, vecteurs,
                 use_basemap, over, vectors_subsampling, myplot_args, monzoom = None):
    """
    Generic method for plotting (ok for plot, plot_both, overlay, but not for diff).
    """

    try:
        vectorcolor = myplot_args["vectorcolor"]
        myplot_args.pop("vectorcolor", None)
    except:
        vectorcolor = "black"

    print("MakeMyPlot starts ")

    # On enleve des options en cas de tracé de vecteur
    myplot_args_vect = copy.copy(myplot_args)
    myplot_args_vect.pop("pointsize", None)
    myplot_args_vect.pop("colormap", None)
    #myplot_args_vect.pop("levelsnumber", None)
    #myplot_args_vect.pop("graphicmode", None)
    myplot_args_vect.pop("minmax", None)
    myplot_args_vect.pop("center_cmap_on_0", None)
    #myplot_args_vect.pop("departments", None)
    myplot_args_vect.pop("contourcolor", None)
    myplot_args_vect.pop("colorbar", None)

    if FF[cle] or vecteurs[cle]:
        field_v = resource.readfield(champ_v[cle])
        if field_v.spectral:
            field_v.sp2gp()
        if monzoom is not None:
            #pass
            #field = field.extract_zoom(monzoom)
            field_v = field_v.extract_zoom(monzoom)
        vectwind = epygram.fields.make_vector_field(field,
                                                    field_v)
        if FF[cle]:
            FF_field = vectwind.to_module()
            myplot = FF_field.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                        use_basemap=use_basemap,
                                        over=over,
                                        **myplot_args)
            if vecteurs[cle]:
                myplot = vectwind.plotfield(over=myplot,
                                            title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                            use_basemap=use_basemap,
                                            subsampling=vectors_subsampling[cle],
                                            plot_module=False,
                                            symbol_options={'color': vectorcolor},
                                            **myplot_args_vect)
            del field_v
            del vectwind
        elif vecteurs[cle]:
            myplot = vectwind.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                        use_basemap=use_basemap,
                                        over=over,
                                        subsampling=vectors_subsampling[cle],
                                        plot_module=False,
                                        symbol_options={'color': vectorcolor},
                                        **myplot_args_vect)
            del field_v
            del vectwind
    else:
        myplot = field.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                 use_basemap=use_basemap,
                                 over=over,
                                 **myplot_args)

    return myplot
'''

def check_for_operation(ope, field):
    try:
        operation_arg1 = ope[0].encode()
        print("operation ", operation_arg1)
        try:
            operation_arg2 = float(ope[1])
            print("operation2 ", operation_arg2)
            field.operation(operation_arg1, operation_arg2)
        except:
            field.operation(operation_arg1)
    except:
        print("no operation...")

    return field


#########################
# Environment functions #######################################################
#########################
def func_open_browser(url, delay=0.):
    """Open a browser tab at **url**, after the optional **delay**."""
    import webbrowser
    import time
    if delay > 0.:
        time.sleep(delay)
    webbrowser.open(url)


def clean_workdir(wkdir):
    """Cleaning of old figures and hardlinks."""
    import shutil
    shutil.rmtree(wkdir, ignore_errors=True)


def init_workdir(wkdir):
    """Set up the working directory."""
    import os
    clean_workdir(wkdir)
    if not os.path.exists(wkdir):
        os.makedirs(wkdir)


#####################
# Enhanced features ###########################################################
#####################
class PortApplication(web.application):
    """Web server on specific port."""
    def run(self, port=8080, *middleware):
        func = self.wsgifunc(*middleware)
        return web.httpserver.runsimple(func, ('0.0.0.0', port))
