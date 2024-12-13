#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2016-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

import os
import json
import uuid
import web

import matplotlib
matplotlib.use("Agg")
#from mpl_toolkits.basemap import Basemap

#from . import web

from footprints.util import rangex
import epygram

from . import epyweb_workdir, vortex_cache, render, basemap_pickle_path, all_fatal_exceptions
from . import util
from .util import getAjaxArg, getAjaxArgSmart, check_for_operation, get_common_args, make_my_plot2

__all__ = ['Epyweb',
           'GetCacheSize', 'GetGeometries', 'GetFieldsAsJSON', 'GetLocalFile',
           'GetFile', 'GetMinMax', 'GetDomain', 'GetPNG',
           'MyPlot', 'MyPlotOverlay', 'MyPlotDiff']


####################
# urls <=> classes ############################################################
####################
class Epyweb(object):
    def GET(self):
        return render.epyweb()


class GetCacheSize(object):
    """Compute (linux only!) and return size of vortex cache"""
    def POST(self):
        try:
            cacheSize = os.popen("du -kshx " + vortex_cache).read()
            return json.dumps(cacheSize)
        except Exception:
            if all_fatal_exceptions:
                raise
            return "Error in cache size retrieval"


class GetGeometries(object):
    """Return the list of existing geometries"""
    def POST(self):
        import usevortex
        toremove_geoms = ['assmp1', 'assmp1sp',
                          'assmp2', 'assmp2sp',
                          'assms1', 'assms1sp',
                          'assms2', 'assms2sp']
        geoms = usevortex.list_vortex_geometries()
        for g in toremove_geoms:
            geoms.remove(g)
        return json.dumps(geoms)


class GetFieldsAsJSON(object):
    """List and returns fields of selected file"""
    def POST(self):
        try:
            fichier = getAjaxArg('file')
            print("fichier = ", fichier)
            if (os.path.isfile(fichier)):
                return json.dumps(util.whichFields(fichier))
            else:
                print("File does not exist => exit")
        except Exception:
            if all_fatal_exceptions:
                raise
            print("Erreur getfieldsasjson")
            return "Erreur getfieldsasjson"


class GetLocalFile(object):
    """Retrieve selected file(s) from path"""

    def POST(self):
        try:
            reponse = {}
            data = web.data()
            args = json.loads(data)
            if not os.path.exists(args['file_path']):
                raise ValueError('file does not exist')
            reponse['localpath'] = args['file_path']
            return reponse
        except ValueError:
            if all_fatal_exceptions:
                raise
            raise Exception('getFile error')


class GetFile(object):
    """Retrieved selected file(s) with usevortex"""

    def POST(self):
        import usevortex
        try:
            reponse = {}
            data = web.data()
            vortexArgs = json.loads(data)

            # Patch pour pbs unicode + utilisation de rangex
            if 'date' in vortexArgs:
                vortexArgs['date'] = util.datex(vortexArgs['date'].encode())
            if 'term' in vortexArgs:
                vortexArgs['term'] = rangex(vortexArgs['term']) #.encode())
            if 'month' in vortexArgs:
                vortexArgs['month'] = rangex(vortexArgs['month'].encode())
            if 'member' in vortexArgs:
                memberstr = '_mb[member]'
                vortexArgs['member'] = rangex(vortexArgs['member']) #.encode())
            else:
                memberstr = ''

            # On ajoute POUR L'INSTANT et par défaut origin=hst (pour les gridpoints)
            vortexArgs['origin'] = 'hst'
            # On enlève le mode demandé du dictionnaire d'arguments : description, existence, get
            mode = vortexArgs['request_mode']
            del vortexArgs['request_mode']

            # Utile pour garder trace du fichier, A ou B, d'origine
            try:
                fromid = vortexArgs['fromid']
                del vortexArgs['fromid']
            except Exception:
                fromid = "A"
            
            # Nom local; cannot use shouldfly, for directory/rights reasons
            localname = os.path.join(epyweb_workdir,
                                     '.'.join(["[date::ymdh]",
                                               "[term]",
                                               memberstr,
                                               fromid,
                                               str(uuid.uuid4())
                                               ]))
            vortexArgs.update(local=localname)
                        
            # Test de complétude de la description
            ressources = usevortex.get_resources(getmode='check',
                                                 **vortexArgs)
            reponse['description'] = [str(ressources)]

            # Si oui allons plus loin
            if (ressources and mode != 'description'):
                # Chemin
                ressources = usevortex.get_resources(getmode='exist',
                                                     **vortexArgs)
                reponse['remotepath'] = '\n'.join([str(exist) + ':\n' + loc.replace(';', '\n') for (loc, exist) in ressources])  # [m for m in ressources]
                if mode == 'existence':
                    # Existence physique : tableau de True False
                    ressources = usevortex.get_resources(getmode='exist',
                                                         **vortexArgs)
                    reponse['existence'] = all([False not in m for m in ressources])
                if mode == 'get':
                    # Rapatriement
                    ressources = usevortex.get_resources(getmode='fetch',
                                                         **vortexArgs)
                    reponse['localpath'] = [str(m) for m in ressources]  # str(m[0])
            return json.dumps(reponse)
        except ValueError:
            if all_fatal_exceptions:
                raise
            raise Exception('getFile error')


class GetMinMax(object):
    """Compute and return min and max of FIRST field"""

    def POST(self):
        try:
            fichier = getAjaxArg('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArg('field_v')
            subzone = getAjaxArg('subzone')
            FF = getAjaxArg('FF')
            ope = getAjaxArg('operation')
            # string vs unicode problems...
            try:
                champ['typeOfLevel'] = champ['typeOfLevel'].encode()
                champ = {str(k):champ[k] for k in champ.keys()}
                champ_v['typeOfLevel'] = champ['typeOfLevel'].encode()
                champ_v = {str(k):champ[k] for k in champ_v.keys()}
            except Exception:
                print("Erreur unicode")

            resource = epygram.formats.resource(fichier, 'r')
            stats = {}
            field = resource.readfield(champ)
            if not field.geometry.grid.get('LAMzone', False):
                subzone = None

            if field.spectral:
                field.sp2gp()
            field = check_for_operation(ope, field)

            if FF:
                field_v = resource.readfield(champ_v)
                if field_v.spectral:
                    field_v.sp2gp()
                vectwind = epygram.fields.make_vector_field(field, field_v)
                FF_field = vectwind.to_module()
                stats['min'] = FF_field.min(subzone=subzone)
                stats['max'] = FF_field.max(subzone=subzone)
                del field_v
                del FF_field
            else:
                stats['min'] = field.min(subzone=subzone)
                stats['max'] = field.max(subzone=subzone)
            resource.close()

            return json.dumps(stats)

        except Exception as ex:
            if all_fatal_exceptions:
                raise
            print(ex.__str__())


class GetDomain(object):
    """Compute and return domain caracteristics"""

    def POST(self):
        try:
            fichier = getAjaxArg('file')
            champ = getAjaxArg('champ')
            resource = epygram.formats.resource(fichier, 'r')

            # On prend la géométrie du 1er champ => compatibilité FA / GRIB
            firstfield = resource.readfield(champ)
            if firstfield.geometry.rectangular_grid:
                (llcrnrlon, llcrnrlat) = firstfield.geometry.gimme_corners_ll()['ll']
                (urcrnrlon, urcrnrlat) = firstfield.geometry.gimme_corners_ll()['ur']
                (ulcrnrlon, ulcrnrlat) = firstfield.geometry.gimme_corners_ll()['ul']
                (lrcrnrlon, lrcrnrlat) = firstfield.geometry.gimme_corners_ll()['lr']
            else:
                (llcrnrlon, llcrnrlat) = (-180, -90)
                (urcrnrlon, urcrnrlat) = (180, 90)

            monzoom = {'lonmin': min(llcrnrlon, ulcrnrlon),
                       'lonmax': max(lrcrnrlon, urcrnrlon),
                       'latmin': min(llcrnrlat, lrcrnrlat),
                       'latmax': max(urcrnrlat, ulcrnrlat),
                       }

            resource.close()
            del resource
            del firstfield

            return json.dumps(monzoom)

        except Exception as ex:
            if all_fatal_exceptions:
                raise
            print(ex.__str__())


class MyPlot(object):
    """Plot of a parameter for 1 or several file(s)"""

    def POST(self):
        try:
            import matplotlib.pyplot as plt
            print('Start MyPlot')

            # For unik name
            local_uuid = str(uuid.uuid4())

            files = getAjaxArgSmart('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            # string vs unicode problems...
            for cle, val in champ.items():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            for cle, val in champ_v.items():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            figax = (None, None)  # New figure (= no overlay)
            basemap_pickle_name = getAjaxArgSmart('basemap_pickle_name')

            monzoom = getAjaxArgSmart('monzoom')

            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            vectors_subsampling = getAjaxArgSmart('vectors_subsampling')

            getcode = getAjaxArgSmart('getcode')
            dpi = getAjaxArgSmart('dpi')

            current_basemap = None
            if os.path.exists(os.path.join(epyweb_workdir, basemap_pickle_name)):
                _current_basemap = pickle.load(open(os.path.join(epyweb_workdir,
                                                                 basemap_pickle_name), 'r'))
                if isinstance(_current_basemap, Basemap):
                    current_basemap = _current_basemap

            # Attention, decuml marche bien entre échéances mais actif aussi entre dates !!
            # Pas (encore ?) pour OVERLAY et DIFF
            decumul = getAjaxArgSmart('decumul')
            operation = getAjaxArgSmart('operation')

            # Algo : au 1er passage on ne trace rien, juste mis en mémoire du champ
            # ensuite on trace champLu - champAvant, si champAvant n'existe pas (cas RR @0h) -> champLu only
            out = {}

            # for fichier in files:
            for cle, val in files.items():
                indiceDecumul = 0
                liste_tmp = []
                myplot_args = get_common_args(cle)
                # Garde fou pour decumul
                if FF[cle] or vecteurs[cle]:
                    decumul[cle] = False
                    print("decumul forced to ", decumul)
                try:
                    myplot_args["meridians"] = myplot_args["meridians"].encode()
                    myplot_args["parallels"] = myplot_args["parallels"].encode()
                except Exception:
                    print("Warning unicode")

                Lexception = False
                for fichier in val:
                    resource = epygram.formats.resource(fichier, 'r')
                    print("FICHIER : " + fichier)
                    try:
                        field = resource.readfield(champ[cle])
                        if field.spectral:
                            field.sp2gp()
                        if cle in operation:
                            field = check_for_operation(operation[cle], field)
                        if decumul[cle] is True:
                            if indiceDecumul == 0:
                                fieldDecumul = field
                                indiceDecumul = +1
                                continue
                    except Exception:  # Cas des RR @0h : param n'existe pas
                        indiceDecumul = +1
                        Lexception = True  # A garder pour reinitialiser le cumul
                        continue

                    if decumul[cle]:
                            waitforme = field
                            try:  # Cas des RR @0h : param n'existe pas
                                validity = field.validity
                                fid = field.fid
                                if Lexception is False:  # cas normal
                                    field = field - fieldDecumul
                                else:  # cas juste après un champ inexistant : field=field et fin de l'exception
                                    Lexception = False
                                field.validity = validity
                                field.fid = fid
                            except Exception:
                                pass
                            fieldDecumul = waitforme

                    if current_basemap is None:
                        print('Actually niou pickle !')
                        if monzoom is not None:
                            pass
                            field = field.extract_zoom(monzoom)
                        #current_basemap = field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                        #                                                  subzone=myplot_args["subzone"],
                        #                                                  specificproj=myplot_args["specificproj"])
                        """
                        current_basemap = field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                                                                      subzone=myplot_args["subzone"],
                                                                      specificproj=myplot_args["specificproj"],
                                                                      zoom=monzoom)
                        """
                        #pickle.dump(current_basemap, open(os.path.join(epyweb_workdir,
                        #                                               basemap_pickle_name), 'w'))
                        # On réutilise le cache

                    #myplot = make_my_plot(resource, field, cle, champ, champ_v, FF, vecteurs,
                    #                      current_basemap, figax, vectors_subsampling, myplot_args, monzoom=monzoom)
                    myplot = make_my_plot2(resource, field, cle, champ, champ_v, FF, vecteurs,
                                          figax, vectors_subsampling, myplot_args, monzoom=monzoom)


                    if getcode:
                        pass
                        #zecode = util.print_code(myplot_args, current_basemap, figax)
                        #print("*** Arguments du plot ***\n" + zecode)

                    # Utilisation d'un nom unique par image, dans un répertoire fixe
                    try:
                        myunikname = os.path.basename(fichier) + "." + ".".join("=".join((str(k), str(v))) for k, v in champ[cle].items())
                    except Exception:
                        myunikname = os.path.basename(fichier) + "." + str(champ[cle].replace(' ', '_'))

                    # On rajoute un petit uuid en cas de rafraichissement d'image
                    myunikfile = os.path.join(epyweb_workdir,
                                              myunikname + "." + local_uuid + '.png')
                    print("Saving figure ", myunikfile)
                    myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                    # memory management
                    resource.close()
                    del resource
                    del field
                    plt.close(myplot[0])
                    del myplot

                    # SLIDE IMAGE STYLE
                    liste_tmp.append('/getPNG/' + os.path.basename(myunikfile))

                out[cle] = liste_tmp

            out2 = []

            del current_basemap

            # On alterne A et B si besoin (cas de plot_both)
            if len(out) == 1:
                out2 = out["A"]
            else:
                for idx, val in enumerate(out["A"]):
                    out2.append(val)
                    out2.append(out["B"][idx])

            print("End Plot")

            return json.dumps(out2)

        except Exception:
            raise
            print("Erreur 3615")


class MyPlotOverlay(object):
    """Plot an overlay of 2 parameters across 1 or several file(s)"""

    def POST(self):
        try:
            import matplotlib.pyplot as plt
            print('Start MyPlot overlay')
            local_uuid = str(uuid.uuid4())

            files = getAjaxArgSmart('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            # string vs unicode problems...
            for cle, val in champ.items():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            for cle, val in champ_v.items():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            figax = (None, None)
            basemap_pickle_name = getAjaxArgSmart('basemap_pickle_name')

            monzoom = getAjaxArgSmart('monzoom')

            # 2adapt
            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            vectors_subsampling = getAjaxArgSmart('vectors_subsampling')

            dpi = getAjaxArgSmart('dpi')

            operation = getAjaxArgSmart('operation')

            #current_basemap = None
            #if os.path.exists(os.path.join(epyweb_workdir, basemap_pickle_name)):
            #    _current_basemap = pickle.load(open(os.path.join(epyweb_workdir,
            #                                                     basemap_pickle_name), 'r'))
            #    if isinstance(_current_basemap, Basemap):
            #        current_basemap = _current_basemap

            # loop sur files["A"] puis concordance avec file["B"]
            liste_tmp = []
            for indice, fichier in enumerate(files["A"]):

                # 1st layer
                myplot_args = get_common_args("A")
                resource = epygram.formats.resource(fichier, 'r')
                field = resource.readfield(champ["A"])

                if field.spectral:
                        field.sp2gp()
                if "A" in operation:
                    field = check_for_operation(operation["A"], field)

                '''
                if current_basemap is None:
                    print('Actually niou pickle !')
                    current_basemap = field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                                                                  subzone=myplot_args["subzone"],
                                                                  specificproj=myplot_args["specificproj"],
                                                                  zoom=monzoom)
                    pickle.dump(current_basemap, open(basemap_pickle_path, 'w'))
                    # On ne le calcule que pour la 1ere itération de la boucle

                '''
                myplot_1 = make_my_plot2(resource, field, "A", champ, champ_v, FF, vecteurs,
                                        figax, vectors_subsampling, myplot_args, monzoom=monzoom)

                resource.close()

                # 2nd layer
                myplot_args = get_common_args("B")
                # on met la légende à gauche pour la champ B pour ne pas enpiéter sur la légende de A
                myplot_args["colorbar"] = "left"
                resource = epygram.formats.resource(files["B"][indice], 'r')
                field = resource.readfield(champ["B"])

                if field.spectral:
                        field.sp2gp()
                if "B" in operation:
                    field = check_for_operation(operation["B"], field)
    
                myplot = make_my_plot2(resource, field, "B", champ, champ_v, FF, vecteurs,
                                      None, myplot_1, vectors_subsampling, myplot_args, monzoom=monzoom)

                try:
                    myunikname = os.path.basename(fichier) + "." + ".".join("=".join((str(k), str(v))) for k, v in champ[cle].items())
                except Exception:
                    myunikname = str(uuid.uuid4())

                # On rajoute un petit uuid en cas de rafraichissement d'image
                myunikfile = os.path.join(epyweb_workdir,
                                          myunikname + "." + local_uuid + '.png')
                myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                # SLIDE IMAGE STYLE
                # liste_tmp.append(myunikfile)
                liste_tmp.append('/getPNG/' + os.path.basename(myunikfile))
                resource.close()
                plt.close(myplot[0])
                del myplot
                del field

            out2 = liste_tmp
            current_basemap = None
            return json.dumps(out2)

        except Exception:
            raise
            print("Erreur 3614")


class MyPlotDiff(object):
    """Plot of a difference between files """

    def POST(self):
        try:
            print('Start MyPlot')
            import matplotlib.pyplot as plt
            local_uuid = str(uuid.uuid4())

            # Les listes de fichiers A et B doivent avoir la même dimension...
            filesA = getAjaxArgSmart('fileA')
            filesB = getAjaxArgSmart('fileB')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            # string vs unicode problems...
            for cle, val in champ.items():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            for cle, val in champ_v.items():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except Exception:
                    print("Warning unicode")

            # For common arguments
            myplot_args = get_common_args("A")

            figax = (None, None)
            basemap_pickle_name = getAjaxArgSmart('basemap_pickle_name')

            monzoom = getAjaxArgSmart('monzoom')
            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            getcode = getAjaxArgSmart('getcode')
            dpi = getAjaxArgSmart('dpi')

            operation = getAjaxArgSmart('operation')

            current_basemap = None
            if os.path.exists(os.path.join(epyweb_workdir, basemap_pickle_name)):
                _current_basemap = pickle.load(open(os.path.join(epyweb_workdir,
                                                                 basemap_pickle_name), 'r'))
                if isinstance(_current_basemap, Basemap):
                    current_basemap = _current_basemap

            out2 = []

            for indice, fichier in enumerate(filesA):
                resourceA = epygram.formats.resource(fichier, 'r')
                fieldA = resourceA.readfield(champ["A"])
                if (filesB[indice] != fichier):
                    resourceB = epygram.formats.resource(filesB[indice], 'r')
                    fieldB = resourceB.readfield(champ["B"])
                else:
                    fieldB = resourceA.readfield(champ["B"])
                if fieldA.spectral:
                    fieldA.sp2gp()
                if "A" in operation:
                    fieldA = check_for_operation(operation["A"], fieldA)
                if fieldB.spectral:
                    fieldB.sp2gp()
                if "B" in operation:
                    fieldB = check_for_operation(operation["B"], fieldB)
                field = fieldB - fieldA

                '''
                if current_basemap is None:
                    print('Actually niou pickle !')
                    current_basemap = field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                                                                  subzone=myplot_args["subzone"],
                                                                  specificproj=myplot_args["specificproj"],
                                                                  zoom=monzoom)
                    pickle.dump(current_basemap, open(basemap_pickle_path, 'w'))
                    # On ne le calcule que pour la 1ere itération de la boucle
                '''
                
                if (FF["A"] and FF["B"]) or (vecteurs["A"] and vecteurs["B"]):
                    fieldA_v = resourceA.readfield(champ_v["A"])
                    fieldB_v = resourceB.readfield(champ_v["B"])

                    if fieldA_v.spectral:
                        fieldA_v.sp2gp()
                    if fieldB_v.spectral:
                        fieldB_v.sp2gp()

                    vectwindA = epygram.fields.make_vector_field(fieldA, fieldA_v)
                    FF_fieldA = vectwindA.to_module()

                    vectwindB = epygram.fields.make_vector_field(fieldB, fieldB_v)
                    FF_fieldB = vectwindB.to_module()

                    if FF["A"] and FF["B"]:
                        FF_field = FF_fieldB - FF_fieldA
                        myplot = FF_field.plotfield(**myplot_args)
                        #title=str(champ["B"]) + ' - \n' + str(champ["A"]) + '\n' + str(fieldB.validity.get()) + "-" + str(fieldA.validity.get()),                                                 **myplot_args)
                        del FF_field
                        del FF_fieldA
                        del FF_fieldB

                        if vecteurs["A"] and vecteurs["B"]:
                            print("no vector difference !")
                    elif vecteurs["A"] and vecteurs["B"]:
                        del FF_fieldA
                        del FF_fieldB
                        print("no vector difference !")
                else:
                    # La méthode générique MakMyPlot n'est pas appelable ici -> duplication légère
                    #myplot_args.pop("vectorcolor", None)

                    #myplot = field.plotfield(title=str(champ["B"]) + ' - \n' + str(champ["A"]) + '\n' + str(fieldB.validity.get()) + "-" + str(fieldA.validity.get()),
                    #                         use_basemap=current_basemap,
                    #                         **myplot_args)
                    #TMP GF DEBUG
                    myplot_args.pop("pointsize", None)
                    myplot_args.pop("gisquality", None)
                    myplot_args.pop("specificproj", None)
                    myplot_args.pop("drawrivers", None)
                    myplot_args.pop("bluemarble", None)
                    #print(myplot_args)
                    myplot = field.cartoplot(title=str(champ[cle]) + '\n' + "Difference manuelle",**myplot_args)


                # Utilisation d'un nom unique par image, dans un répertoire fixe
                try:
                    myunikname = os.path.basename(fichier) + "." + ".".join("=".join((str(k), str(v))) for k, v in champ[cle].items())
                except Exception:
                    myunikname = str(uuid.uuid4())
                # On rajoute un petit uuid en cas de rafraichissement d'image
                myunikfile = os.path.join(epyweb_workdir,
                                          myunikname + "." + local_uuid + '.png')
                myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                # SLIDE IMAGE STYLE
                out2.append('/getPNG/' + os.path.basename(myunikfile))

                print("Closing figure...")
                plt.close(myplot[0])
                del myplot
                del field
                del fieldA
                del fieldB

                resourceA.close()
                try:
                    resourceB.close()
                except Exception:
                    pass

            del current_basemap

            return json.dumps(out2)
        except Exception:
            raise
            print("Erreur 3615 diff")


class GetPNG(object):
    """Retrieve figure file."""
    def GET(self, png):
        web.header('Content-type', 'image/png')
        with open(os.path.join(epyweb_workdir, png), 'rb') as f:
            return f.read()
