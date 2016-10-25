#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime, timedelta
import cPickle
import uuid
import glob
import copy
import web
import shutil

import matplotlib
matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap

from footprints.util import rangex
from opinel import interrupt

import epygram
from epygram import epylog
import usevortex
print 'Epygram version:', epygram.__version__





# Vortex cache
location_of_vortex_cache = 'MTOOLDIR'
vortex_cache_dir = os.getenv(location_of_vortex_cache)
if not vortex_cache_dir:
    location_of_vortex_cache = 'FTDIR'
    vortex_cache_dir = os.getenv(location_of_vortex_cache)
    if not vortex_cache_dir:
        location_of_vortex_cache = 'WORKDIR'
        vortex_cache_dir = os.getenv(location_of_vortex_cache)
        if not vortex_cache_dir:
            location_of_vortex_cache = 'TMPDIR'
            vortex_cache_dir = os.getenv(location_of_vortex_cache)
if location_of_vortex_cache == 'MTOOLDIR':
    vortex_cache = os.path.join(vortex_cache_dir, 'cache')
elif location_of_vortex_cache in ('FTDIR', 'WORKDIR', 'TMPDIR'):
    vortex_cache = os.path.join(vortex_cache_dir, 'mtool', 'cache')
if not vortex_cache_dir:
    raise ValueError('no rootdir has been defined for the Vortex cache: please export $MTOOLDIR.')
if not os.path.exists(vortex_cache):
    os.makedirs(vortex_cache)

# Epyweb workdir for tmp files (basemap pickle, resources hardlinks and figures)
epyweb_workdir = os.path.join(vortex_cache_dir, 'epyweb')
basemap_pickle_path = os.path.join(epyweb_workdir, 'basemap.cPickle')

#############
# Web stuff #
#############
urls = (
    '/',
    'index',
    '/faplot',
    'faplot',
    '/epyweb',
    'epyweb',
    '/myplot',
    'myplot',
    '/myplot_overlay',
    'myplot_overlay',
    '/myplot_diff',
    'myplot_diff',
    '/getfieldsasjson',
    'getfieldsasjson',
    '/getminmax',
    'getminmax',
    '/getdomain',
    'getdomain',
    '/getFile',
    'getFile',
    '/getminmaxasjson',
    'getminmaxasjson',
    '/getCacheSize',
    'getCacheSize',
    '/GetPNG/(.+)',
    'GetPNG'
    )

render = web.template.render('templates', base='base')


class getCacheSize:
    """Compute (linux only!) and return size of vortex cache"""
    def POST(self):
        try:
            cacheSize = os.popen("du -kshx " + vortex_cache).read()
            return json.dumps(cacheSize)
        except:
            return "Error in cache size retrieval"


class index:
    def GET(self):
        return render.index()


class faplot:
    def GET(self):
        return render.faplot()


class epyweb:
    def GET(self):
        return render.epyweb()


class getfieldsasjson:
    """List and returns fields of selected file"""
    def POST(self):
        try:
            fichier = getAjaxArg('file')
            print("fichier = ", fichier)
            if (os.path.isfile(fichier)):
                return json.dumps(whichFields(fichier))
            else:
                print("File does not exist => exit")
        except Exception:
            print "Erreur getfieldsasjson"
            return "Erreur getfieldsasjson"


def whichFields(fichier):
    """List and returns fields of selected file"""

    try:
        resource = epygram.formats.resource(fichier, 'r')
        listoffields = resource.listfields()
        resource.close()
        return listoffields
    except ValueError:
        raise Exception('whichField error')


class getFile:
    """Retrieved selected file(s) with usevortex"""

    def POST(self):
        try:
            reponse = {}
            data = web.data()
            vortexArgs = json.loads(data)

            # Patch pour pbs unicode + utilisation de rangex
            if 'date' in vortexArgs:
                vortexArgs['date'] = datex(vortexArgs['date'].encode())
            if 'term' in vortexArgs:
                vortexArgs['term'] = rangex(vortexArgs['term'].encode())
            if 'month' in vortexArgs:
                vortexArgs['month'] = rangex(vortexArgs['month'].encode())

            # On ajoute POUR L'INSTANT et par défaut origin=hst (pour les gridpoints)
            vortexArgs['origin'] = 'hst'
            # On enlève le mode demandé du dictionnaire d'arguments : description, existence, get
            mode = vortexArgs['request_mode']
            del vortexArgs['request_mode']

            #Utile pour garder trace du fichier, A ou B, d'origine
            try:
                fromid = vortexArgs['fromid']
                del vortexArgs['fromid']
            except Exception:
                fromid = "A"

            # Test de complétude de la description
            ressources = usevortex.get_resources(getmode='check',
                                                 **vortexArgs)
            reponse['description'] = [str(ressources)]

            # Si oui allons plus loin
            if (ressources and mode != 'description'):
                # Chemin
                ressources = usevortex.get_resources(getmode='locate',
                                                     **vortexArgs)
                reponse['remotepath'] = ressources  #[m for m in ressources]
                if mode == 'existence':
                    # Existence physique : tableau de True False
                    ressources = usevortex.get_resources(getmode='exist',
                                                         **vortexArgs)
                    reponse['existence'] = [False not in m for m in ressources]
                if mode == 'get':
                    # Rapatriement + nom local
                    ressources = usevortex.get_resources(getmode='fetch',
                                                         local=os.path.join(epyweb_workdir,
                                                                            '_'.join(["[date::ymdh]",
                                                                                      "[term]",
                                                                                      fromid,
                                                                                      str(uuid.uuid4())
                                                                                      ])),
                                                         **vortexArgs)
                    reponse['localpath'] = [str(m) for m in ressources]  #str(m[0])
            return json.dumps(reponse)
        except ValueError:
            raise Exception('getFile error')


class getminmax:
    """Compute and return min and max of FIRST field"""

    def POST(self):
        try:
            fichier = getAjaxArg('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArg('field_v')
            #subzone = getAjaxArg('subzone')  #TODO: ajouter subzone dans l'appel à stats() ?
            FF = getAjaxArg('FF')
            ope = getAjaxArg('operation')
            #string vs unicode problems...
            try:
                champ['typeOfLevel'] = champ['typeOfLevel'].encode()
                champ = {str(k):champ[k] for k in champ.keys()}
                champ_v['typeOfLevel'] = champ['typeOfLevel'].encode()
                champ_v = {str(k):champ[k] for k in champ_v.keys()}
            except:
                print("Erreur unicode")

            resource = epygram.formats.resource(fichier, 'r')
            stats = {}
            field = resource.readfield(champ)

            if field.spectral:
                field.sp2gp()
            field = CheckForOperation(ope, field)

            if FF:
                field_v = resource.readfield(champ_v)
                if field_v.spectral:
                    field_v.sp2gp()
                vectwind = epygram.fields.make_vector_field(field, field_v)
                FF_field = vectwind.to_module()
                stats['min'] = FF_field.stats()['min']
                stats['max'] = FF_field.stats()['max']
                del field_v
                del FF_field
            else:
                stats['min'] = field.stats()['min']
                stats['max'] = field.stats()['max']
            resource.close()

            return json.dumps(stats)

        except Exception, ex:
            print ex.__str__()


class getdomain:
    """Compute and return domain caracteristics"""

    def POST(self):
        try:
            fichier = getAjaxArg('file')
            resource = epygram.formats.resource(fichier, 'r')

            #On prend la géométrie du 1er champ => compatibilité FA / GRIB
            firstfield = resource.readfield(resource.listfields()[0])
            if firstfield.geometry.rectangular_grid:
                (llcrnrlon, llcrnrlat) = firstfield.geometry.gimme_corners_ll()['ll']
                (urcrnrlon, urcrnrlat) = firstfield.geometry.gimme_corners_ll()['ur']
                (ulcrnrlon, ulcrnrlat) = firstfield.geometry.gimme_corners_ll()['ul']
                (lrcrnrlon, lrcrnrlat) = firstfield.geometry.gimme_corners_ll()['lr']
            else:
                (llcrnrlon, llcrnrlat) = (-180, -90)
                (urcrnrlon, urcrnrlat) = (180, 90)

            monzoom = {
                       'lonmin': min(llcrnrlon, ulcrnrlon),
                       'lonmax': max(lrcrnrlon, urcrnrlon),
                       'latmin': min(llcrnrlat, lrcrnrlat),
                       'latmax': max(urcrnrlat, ulcrnrlat),
                       }

            resource.close()
            del resource
            del firstfield

            return json.dumps(monzoom)

        except Exception, ex:
            print ex.__str__()


class myplot:
    """Plot of a parameter for 1 or several file(s)"""

    def POST(self):
        try:
            import matplotlib.pyplot as plt
            print 'Start MyPlot'

            #For unik name
            local_uuid = str(uuid.uuid4())

            files = getAjaxArgSmart('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            #string vs unicode problems...
            for cle, val in champ.iteritems():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except Exception:
                    print ("Warning unicode")

            for cle, val in champ_v.iteritems():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except Exception:
                    print ("Warning unicode")

            existingfigure = (None, None)  #New figure (= no overlay)
            new_pickle = getAjaxArgSmart('new_pickle')

            monzoom = getAjaxArgSmart('monzoom')

            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            vectors_subsampling = getAjaxArgSmart('vectors_subsampling')

            getcode = getAjaxArgSmart('getcode')
            dpi = getAjaxArgSmart('dpi')

            try:
                existingbasemap = cPickle.load(open(basemap_pickle_path, 'r'))
                if not isinstance(existingbasemap, Basemap):
                    raise Exception('no valid Basemap in pickle.')
            except Exception:
                existingbasemap = None

            #Attention, decuml marche bien entre échéances mais actif aussi entre dates !!
            #Pas (encore ?) pour OVERLAY et DIFF
            decumul = getAjaxArgSmart('decumul')
            operation = getAjaxArgSmart('operation')

            # Algo : au 1er passage on ne trace rien, juste mis en mémoire du champ
            # ensuite on trace champLu - champAvant, si champAvant n'existe pas (cas RR @0h) -> champLu only
            out = {}

            #for fichier in files:
            for cle, val in files.iteritems():
                indiceDecumul = 0
                liste_tmp = []
                myplot_args = get_common_args(cle)
                #Garde fou pour decumul
                if FF[cle] or vecteurs[cle]:
                    decumul[cle] = False
                    print("decumul forced to ", decumul)
                try:
                    myplot_args["meridians"] = myplot_args["meridians"].encode()
                    myplot_args["parallels"] = myplot_args["parallels"].encode()
                except Exception:
                    print("Warning unicode")

                for fichier in val:
                    resource = epygram.formats.resource(fichier, 'r')
                    print("FICHIER : " + fichier)
                    try:
                        field = resource.readfield(champ[cle])
                        if field.spectral:
                            field.sp2gp()
                        if cle in operation:
                            field = CheckForOperation(operation[cle], field)
                        if decumul[cle] == True:
                            if indiceDecumul == 0:
                                fieldDecumul = field
                                indiceDecumul = +1
                                continue
                    except Exception:  #Cas des RR @0h : param n'existe pas
                        indiceDecumul = +1
                        continue

                    if decumul[cle] == True:
                            waitforme = field
                            try:  #Cas des RR @0h : param n'erxiste pas
                                field = field - fieldDecumul
                            except Exception:
                                pass
                            fieldDecumul = waitforme

                    if existingbasemap == None or new_pickle == True:
                        print 'Actually niou pickle !'
                        existingbasemap = \
                            field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                                subzone=myplot_args["subzone"], specificproj=myplot_args["specificproj"],
                                zoom=monzoom)
                        cPickle.dump(existingbasemap, open(basemap_pickle_path, 'w'))
                        #On réutilise le cache
                        new_pickle = False

                    myplot = MakeMyPlot(resource, field, cle, champ, champ_v, FF, vecteurs, existingbasemap, existingfigure, vectors_subsampling, myplot_args)

                    if getcode:
                        zecode = print_code(myplot_args, existingbasemap, existingfigure)
                        print("*** Arguments du plot ***\n" + zecode)

                    # Utilisation d'un nom unique par image, dans un répertoire fixe
                    try:
                        myunikname = os.path.basename(fichier) + "_" + "_".join("=".join((str(k), str(v))) for k, v in champ[cle].iteritems())
                    except:
                        myunikname = os.path.basename(fichier) + "_" + str(champ[cle])

                    #On rajoute un petit uuid en cas de rafraichissement d'image
                    myunikfile = os.path.join(epyweb_workdir,
                                              myunikname + "_" + local_uuid + '.png')
                    print("Saving figure ", myunikfile)
                    myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                    #memory management
                    resource.close()
                    del resource
                    del field
                    plt.close(myplot[0])
                    del myplot

                    #SLIDE IMAGE STYLE
                    liste_tmp.append('/GetPNG/' + os.path.basename(myunikfile))

                out[cle] = liste_tmp

            out2 = []

            del existingbasemap

            #On alterne A et B si besoin (cas de plot_both)
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


class myplot_overlay:
    """Plot an overlay of 2 parameters across 1 or several file(s)"""

    def POST(self):
        try:
            import matplotlib.pyplot as plt
            print 'Start MyPlot overlay'
            local_uuid = str(uuid.uuid4())

            files = getAjaxArgSmart('file')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            #string vs unicode problems...
            for cle, val in champ.iteritems():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except:
                    print ("Warning unicode")

            for cle, val in champ_v.iteritems():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except:
                    print ("Warning unicode")

            existingfigure = (None, None)
            new_pickle = getAjaxArgSmart('new_pickle')

            monzoom = getAjaxArgSmart('monzoom')

            #2adapt
            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            vectors_subsampling = getAjaxArgSmart('vectors_subsampling')

            #getcode = getAjaxArgSmart('getcode')  #TODO: useless ?
            dpi = getAjaxArgSmart('dpi')

            operation = getAjaxArgSmart('operation')

            try:
                existingbasemap = cPickle.load(open(basemap_pickle_path, 'r'))
                if not isinstance(existingbasemap, Basemap):
                    raise Exception('no valid Basemap in pickle.')
            except Exception:
                existingbasemap = None

            #out = {}  #TODO: cleanme ?

            #loop sur files["A"] puis concordance avec file["B"]
            liste_tmp = []
            for indice, fichier in enumerate(files["A"]):

                #1st layer
                myplot_args = get_common_args("A")
                resource = epygram.formats.resource(fichier, 'r')
                field = resource.readfield(champ["A"])

                if field.spectral:
                        field.sp2gp()
                if "A" in operation:
                    field = CheckForOperation(operation["A"], field)

                if existingbasemap == None or new_pickle == True:
                        print 'Actually niou pickle !'
                        existingbasemap = \
                            field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                                subzone=myplot_args["subzone"], specificproj=myplot_args["specificproj"],
                                zoom=monzoom)
                        cPickle.dump(existingbasemap, open(basemap_pickle_path, 'w'))
                        #On ne le calcule que pour la 1ere itération de la boucle
                        new_pickle = False

                myplot_1 = MakeMyPlot(resource, field, "A", champ, champ_v, FF, vecteurs, existingbasemap, existingfigure, vectors_subsampling, myplot_args)

                resource.close()

                #2nd layer
                myplot_args = get_common_args("B")
                #on met la légende à gauche pour la champ B pour ne pas enpiéter sur la légende de A
                myplot_args["colorbar"] = "left"
                resource = epygram.formats.resource(files["B"][indice], 'r')
                field = resource.readfield(champ["B"])

                if field.spectral:
                        field.sp2gp()
                if "B" in operation:
                    field = CheckForOperation(operation["B"], field)

                myplot = MakeMyPlot(resource, field, "B", champ, champ_v, FF, vecteurs, None, myplot_1, vectors_subsampling, myplot_args)

                try:
                    myunikname = os.path.basename(fichier) + "_" + "_".join("=".join((str(k), str(v))) for k, v in champ[cle].iteritems())
                except Exception:
                    myunikname = str(uuid.uuid4())


                #On rajoute un petit uuid en cas de rafraichissement d'image
                myunikfile = os.path.join(epyweb_workdir,
                                          myunikname + "_" + local_uuid + '.png')
                myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                #SLIDE IMAGE STYLE
                #liste_tmp.append(myunikfile)
                liste_tmp.append('/GetPNG/' + os.path.basename(myunikfile))
                resource.close()
                plt.close(myplot[0])
                del myplot
                del field

            out2 = liste_tmp
            existingbasemap = None
            return json.dumps(out2)

        except Exception:
            raise
            print("Erreur 3614")


class myplot_diff:
    """Plot of a difference between files """

    def POST(self):
        try:
            print 'Start MyPlot'
            import matplotlib.pyplot as plt
            local_uuid = str(uuid.uuid4())

            #Les listes de fichiers A et B doivent avoir la même dimension...
            filesA = getAjaxArgSmart('fileA')
            filesB = getAjaxArgSmart('fileB')
            champ = getAjaxArgSmart('field')
            champ_v = getAjaxArgSmart('field_v')

            #string vs unicode problems...
            for cle, val in champ.iteritems():
                try:
                    champ[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ[cle] = {str(k):champ[cle][k] for k in val.keys()}
                except Exception:
                    print ("Warning unicode")

            for cle, val in champ_v.iteritems():
                try:
                    champ_v[cle]['typeOfLevel'] = val['typeOfLevel'].encode()
                    champ_v[cle] = {str(k):champ_v[cle][k] for k in val.keys()}
                except Exception:
                    print ("Warning unicode")

            #For common arguments
            myplot_args = get_common_args("A")

            existingfigure = (None, None)
            new_pickle = getAjaxArgSmart('new_pickle')

            monzoom = getAjaxArgSmart('monzoom')
            FF = getAjaxArgSmart('FF')
            vecteurs = getAjaxArgSmart('vecteurs')
            #vectors_subsampling = getAjaxArgSmart('vectors_subsampling') #TODO: useless ou à rajouter ?
            getcode = getAjaxArgSmart('getcode')
            dpi = getAjaxArgSmart('dpi')

            operation = getAjaxArgSmart('operation')

            try:
                existingbasemap = cPickle.load(open(basemap_pickle_path, 'r'))
                if not isinstance(existingbasemap, Basemap):
                    raise Exception('no valid Basemap in pickle.')
            except Exception:
                existingbasemap = None

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
                    fieldA = CheckForOperation(operation["A"], fieldA)
                if fieldB.spectral:
                    fieldB.sp2gp()
                if "B" in operation:
                    fieldB = CheckForOperation(operation["B"], fieldB)
                field = fieldB - fieldA

                if existingbasemap == None or new_pickle == True:
                    print 'Actually niou pickle !'
                    existingbasemap = field.geometry.make_basemap(gisquality=myplot_args["gisquality"],
                        subzone=myplot_args["subzone"], specificproj=myplot_args["specificproj"], zoom=monzoom)
                    cPickle.dump(existingbasemap, open(basemap_pickle_path, 'w'))
                    #On ne le calcule que pour la 1ere itération de la boucle
                    new_pickle = False

                if (FF["A"] and FF["B"]) or (vecteurs["A"] and vecteurs["B"]):
                    fieldA_v = resourceA.readfield(champ_v["A"])
                    fieldB_v = resourceB.readfield(champ_v["B"])

                    if fieldA_v.spectral:
                        fieldA_v.sp2gp()
                    if fieldB_v.spectral:
                        fieldB_v.sp2gp()
                    #field_v = fieldB_v - fieldA_v  #TODO: cleanme ?

                    vectwindA = epygram.fields.make_vector_field(fieldA, fieldA_v)
                    FF_fieldA = vectwindA.to_module()

                    vectwindB = epygram.fields.make_vector_field(fieldB, fieldB_v)
                    FF_fieldB = vectwindB.to_module()

                    if FF["A"] and FF["B"]:
                        FF_field = FF_fieldB - FF_fieldA
                        myplot = FF_field.plotfield(title=str(champ["B"]) + ' - \n' + str(champ["A"]) + '\n' + str(fieldB.validity.get()) + "-" + str(fieldA.validity.get()),
                                                    use_basemap=existingbasemap,
                                                    **myplot_args)
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
                    #La méthode générique MakMyPlot n'est pas appelable ici -> duplication légère
                    myplot_args.pop("vectorcolor", None)

                    myplot = field.plotfield(title=str(champ["B"]) + ' - \n' + str(champ["A"]) + '\n' + str(fieldB.validity.get()) + "-" + str(fieldA.validity.get()),
                                             use_basemap=existingbasemap,
                                             **myplot_args)

                if getcode:
                    zecode = print_code(myplot_args, existingbasemap, existingfigure)
                    print("*** Arguments du plot ***\n" + zecode)

                # Utilisation d'un nom unique par image, dans un répertoire fixe
                try:
                    myunikname = os.path.basename(fichier) + "_" + "_".join("=".join((str(k), str(v))) for k, v in champ[cle].iteritems())
                except Exception:
                    myunikname = str(uuid.uuid4())
                #On rajoute un petit uuid en cas de rafraichissement d'image
                myunikfile = os.path.join(epyweb_workdir,
                                          myunikname + "_" + local_uuid + '.png')
                myplot[0].savefig(myunikfile, dpi=dpi, bbox_inches='tight')

                #SLIDE IMAGE STYLE
                out2.append(myunikfile)
                out2.append('/GetPNG/' + os.path.basename(myunikfile))

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

            del existingbasemap

            return json.dumps(out2)
        except:
            raise
            print("Erreur 3615 diff")


def whichMinMax(fichier, champ):
    try:
        #subzone = None  #TODO: init + appel à stats() ?
        resource = epygram.formats.resource(fichier, 'r')
        field = resource.readfield(champ)
        resource.close()
        if field.spectral:
            field.sp2gp()

        return field.stats()

    except Exception, ex:
        print ex.__str__()


def getAjaxArg(sArg, sDefault=''):
    """Picks out and returns a single value, regardless of GET or POST."""

    try:
        data = web.data()
        #print ('DATA :', data)
        dic = None
        if data:
            dic = json.loads(data)
        else:
            # maybe it was a GET?  check web.input()
            dic = dict(web.input())
        if dic:
            # print("DIC : ",dic)
            # print("DIC[1] : ",dic[sArg])
            if dic.has_key(sArg):
                if dic[sArg]:
                    return dic[sArg]
                else:
                    return sDefault
            else:
                return sDefault
        else:
            return sDefault
    except ValueError:
        raise Exception('getAjaxArg - no JSON arguments to decode. This method required a POST with JSON arguments.'
                        )


def getAjaxArgSmart(sArg, sDefault=''):
    """Picks out and returns a single value, regardless of GET or POST.
    Convert "" to None and true / false strings to True False, and string with starting { to dict 
    USELESS when correct use of Javascript is done !!!!!!!!!!!!!!!!!!!!!!"""

    try:
        data = web.data()
        dic = None
        if data:
            dic = json.loads(data)
        else:
            # Trick for grib
            # dic = json.loads(data, separators=(',',':'))
            # maybe it was a GET?  check web.input()
            dic = dict(web.input())
        if dic:
            if dic.has_key(sArg):
                if dic[sArg] == '':
                    return None
                elif dic[sArg] == 'false':
                    # print(sArg,": None")
                    return False
                elif dic[sArg] == 'true':
                    # print(sArg,": False")
                    return True
                elif type(dic[sArg]) == type(str()) and dic[sArg][0] \
                    == '{':
                    # print(sArg,": True")
                    #print 'on entre'
                    return json.loads(dic[sArg])
                    #print type(dic[sArg])
                else:
                    return dic[sArg]
            else:
                return sDefault
        else:
            return sDefault
    except ValueError:
        raise Exception('getAjaxArg - no JSON arguments to decode. This method required a POST with JSON arguments.'
                        )


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

    except Exception, ex:
        print "Error print_code !"
        print ex.__str__()


def get_common_args(fileid):

            myplot_args = {}
            mini = getAjaxArgSmart('min')
            maxi = getAjaxArgSmart('max')

            myplot_args["minmax"] = (mini[fileid], maxi[fileid])
            myplot_args["levelsnumber"] = getAjaxArgSmart('levelsnumber')[fileid]  # nombre/slidebar ?
            myplot_args["colormap"] = getAjaxArg('colormap')[fileid]  # un menu déroulant avec des mini-images de chaque colormap ?
            myplot_args["graphicmode"] = getAjaxArgSmart('graphicmode')[fileid]  # cases radiobutton [colorshades,contourlines,points]
            mypointsize = getAjaxArgSmart('pointsize')[fileid]
            if (mypointsize != ""):
                myplot_args["pointsize"] = mypointsize  # nombre/slidebar ?
            myplot_args["contourcolor"] = getAjaxArgSmart('contourcolor')[fileid]
            myplot_args["vectorcolor"] = getAjaxArgSmart('vectorcolor')[fileid]
            myplot_args["subzone"] = getAjaxArgSmart('subzone')[fileid]  # cases radiobutton [C,CI,CIE]
            myplot_args["gisquality"] = getAjaxArgSmart('gisquality')  # pour régler la finesse des traits de côte : cases radiobutton [c,l,i,h,f]
            myplot_args["specificproj"] = getAjaxArgSmart('specificproj')  # pour utiliser une projection de la carte particulière [kav7,cyl,ortho,nsperXXXX]
            myplot_args["center_cmap_on_0"] = getAjaxArgSmart('center_cmap_on_0')  # checkbox
            myplot_args["drawrivers"] = getAjaxArgSmart('drawrivers')  # checkbox
            myplot_args["meridians"] = getAjaxArgSmart('meridians')  # nombre/slidebar
            myplot_args["parallels"] = getAjaxArgSmart('parallels')
            myplot_args["bluemarble"] = getAjaxArgSmart('bluemarble')
            #myplot_args["minmax_in_title"] = True  # getAjaxArg('minmax_in_title') # à ignorer ?
            myplot_args["departments"] = getAjaxArgSmart('departments')  # checkbox

            return myplot_args


def datex(start, end=None, step=None):
    """
    Extended date expansion : YYYYMMDDHH-YYYYMMDDHH-HH
    """
    rangevalues = list()
    arguments = start.split('-')
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
        print ("Uncorrect date range")

    start_date = arg2date(start_arg)
    end_date = arg2date(end_arg)
    delta = timedelta(days=0, seconds=int(delta_arg) * 3600)

    d = start_date
    while d <= end_date:
        rangevalues.append(d.strftime("%Y%m%d%H"))
        d += delta

    return (rangevalues)


def arg2date(myarg):
    out = datetime(int(myarg[0:4]), int(myarg[4:6]), int(myarg[6:8]), int(myarg[8:10]))
    return out


def MakeMyPlot(resource, field, cle, champ, champ_v, FF, vecteurs, existingbasemap, existingfigure, vectors_subsampling, myplot_args):
    """Generic method for plotting (ok for plot, plot_both, overlay, but not for diff"""

    try:
        vectorcolor = myplot_args["vectorcolor"]
        myplot_args.pop("vectorcolor", None)
    except:
        vectorcolor = "black"

    print("MakeMyPlot starts ")

    #On enleve des options en cas de tracé de vecteur
    myplot_args_vect = copy.copy(myplot_args)
    myplot_args_vect.pop("pointsize", None)
    myplot_args_vect.pop("colormap", None)
    myplot_args_vect.pop("levelsnumber", None)
    myplot_args_vect.pop("graphicmode", None)
    myplot_args_vect.pop("minmax", None)
    myplot_args_vect.pop("center_cmap_on_0", None)
    myplot_args_vect.pop("departments", None)
    myplot_args_vect.pop("contourcolor", None)
    myplot_args_vect.pop("colorbar", None)

    if FF[cle] or vecteurs[cle]:
        field_v = resource.readfield(champ_v[cle])
        if field_v.spectral:
            field_v.sp2gp()
        vectwind = epygram.fields.make_vector_field(field,
                                                    field_v)

        if FF[cle]:
            FF_field = vectwind.to_module()
            myplot = FF_field.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                        use_basemap=existingbasemap,
                                        over=existingfigure,
                                        **myplot_args
                                        )
            if vecteurs[cle]:
                myplot = vectwind.plotfield(over=myplot,
                                            title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                            use_basemap=existingbasemap,
                                            subsampling=vectors_subsampling[cle],
                                            plot_module=False,
                                            symbol_options={'color': vectorcolor},
                                            **myplot_args_vect
                                            )
            del field_v
            del vectwind
        elif vecteurs[cle]:
                myplot = vectwind.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                            use_basemap=existingbasemap,
                                            over=existingfigure,
                                            subsampling=vectors_subsampling[cle],
                                            plot_module=False,
                                            symbol_options={'color': vectorcolor},
                                            **myplot_args_vect
                                            )
                del field_v
                del vectwind
    else:
                myplot = field.plotfield(title=str(champ[cle]) + '\n' + str(field.validity.get()),
                                         use_basemap=existingbasemap,
                                         over=existingfigure,
                                         **myplot_args
                                         )

    return myplot


def CheckForOperation(ope, field):
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


class GetPNG(object):
    def GET(self, png):
        web.header('Content-type', 'image/png')
        with open(os.path.join(epyweb_workdir, png), 'rb') as f:
            return f.read()


#


def func_open_browser(url=os.getcwd(), delay=0.):
    import webbrowser
    import time
    if delay > 0.:
        time.sleep(delay)
    webbrowser.open(url)


def clean_workdir():
    """Cleaning of old figures and hardlinks."""
    shutil.rmtree(epyweb_workdir, ignore_errors=True)


def init_workdir():
    """Set up the working directory."""

    clean_workdir()
    if not os.path.exists(epyweb_workdir):
        os.makedirs(epyweb_workdir)


def main(open_browser=False,
         verbose=True):

    init_workdir()

    epygram.init_env()
    epylog.setLevel('WARNING')
    if verbose:
        epylog.setLevel('INFO')

    # to avoid any path issues, "cd" to the web root. #FIXME: ? needed for the templates
    web_root = os.path.abspath(os.path.dirname(__file__))
    if os.getcwd() != web_root:
        os.chdir(web_root)

    epyweb_url = 'http://' + os.getenv('HOSTNAME') + ':8080/epyweb'
    print "====================="
    print '*epyweb* Interface =>', epyweb_url
    print '*epyweb* Workdir   =>', epyweb_workdir
    print '*vortex* Cache     =>', vortex_cache
    print '(based on $' + location_of_vortex_cache + '=' + vortex_cache_dir + ')'
    print '(the *vortex* cache location is accessible by priority order through:'
    print '$MTOOLDIR, $FTDIR, $WORKDIR, $TMPDIR'
    print "====================="
    if location_of_vortex_cache == 'TMPDIR':
        epylog.warning(' '.join(['the use of $TMPDIR as rootdir for the Vortex',
                                 'cache is hazardous. You should define a',
                                 'better rootdir using $MTOOLDIR.']))

    if open_browser:
        import threading
        t = threading.Thread(target=func_open_browser, kwargs={'url':epyweb_url,
                                                               'delay':1.})
        t.start()
    app = web.application(urls, globals())
    try:
        app.run()
    finally:
        clean_workdir()


if __name__ == '__main__':
    main()
