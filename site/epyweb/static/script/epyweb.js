// Copyright (c) Météo France (2016-)
// This software is governed by the CeCILL-C license under French law.
// http://www.cecill.info

//Descripteurs des file/fields A et B
//Clés : vortex_path, grib_param, format, local_path
Adesc = {}
Bdesc = {}

Adesc["vortex_path"] = "undefined";
Adesc["grib_param"] = [];
Adesc["format"] = "unknown";
Adesc["local_path"] = "";
Adesc["suffixe"] = "" ;
Adesc["graphicmode"] = "colorshades";

Bdesc["vortex_path"] = "undefined";
Bdesc["grib_param"] = [];
Bdesc["format"] = "unknown";
Bdesc["local_path"] = "";
Bdesc["suffixe"] = "_cloned" ;
Bdesc["graphicmode"] = "colorshades";

//Counter for unique slideshow id
var slidenumber = 0;

//On force le recalcul du basemap à chaque rechargement de la page
var basemap_pickle_name = generateUUID()+"_bm.cPickle";

//VORTEX SUGGESTIONS 
var vortex_autocomplete = {} ;
vortex_autocomplete["vapp"] = ["arome", "aladin", "arpege"] ;
vortex_autocomplete["vconf"] = ["3dvarfr","pefrance", "pifrance", "ifsfr",
                                "indien", "polynesie", "caledonie", "antilles","guyane",
                                "4dvarfr", "pearp"] ;
vortex_autocomplete["model"] = ["arome", "aladin","surfex", "arpege"] ;
vortex_autocomplete["kind"] = ["historic", "gridpoint", "ic"] ;
vortex_autocomplete["cutoff"] = ["assim", "production"] ;
vortex_autocomplete["nativefmt"]  = ["fa", "grib"] ;
vortex_autocomplete["namespace"]  = ["oper.multi.fr", "olive.multi.fr", "vortex.multi.fr", "vortex.cache.fr"] ;
vortex_autocomplete["suite"] = ["oper", "dble"] ;
vortex_autocomplete["experiment"] = ["oper", "dble"] ;

//Pas que Vortex...
vortex_autocomplete["operation"] = ["-,273.15", "exp","/,9.81","/,100"] ;

//Parameter presets

var param_presets ={
    hu2m: {
        field: "1 / 2 metre relative humidity / 1 / 2r / 0",
        TypeOfLevel: "103 /  None",
        Level: 2,
        FF: false,
        vecteurs: false,
    },
    ff10m: {
        field: "2 / 10 metre U wind component / 2 / 10u / 0",
        field_v: "3 / 10 metre V wind component / 2 / 10v / 0",
        TypeOfLevel: "103 /  None",
        Level: 10,
        FF: true,
        vecteurs: false,
    },
    isp: {
        field: "7 / unknown / 5 / unknown / 0",
        TypeOfLevel: "",
        Level: "",
        FF: false,
        vecteurs: false,
    },    
    none: {
        field: "",
        field_v: "",
        TypeOfLevel: "",
        Level: "",
        FF: false,
        vecteurs: false,
    },    
    rr: {
        field: "65 / Rain precipitation rate / 1 / rprate / 0",
        TypeOfLevel: "1 /  None",
        Level: "0",
        FF: false,
        vecteurs: false,
    },
    tpw850: {
        field: "3 / Pseudo-adiabatic potential temperature / 0 / papt / 0",
        TypeOfLevel: "100 /  None",
        Level: 850,
        FF: false,
        vecteurs: false,
    },
    A2B: { //Used for copy between A and B
        field: null,
        field_v: null,
        TypeOfLevel: null,
        Level: null,
        FF: null,
        vecteurs: null,
    },
    
    
    

}

//VORTEX presets
//Definitions of resources
var vortex_presets_resources = {
        ClearAll: {
            vapp: "",
            vconf: "",
            model: "",
            geometry: "",
            kind: "",
            nativefmt: "",
            date: "",
            cutoff: "",
            term: "",
            member: ""
            },
ArpegeForecastGrib: {
            vapp: "arpege",
            vconf: "4dvarfr",
            model: "arpege",
            geometry: "GLOB025",
            cutoff: "production",
            kind: "gridpoint",
            nativefmt: "grib",
            member: ""
        },
ArpegeForecastFA: {
            vapp: "arpege",
            vconf: "4dvarfr",
            model: "arpege",
            geometry: "global1798",
            cutoff: "production",
            kind: "historic",
            nativefmt: "fa",
            member: ""
        },
PearpGrib: {
            vapp: "arpege",
            vconf: "pearp",
            model: "arpege",
            geometry: "GLOB025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production",
            member: "1"
        },
AromeFranceForecastGrib: {
            vapp: "arome",
            vconf: "3dvarfr",
            model: "arome",
            geometry: "EURW1S40",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production",
            member: ""
        },
AromeFranceForecastFA: {
            vapp: "arome",
            vconf: "3dvarfr",
            model: "arome",
            geometry: "franmgsp",
            kind: "historic",
            cutoff: "production",
            nativefmt: "fa",
            member: ""
        },
AromeIfsForecastGrib: {
            vapp: "arome",
            vconf: "ifsfr",
            model: "arome",
            geometry: "EURW1S40",
            kind: "gridpoint",
            nativefmt: "grib",
            member: ""
        },        
PearoGrib: {
            vapp: "arome",
            vconf: "pefrance",
            model: "arome",
            geometry: "EURW1S40",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production",
            member: "1"
        },
AromeAntillesForecastGrib: {
            vapp: "arome",
            vconf: "antilles",
            model: "arome",
            geometry: "CARAIB0025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production"
        },

AromeCaledonieForecastGrib: {
            vapp: "arome",
            vconf: "caledonie",
            model: "arome",
            geometry: "NCALED0025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production"
        },
AromePolynesieForecastGrib: {
            vapp: "arome",
            vconf: "polynesie",
            model: "arome",
            geometry: "POLYN0025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production"
        },
AromeIndienForecastGrib: {
            vapp: "arome",
            vconf: "indien",
            model: "arome",
            geometry: "INDIEN0025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production"
        },
AromeGuyaneForecastGrib: {
            vapp: "arome",
            vconf: "guyane",
            model: "arome",
            geometry: "GUYANE0025",
            kind: "gridpoint",
            nativefmt: "grib",
            cutoff: "production"
        },

    
};

//Definitions of providers
    var vortex_presets_providers = {
        ClearAll: {
            namespace: "",
            suite: "",
            experiment: "",
            block: ""
        },
        DsiShell: {
            namespace: "oper.multi.fr",
            suite: "oper",
            experiment: "",
            block: "forecast"
        },
        DsiVortex: {
            namespace: "vortex.multi.fr",
            suite: "",
            experiment: "oper",
            block: "forecast"
        },
        OlivePerl: {
            namespace: "olive.multi.fr",
            suite: "",
            experiment: "XPID",
            block: ""
        },
        OliveVortex: {
            namespace: "vortex.multi.fr",
            suite: "",
            experiment: "XPID",
            block: "forecast"
        }
};


var zoom_presets_zones = {
        ClearAll: {
            north: 90,
            east: 180,
            west: -180,
            south: -90
        },
        Antilles: {
            north: 18,
            east: -59,
            west: -63,
            south: 13
        },
        Reunion: {
            north: -20.3,
            east: 56.5,
            west: 54.5,
            south: -21.8
        },
        Tahiti: {
            north: -16.8,
            east: -148.8,
            west: -150.3,
            south: -18.5
        },
        NlleCaledonie: {
            north: -19.8,
            east: 168.3,
            west: 163.5,
            south: -22.8
        },
        Guyane: {
            north: 6,
            east: -51.3,
            west: -54.9,
            south: 2
        },
};

$(document).ready(function() {
            /*
             *g
             */
	getGeometries("");
	var currentday = $.datepicker.formatDate('yymmdd', new Date());
    var prevDiv = $("#graphicplot0") //on garde l'historique du plot précédent pour décaler les dialogues contenant les figures

    //Variable with the zoom spinners characteristics
    var zoom_spinner = {};
    //Min and max range
    zoom_desc_init = {};
    zoom_desc_init["north"] = 90;
    zoom_desc_init["south"] = 90;
    zoom_desc_init["east"] = 180;
    zoom_desc_init["west"] = 180;

    //Make zoom inputs as spinner with min/max values....
    $.each(zoom_desc_init, function(key, value) {
        zoom_spinner[key] = $("#zoom_" + key).spinner();
        zoom_spinner[key].spinner("enable");
        zoom_spinner[key].spinner("option", "max", value);
        zoom_spinner[key].spinner("option", "min", -value);
        zoom_spinner[key].spinner("option", "numberFormat", "n");
        zoom_spinner[key].spinner("option", "step", 0.1);
    });    
    
    var waiting_text = '<img id="img-wait-spinner" src="static/loader.gif" alt="Loading"/></img>';

    ///////// V O R T E X  I N I T ///////// 
    
    $( "#filenfield_tabs" ).tabs();
    $( "#modetabs" ).tabs();
    
    //Initialisations des presets
    VortexPresetInit("");
    
    //Initialisations des autocomplete
    VortexAutoCompleteInit("");

    //VORTEX initialisations 
    $('#checkexistence').prop("disabled", true);
    $('#getfile').prop("disabled", true);
    $('#checkdescription').prop("disabled", false);
    $("#vortex_preset_resources").val('ClearAll');
    $("#vortex_preset_providers").val('ClearAll');
    


    ///////// O T H E R  I N I T S ///////// 

    $("#cloneFileField").prop("disabled", false);
    
    $(".gribonly").prop("disabled", true);
     


    //Initialisation ici

    //Accordion defauts, and fine tune according to "collapse" class
    $(function() {
        $(".accordion").accordion({
            collapsible: true,
            heightStyle: "content"
        });
        $(".accordion.collapse").accordion({
            active: false,
        });
    });
    $(".resizable").resizable({});

    //Misc initialisations
    $(function() {
        $("#tabs_plot").tabs();
        $("#tabs_plot").tabs("disable", "#vectors");
    });

    
    $( function() { 
        $( document ).tooltip(
            {
            show: {
                //effect: "blind",
                delay: 1000
            },
            position: {
                my: "right",
                at: "left bottom+55",
                of: "#fortooltip"
            }
            });
            } );

    //Colormap list
    $.widget("custom.iconselectmenu", $.ui.selectmenu, {
        _renderItem: function(ul, item) {
            var li = $("<li>", {
                text: item.label
            });
            if (item.disabled) {
                li.addClass("ui-state-disabled");
            }
            $("<span>", {
                    style: item.element.attr("data-style"),
                    "class": "ui-icon " + item.element.attr("data-class")
                })
                .appendTo(li);
            return li.appendTo(ul);
        }
    });
    
    InitColorMapList("","");

    //Initialisations of
    // #graphicmode, #operation, #gisquality, #subzone et #specificproj
    $(".selectmenu").val('');
    $(".selectmenu").selectmenu();

    //On rend draggable les fieldsets 
    //$( ".important" ).draggable();
    //$( ".option" ).draggable();

    //Clear previous inputs
    $(':input').not(':button, :submit, :reset, :hidden').val(''); //, :checkbox, :radio').val('');

    //Initialisation diverses et variées
    $(".initchecktrue").prop("checked", true);
    $(".initcheckfalse").prop("checked", false);
    $("#field2").prop("disabled", true);
    $("#file2").prop("disabled", true);
    $("#date").val(currentday + "00");

    //Gestion dynamique du champ (field) ou des champs (U et V)
    $("#field_v").prop("disabled", false);
    ManageUandV("");

    $("#FF").click(function() {
        ManageUandV("");
    });
    $("#vecteurs").click(function() {
        ManageUandV("");
        if (this.checked) {
            $("#tabs_plot").tabs("enable", "#vectors");
        } else {
            $("#tabs_plot").tabs("disable", "#vectors");
        }
    });

    var spinner_sub = $("#vectors_subsampling").spinner();
    spinner_sub.spinner("enable");
    spinner_sub.spinner("value", 10);
    spinner_sub.spinner("option", "max", 100);
    spinner_sub.spinner("option", "min", 1);


    //Initialisation du slider min/max
    changeSlider(0, 1,"");
    smart_slider("");
    //$("#range1").val($("#slider-range").slider("values", 0) + " to " + $("#slider-range").slider("values", 1));


    //BlueMarble stuff
    
    $( "#bluemarble" ).slider({
      value:0,
      min: 0,
      max: 1.05,
      step: 0.1,
      slide: function( event, ui ) {
        $( "#bluemarble_value" ).val( ui.value );
      }
    });
    $( "#bluemarble_value" ).val( $( "#bluemarble" ).slider( "value" ) );
    
    //Prise en compte des changements des min/max entrés au clavier (avec réattribution des min/max possibles du slider)
    //Et désactivation dans ce cas de l'autoupdate
    //smart_slider("suffixe")

    //Désactivation de l'autoupdate si on change les curseurs de place
    //Ainsi tous les cas sont traités sauf MàJ auto qui ne désactive pas bien sur
    $('#slider-range').on('slidestop', function() {
        desactivateAutoUpdate("");
    });

    //Initialisation (+ style) des boutons
    $(".boutonlike").buttonset();


    //On désactive le plot simple ; sera activé après un getfile
    $("#getplot").prop("disabled", true);
    //Les autres sont invisibles
    $("#getplotboth").prop("disabled", true);
    $("#getplotboth").hide()
    $("#overlay").prop("disabled", true);
    $("#overlay").hide();
    $("#difference").prop("disabled", true);
    $("#difference").hide();

    //Spinners initialisations
    var spinner_lvl1 = $("#levelsnumber1").spinner();
    spinner_lvl1.spinner("enable");
    spinner_lvl1.spinner("value", 50);

    var spinner_lvl2 = $("#levelsnumber2").spinner();
    spinner_lvl2.spinner("enable");
    spinner_lvl2.spinner("value", 50);
    spinner_lvl2.spinner("disable");

    var spinner_pt = $("#pointsize").spinner();
    spinner_pt.spinner("enable");
    spinner_pt.spinner("value", 20);
    //spinner_pt.spinner( "disable" );

    var spinner_ll = $("#llgridstep").spinner();
    spinner_ll.spinner("enable");

    var spinner_dpi = $("#dpi").spinner();
    spinner_dpi.spinner("enable");
    spinner_dpi.spinner("value", 200);

    //Garde fou pour lat lon gridstep (valeur entre 0 et 180)
    $(function() {
        $("#llgridstep").spinner({
            spin: function(event, ui) {
                if (ui.value < 0) {
                    $(this).spinner("value", 180);
                    return false;
                } else if (ui.value > 180) {
                    $(this).spinner("value", 0);
                    return false;
                }
            }
        });
    });



    LatLon_sets = ["none","1","auto","tropics"] ;
    MakeAutoComplete("#meridians",LatLon_sets)
    MakeAutoComplete("#parallels",LatLon_sets)


    //Labels pré-initialisés pour éviter un changement lors du 1er clic
      //Partie render
    setup_render_labels("");
    
      //Autres
    $("#drawdepartments_label").text("Draw departments is off");
    $("#center_cmap_on_0_label").text("Center on 0 is off");
    $("#drawrivers_label").text("Draw rivers is off");
    $("#active_zoom").prop('checked', false).button('refresh');
    $("#active_zoom_label").text("Zooming is off");


    
    $("#center_cmap_on_0").click(function() {
        $("#center_cmap_on_0_label").text(this.checked ? "Center on 0 is on " : "Center on 0 is off");
    });

    $("#drawdepartments").click(function() {
        $("#drawdepartments_label").text(this.checked ? "Draw departments is on " : "Draw departments is off");
    });

    $("#drawrivers").click(function() {
        $("#drawrivers_label").text(this.checked ? "Draw rivers is on " : "Draw rivers is off");
    });

    $("#active_zoom").click(function() {
        $("#active_zoom_label").text(this.checked ? "Zooming is on" : "Zooming is off");
    });
    
    
    


    //Gestion du pickle en fonction des choix de tracé :
    //basemap_pickle_name = generateUUID()+"_bm.cPickle"; quand il faut le retracer
    //pour le zoom, on en profite pour l'activer automatiquement
    $("#gisquality").selectmenu({
        change: function() {
            basemap_pickle_name = generateUUID()+"_bm.cPickle";
        }
    });
    $("#subzone").selectmenu({
        change: function() {
            basemap_pickle_name = generateUUID()+"_bm.cPickle";
        }
    });
    $("#specificproj").selectmenu({
        change: function() {
            basemap_pickle_name = generateUUID()+"_bm.cPickle";
        }
    });
    $("#zoom_east").on("focus", function() {
        basemap_pickle_name = generateUUID()+"_bm.cPickle";
        activateZoom();
    });
    $("#zoom_north").on("focus", function() { //old : spinner({change:function(){
        basemap_pickle_name = generateUUID()+"_bm.cPickle";
        activateZoom();
    });
    $("#zoom_west").on("focus", function() {
        basemap_pickle_name = generateUUID()+"_bm.cPickle";
        activateZoom();
    });
    $("#zoom_south").on("focus", function() {
        basemap_pickle_name = generateUUID()+"_bm.cPickle";
        activateZoom();
    });
    $("#active_zoom").change(function() {
        basemap_pickle_name = generateUUID()+"_bm.cPickle";
    });
    
    //Si on modifie la description on met en evidence l'incoherence entre la description et ce qui est en mémoire
    $(".vortex").change(function() {
    	$('#vortex_existence span').text("Unchecked");;
    	$('#vortex_get span').text("Current description not fetched");
    });
    
    
    
    
    //Zones pre-configurees pour zoom
    
    $("#zoom_preset_zones").selectmenu({
        change: function(event, data) {
            for (var key in zoom_presets_zones[data.item.value]) {
                //alert(key);
                zoom_spinner[key].spinner("value", zoom_presets_zones[data.item.value][key]);
            }
             basemap_pickle_name = generateUUID()+"_bm.cPickle";
             activateZoom();
        }
    });


    //Quand un onglet est activé : 
    //1. on garde son Id pour savoir le type de plot voulu (sauf si vectors => colorshades
    //2. on déplace les élements persistants d'un tab à l'autre (#move_id -> moving_id_content_nb)
    activate_tabs_plot(Adesc);
    
    $("#closeall").click(function() {
        $(".ui-dialog-content").dialog("close");
    })
    
    $("#resetall").click(function() {
        $("#getplot").prop("disabled", false);
        $("#getplotboth").prop("disabled", false);
        $("#overlay").prop("disabled", false);
        $("#difference").prop("disabled", false);
        $("#getfile").prop("disabled", false);
        $("#getfile_cloned").prop("disabled", false);
        })

    Param_Presets_Init("")
    

$("#cachesize").click(function() {
    updateCacheSize();
});
    

    ///////// L E T ' S  P L O T ! ///////// 

    //Création d'un graphique (initié par un clic) : calcul, puis affichage du code html de retour dans une boite de dialogue
    //2DO: on plot le graphique correspondant à l'onglet actif !
$("#getplot").click(function() {
    
	console.log("On démarre le plot");
	
	cleantTooltips();
        
	//Customisation en fonction de l'onglet actif
    var currentTab = $("#filenfield_tabs").tabs('option', 'active');
	//alert(currentTab);
    var activeDesc = Adesc ;
    var activeSuffixe = "" ;
    if (currentTab == 1) {
    	activeDesc = Bdesc ; 
    	activeSuffixe = "_cloned" ;
    	} 
        
    	
    	
        //Validations de base pour éviter des erreurs
        if ($("#field").hasClass('fieldnothere')) {  
            alert("Enter a -relevant- field...");
        } else {
            //Mise en attente pour signifier qu'un plot est en création
            //et désactivation du bouton plot pour éviter tout conflit
            $("#getplot").html(waiting_text);
            $("#getplot").prop("disabled", true);

            //Pour entête boite de dialogue
            //label = $("#field").val();
            label = $("#experiment" + activeSuffixe).val() + $("#suite" + activeSuffixe).val() ;

            //Si demandé recherche des min/max
            if ($("#autoupdate" + activeSuffixe).prop('checked') == true) {
                var args_minmax_json = JSON.stringify(giveme_args_minmax(activeDesc));
                AutoUpdateMinMax(args_minmax_json, activeSuffixe)
            } else {};
            //On met le zoom à jour si besoin (cas d'un fichier rapatrié avec zoom activé, et qu'on veut avoir les coordonnees du nouveau fichier)
            //FBI !! GetCoordinates(zoom_spinner)

            //Construction des variables à passer
            var args_plot = {};
            
            args_plot["file"] = {};
            args_plot["field"] = {};
            args_plot["field_v"] = {};
            
            //On utilise .sort() pour avoir les figures dans l'ordre : toutes les échéances d'un réseau se suivent
            args_plot["file"]["A"] = activeDesc["local_path"].sort(); 

             //Gestion de la lecture des paramètres selon le format de fichier
            //Cas d'un fichier grib1
            if (activeDesc["format"] == "grib1") {
                args_plot["field"]["A"] = giveme_grib1_handgrip(activeSuffixe) ;
                //Pour la composante V, tout idem sauf indicatorOfParameter
                args_plot["field_v"]["A"] = giveme_grib1_handgrip(activeSuffixe) ;
                try {
                    args_plot["field_v"]["A"]["indicatorOfParameter"] = parseInt(GetGribKey("#field_v" + activeSuffixe,0));                    
                }
                catch (err) {}
            }

            //Cas d'un fichier grib2
            else if (activeDesc["format"] == "grib2") {
                args_plot["field"]["A"] = giveme_grib2_handgrip(activeSuffixe) ;
                //Pour la composante V, tout idem sauf indicatorOfParameter
                args_plot["field_v"]["A"] = giveme_grib2_handgrip(activeSuffixe) ;
                try {
                    args_plot["field_v"]["A"]["parameterNumber"] = parseInt(GetGribKey("#field_v" + activeSuffixe,0));                    
                }
                catch (err) {}
            }

            else {
                args_plot["field"]["A"] = $("#field" + activeSuffixe).val();
                args_plot["field_v"]["A"] = $("#field_v" + activeSuffixe).val();
            }      
            
            //Common args between A and B
            args_plot = get_common_args(args_plot);
            //Specific args between A and B
            args_plot = get_specific_args(args_plot, activeSuffixe,"A");

            args_plot["graphicmode"] = {}
            args_plot["graphicmode"]["A"] = activeDesc["graphicmode"];
            args_plot["basemap_pickle_name"] = basemap_pickle_name;
            
            var args_plot_json = JSON.stringify(args_plot);
            //console.log(args_plot);
            
            // Plot demandé
            $.ajax({
                type: "POST",
                async: true,
                url: "/myplot",
                data: args_plot_json,
                contentType: "application/json; charset=utf-8",
                dataType: "json",
                success: function(response) {

                    var newDiv = $(document.createElement('div'));

                    plot2html(newDiv, prevDiv, response, 1)

                    $("#getplot").html('Plot');
                    $("#getplot").prop("disabled", false);
                    prevDiv = newDiv

                },
                error: function(response) {
                    alert(response.responseText);
                    $("#getplot").html('Plot');
                    $("#getplot").prop("disabled", false);
                }
            });
        }
    });


// O V E R L A Y
$("#overlay").click(function() {
    
    cleantTooltips();
    
    console.log("On démarre le plot overlay");

        //Validations de base pour éviter des erreurs
        if ($("#field").hasClass('fieldnothere')) {  
            alert("Enter a -relevant- field...");
        } else {
            //Mise en attente pour signifier qu'un plot est en création
            //et désactivation du bouton plot pour éviter tout conflit
            $("#overlay").html(waiting_text);
            $("#overlay").prop("disabled", true);

            //Pour entête boite de dialogue
            //label = "Overlay " + $("#field").val(); 
            label = $("#experiment_cloned").val()+ $("#suite_cloned").val() + " on top of " + $("#experiment").val()+ $("#suite").val();

            //Si demandé recherche des min/max
            if ($("#autoupdate").prop('checked') == true) {
                var args_minmax_json = JSON.stringify(giveme_args_minmax(Adesc));
                AutoUpdateMinMax(args_minmax_json,"")
            }
            if ($("#autoupdate_cloned").prop('checked') == true) {
                args_minmax_json = JSON.stringify(giveme_args_minmax(Bdesc));
                AutoUpdateMinMax(args_minmax_json,"_cloned")
            } else {};

            //Construction des variables à passer
            var args_plot_both = {};
            
            //Dict pour les chemins
            args_plot_both["file"] = {}
            args_plot_both["file"]["A"] = Adesc["local_path"].sort() ;
            args_plot_both["file"]["B"] = Bdesc["local_path"].sort() ;

            //$.each(Adesc["local_path"], function( index, value ) {
            //        	args_plot_both["file"].push(Adesc["local_path"][index]);
            //        	args_plot_both["file"].push(Bdesc["local_path"][index]);
            //        });

             //Gestion de la lecture des paramètres selon le format de fichier
             //On différencie les paramètres de A et de B
            //Cas d'un fichier grib
            args_plot_both["field"] = {}
            args_plot_both["field_v"] = {}
            
            //On récupère field et field_v, pour les 2 descripteurs A et B
            gime_fields_args(args_plot_both,Adesc,"A","")
            gime_fields_args(args_plot_both,Bdesc,"B","_cloned")
            
            //Common args between A and B
            args_plot_both = get_common_args(args_plot_both);

            //Specific args between A and B
            //Ordre important pour le moment, car il existe des champs None pour "B"
            //Not anymore
            args_plot_both = get_specific_args(args_plot_both,"_cloned","B");
            args_plot_both = get_specific_args(args_plot_both,"","A");
            

            args_plot_both["graphicmode"] = {}
            args_plot_both["graphicmode"]["A"] = Adesc["graphicmode"] ;
            args_plot_both["graphicmode"]["B"] = Bdesc["graphicmode"] ;
            args_plot_both["basemap_pickle_name"] = basemap_pickle_name;
            
            
            var args_plot_both_json = JSON.stringify(args_plot_both);
            
            // Plot demandé
            $.ajax({
                type: "POST",
                async: true,
                url: "/myplot_overlay",
                data: args_plot_both_json,
                contentType: "application/json; charset=utf-8",
                dataType: "json",
                    success: function(response) {
                        
                    var newDiv = $(document.createElement('div'));
                    
                    plot2html(newDiv,prevDiv,response,1)
                    
                    $("#overlay").html('Overlay A and B');
                    $("#overlay").prop("disabled", false);
                    prevDiv = newDiv
                
                },
                error: function(response) {
                    alert(response.responseText);
                    $("#overlay").html('Overlay A and B');
                    $("#overlay").prop("disabled", false);
                }
            });

            //Pas besoin d'un nouveau pickle sauf exceptions listées plus haut
            //new_pickle = "false";
        }
        
})

// D I F F E R E N C E 
$("#difference").click(function() {
    
    cleantTooltips();
    	
    	console.log("On démarre la diff");

        //Validations de base pour éviter des erreurs
    	if ($("#field").hasClass('fieldnothere')) {  
            alert("Enter a -relevant- field...");
        } else {
            //Mise en attente pour signifier qu'un plot est en création
            //et désactivation du bouton plot pour éviter tout conflit
            $("#difference").html(waiting_text);
            $("#difference").prop("disabled", true);

            //Pour entête boite de dialogue
            //label = "Difference " + $("#field").val();
            label = $("#experiment_cloned").val()+ $("#suite_cloned").val() + " - " + $("#experiment").val()+ $("#suite").val();

            //Si demandé recherche des min/max
            if ($("#autoupdate").prop('checked') == true) {
                var args_minmax_json = JSON.stringify(giveme_args_minmax(Adesc));
                AutoUpdateMinMax(args_minmax_json,"")
            } else {};

            //Construction des variables à passer
            var args_plot_diff = {};
            
            args_plot_diff["fileA"] = Adesc["local_path"].sort(); 
            args_plot_diff["fileB"] = Bdesc["local_path"].sort();

            //Gestion de la lecture des paramètres selon le format de fichier
            //On récupère field et field_v, pour les 2 descripteurs A et B
            args_plot_diff["field"] = {}
            args_plot_diff["field_v"] = {}
            gime_fields_args(args_plot_diff,Adesc,"A","")
            gime_fields_args(args_plot_diff,Bdesc,"B","_cloned")
            
            //Common args between A and B
            args_plot_diff = get_common_args(args_plot_diff);

            //Specific args between A and B

            //Ordre important pour le moment, car il existe des champs None pour "B"
            //Not anymore...
            args_plot_diff = get_specific_args(args_plot_diff,"_cloned","B");
            args_plot_diff = get_specific_args(args_plot_diff,"","A");
            

            args_plot_diff["graphicmode"] = {}
            args_plot_diff["graphicmode"]["A"] = Adesc["graphicmode"] ;
            args_plot_diff["basemap_pickle_name"] = basemap_pickle_name;
            
            
            var args_plot_diff_json = JSON.stringify(args_plot_diff);
            
            // Plot demandé
            $.ajax({
                type: "POST",
                async: true,
                url: "/myplot_diff",
                data: args_plot_diff_json,
                contentType: "application/json; charset=utf-8",
                dataType: "json",
                    success: function(response) {
                    var newDiv = $(document.createElement('div'));
                    
                    plot2html(newDiv,prevDiv,response,1)
                    
                    $("#difference").html('Difference B-A');
                    $("#difference").prop("disabled", false);
                    prevDiv = newDiv
                
                },
                error: function(response) {
                    alert(response.responseText);
                    $("#difference").html('Difference B-A');
                    $("#difference").prop("disabled", false);
                }
            });

            //Pas besoin d'un nouveau pickle sauf exceptions listées plus haut
            //new_pickle = "false";
        }
    });


// P L O T  A  A N D  B
$("#getplotboth").click(function() {
    
    cleantTooltips();
    	
    	console.log("On démarre le plot both");

        //Validations de base pour éviter des erreurs
    	if ($("#field").hasClass('fieldnothere')) {  
            alert("Enter a -relevant- field...");
        } else {
            //Mise en attente pour signifier qu'un plot est en création
            //et désactivation du bouton plot pour éviter tout conflit
            $("#getplotboth").html(waiting_text);
            $("#getplotboth").prop("disabled", true);

            //Pour entête boite de dialogue
            label = $("#experiment").val()+ $("#suite").val() + " (left) vs " + $("#experiment_cloned").val()+ $("#suite_cloned").val() + " (right)"; 

            //Si demandé recherche des min/max
            if ($("#autoupdate").prop('checked') == true) {
                var args_minmax_json = JSON.stringify(giveme_args_minmax(Adesc));
                AutoUpdateMinMax(args_minmax_json,"")
            }
            if ($("#autoupdate_cloned").prop('checked') == true) {
                args_minmax_json = JSON.stringify(giveme_args_minmax(Bdesc));
                AutoUpdateMinMax(args_minmax_json,"_cloned")
            } else {};

            //Construction des variables à passer
            var args_plot_both = {};
            
            //Dict pour les chemins
            args_plot_both["file"] = {}
            args_plot_both["file"]["A"] = Adesc["local_path"].sort() ;
            args_plot_both["file"]["B"] = Bdesc["local_path"].sort() ;

            //$.each(Adesc["local_path"], function( index, value ) {
            //        	args_plot_both["file"].push(Adesc["local_path"][index]);
            //        	args_plot_both["file"].push(Bdesc["local_path"][index]);
            //        });

             //Gestion de la lecture des paramètres selon le format de fichier
             //On différencie les paramètres de A et de B
            //Cas d'un fichier grib
            args_plot_both["field"] = {}
            args_plot_both["field_v"] = {}
            
            //On récupère field et field_v, pour les 2 descripteurs A et B
            gime_fields_args(args_plot_both,Adesc,"A","")
            gime_fields_args(args_plot_both,Bdesc,"B","_cloned")
            
            //Common args between A and B
            args_plot_both = get_common_args(args_plot_both);

            //Specific args between A and B
            //Ordre important pour le moment, car il existe des champs None pour "B"
            //Not anymore
            args_plot_both = get_specific_args(args_plot_both,"_cloned","B");
            args_plot_both = get_specific_args(args_plot_both,"","A");
            

            args_plot_both["graphicmode"] = {}
            args_plot_both["graphicmode"]["A"] = Adesc["graphicmode"] ;
            args_plot_both["graphicmode"]["B"] = Bdesc["graphicmode"] ;
            args_plot_both["basemap_pickle_name"] = basemap_pickle_name;
            
            
            var args_plot_both_json = JSON.stringify(args_plot_both);
            
            // Plot demandé
            $.ajax({
                type: "POST",
                async: true,
                url: "/myplot",
                data: args_plot_both_json,
                contentType: "application/json; charset=utf-8",
                dataType: "json",
                    success: function(response) {
                        
                    var newDiv = $(document.createElement('div'));
                    
                    plot2html(newDiv,prevDiv,response,2)
                    
                    $("#getplotboth").html('Plot A and B');
                    $("#getplotboth").prop("disabled", false);
                    prevDiv = newDiv
                
                },
                error: function(response) {
                    alert(response.responseText);
                    $("#getplotboth").html('Plot A and B');
                    $("#getplotboth").prop("disabled", false);
                }
            });

            //Pas besoin d'un nouveau pickle sauf exceptions listées plus haut
            //new_pickle = "false";
        }
    });




    ///////// L E T ' S  G R I B ! ///////// 

    //On ajuste l'autocomplétion à la saisie (utile pour fichier grib)
    //$("#range1").change(function() {

    $(".smart_grib").change(function() {
        if (Adesc["format"] == "grib1") {SmartGribSelect(Adesc);}
        if (Adesc["format"] == "grib2") {SmartGrib2Select(Adesc);}
    });

    $(".smart_grib").click(function() {
        if (Adesc["format"] == "grib1") {SmartGribSelect(Adesc);}
        if (Adesc["format"] == "grib2") {SmartGrib2Select(Adesc);}
    });


    ///////// L E T ' S  V O R T E X ! ///////// 

    $("#file_path").button({icons: {primary: 'ui-icon-info'}});
    //$( "#file_path" ).button( "refresh" );
    
    
    
    
    $("#file_path")
       //.button()
       .click(function() {
        alert(Adesc["vortex_path"]);
    });
    
    //On vérifie la complétude de la description à chaque changement de valeur
    $(".vortex").change(function() {
        CheckVortexDescription("");
    });

    //On demande explicitement la vérification de la description
    $("#checkdescription")
       //.button()
       .click(function() {
        CheckVortexDescription("");
    });

    //Cas où on demande la vérification d'existence
    $("#checkexistence")
       //.button()
       .click(function() {
           
        cleantTooltips();   
        
        var args_vortex_check_json = CreateVortexArgs("existence","");

        $('#vortex_existence span').text("Checking...");

        $.ajax({
            type: "POST",
            async: true,
            url: "/getFile",
            data: args_vortex_check_json,
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(vortexAnswer) {
                Adesc["vortex_path"] = vortexAnswer["remotepath"];
                $('#vortex_existence span').text(vortexAnswer["existence"]);
                $("#vortex_existence").effect("pulsate", 500);
            },
            error: function(response) {
                $('#vortex_existence span').text("Error");
            },
        });
    });

    //Cas où on demande le rapatriement 
$("#getfile")
       //.button()
       .click(function() {
           
        cleantTooltips();   

        $('#vortex_get span').text("Transfering data...");
        
        var args_vortex_get_json = CreateVortexArgs("get","");
       
        $.ajax({
            type: "POST",
            async: true,
            url: "/getFile",
            data: args_vortex_get_json,
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(vortexAnswer) {
                
                Adesc["local_path"] = vortexAnswer["localpath"];
        
                //UpdateFields(Adesc["local_path"][0],zoom_spinner);
                UpdateFields(Adesc,zoom_spinner); 

                $('#vortex_get span').text("Current description fetched!");
                
                $("#getplot").prop("disabled", false);
                $("#getplotboth").prop("disabled", false);
                $("#overlay").prop("disabled", false);
                $("#difference").prop("disabled", false);
                
                $("#accordion_vortex").accordion("option", "active", 1);
                
                makeTabTitle(args_vortex_get_json,"fileA")
                //var backToJsonObj = JSON.parse(args_vortex_get_json);
                //var jsonPretty = JSON.stringify(backToJsonObj, ["model","vconf","geometry","date","cutoff","term","expériment","block","suite"], 4);
                //alert(jsonPretty);
                //$('#fileA').prop('title', "Currently in memory : " + jsonPretty.replace(/\"/g,"") );
                //alert(JSON.stringify(args_vortex_get_json,null,2))


            },
            error: function(response) {
                $('#vortex_get span').text("Transfer error!");
            }
        });
        
    });
$('#filesystemOK').prop("disabled", false);
$("#filesystemOK")
       //.button()
       .click(function() {
       Adesc["local_path"] = $("#filesystem_path").val().split('\n');
       //alert(Adesc["local_path"][1]);
       UpdateFields(Adesc,zoom_spinner);         
        $("#getplot").prop("disabled", false);
        $("#getplotboth").prop("disabled", false);
        $("#overlay").prop("disabled", false);
        $("#difference").prop("disabled", false);
        
        $("#accordion_vortex").accordion("option", "active", 1);
        
        makeTabTitleFS(Adesc["local_path"],"fileA")
	});
	
	
	
	
    /////// C L O N E  E T  A L 

$("#cloneFileField").click(function () {
    
    cleantTooltips();
        
         //Bouton à usage unique !
         $("#cloneFileField").prop("disabled", true);
         $("#cloneFileField").hide();

         //Valeur de la colormap
         currentColorMap = $("#ColorMapList").val()
         
         //On arrête les select et tabs avant clonage
         $("#vortex_preset_resources" ).selectmenu( "destroy" );
         $("#vortex_preset_providers" ).selectmenu( "destroy" );
         $("#ColorMapList" ).iconselectmenu( "destroy" );
         $("#tabs_plot" ).tabs( "destroy" );
         $("#modetabs" ).tabs( "destroy" );
         
         //Idem pour les boutons
         $(".boutonlike").buttonset("destroy");
         
         //New tab
         var cloneTabVortex =  $('#filenfield_tab1').clone(false);//.off();
         cloneTabVortex.prop('id', 'filenfield_tab2' );
        
        //On renomme tous les id,car elles doivent rester uniques
        cloneTabVortex.find("*").each(function(index, element) { // And all inner elements.
            if (element.id) {
                element.id = element.id + "_cloned";
            }
            if (element.className.search("vortex") > -1) {
                element.className = "vortexcloned";
            }
            if (element.className.search("smart_grib") > -1) {
                //element.className = element.className + " smart_grib_cloned";
                element.className = "smart_grib_cloned";
            }
        });         
        //On rajoute le tab (corps + entete)
        cloneTabVortex.insertAfter("#filenfield_tab1");
        
        var fileAinMem = $('#fileA').attr("title") ;
        //alert($('#fileA').attr("title"));
        $(".roottabs").append('<li><a id="fileB" title="' + fileAinMem + '" href="#filenfield_tab2"><span class="filenfield_tab_list">___ B ___</span></a></li>');
        
        
    	 
    	//On actualise les tabs et on différencie leur classe pour un affichage plus coloré
    	$( "#filenfield_tabs" ).tabs("refresh");
    	$('#filenfield_tabs .ui-tabs-nav a[href="#filenfield_tab1"], #filenfield_tab1').addClass('original');
    	$('#filenfield_tabs .ui-tabs-nav a[href="#filenfield_tab2"], #filenfield_tab2').addClass('cloned');
    	$('#accordion_render_cloned').addClass('cloned');
    	
    	
    	
    	 
        $("#accordion_vortex_cloned").accordion({
                active: false,
                collapsible: true,
                heightStyle: "content"
            });
        $("#accordion_actions_cloned").accordion({
                collapsible: true,
                active: false,
                heightStyle: "content"
            });
        $("#accordion_render_cloned").accordion({
                collapsible: true,
                active: false,
                heightStyle: "content"
            });
            
        
        //On rend visibles les boutons
        //$("#getplotboth").prop("disabled", false);
        $("#getplotboth").show()
        //$("#overlay").prop("disabled", false);
        $("#overlay").show()
        //$("#difference").prop("disabled", false);
        $("#difference").show();
        
        VortexAutoCompleteInit("_cloned");
        
        //On réactive les select 
        VortexPresetInit("");
        VortexPresetInit("_cloned");

        //Auto complétion des geometries
        getGeometries("_cloned");
        
        
        //Gestion auto completion grib
        $(".smart_grib_cloned").change(function() {
            if (Bdesc["format"] == "grib1") {
                SmartGribSelect(Bdesc);
            }
            if (Bdesc["format"] == "grib2") {
                SmartGrib2Select(Bdesc);
            }
        });

        $(".smart_grib").click(function() {
            if (Bdesc["format"] == "grib1") {
                SmartGribSelect(Bdesc);
            }
            if (Bdesc["format"] == "grib2") {
                SmartGrib2Select(Bdesc);
            }
        });
        
        
        //Gestion des buttons


        //Gestion dynamique du champ (field) ou des champs (U et V)
        $("#field_v_cloned").prop("disabled", false);
        ManageUandV("_cloned");

        $("#FF_cloned").click(function() {
            ManageUandV("_cloned");
        });
        $("#vecteurs_cloned").click(function() {
            ManageUandV("_cloned");
            if (this.checked) {
                $("#tabs_plot_cloned").tabs("enable", "#vectors");
            } else {
                $("#tabs_plot_cloned").tabs("disable", "#vectors");
            }
        });
        //Attention aux attributs for des clone !
        $("#FF_label_cloned").attr('for','FF_cloned');
        $("#vecteurs_label_cloned").attr('for','vecteurs_cloned');
        $("#update_label_cloned").attr('for','autoupdate_cloned');
        $("#reverse_colormap_label_cloned").attr('for','reverse_colormap_cloned');
        $("#decumul_label_cloned").attr('for','decumul_cloned');
        //On réactive le style des boutons
        $(".boutonlike").buttonset();
        
        
        $("#file_path_cloned").button({
            icons: {primary: 'ui-icon-info'}
        });

        $("#file_path_cloned").click(function() {
            alert(Bdesc["vortex_path"]);
        });
        
        //On vérifie la complétude de la description à chaque changement de valeur
        $(".vortex_cloned").change(function() {
            CheckVortexDescription("_cloned");
        });

        //On demande explicitement la vérification de la description
        $("#checkdescription_cloned")
           //.button()
           .click(function() {
            CheckVortexDescription("_cloned");
        });

        //Cas où on demande la vérification d'existence
        $("#checkexistence_cloned")
           //.button()
           .click(function() {

                    var args_vortex_check_json_cloned = CreateVortexArgs("existence", "cloned");

                    $('#vortex_existence_cloned span').text("Checking...");

                    $.ajax({
                        type: "POST",
                        async: true,
                        url: "/getFile",
                        data: args_vortex_check_json_cloned,
                        contentType: "application/json; charset=utf-8",
                        dataType: "json",
                        success: function(vortexAnswer) {
                            Bdesc["vortex_path"] = vortexAnswer["remotepath"];
                            $('#vortex_existence_cloned span').text(vortexAnswer["existence"]);
                            $("#vortex_existence_cloned").effect("pulsate", 500);
                        },
                        error: function(response) {
                            $('#vortex_existence_cloned span').text("Error");
                        },
                    });
        })
        
        $("#getfile_cloned")
         //.button()
            
         .click(function() {
        
          $('#vortex_get_cloned span').text("Transfering data...");
        
          var args_vortex_get_cloned_json = CreateVortexArgs("get","cloned");

          $.ajax({
            type: "POST",
            async: true,
            url: "/getFile",
            data: args_vortex_get_cloned_json,
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(vortexAnswer) {
                Bdesc["local_path"] = vortexAnswer["localpath"];
                UpdateFields(Bdesc,zoom_spinner); 
                $('#vortex_get_cloned span').text("Current description fetched!");
                //NIOU 2 TEST
                $("#accordion_vortex_cloned").accordion("option", "active", 1);
                makeTabTitle(args_vortex_get_cloned_json,"fileB")
            },
            error: function(response) {
                $('#vortex_get_cloned span').text("Transfer error!");
            }
          });
        });
        
        //On initialise le descriptif du 2e fichier avec le 1er (évite de rapatrier par exemple)
        //puis on actualise l'auto complétion (comme après rapatriement de fichier)
        Bdesc["vortex_path"] = Adesc["vortex_path"];
        Bdesc["grib_param"] = Adesc["grib_param"];
        Bdesc["format"] = Adesc["format"];
        Bdesc["local_path"] = Adesc["local_path"];
        UpdateFields(Bdesc,zoom_spinner);
        
        //On s'occupe de la partie render
        //ON réactive les colormap list

        var listTabs = $( '#filenfield_tab2 .hardtobecloned' );
        var listTabsItems = listTabs.find('a');
        listTabsItems.each( function(i) {
            var original = $(this).attr("href") ;
            $(this).attr("href",original + "_cloned");
            //$(this).attr("id",original + "_cloned");
        })
        
        
        
        InitColorMapList("",currentColorMap);
        InitColorMapList("_cloned",currentColorMap);
        setup_render_labels("_cloned");
        
        /*
        changeSlider(minmax["min"], minmax["max",]"_cloned")
        changeSlider(minmax["min"], minmax["max"],suffixe);
        $("#range1" + suffixe).val($("#slider-range" + suffixe).slider("values", 0) + " to " + $("#slider-range"+ suffixe).slider("values", 1));
        */
        var zemini =  $( "#slider-range" ).slider( "values")[0]
        var zemaxi = $( "#slider-range" ).slider( "values")[1]
        //alert(zemini,zemaxi);
        changeSlider(zemini,zemaxi,"_cloned")
        smart_slider("_cloned")
        $('#slider-range_cloned').on('slidestop', function() {
            desactivateAutoUpdate("_cloned");
            });
            
            
            
        Param_Presets_Init("_cloned")

        
        $('#filesystemOK_cloned').prop("disabled", false);
		$("#filesystemOK_cloned")
		       //.button()
		       .click(function() {
		       Bdesc["local_path"] = $("#filesystem_path_cloned").val().split('\n');
		       UpdateFields(Bdesc,zoom_spinner);         
		        $("#getplot").prop("disabled", false);
		        $("#getplotboth").prop("disabled", false);
		        $("#overlay").prop("disabled", false);
		        $("#difference").prop("disabled", false);
		        $("#accordion_vortex_cloned").accordion("option", "active", 1);
		        makeTabTitleFS(Bdesc["local_path"],"fileB")
			});
        
        
        /*
        $(".hardtobecloned").each(function() {
            //$(this).find('li').each(function(){
            alert("begin")
            var $this = $(this);
            $this.children("a").each(function(i) {
                //$this; // parent li
                //this; // child li
                alert(i,this.attr("href"))
                //var value = $(this).attr('href');
                //alert(value);
                //$(this).attr('href', value + "_cloned");
            });
            });
        */
        

        //Cloner les tabs c'est vraiment une plaie
         
        $( "#modetabs" ).tabs();
         
        listTabs2 = $( '#filenfield_tab2 .hardtobecloned2' );
        listTabsItems2 = listTabs2.find('a');
        listTabsItems2.each( function(i) {
            var original = $(this).attr("href") ;
            $(this).attr("href",original + "_cloned");
            //$(this).attr("id",original + "_cloned");
        })
        
         $( "#modetabs_cloned" ).tabs();
        
        
        //On créé la mécanique de chaque tab
        activate_tabs_plot(Adesc);
        activate_tabs_plot(Bdesc);
        
        //On focus sur l'onglet créé
        //New style
        $("#filenfield_tabs").tabs( "option", "active", 1);
        //Old style
        //$("#filenfield_tabs").tabs( "select", 1);
        
    });

    /////// E N D  O F  C L O N E  E T  A L 





});
// END OF DOCUMENT READY {}


//Gestion des images
function MySlideShow(myclass,nb) {
    $('.'+myclass).slick({
      arrows: false,
      dots: true,
      slidesToShow: nb,
      slidesToScroll: nb,
      infinite: true,
      speed: 0,
      cssEase: 'linear',
      autoplaySpeed: 1000,
      afterChange: function (slickSlider, i) {
                        //remove all active class
                        $('.' + myclass + '_nav .slick-slide').removeClass('slick-active');
                        //set active class for current slide
                        $('.' + myclass + '_nav .slick-slide').eq(i).addClass('slick-active');
                    }
      });
      
      $('.slider-nav').slick({
      	slidesToShow: 3,
      	slidesToScroll: 1,
      	asNavFor: '.slider-for',
      	dots: true,
      	centerMode: true,
      	focusOnSelect: true
      	
      });
      
      
      
             
$('.slideshow').slick({
	infinite: true,
});	
}

/*
$('.pause').on('click', function() {
    $('.slider')
        .slick('slickPause')
        .slick('slickSetOption', 'pauseOnDotsHover', false);
});

$('.play').on('click', function() {
    $('.slider')
        .slick('slickPlay')
        .slick('slickSetOption', 'pauseOnDotsHover', true);
});
*/







function disableAccordion() {
    //Dynamically disables a part of an accordion (class = optiondisable)
    // Add the class ui-state-disabled to the headers that you want disabled
    $(".optiondisable").addClass("ui-state-disabled");
    // Now the hack to implement the disabling functionnality
    var myaccordion = $("#accordion_plot").data("uiAccordion");
    myaccordion._std_clickHandler = myaccordion._clickHandler;
    myaccordion._clickHandler = function(event, target) {
        var clicked = $(event.currentTarget || target);
        if (!clicked.hasClass("ui-state-disabled")) this._std_clickHandler(event, target);
    };
};

function activateZoom() {
    //$("#getplot").toggle( "shake",2000 );
    if (active_zoom.checked == false) {
        $("#active_zoom").prop('checked', true).button('refresh');
        $("#active_zoom_label").text("Zooming is on");
        //$("#active_zoom_label").effect( "highlight", {color: 'red'}, 500 );
        $("#active_zoom_label").effect("pulsate", 500);
        //$("#active_zoom_label").fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);
        //$("#getplot").toggle( "highlight" );
    }
};

function desactivateAutoUpdate(suffixe) {
    //if (autoupdate.checked == true) {
    if ($("#autoupdate" + suffixe).prop('checked') == true) {    
        $("#autoupdate" + suffixe).prop('checked', false).button('refresh');
        $("#update_label" + suffixe).text("AutoUpdate is off");
        //$("#update_label").effect( "highlight", {color: 'red'}, 500 );
        $("#update_label" + suffixe).effect("pulsate", 500);
    }
};

function callback() {
    setTimeout(function() {
        $("#effect:visible").removeAttr("style").fadeOut();
    }, 1000);
};

function changeSlider(zemin, zemax,suffixe) {

    //On fait un slider step intelligent en fonction de l'ordre de grandeur du champ
    //Le step = 1/100 de l'amplitude du champ
    delta = (zemax - zemin) 
    delta_10power = delta.toExponential().split("e")[1]
    delta_slider = Math.pow(10, delta_10power) / 100.
    
    $("#slider-range" + suffixe).slider({
        orientation: "horizontal",
        range: true,
        //step: (zemax - zemin) / 100.,
        step: delta_slider,
        min: zemin - Math.abs(0.2 * zemin),
        max: zemax + Math.abs(0.2 * zemax),
        values: [zemin, zemax],
        slide: function(event, ui) {
            $("#range1" + suffixe).val(ui.values[0] + " to " + ui.values[1]);
        }
    });
    //$( "#range1" ).val( $( "#slider-range" ).slider( "values", 0 ) + " - " + $( "#slider-range" ).slider( "values", 1 ) );

}

function ManageUandV(suffixe) {
    //Gestion du passage entre un seul "field" et l'apparition de 2 "U et V"
    //selon l'activation des boutons FF et vecteurs

    if ($("#FF" + suffixe).is(':checked') || $("#vecteurs" + suffixe).is(':checked')) {
        $("#field_v" + suffixe).show();
        $("#field_v_label" + suffixe).show();
        //$( "#field_v" ).prop("disabled", false );
        $("#field_v_label" + suffixe).text("V :");
        $("#field_label" + suffixe).text("U :");
    } else {
        $("#field_v" + suffixe).hide();
        $("#field_v_label" + suffixe).hide();
        //$( "#field_v" ).prop("disabled", true );
        $("#field_label" + suffixe).text("Field :");
    }
};




function UpdateCheckGet() {

    //On désactive le check et get si la description n'est pas complete

    //Upalert($('#vortex_description span').text()[0])

    //Il n'y a aucun False dans la description des ressources
    if ($('#vortex_description span').text().match('False') == null) {
        $('#checkexistence').prop("disabled", false);
        $('#getfile').prop("disabled", false);
        $("#checkexistence").effect("pulsate", 500);
        $("#getfile").effect("pulsate", 500);
        $("#file_path").prop("disabled", false);
        $("#file_path").effect("pulsate", 500);

    } else {
        $('#checkexistence').prop("disabled", true);
        $('#getfile').prop("disabled", true);

    }
};


//passer A ou B
function CheckVortexDescription(suffixe) {

        var args_vortex_json = CreateVortexArgs("description",suffixe);

        $.ajax({
            type: "POST",
            async: true,
            url: "/getFile",
            data: args_vortex_json,
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(vortexAnswer) {
                //On garde la valeur precedente en memoire
                previous_desc = $('#vortex_description' + suffixe + ' span').text()
                
                //A priori la completude est identique pour tous les éléments => on n'affiche que le premier (vs l'existenc)
                $('#vortex_description' + suffixe + ' span').text(vortexAnswer["description"][0]);
                //On renseigne le bouton de localisation de manière moche
                if (suffixe == "") { Adesc["vortex_path"] = vortexAnswer["remotepath"];}
                else { Bdesc["vortex_path"] = vortexAnswer["remotepath"];}

                //Blink and update only if change
                if (vortexAnswer["description"] != previous_desc) {
                    $("#vortex_description" + suffixe).effect("pulsate", 500);
                    UpdateCheckGet();
                }
            },
            error: function(response) {
                console.log(response.responseText);
            },
        });
    };
    
    
function updateCacheSize() {
    
        //var nada = null;
        var args_vortex_json = CreateVortexArgs("description","");
        

        $.ajax({
            type: "POST",
            async: true,
            url: "/getCacheSize",
            data: args_vortex_json,
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(myAnswer) {
                //On garde la valeur precedente en memoire
                //previous_desc = $('#vortex_description' + suffixe + ' span').text()
                //alert(myAnswer);
                $('#cachesize span').text(myAnswer);
               },
            error: function(response) {
                console.log("Erreur !!" + response.responseText);
            },
        });
    };

function getGeometries(suffixe) {

        $.ajax({
            type: "POST",
            async: true,
            url: "/getGeometries",
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(myAnswer) {
                //On garde la valeur precedente en memoire
                //previous_desc = $('#vortex_description' + suffixe + ' span').text()
                //alert(myAnswer[1]);
                MakeAutoComplete("#geometry"+suffixe, myAnswer);
               },
            error: function(response) {
                console.log("Erreur !!" + response.responseText);
            },
        });
    };    

//Passer A ou B ["suffixe"] ?
function CreateVortexArgs(mymode,suffixe) {

        //Seuls les champs non vides de class vortex+suffixe seront passes en arguement
        var args_vortex_create = {};
        //if (suffixe != "") {suffixe = "." + suffixe;}
        
        $('.vortex'+suffixe).each(function(i, obj) {
            if ($(this).val() != "") {
                ZeCleanID = $(this).attr("id").replace('_cloned','') //car l'id est passé comme clé au dict pour vortex : vapp, etc.
                args_vortex_create[ZeCleanID] = $(this).val();
            }
        });

        //On rajoute le mode d'exécution demandé
        args_vortex_create["request_mode"] = mymode
        //On rajoute le fichier d'origine
        if (suffixe == "") {
            args_vortex_create["fromid"] = "A";
        }
        else {
            args_vortex_create["fromid"] = "B";
        }

        var args_vortex_create_json = JSON.stringify(args_vortex_create);
        
        return args_vortex_create_json ;

    };



function MakeAutoComplete(id,liste) {

            //try { console.log("Destroy!!");  $("#field" + oneDesc["suffixe"]).autocomplete("close");}
            //catch (error) {console.log(error)}
            
            $(id).autocomplete({
                source: liste,
                minLength: 0,
                delay: 0,
            });
            
            //Autocomplete list is refresh onFocus and onClick
            $(id).focus(function() {
            	$(this).autocomplete('search', $(this).val());
            });
            $(id).click(function() {
              $(this).autocomplete('search', $(this).val());
            });
            
}

//Mise à jour auto de la liste des champs quand on change le fichier d'entrée  
function UpdateFields(oneDesc,zoom_spinner) {
    // Arguments : chemin absolu du fichier
    var args_fields = {};
    args_fields["file"] = oneDesc["local_path"][oneDesc["local_path"].length-1]; 
    
    //Cas du clonage avant remplissage des champs : on vérifie que le fichier existe bien! 
    
    var args_fields_json = JSON.stringify(args_fields);

    //Appel pour obtenir la liste des champs en JSON, utilisée dans une boite d'autocomplétion
    $.ajax({
        type: "POST",
        async: true,
        url: "/getfieldsasjson",
        data: args_fields_json,
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        success: function(fields) {

            //Reco empirique de FA vs GRIB1 vs GRIB2 -> attribution d'un type et (dés)activation des champs supplémentaires
            if (fields[0]["typeOfLevel"] != undefined ) {
                oneDesc["format"] = "grib1";
                $(".gribonly").prop("disabled", false);
                $("#TypeOfLevel" + oneDesc["suffixe"]).prop("disabled", false);
                $("#Level" + oneDesc["suffixe"]).prop("disabled", false);
            }
            if (fields[0]["discipline"] != undefined ) {
                oneDesc["format"] = "grib2";
                $(".gribonly").prop("disabled", false);
                $("#TypeOfLevel" + oneDesc["suffixe"]).prop("disabled", false);
                $("#Level" + oneDesc["suffixe"]).prop("disabled", false);
            }
            if (fields[0]["typeOfLevel"] == undefined && fields[0]["discipline"] == undefined) {
                oneDesc["format"] = "FA";
                $(".gribonly").prop("disabled", true);
                $("#TypeOfLevel" + oneDesc["suffixe"]).prop("disabled", true);
                $("#Level" + oneDesc["suffixe"]).prop("disabled", true);
            }

            if (oneDesc["format"] == "grib1") {
                oneDesc["grib_param"] = fields;
                SmartGribSelect(oneDesc);
            }
            else if (oneDesc["format"] == "grib2") {
                oneDesc["grib_param"] = fields;
                SmartGrib2Select(oneDesc);
            }
            else if (oneDesc["format"] == "FA") {
                var ListeParam = [];
                $.each(fields, function(index, field) {
                    ListeParam.push(field);
                })
                MakeAutoComplete("#field" + oneDesc["suffixe"],ListeParam);
                MakeAutoComplete("#field_v" + oneDesc["suffixe"],ListeParam);
            }
            $("#accordion_actions").accordion("option", "active", 0);
            
            //On met la lecture du domaine à ce niveau pour éviter les conflits de lecture simultané de FA
            //On ne met pas à jour le zoom pour le fichier B (= cloned)
            if (oneDesc["suffixe"] == "") {
                GetCoordinates(zoom_spinner,fields[0]);
            }
            
        },
        error: function(response) {
            console.log(response.responseText + " File is probably non existent ");
        }
    });
}


function GetGribKey(ID,number) {
    //En accord avec le remplissage, on retourne la n-ème valeur du champ ID (n=0, 1, ..)
    champs = $(ID).val().split('/') ;
    //console.log(champs)
    try {
        return champs[number].trim();
    }
    catch (err) {
    	return $(ID).val().split('/')[0].trim();
    }
}


//Gives unique description of parameter (grib1)
// with only "indicatorOfParameter", "table2Version","level","typeOfLevel"
function giveme_grib1_handgrip(suffixe) {
    
    handgrip = {}
    handgrip["indicatorOfParameter"] = parseInt(GetGribKey("#field"+suffixe,0));
    handgrip["table2Version"] = parseInt(GetGribKey("#field"+suffixe,1));
    handgrip["level"] = parseInt(GetGribKey("#Level"+suffixe,0));
    handgrip["typeOfLevel"] = GetGribKey("#TypeOfLevel"+suffixe,0);
    //console.log ("handgrip min max :" + handgrip);
    return handgrip;
    
}

//Gives unique description of parameter (grib2)
// with only 'parameterNumber','parameterCategory','discipline','typeOfFirstFixedSurface'[-typeOfStatisticalProcessing],[ 'level' OR  'scaledValueOfCentralWaveNumber']
//Web interface:
//     #field           [ 'parameterNumber','name','parameterCategory',"shortname",'discipline' ];                    ]
//     #TypeOfLevel     [ 'typeOfFirstFixedSurface' OR typeOfStatisticalProcessing] 
//     #Level           [ 'level' OR  'scaledValueOfCentralWaveNumber']
function giveme_grib2_handgrip(suffixe) {
    
    handgrip = {}
    handgrip["parameterNumber"] = parseInt(GetGribKey("#field"+suffixe,0));
    handgrip["parameterCategory"] = parseInt(GetGribKey("#field"+suffixe,2));
    handgrip["discipline"] = parseInt(GetGribKey("#field"+suffixe,4));
    //ISP and PTYPE (bad) hacks
    if (handgrip["parameterNumber"] == 7) {
      handgrip["scaledValueOfCentralWaveNumber"] = parseInt(GetGribKey("#Level"+suffixe,0));
    }
    else if (handgrip["parameterNumber"] == 19) {
      handgrip["typeOfStatisticalProcessing"] = parseInt(GetGribKey("#TypeOfLevel"+suffixe,1));
      handgrip["level"] = parseInt(GetGribKey("#Level"+suffixe,0));
    }
    else {
      handgrip["level"] = parseInt(GetGribKey("#Level"+suffixe,0));
      handgrip["typeOfFirstFixedSurface"] = parseInt(GetGribKey("#TypeOfLevel"+suffixe,0));
    }
    
    return handgrip;
    
}


//dict A ou B (avec suffixe)
//from epygram list of field to text boxes...
//GRIB1 only !
function SmartGribSelect(oneDesc) {
    
            var ListeParam = [];
            var ListeTypeLevel = [];
            var ListeLevel = [];
            
            //On récupère la 1ère valeur de chaque champ
            var user_Param = GetGribKey("#field" + oneDesc["suffixe"],0);
            var user_table2 = GetGribKey("#field" + oneDesc["suffixe"],1);
            var user_TypeOfLevel = GetGribKey("#TypeOfLevel" + oneDesc["suffixe"],0);
            //.val().split('/')[0].trim();
            var user_Level = GetGribKey("#Level" + oneDesc["suffixe"],0);

            //console.log(user_Param,user_TypeOfLevel,user_Level,user_table2,field["table2Version"])

                //On parcourt la liste des clés du grib et on filtre 
                //si valeurs déjà saisies dans les autres champs -> on filtre avec
                //sinon on prend tout
                //Exemple pour le paramètre : conditions sur le type de niveau et le niveau
                $.each(oneDesc["grib_param"], function(index, field) {
                    if ((user_TypeOfLevel == "" || user_TypeOfLevel == field["typeOfLevel"]) &&
                    (user_Level == "" || user_Level == field["level"])) {
                    description = field["indicatorOfParameter"] + " / " + field["table2Version"] + " / " + field["shortName"] + " / " + field["name"];
                        //On n'ajoute que si pas déjà present
                        if (jQuery.inArray(description,ListeParam) == -1) {
                            ListeParam.push(description);
                        }
                    }
                    
                //TYPE OF LEVEL : idem
                //indicatorOfParameter
                    if ((user_Param == "" || user_Param == field["indicatorOfParameter"]) &&
                    (user_table2 == "" || user_table2 == field["table2Version"]) &&
                    (user_Level == "" || user_Level == field["level"])) {
                        description = field["typeOfLevel"] + " / " + field["indicatorOfTypeOfLevel"];
                        if (jQuery.inArray(description,ListeTypeLevel) == -1) {
                            ListeTypeLevel.push(description);
                        }    
                    }
                    
               //LEVEL
                    if ((user_TypeOfLevel == "" || user_TypeOfLevel == field["typeOfLevel"]) &&
                    (user_Param == "" || user_Param == field["indicatorOfParameter"])) {
                        description = field["level"].toString();}
                        if (jQuery.inArray(description,ListeLevel) == -1) {
                            ListeLevel.push(description);
                        }
                   });

                MakeAutoComplete("#field" + oneDesc["suffixe"],ListeParam);
                MakeAutoComplete("#field_v" + oneDesc["suffixe"],ListeParam);
                MakeAutoComplete("#TypeOfLevel" + oneDesc["suffixe"],ListeTypeLevel);
                MakeAutoComplete("#Level" + oneDesc["suffixe"],ListeLevel);
	
}


//We gather them this way:
//          [ 'parameterNumber','name','parameterCategory',"shortname",'discipline' ];                    ]
//          [ 'typeOfFirstFixedSurface' - 'typeOfStatisticalProcessing' ] 
//          [ 'level' OR  'scaledValueOfCentralWaveNumber']
//to produce a handgrip for epygram plot with only:
//          'parameterNumber','parameterCategory','discipline','typeOfFirstFixedSurface', ['level' OR  'scaledValueOfCentralWaveNumber']

function SmartGrib2Select(oneDesc) {
    
            var ListeParam = [];
            var ListeTypeLevel = [];
            var ListeLevel = [];
            
            console.log("SmartGrib2Select")
            
            MakeAutoComplete("#field_v" + oneDesc["suffixe"],"");
            MakeAutoComplete("#TypeOfLevel" + oneDesc["suffixe"],"");
            MakeAutoComplete("#Level" + oneDesc["suffixe"],"");
            
            var subSamples = oneDesc["grib_param"]
            
            //Pré filtrage
            //On récupère la 1ère valeur de chaque champ et on filtre avec si non nulle
            var user_Param = GetGribKey("#field" + oneDesc["suffixe"],0);
            if (user_Param !== "") {
              var test = subSamples.filter(function(d) { return d.parameterNumber == user_Param; })
              subSamples = test;
            }
            var user_table2 = GetGribKey("#field" + oneDesc["suffixe"],2);
            if (user_table2 !== "") {
              var test = subSamples.filter(function(d) { return d.parameterCategory == user_table2; })
              subSamples = test;
            }
            var user_TypeOfLevel = GetGribKey("#TypeOfLevel" + oneDesc["suffixe"],0);
            if (user_TypeOfLevel !== "" && user_TypeOfLevel !== "undefined" ) {
              var test = subSamples.filter(function(d) { return d.typeOfFirstFixedSurface == user_TypeOfLevel; })
              subSamples = test;
            }
            var user_Level = GetGribKey("#Level" + oneDesc["suffixe"],0);
            if (user_Level !== "") {
              test = subSamples.filter(function(d) { return d.level == user_Level || d.scaledValueOfCentralWaveNumber == user_Level; })
              subSamples = test;
            }
            
              //On parcourt la liste des clés du grib et on remplit l'autocomplétion (avec le préfiltre cela est très rapide)
              //si valeurs déjà saisies dans les autres champs -> on filtre avec
              //sinon on prend tout
              //Exemple pour le paramètre : conditions sur le type de niveau et le niveau
              //$.each(oneDesc["grib_param"], function(index, field) {
              $.each(subSamples, function(index, field) {  
                  
                    //Dealing with (missing) level in ISP case... and add typeOfFirstFixedSurface for ptype
                    try {
                          levelasgrib1 = field["level"].toString();
                        }
                        catch(error) {
                          console.log("ISP initialization " + field["scaledValueOfCentralWaveNumber"].toString())
                          levelasgrib1 = field["scaledValueOfCentralWaveNumber"].toString();                          
                        }
                    
                    
                  //TYPE OF PARAM
                    if ((user_TypeOfLevel == "" || user_TypeOfLevel == "undefined" || user_TypeOfLevel == field["typeOfFirstFixedSurface"]) &&
                    (user_Level == "" || user_Level == levelasgrib1)) {
                    description = field["parameterNumber"] + " / " + field["name"] + " / " + field["parameterCategory"] + " / " + field["shortName"] + " / " + field["discipline"];
                        //On n'ajoute que si pas déjà present
                        if (jQuery.inArray(description,ListeParam) == -1) {
                            ListeParam.push(description);
                        }
                    }
   
                //TYPE OF LEVEL : idem
                    if ((user_Param == "" || user_Param == field["parameterNumber"]) &&
                    (user_table2 == "" || user_table2 == field["parameterCategory"]) &&
                    (user_Level == "" || user_Level == levelasgrib1)) {
                        try {
                            description = field["typeOfFirstFixedSurface"] + " / " + field["typeOfStatisticalProcessing"] 
                        }
                        catch(error) {
                          description = field["typeOfFirstFixedSurface"] + " / " + " Nada"
                        }
                        if (jQuery.inArray(description,ListeTypeLevel) == -1) {
                            ListeTypeLevel.push(description);
                        }    
                    }
                    

               //LEVEL
                    if ((user_TypeOfLevel == "" || user_TypeOfLevel == "undefined" || user_TypeOfLevel == field["typeOfFirstFixedSurface"]) &&
                    (user_Param == "" || user_Param == field["parameterNumber"])) {
                        try {
                          description = field["level"].toString();
                        }
                        catch(error) {
                          console.log("ISP level " + field["scaledValueOfCentralWaveNumber"].toString())
                          description = field["scaledValueOfCentralWaveNumber"].toString();                          
                        }
                        
                        if (jQuery.inArray(description,ListeLevel) == -1) {
                            ListeLevel.push(description);
                        }
                    }
                   });

                MakeAutoComplete("#field" + oneDesc["suffixe"],ListeParam);
                MakeAutoComplete("#field_v" + oneDesc["suffixe"],ListeParam);
                MakeAutoComplete("#TypeOfLevel" + oneDesc["suffixe"],ListeTypeLevel);
                MakeAutoComplete("#Level" + oneDesc["suffixe"],ListeLevel);  
}












function Param_Presets_Init(suffixe) {
    
    /*
    for (var zeparam in param_presets) {
        $("#param_" + zeparam + suffixe).click(function() {
            alert(zeparam);
            alert(param_presets[zeparam]["field"]);
        //ParameterPreset(key,suffixe);
        for (var key in param_presets[param]) {
            
            $("#" + key  + suffixe).val(param_presets[param][key]);
        }
    
        //$("#field").val(param_presets[param]["field"]);
        //$("#TypeOfLevel").val("heightAboveGround / 105");
        //$("#Level").val("2");
        $("#FF"+suffixe).prop("checked", param_presets[param]["FF"]);
        $("#FF"+suffixe).trigger('click');

        $("#Vecteurs"+suffixe).prop("checked", param_presets[param]["Vecteurs"]);
        ManageUandV(suffixe);
    })
    }
    */
    $("#param_hu2m"+suffixe).click(function() {
        ParameterPreset("hu2m",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
    })
    
    $("#param_rr"+suffixe).click(function() {
        ParameterPreset("rr",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
    })
    
    $("#param_none"+suffixe).click(function() {
        ParameterPreset("none",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
    })
    
    $("#param_ff10m"+suffixe).click(function() {
        ParameterPreset("ff10m",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
    })
    
    $("#param_isp"+suffixe).click(function() {
        ParameterPreset("isp",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
        })
    
    $("#param_tpw850"+suffixe).click(function() {
        ParameterPreset("tpw850",suffixe);
        if (suffixe == "") {SmartGribSelect(Adesc);}
        else { SmartGribSelect(Bdesc); }
        })
        
    $("#param_A2B"+suffixe).click(function() {
        ParameterPreset("A2B",suffixe);
        //Trick !
        if (suffixe == "") {SmartGribSelect(Bdesc);}
        else { SmartGribSelect(Adesc); }
        })
    
    }        

function ParameterPreset(param,suffixe) {
    
    //Recopie des paramètres de A vers B et réciproquement
    
    var suffixe2 = suffixe //suffixe2 : là ou on initialise les valeurs, different de suffixe seulement en cas de copie
    //Cas de la copie
    if (param == "A2B") {
        if (suffixe == "") {suffixe2 = "_cloned";}
        if (suffixe == "_cloned") {suffixe2 = "";}
    }
    for (var key in param_presets[param]) {
        //Cas de la copie
        if (param == "A2B") {
            param_presets[param][key] = $("#" + key  + suffixe).val();
            //alert (key + param_presets[param][key]);
        }
        if (param_presets[param][key] != "" || param == "none") {
            $("#" + key  + suffixe2).val(param_presets[param][key]);
        }
    }
        
        //Cohérence des boutons FF et vecteurs
        if ($("#FF"+suffixe2).prop("checked") != param_presets[param]["FF"]) {
        	$("#FF"+suffixe2).trigger('click');
        }
        if ($("#Vecteurs"+suffixe2).prop("checked") != param_presets[param]["Vecteurs"]) {
        	$("#Vecteurs"+suffixe2).trigger('click');
        }
        
        //$("#Vecteurs"+suffixe2).prop("checked", param_presets[param]["Vecteurs"]);
        ManageUandV(suffixe2);


}

//Passer A ou B["suffixe"]            
function VortexAutoCompleteInit(suffixe) {
        
        $.each(vortex_autocomplete, function(key, val) {
            MakeAutoComplete("#" + key + suffixe,val);
            //$("#" + key + suffixe).autocomplete({
            //    source: vortex_autocomplete[key]
        });
}            

//Idem
function VortexPresetInit(suffixe) {

    $("#vortex_preset_resources" + suffixe).selectmenu({
        change: function(event, data) {
            for (var key in vortex_presets_resources[data.item.value]) {
                $("#" + key  + suffixe).val(vortex_presets_resources[data.item.value][key]);
            }
            //Bogus trigger for lauching description check just once
            $("#vapp" + suffixe).trigger('change');
        }
    });

    $("#vortex_preset_providers" + suffixe).selectmenu({
        change: function(event, data) {
            for (var key in vortex_presets_providers[data.item.value]) {
                $("#" + key + suffixe).val(vortex_presets_providers[data.item.value][key]);
            }
            //Bogus trigger for lauching description check just once
            $("#vapp" + suffixe).trigger('change');
        }
    });

}


function AutoUpdateMinMax(args_minmax_json,suffixe) {
//Calcul du min max et création d'un slider ad hoc		
            //Mise à jour du min/max si les valeurs ne sont pas figées
                $.ajax({
                    type: "POST",
                    async: false,
                    url: "/getminmax",
                    data: args_minmax_json,
                    contentType: "application/json; charset=utf-8",
                    dataType: "json",
                    success: function(minmax) {
                        changeSlider(minmax["min"], minmax["max"],suffixe);
                        $("#range1" + suffixe).val($("#slider-range" + suffixe).slider("values", 0) + " to " + $("#slider-range"+ suffixe).slider("values", 1));
                    },
                    error: function(response) {
                        alert(response.responseText);
                    }
                });

}

function giveme_subzone(suffixe){
    
    var subzone = "CI";
    if ($("#subzone").val() == "C") {
        subzone = "C";
        }
    if ($("#subzone").val() == "CIE") {
        subzone = null;           } 
    //Cas des gribs : subzone prohibé !
    if (suffixe == "" & (Adesc["format"] == "grib1" || Adesc["format"] == "grib2")) {
        subzone = null ;            }
    if (suffixe == "_cloned" & (Bdesc["format"] == "grib1" || Bdesc["format"] == "grib2")) {
        subzone = null ;            }
    return subzone;
}


function giveme_args_minmax(Onedesc) {
            
            
            
            var args_minmax_tmp = {};
            
            var subzone = giveme_subzone(Onedesc["suffixe"]);
            
            if (Onedesc["format"] == "grib1") {
                args_minmax_tmp["field"] = giveme_grib1_handgrip(Onedesc["suffixe"]);
                //Pour la composante V, tout idem sauf indicatorOfParameter
                args_minmax_tmp["field_v"] = giveme_grib1_handgrip(Onedesc["suffixe"]) ;
                try {
                    args_minmax_tmp["field_v"]["indicatorOfParameter"] = parseInt(GetGribKey("#field_v",0));                    
                }
                catch (err) {}
            }
            
            else if (Onedesc["format"] == "grib2") {
                args_minmax_tmp["field"] = giveme_grib2_handgrip(Onedesc["suffixe"]);
                //Pour la composante V, tout idem sauf indicatorOfParameter
                args_minmax_tmp["field_v"] = giveme_grib2_handgrip(Onedesc["suffixe"]) ;
                try {
                    args_minmax_tmp["field_v"]["parameterNumber"] = parseInt(GetGribKey("#field_v",0));                    
                }
                catch (err) {}
            }
            
            else {
                args_minmax_tmp["field"] = $("#field"+Onedesc["suffixe"]).val();
                args_minmax_tmp["field_v"] = $("#field_v"+Onedesc["suffixe"]).val();
            }
            
            if ($("#operation" + Onedesc["suffixe"]).val() != "") {
                args_minmax_tmp["operation"] = $("#operation" + Onedesc["suffixe"]).val().split(/,/);
            }

            //JUST THE LAST ONE FOR THE MOMENT
            args_minmax_tmp["file"] = Onedesc["local_path"][Onedesc["local_path"].length-1]; 
            args_minmax_tmp["subzone"] = subzone;
            args_minmax_tmp["FF"] = $("#FF"+Onedesc["suffixe"]).is(':checked');
            
            
            return args_minmax_tmp;
            
}


function get_common_args(args_plot) {

            //Gestion à la main des défauts car pb d'affichage quand le défaut est choisi dans le code html
            var specificproj = "";
            if ($("#specificproj").val() != "") {
                specificproj = $("#specificproj").val();
            }
            var llgridstep = ""; //syntaxe plus ?? importante pour faire passer un string vide au lieu d'un réel de manière transparente
            if ($("#llgridstep").val() != "") {
                llgridstep = parseFloat($("#llgridstep").val());
            }
            var meridians = "auto"; //syntaxe plus ?? importante pour faire passer un string vide au lieu d'un réel de manière transparente
            if ($("#meridians").val() != "") {
                if ($("#meridians").val() == "none") {
                    meridians = null ;
                }
                else {
                    meridians = $("#meridians").val();
                }
            }
            var parallels = "auto"; //syntaxe plus ?? importante pour faire passer un string vide au lieu d'un réel de manière transparente
            if ($("#parallels").val() != "") {
                if ($("#parallels").val() == "none") {
                    parallels = null ;
                }
                else {
                    parallels = $("#parallels").val();
                }
                
            }
            
            var gisquality = "i";
            if ($("#gisquality").val() != "") {
                gisquality = $("#gisquality").val();
            }

            var monzoom = "";
            if ($("#active_zoom").is(':checked')) {
                monzoom = {};
                monzoom["lonmin"] = parseFloat($("#zoom_west").val());
                monzoom["lonmax"] = parseFloat($("#zoom_east").val());
                monzoom["latmin"] = parseFloat($("#zoom_south").val());
                monzoom["latmax"] = parseFloat($("#zoom_north").val());
            }

            //A and B common
            //args_plot["subzone"] = giveme_subzone();
            args_plot["gisquality"] = gisquality;
            args_plot["specificproj"] = specificproj;
            //args_plot["center_cmap_on_0"] = $("#center_cmap_on_0").val();
            args_plot["drawrivers"] = ($("#drawrivers").is(':checked'));
            args_plot["monzoom"] = monzoom;
            //args_plot["llgridstep"] = llgridstep;
            args_plot["meridians"] = meridians;
            args_plot["parallels"] = parallels;
            args_plot["departments"] = ($("#drawdepartments").is(':checked'));
            args_plot["dpi"] = parseInt($("#dpi").val());
            args_plot["getcode"] = "True" //on force par défaut à sortir les arguments de plotfield $("#getcode").is(':checked') ;
            args_plot["bluemarble"] = parseFloat($("#bluemarble").slider("value"));
            //args_plot["contourcolor"] = $("#contourcolor").val();

            

return args_plot;

}


function get_specific_args(my_args_plot,suffixe,filename) {

            //Initialisations des dics non encore utilisées
            var toCheck = ["FF","min","max","levelsnumber","colormap","pointsize","vectors_subsampling","contourcolor","vectorcolor","vecteurs","decumul","subzone","operation"]
            for(var i=0, len=toCheck.length; i < len; i++){
                if (!(toCheck[i] in my_args_plot)) {
                my_args_plot[toCheck[i]]={};
            }
            }
            
            var colormap = $("#ColorMapList" + suffixe).val()
            if (colormap == "") {
                colormap = "jet";
            }
            if ($("#reverse_colormap" + suffixe).is(':checked')) {
                colormap = colormap + "_r";
            }

            my_args_plot["min"][filename] = parseFloat($("#slider-range" + suffixe).slider("values", 0));
            my_args_plot["max"][filename] = parseFloat($("#slider-range" + suffixe).slider("values", 1));
            my_args_plot["levelsnumber"][filename] = parseInt($("#levelsnumber1" + suffixe).val());
            my_args_plot["colormap"][filename] = colormap;
            my_args_plot["pointsize"][filename] = parseInt($("#pointsize" + suffixe).val());
            my_args_plot["vectors_subsampling"][filename] = parseInt($("#vectors_subsampling" + suffixe).val());
            my_args_plot["contourcolor"][filename] = $("#contourcolor" + suffixe).val();
            my_args_plot["vectorcolor"][filename] = $("#vectorcolor" + suffixe).val();
            my_args_plot["subzone"][filename] = giveme_subzone(suffixe);
            
            //On n'envoie l'opération que si non nulle
            if ($("#operation" + suffixe).val() != "") {
                my_args_plot["operation"][filename] = $("#operation" + suffixe).val().split(/,/);
            }
            //alert(my_args_plot["operator"][filename].split(/,/));
            
            
            
            //DONE
            
            //if (!("FF" in my_args_plot)) {
            //    my_args_plot["FF"]={};
            //}
            my_args_plot["FF"][filename] = $("#FF" + suffixe).is(':checked');
            
            
            //my_args_plot["FF"][filename] = $("#FF" + suffixe).is(':checked');
            my_args_plot["vecteurs"][filename] = $("#vecteurs" + suffixe).is(':checked');
            
            my_args_plot["decumul"][filename] = ($("#decumul" + suffixe).is(':checked'));
            
         return my_args_plot;

}   

function GetCoordinates (zoom_spinner,champ) {
                //Appel pour obtenir les coordonnées du domaine pour préremplir le zoom
                //Uniquement si le zoom n'est pas déjà activé... (pour ne pas perdre les valeurs
                //en changeant de fichier
                var args_domain = {};
                args_domain["file"] = Adesc["local_path"][0];
                args_domain["champ"] = champ;
                var args_domain_json = JSON.stringify(args_domain);

                if ($("#active_zoom").is(':checked')) {
                    //new_pickle = false;
                } else {
                    basemap_pickle_name = generateUUID()+"_bm.cPickle";
                    $.ajax({
                        type: "POST",
                        async: true,
                        url: "/getdomain",
                        data: args_domain_json,
                        contentType: "application/json; charset=utf-8",
                        dataType: "json",
                        success: function(domaine) {
                            zoom_spinner["west"].spinner("value", domaine['lonmin']);
                            zoom_spinner["east"].spinner("value", domaine['lonmax']);
                            zoom_spinner["south"].spinner("value", domaine['latmin']);
                            zoom_spinner["north"].spinner("value", domaine['latmax']);

                        },
                        error: function(response) {
                            alert("Domain inquiry pb !");
                            console.log(response.responseText+"Domain inquiry pb !");
                        }
                    });
                };
}   

function plot2html(targetDiv, previousDiv, response, number) {
//function plot2html(response) {

    //From a targetDiv, a previousDiv, a list of figure and a number of figures to show at the same time
    //render the html code for a slideshow
    //and presents it inside a dialog box
    //that is slightly offset from the previousDiv
    
    mywidth = 700;
    if (number == 2) { mywidth = 1000 ; }

    var liste = "";
    slidenumber = slidenumber + 1;
    var slideshow_class = "slideshow" + slidenumber.toString();

    $.each(response, function(index, value) {
        liste = liste + '\n<div>\n<a href="' + value + '" target=_blank><img style="max-width:95%;" src=' + value + '></a>\n</div>\n';
    });

    var htmlCode = '<div class="' + slideshow_class + '">\n' + liste + "</div>"
    
    var htmlCode_nav = '<div class="' + slideshow_class + '_nav">\n' + liste + "</div>"
    
    var autoplay_html = '<div class="lower"> <button class="prev ' + slideshow_class + '">Prev</button> <button class="play ' + slideshow_class + '">Play</button> <button class="pause ' + slideshow_class + '">Pause</button> <button class="next ' + slideshow_class + '">Next</button> </div>'

    targetDiv.html(htmlCode + autoplay_html); //htmlCode_nav
    
    $('#graphicplot0').append(targetDiv);

    targetDiv.dialog({
        width: mywidth,
        height: "auto",
        title: label,
        position: {
            my: "right top",
            at: "right+10 top-10",
            of: previousDiv
        }
    });
    targetDiv.dialogExtend({
        "maximizable" : true,
        "minimizable" : true,
        "dblclick" : "collapse",
        "icons" : { "maximize" : "ui-icon-arrow-4-diag" }
      });

    MySlideShow(slideshow_class, number);
    
    //set active class to first slide
    /*
    $('.' + slideshow_class + '_nav .slick-slide').eq(0).addClass('slick-active');
    
    $('.' + slideshow_class + '_nav').slick({
    slidesToShow: 4,
    slidesToScroll: 1,
    asNavFor: '.' + slideshow_class,
    autoplay: false,
    dots: true,
    centerMode: true,
    focusOnSelect: true,
    vertical: false
        
    });
    
    $('.' + slideshow_class + '_nav').on('mouseenter', '.slick-slide', function (e) {
	var $currTarget = $(e.currentTarget), 
    	index = $currTarget.data('slick-index'),
        slickObj = $('.' + slideshow_class).slick('getSlick');
    
    slickObj.slickGoTo(index);
    });
    */
    
    $('.pause.' + slideshow_class).on('click', function() {
    $("." + slideshow_class)
        .slick('slickPause')
        .slick('slickSetOption', 'pauseOnDotsHover', false);
        
    });
    
    $('.play.' + slideshow_class).on('click', function() {
    $('.' + slideshow_class)
        .slick('slickPlay')
        .slick('slickSetOption', 'pauseOnDotsHover', true);
        
    });
    
    $('.next.' + slideshow_class).on('click', function() {
    $("." + slideshow_class)
        .slick('slickNext')

    });
    
    $('.prev.' + slideshow_class).on('click', function() {
    $("." + slideshow_class)
        .slick('slickPrev')

    });
    


    //return (htmlCode,slideshow_class)

}

function gime_fields_args(ze_args_plot,zeDesc,ze_id,suffixe) {
            
            if (zeDesc["format"] == "grib1") {
                ze_args_plot["field"][ze_id] = giveme_grib1_handgrip(suffixe) ;
                ze_args_plot["field_v"][ze_id] = giveme_grib1_handgrip(suffixe) ;
                try {
                    ze_args_plot["field_v"][ze_id]["indicatorOfParameter"] = parseInt(GetGribKey("#field_v" + suffixe,0));                    
                }
                catch (err) {}
            }
            else if (zeDesc["format"] == "grib2") {
                ze_args_plot["field"][ze_id] = giveme_grib2_handgrip(suffixe) ;
                ze_args_plot["field_v"][ze_id] = giveme_grib2_handgrip(suffixe) ;
                try {
                    ze_args_plot["field_v"][ze_id]["parameterNumber"] = parseInt(GetGribKey("#field_v" + suffixe,0));                    
                }
                catch (err) {}
            }
            else {
                ze_args_plot["field"][ze_id] = $("#field" + suffixe).val();
                ze_args_plot["field_v"][ze_id] = $("#field_v" + suffixe).val();
            }

}

function InitColorMapList(suffixe,current) {
    $("#ColorMapList" + suffixe).val(current);
    $("#ColorMapList"+ suffixe)
        .iconselectmenu()
        .iconselectmenu("menuWidget")
        .addClass("ui-menu-icons avatar overflow");
        }
        
        
function setup_render_labels(suffixe) {
	//Les défauts sont geres directment dans le code html, evite des incoherences lors du clonage
    //$("#decumul_label" + suffixe).text("Decumul is off ");
    //$("#update_label" + suffixe).text("AutoUpdate is on ");
    //$("#reverse_colormap_label" + suffixe).text("Reversed colormap is off");
    $("#autoupdate" + suffixe).click(function() {
        $("#update_label" + suffixe).text(this.checked ? "AutoUpdate is on " : "AutoUpdate is off");
    });
    $("#reverse_colormap" + suffixe).click(function() {
        $("#reverse_colormap_label" + suffixe).text(this.checked ? "Reversed colormap is on " : "Reversed colormap is off");
    });
    $("#decumul" + suffixe).click(function() {
        $("#decumul_label" + suffixe).text(this.checked ? "Decumul is on " : "Decumul is off");
    });
    }
    
function activate_tabs_plot(zeDesc) {
    var suffixe = zeDesc["suffixe"]
    $("#tabs_plot" + suffixe).tabs({
        beforeActivate: function(event, ui) {
            zeDesc["graphicmode"] = ui.newPanel.attr('id').replace('_cloned','');
            
            if (zeDesc["graphicmode"] == "points") {
                $("#move_col"+ suffixe).detach().appendTo($("#moving_color_content_2"+ suffixe));
            } else if (zeDesc["graphicmode"] == "contourlines") {
                $("#move_nb"+ suffixe).detach().appendTo($("#moving_nb_content_2"+ suffixe));
            } else if (zeDesc["graphicmode"] == "colorshades") {
                $("#move_col"+ suffixe).detach().appendTo($("#moving_color_content_1"+ suffixe));
                $("#move_nb"+ suffixe).detach().appendTo($("#moving_nb_content_1"+ suffixe));
            } else if (zeDesc["graphicmode"] == "vectors") {
                zeDesc["graphicmode"] = "colorshades"; //si l'onglet sélectionné est le vecteur, on force à colorshade pour FF                   
            }
        }
    });
    }

function smart_slider(suffixe) {

            $("#range1" + suffixe).change(function() {
                zemin = parseFloat($(this).val().split(' to ')[0]);
                zemax = parseFloat($(this).val().split(' to ')[1]);
                //alert(zemax);

                changeSlider(zemin, zemax, suffixe)
                $("#range1" + suffixe).val($("#slider-range" + suffixe).slider("values", 0) + " to " + $("#slider-range" + suffixe).slider("values", 1));
                desactivateAutoUpdate(suffixe);
            })
}


function cleantTooltips() {
    $('.ui-tooltip-content').remove();
        $('.ui-tooltip').remove();
};
    

function makeTabTitle(args_json,objId) {
    
    //Put ressources description on tab title for tooltip
    var backToJsonObj = JSON.parse(args_json);
    var jsonPretty = JSON.stringify(backToJsonObj, ["model","vconf","geometry","date","cutoff","term","expériment","block","suite"], 4);
    //alert(jsonPretty);
    $('#'+objId).prop('title', "Currently in memory : " + jsonPretty.replace(/\"/g,"") );
    //alert(JSON.stringify(args_vortex_get_json,null,2))
                
            }
            
function makeTabTitleFS(fileslist,objId) {
    
    $('#'+objId).prop('title', "Currently in memory : " + JSON.stringify(fileslist));
    //alert(JSON.stringify(args_vortex_get_json,null,2))
                
            }

function generateUUID() {
    var d = new Date().getTime();
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        var r = (d + Math.random()*16)%16 | 0;
        d = Math.floor(d/16);
        return (c=='x' ? r : (r&0x3|0x8)).toString(16);
    });
    return uuid;
};
