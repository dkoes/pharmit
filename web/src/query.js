/*
 * Pharmit Web Client
 * Copyright 2015 David R Koes and University of Pittsburgh
 *  The JavaScript code in this page is free software: you can
    redistribute it and/or modify it under the terms of the GNU
    General Public License (GNU GPL) as published by the Free Software
    Foundation, either version 2 of the License, or (at your option)
    any later version.  The code is distributed WITHOUT ANY WARRANTY;
    without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU GPL for more details.
 */
/*
	query.js
	Left div that manages the query information.
 */

var Pharmit = Pharmit || {};

Pharmit.Query = (function() {

    var shapeMode = "filter";
    var defaultFeature = {name:"Hydrophobic",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
    var defaultInShapeFeature = {name:"InclusionSphere",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
    var defaultExShapeFeature = {name:"ExclusionSphere",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
    var pharmaGistRegEx = /@<TRIPOS>MOLECULE[^@]*?@<TRIPOS>ATOM\n(\s*\d+\s*(ACC|DON|CAT|ANI|HYD|AR).*)*\n@<TRIPOS>BOND\n/g;
    var privatedialog = null;
    var getQueryObj = null;
    var endsWith = function(str, suffix) {
        return str.indexOf(suffix, str.length - suffix.length) !== -1;
    };

    function Query(element, viewer, results) {
        //private variables and functions
        var querydiv = $('<div>').addClass('pharmit_query pharmit_overlay');
        var body = null;
        var features = null;
        var inshapefeatures = null;
        var exshapefeatures = null;

        var featuregroup = null;
        var receptorData = null;
        var receptorName = null; //filename (need ext)
        var receptorKey = null; //md5 key to avoid transfering full structure
        var ligandData = null;
        var ligandName = null;

        var doSearch = function() {
            var qobj = getQueryObj(true);

            if(shapeMode == 'search') {
                ga('send','event','query','shape',qobj.subset);
                results.shquery(qobj, receptorData);
            } else {
                //results manages queries
                ga('send','event','query','pharmacophore',qobj.subset);
                results.phquery(qobj, receptorData);
            }
        };

        //boiler plate for asynchronously extracting text from a file input
        var readText = function(input,func) {
            if(input.files.length > 0) {
                var file = input.files[0];
                var reader = new FileReader();
                reader.onload = function(evt) {
                    func(evt.target.result,file.name);
                };
                reader.readAsText(file);
                $(input).val('');
            }
        };


        var setFeaturesHelper = function(featurearray, features, type) {
            var container = features.parent();
            features.detach();
            //replace features
            var old = features.children();
            $.each(old, function(i, fdiv) {
                fdiv.feature.deleteFeature();
            });

            features.empty();
            if(featurearray) {
                $.each(featurearray, function(i, pt) {
                    new Feature(viewer, features, pt, type);
                });
            }
            features.accordion("option","active",false);
            features.accordion("refresh");

            container.prepend(features); 
        };

        //take an array of pharmacophore features (query.points) and
        //put them in the query view
        var setFeatures = function(featurearray) {

            if(!featurearray) return;
            var start = new Date().getTime();

            viewer.disableRendering();
            //while we're removing/adding bunches of features, don't bother rendering until the end

            var phfeatures = [];
            var inspheres = [];
            var exspheres = [];

            //split up the featurearray by feature type so shape features go in the right place

            for(var i = 0, n = featurearray.length; i < n; i++) {
                var f = featurearray[i];
                if(f.name == "InclusionSphere") {
                    inspheres.push(f);
                }
                else if(f.name == "ExclusionSphere") {
                    exspheres.push(f);
                }
                else {
                    phfeatures.push(f);
                }
            }

            setFeaturesHelper(phfeatures, features);
            setFeaturesHelper(inspheres, inshapefeatures, Feature.INCLUSIVESHAPE);
            setFeaturesHelper(exspheres, exshapefeatures, Feature.EXCLUSIVESHAPE);
            $('#exselect').change();

            viewer.enableRendering();
            var end = new Date().getTime();
            var time = end - start;
            //console.log('setFeatures time: ' + time);
        };

        //query server to get pharmacophore
        //result replaces any existing featuers
        var loadFeatures = function(data, lname) {

            ligandData = null;
            ligandName = null;
            var postData = {
                    cmd: 'getpharma',
                    ligand: data,
                    ligandname: lname,
            };

            if(receptorName) {
                if(receptorKey) { //most likely
                    postData.reckey = receptorKey;
                    postData.recname = receptorName;
                } else {
                    postData.receptor = receptorData;
                    postData.recname = receptorName;	
                }
            }

            $.post(Pharmit.server, postData, null, 'json').done(function(ret) {
                if(ret.status) { //success
                    setFeatures(ret.points);					
                    if(ret.mol) {
                        //this was molecular data, save it
                        ligandName = lname;
                        ligandData = data;

                        //pharmagist files embed the pharmacophore with the molecule, need to remove it
                        if(endsWith(lname,"mol2") && pharmaGistRegEx.test(data)) {
                            data = data.replace(pharmaGistRegEx, '');
                        }

                        viewer.setLigand(data, lname);						
                    }
                    else {
                        viewer.setView(); //orient on features in absence of ligand
                    }

                } else {
                    alert("Error: "+ret.msg);
                }
            }).fail(function(obj,err) {
                if(err == "parsererror") {
                    alert("Didn't understand the the server resonse.  If you can reproduce this error, please "+
                            "send your input files to "+Pharmit.email+".");
                }
                else {
                    alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
                }
            });

        };


        //set receptor variables, show receptor in viewer,
        //and register receptor with server
        var loadReceptor = function(data, fname) {

            if(!data || !fname)
                return;

            receptorData = data;
            receptorName = fname;
            viewer.setReceptor(data, fname);

            //calculate md5 of receptor
            receptorKey = null; //but don't set it until we get ack from server
            var rKey = Pharmit.hex_md5(receptorData);	

            $.post( Pharmit.server, { 
                cmd: "setreceptor",
                key: rKey,
                receptor: receptorData
            }).done(function() {
                receptorKey = rKey;
            }); //key setting isn't critical, so skip the the fail handler
        };

        //given ligand and (optional) receptor data, load the structures and compute interaction features
        //this creates a new query
        this.setLigandAndReceptor = function(ligand, ligandName, receptor, receptorName) {
            if(receptor) {
                loadReceptor(receptor, receptorName); 
            }
            loadFeatures(ligand, ligandName);			
        };

        //order features so enabled are on top and within the enabled/disabled
        //categories features are sorted by type
        var sortFeatures = function() {
            var fdivs = features.children().detach();

            fdivs.sort(function(a,b) {
                var x = a.feature.obj;
                var y = b.feature.obj;

                if(x.enabled != y.enabled) {
                    return y.enabled-x.enabled;
                }
                else if(x.name != y.name) {
                    return x.name.localeCompare(y.name);
                }
                return x.radius-y.radius;

            });

            features.append(fdivs);
        };


        var loadSession = this.loadSession = function(data) {

            var query = data; //support passing an object directly
            if(typeof(data) == "string") 
                query = $.parseJSON(data);

            setFeatures(query.points);

            loadReceptor(query.receptor, query.recname);		

            viewer.setReceptor(receptorData, receptorName);

            if(query.sdf) { //backwards compat with zincpharmer
                ligandData = decodeURIComponent(query.sdf);
                //try to guess format
                if(ligandData.match(/^@<TRIPOS>MOLECULE/)) {
                    ligandName = ".mol2";
                } else if(ligandData.match(/^HETATM/) || ligandData.match(/^ATOM/)) {
                    ligandName = ".pdb";
                } else if(ligandData.match(/^.*\n.*\n.\s*(\d+)\s+(\d+)/)){
                    ligandName = ".sdf"; //could look at line 3
                } else {
                    ligandName = ".xyz";
                }
            } else {
                ligandData = query.ligand;
                ligandName = query.ligandFormat;
            }
            viewer.setLigand(ligandData, ligandName);		

            //get named settings, including visualization
            $.each(query, function(key,value) {
                var i = $('input[name='+key+']');
                if(i.length) {
                    i.val(value).change();
                }
                else {
                    i = $('select[name='+key+']');
                    if(i.length) {
                        i.val(value).change();
                    }
                    else {
                        i = $('button[name='+key+']');
                        if(i.length) {
                            i.val(value).change();
                        }	
                    }
                }
            });

            //set radiobuttons
            $.each($('.ui-buttonset',querydiv), function(index, elem) {
                //for every radiobutton set
                var name = elem.id; //id of buttonset is property
                if(query[name]) {
                    var id = '#'+query[name];
                    $(id).prop('checked',true).change();
                }
            });

            viewer.setView(query.view);			

        };

        //return the query object
        //if trimreceptor is true, will delete receptor information if a key is set
        getQueryObj = function(trimreceptor) {

            //get features
            var ret = {};
            ret.points = [];

            $.each(features.children(), function(key, fdiv) {
                ret.points.push(fdiv.feature.obj);
            });
            $.each(inshapefeatures.children(), function(key, fdiv) {
                ret.points.push(fdiv.feature.obj);
            });
            $.each(exshapefeatures.children(), function(key, fdiv) {
                ret.points.push(fdiv.feature.obj);
            });
            //everything with a name is something we want to save

            $.each($('[name]',querydiv), function(index, elem) {
                var name = elem.name;
                if(name) {
                    var val = elem.value;
                    if($.isNumeric(elem.value)) {
                        val = Number(elem.value);
                    }
                    ret[name] = val;
                }
            });

            //radio buttons
            $.each($('.ui-buttonset',querydiv), function(index, elem) {
                //for every radiobutton set
                var name = elem.id; //id of buttonset is property
                if(name && name.length > 0) {
                    var val = $(':checked', elem).attr('id');
                    ret[name] = val;
                }
            });

            //structures
            ret.ligand = ligandData;
            ret.ligandFormat = ligandName;

            if(!receptorKey || !trimreceptor) ret.receptor = receptorData; //send full receptor
            ret.recname = receptorName;
            ret.receptorid = receptorKey;

            ret.view = viewer.getView();
            return ret;
        };

        var saveSession = function() {
            Pharmit.inFormSubmit = true;			

            //IE doesn't support arbitrary data url's so much echo through a server
            //to download a file that is already on the client's machine
            // echo data back as a file to save
            var cmd = Pharmit.server+'?cmd=savedata&type="text%2Fphjson"&fname="pharmit.json"';
            var form = $('<form>', { 'action': cmd, 'method': 'post'});
            var qobj = getQueryObj(false);
            form.append($('<input>', {'name':"data",'type':"hidden",value:JSON.stringify(qobj,null,4)}));
            form.appendTo(document.body);
            form.submit();
            $(form).remove();			

        };

        //escape special characters to avoid injections
        var escHTML = function(str) { return $('<div/>').text(str).html(); };	

        var updateShapeMesh = function(sel) {
            //send query to server to get mesh
            var qobj = getQueryObj(true);
            var kind = sel.meshstyle.kind;

            var postData = {cmd: 'getmesh',
                    type: kind,
                    json: JSON.stringify(qobj)
            };

            if(sel.mesh === null || sel.meshtolerance != sel.val()) {				
                $.post(Pharmit.server, postData, null, 'json').done(function(ret) {

                    if(ret.tolerance == sel.val() || sel.mesh === null) {
                        if(sel.mesh) viewer.removeMesh(sel.mesh);
                        sel.mesh = viewer.addMesh(ret, sel.meshstyle);
                        sel.meshtolerance = ret.tolerance;
                    }

                }
                );
            }
        };

        //setup all the shape query controls
        var createShapeQuery = function(body) {
            //shape
            var shapegroup = $('<div>').appendTo(body);
            $('<div>Shape</div>').appendTo(shapegroup).addClass('pharmit_heading');

            //inclusive shape
            var shapers = $('<div>').appendTo(shapegroup).addClass('pharmit_shapeconstraints');	

            var iheading = $('<h3 id="inshapehead">Inclusive Shape<br></h3>').appendTo(shapers);
            var inclusivediv = $('<div id="pharmit_inclusiveshape">').appendTo(shapers);
            var cell = null;

            var inselectdiv = $('<div>').appendTo(inclusivediv);
            var inselect = $('<select name="inselect" id="inselect">').addClass('pharmit_shapeselector').appendTo(inselectdiv);

            $('<option value="none">None</option>').appendTo(inselect);	
            $('<option value="ligand">Ligand</option>').appendTo(inselect);
            $('<option value="points">Spheres</option>').appendTo(inselect);

            inselect.val("none");
            inselect.selectmenu({
                width: '14em', 
                appendTo: inselectdiv, 
                change: function() { inselect.change();},
                position: {my: "left top", at: "left bottom", collision: "flip"}
            });			

            //none
            var nonediv = $('<div id="none-indiv">').appendTo(inclusivediv);
            $('<div>No inclusive shape will be used to filter pharmacophore matches.</div>').appendTo(nonediv).addClass("pharmit_shapedesc pharmit_shapefiltertext");
            $('<div>No inclusive shape will be used for shape screening.</div>').appendTo(nonediv).addClass("pharmit_shapedesc pharmit_shapesearchtext");


            //ligand
            var liganddiv = $('<div id="ligand-indiv">').appendTo(inclusivediv);
            $('<div>At least one heavy atom in a pharmacophore aligned pose must fall within the inclusive shape.</div>').appendTo(liganddiv).addClass("pharmit_shapedesc pharmit_shapefiltertext");
            $('<div>The entirety of the inclusive shape must be contained within the aligned hit.</div>').appendTo(liganddiv).addClass("pharmit_shapedesc pharmit_shapesearchtext");

            var inltable = $('<table>').appendTo(liganddiv);
            var inlrow = $('<tr>').appendTo(inltable);
            $('<td>').append('<label title="Depth in Angstroms to reduce surface of ligand inclusive shape by." value = "1" for="intolerance">Tolerance:</label>').appendTo(inlrow);
            cell = $('<td>').appendTo(inlrow).addClass('pharmit_shapecell');
            var intol = $('<input id="intolerance" name="intolerance">').appendTo(cell);
            intol.mesh = null;
            intol.meshstyle = {kind: "inclusive"};

            intol.change(function() { updateShapeMesh(intol); });
            intol.spinner({step: 0.5, stop: function() { intol.change();}});
            intol.val(1);

            var instylediv = $('<div id="inshapestyle">').appendTo(liganddiv).addClass('pharmit_meshstylediv');
            $('<input type="radio" id="inshapestyle-hide" name="inshapestyle"><label for="inshapestyle-hide">Hide</label>').appendTo(instylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    intol.meshstyle.hidden = true;
                    viewer.updateMesh(intol.mesh, intol.meshstyle);
                }
                instylediv.buttonset("refresh");
            });
            $('<input type="radio" id="inshapestyle-wire" name="inshapestyle"><label for="inshapestyle-wire">Wire</label>').appendTo(instylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    intol.meshstyle.hidden = false;
                    intol.meshstyle.wireframe = true;
                    viewer.updateMesh(intol.mesh, intol.meshstyle);
                }
                instylediv.buttonset("refresh");
            });
            $('<input type="radio" id="inshapestyle-solid" name="inshapestyle"><label for="inshapestyle-solid">Solid</label>').appendTo(instylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    intol.meshstyle.hidden = false;
                    intol.meshstyle.wireframe = false;
                    viewer.updateMesh(intol.mesh, intol.meshstyle);
                }
                instylediv.buttonset("refresh");
            }).prop("checked",true);
            instylediv.buttonset();


            //interaction points
            var inpoints = $('<div id="points-indiv">').appendTo(inclusivediv);
            $('<div>At least one heavy atom center of a pharmacophore aligned pose must fall within <b>each</b> sphere.</div>').appendTo(inpoints).addClass("pharmit_shapedesc pharmit_shapefiltertext");
            $('<div>The entirety of the interaction point shapes must be contained within the aligned hit.</div>').appendTo(inpoints).addClass("pharmit_shapedesc pharmit_shapesearchtext");

            inshapefeatures = $('<div>').appendTo(inpoints);
            inshapefeatures.accordion({header: "> div > h3", 
                animate: true, 
                active: false,
                collapsible: true,
                heightStyle:'content',
                beforeActivate: function( event, ui ) { 
                    var fdiv = null;

                    //deslect all features
                    var fdivs = inshapefeatures.children();
                    $.each(fdivs, function(key,fdiv) {
                        fdiv.feature.deselectFeature();
                    });
                    if(ui.newHeader.length > 0) { //being activated
                        fdiv = ui.newHeader.parent();
                        fdiv.get(0).feature.selectFeature();
                    }

                }})
                .sortable({ //from jquery ui example
                    axis: "y",
                    handle: "h3",
                    stop: function( event, ui ) {
                        // IE doesn't register the blur when sorting
                        // so trigger focusout handlers to remove .ui-state-focus
                        ui.item.children( "h3" ).triggerHandler( "focusout" );
                        // Refresh accordion to handle new order
                        $( this ).accordion( "refresh" );
                    }
                });			

            var adddiv = $('<div class="pharmit_pointsbuttons">').appendTo(inpoints);
            var ptaddbutton = $('<button>Add</button>').appendTo(adddiv)
            .button({text: true, icons: {secondary: "ui-icon-circle-plus"}})
            .click(function() {new Feature(viewer, inshapefeatures, defaultInShapeFeature, Feature.INCLUSIVESHAPE);}); //feature adds a reference to itself in its container

            //handler for choosing inclusive mode
            inselect.change(function() {
                nonediv.hide();
                liganddiv.hide();
                inpoints.hide();
                if(this.value == 'ligand'){ 
                    $.each(inshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.hideFeature();
                    });
                    liganddiv.show();
                    $('#inshapehead').addClass("pharmit_oblique");
                    intol.change();
                }
                else if(this.value == 'points') {
                    $.each(inshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.unhideFeature();
                    });
                    $('#inshapehead').addClass("pharmit_oblique");
                    inpoints.show();
                    if(intol.mesh) viewer.removeMesh(intol.mesh);
                    intol.mesh = null;
                }
                else {
                    $.each(inshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.hideFeature();
                    });
                    $('#inshapehead').removeClass("pharmit_oblique");
                    nonediv.show();
                    if(intol.mesh) viewer.removeMesh(intol.mesh);
                    intol.mesh = null;
                }
                inselect.selectmenu("refresh");	       
            });
            inselect.change();

            //exclusive shape
            var eheading = $('<h3 id="exshapehead">Exclusive Shape<br></h3>').appendTo(shapers);
            var exclusivediv = $('<div id="pharmit_exclusiveshape">').appendTo(shapers);

            var exselectdiv = $('<div>').appendTo(exclusivediv);
            var exselect = $('<select name="exselect" id="exselect">').addClass('pharmit_shapeselector').appendTo(exselectdiv);

            $('<option value="none">None</option>').appendTo(exselect);	
            $('<option value="receptor">Receptor</option>').appendTo(exselect);
            $('<option value="points">Spheres</option>').appendTo(exselect);

            exselect.val("none");
            exselect.selectmenu({
                width: '11em', 
                appendTo: exselectdiv, 
                change: function() { exselect.change(); },
                position: {my: "left top", at: "left bottom", collision: "flip"}
            });

            var exnonediv = $('<div id="none-exdiv">').appendTo(exclusivediv);
            $('<div>No exclusive shape constraint will be applied.</div>').appendTo(exnonediv).addClass("pharmit_shapedesc");

            var recdiv = $('<div id="receptor-exdiv">').appendTo(exclusivediv);
            var phexdesc = $('<div>').appendTo(recdiv).addClass("pharmit_shapedesc pharmit_shapefiltertext");
            phexdesc.text("Hits that have heavy atom centers within the exclusive shape in their pharmacophore aligned pose will be filtered out.");

            var shexdesc = $('<div>').appendTo(recdiv).addClass("pharmit_shapedesc pharmit_shapesearchtext");
            shexdesc.text("Aligned hits may not overlap any portion of the exclusive shape.");


            var extable = $('<table>').appendTo(recdiv);
            var exrow = $('<tr>').appendTo(extable);
            $('<td>').append('<label title="Depth in Angstroms to reduce surface of receptor exclusive shape by." value = "1" for="extolerance">Tolerance:</label>').appendTo(exrow);
            cell = $('<td>').appendTo(exrow).addClass('pharmit_shapecell');
            var extol = $('<input id="extolerance" name="extolerance">').appendTo(cell).val(1.0);
            extol.mesh = null;
            extol.meshstyle = {kind: "exclusive"};

            extol.change(function() { updateShapeMesh(extol); });
            extol.spinner({step: 0.5, stop: function() { extol.change();}});

            var exstylediv = $('<div id="exshapestyle">').appendTo(recdiv).addClass('pharmit_meshstylediv');
            $('<input type="radio" id="exshapestyle-hide" name="exshapestyle"><label for="exshapestyle-hide">Hide</label>').appendTo(exstylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    extol.meshstyle.hidden = true;
                    viewer.updateMesh(extol.mesh, extol.meshstyle);
                }
                exstylediv.buttonset("refresh");
            });
            $('<input type="radio" id="exshapestyle-wire" name="exshapestyle"><label for="exshapestyle-wire">Wire</label>').appendTo(exstylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    extol.meshstyle.hidden = false;
                    extol.meshstyle.wireframe = true;
                    viewer.updateMesh(extol.mesh, extol.meshstyle);
                }
                exstylediv.buttonset("refresh");
            });
            $('<input type="radio" id="exshapestyle-solid" name="exshapestyle"><label for="exshapestyle-solid">Solid</label>').appendTo(exstylediv)
            .change(function() {
                if($(this).prop("checked")) {
                    extol.meshstyle.hidden = false;
                    extol.meshstyle.wireframe = false;
                    viewer.updateMesh(extol.mesh, extol.meshstyle);
                }
                exstylediv.buttonset("refresh");
            }).prop("checked",true);
            exstylediv.buttonset();



            var expoints = $('<div id="points-exdiv">').appendTo(exclusivediv);
            $('<div>Heavy atom centers of pharmacophore aligned poses may not fall within any exclusion sphere.</div>').appendTo(expoints).addClass("pharmit_shapedesc pharmit_shapefiltertext");
            $('<div>Aligned hits may not overlap any exclusion sphere.</div>').appendTo(expoints).addClass("pharmit_shapedesc pharmit_shapesearchtext");

            exshapefeatures = $('<div>').appendTo(expoints);
            exshapefeatures.accordion({header: "> div > h3", 
                animate: true, 
                active: false,
                collapsible: true,
                heightStyle:'content',
                beforeActivate: function( event, ui ) { 
                    var fdiv = null;

                    //deslect all features
                    var fdivs = exshapefeatures.children();
                    $.each(fdivs, function(key,fdiv) {
                        fdiv.feature.deselectFeature();
                    });
                    if(ui.newHeader.length > 0) { //being activated
                        fdiv = ui.newHeader.parent();
                        fdiv.get(0).feature.selectFeature();
                    }

                }})
                .sortable({ //from jquery ui example
                    axis: "y",
                    handle: "h3",
                    stop: function( event, ui ) {
                        // IE doesn't register the blur when sorting
                        // so trigger focusout handlers to remove .ui-state-focus
                        ui.item.children( "h3" ).triggerHandler( "focusout" );
                        // Refresh accordion to handle new order
                        $( this ).accordion( "refresh" );
                    }
                });			

            adddiv = $('<div class="pharmit_pointsbuttons">').appendTo(expoints);
            ptaddbutton = $('<button>Add</button>').appendTo(adddiv)
            .button({text: true, icons: {secondary: "ui-icon-circle-plus"}})
            .click(function() {new Feature(viewer, exshapefeatures, defaultExShapeFeature, Feature.EXCLUSIVESHAPE);}); //feature adds a reference to itself in its container


            //handler for choosing exclusive mode
            exselect.change(function() {
                if(this.value == "receptor") {
                    $.each(exshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.hideFeature();
                    });
                    exnonediv.hide();
                    expoints.hide();
                    recdiv.show();
                    extol.change(); //mesh
                    $('#exshapehead').addClass("pharmit_oblique");
                } else if(this.value == "points") {
                    $.each(exshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.unhideFeature();
                    });
                    exnonediv.hide();
                    recdiv.hide();
                    expoints.show();
                    $('#exshapehead').addClass("pharmit_oblique");
                    if(extol.mesh) viewer.removeMesh(extol.mesh);
                    extol.mesh = null;
                } 
                else { //none
                    $.each(exshapefeatures.children(), function(key, fdiv) {
                        fdiv.feature.hideFeature();
                    });

                    $('#exshapehead').removeClass("pharmit_oblique");
                    exnonediv.show();
                    expoints.hide();
                    recdiv.hide();
                    if(extol.mesh) viewer.removeMesh(extol.mesh);
                    extol.mesh = null;
                }
                exselect.selectmenu("refresh");	       
            });
            exselect.change();


            shapers.accordion({animate: true, active: false, collapsible: true, heightStyle:'content'});

        };

        //create controls for specifying filters of query
        var createFilterQuery = function(body) {
            //filters
            var filtergroup = $('<div>').appendTo(body);
            $('<div>Filters</div>').appendTo(filtergroup).addClass('pharmit_heading');
            var filters = $('<div>').appendTo(filtergroup);		

            var heading = $('<h3>Hit Reduction<br></h3>').appendTo(filters);

            var hitreductions = $('<div>').addClass("pharmit_hitreduction").appendTo(filters);
            var reducetable = $('<table>').appendTo(hitreductions);

            var setReductionStyle = function() { //change style of headings of filters are specified
                if($('#reduceorienttext').val() !== '' ||
                        $('#reduceconfstext').val() !== '' ||
                        $('#reducehitstext').val() !== '') {
                    heading.addClass('pharmit_filtermodified');
                } else {
                    heading.removeClass('pharmit_filtermodified');
                }
            };


            var row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
            $('<td>').append('<label title="Maximum number of orientations returned for each conformation" value="1" for="reduceorienttext">Max Hits per Conf:</label>').appendTo(row);
            var cell = $('<td>').appendTo(row);
            $('<input id="reduceorienttext" name="max-orient">').appendTo(cell).spinner({min: 0, stop: setReductionStyle}).change(setReductionStyle);

            row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
            $('<td>').append('<label title="Maximum number of conformations returned for each compound" value="1" for="reduceconfstext">Max Hits per Mol:</label>').appendTo(row);
            cell = $('<td>').appendTo(row);
            $('<input id="reduceconfstext" name="reduceConfs">').appendTo(cell).spinner({min: 0, stop: setReductionStyle}).change(setReductionStyle);

            row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
            $('<td>').append('<label title="Maximum number of hits returned" value="1" for="reducehitstext">Max Total Hits:</label>').appendTo(row);
            cell = $('<td>').appendTo(row);
            $('<input id="reducehitstext" name="max-hits">').appendTo(cell).spinner({min: 0, stop: setReductionStyle}).change(setReductionStyle);


            var screenheading = $('<h3>Hit Screening<br></h3>').appendTo(filters);
            var hitscreening = $('<div>').appendTo(filters).addClass('pharmit_hitscreening');
            var screentable = $('<table>').appendTo(hitscreening);

            var setScreensStyle = function() { //change style of headings of filters are specified
                var nothingset = true;
                $('[name]',hitscreening).each( function(index, element) {
                    if(element.value !== '')
                        nothingset = false;
                });
                if(!nothingset)  {
                    screenheading.addClass('pharmit_filtermodified');
                } else {
                    screenheading.removeClass('pharmit_filtermodified');
                }
            };

            var addScreeningRow = function(name, label, title, minval) {
                //boilerplate for screening row
                var row = $('<tr>').addClass('pharmit_filterrow').appendTo(screentable);
                var cell = $('<td>').appendTo(row);
                $('<input name="min'+name+'">').appendTo(cell).spinner({min: minval, stop: setScreensStyle}).change(setScreensStyle);
                $('<td>').appendTo(row).append($('<label title="'+title+'" value="1" >&le;  '+label+' &le;</label>'));
                cell = $('<td>').appendTo(row);
                $('<input name="max'+name+'">').appendTo(cell).spinner({min: minval, stop: setScreensStyle}).change(setScreensStyle);

            };

            addScreeningRow("MolWeight", "MolWeight", "Minimum/maximum molecular weight (weights are approximate)", 0);
            addScreeningRow("rotbonds", "RotBonds", "Minimum/maximum number of rotatable bonds", 0);
            addScreeningRow("logp", "LogP", "Minimum/maximum partition coefficient as calculated by OpenBabel");
            addScreeningRow("psa", "PSA", "Minimum/maximum polar surface area as calculated by OpenBabel");
            addScreeningRow("aromatics", "Aromatics", "Minimum/maximum number of smallest aromatic rings", 0);
            addScreeningRow("hba", "HBA", "Minimum/maximum number of hydrogen bond acceptors",0);
            addScreeningRow("hbd", "HBD", "Minimum/maximum number of hydrogen bond donors",0);

            filters.accordion({animate: true, active: false, collapsible: true, heightStyle:'content'});
        };

        //setup select menu for choosing search mode
        var prependModeSelect = function(header) {
            var shapemodeid = 'ShapeModeSelect';
            var shapemodediv = $('<div>').prependTo(header).addClass('pharmit_shapemodediv');					
            var shapeselect = $('<button name="'+shapemodeid+'" id="'+shapemodeid+' value="filter">Pharmacophore Search -&gt; Shape Filter</button>').addClass('pharmit_styleselector').appendTo(shapemodediv);
            shapeMode = 'filter';
            shapeselect.button();

            //setup button according to itsvalue, do not change value
            shapeselect.change(function() {
                shapeMode = this.value;
                if(shapeMode == 'search') {
                    shapeselect.button("option","label","Shape Search -&gt; Pharmacophore Filter");
                    $('.pharmit_shapefiltertext').hide();
                    $('.pharmit_shapesearchtext').show();
                } else { //filter
                    shapeMode = 'filter';
                    shapeselect.button("option","label","Pharmacophore Search -&gt; Shape Filter");
                    $('.pharmit_shapefiltertext').show();
                    $('.pharmit_shapesearchtext').hide();
                }
                shapeselect.button("refresh");	 
            });
            shapeselect.val("filter");
            shapeselect.change();

            //invert value on click
            shapeselect.click(function() {

                if(this.value == 'filter') {
                    this.value = 'search';
                } else { //filter
                    this.value = 'filter';
                }
                shapeselect.change();      

            });

        };


        //public variables and functions

        var closer = $('<div>').appendTo(querydiv).addClass('pharmit_leftclose');
        var closericon = $('<span>').addClass("ui-icon ui-icon-carat-1-w").appendTo(closer);

        //initialization code
        querydiv.resizable({handles: "e",
            resize: function(event, ui) {
                viewer.setLeft(ui.size.width);
            }
        });
        querydiv.disableSelection();


        closer.click(function() {
            if(closer.hasClass('pharmit_leftisclosed')) {
                closer.removeClass('pharmit_leftisclosed');
                closericon.removeClass('ui-icon-carat-1-e');
                closericon.addClass('ui-icon-carat-1-w');
                var start = querydiv.width();
                querydiv.css('width', ''); //restore stylesheet width	
                var target = querydiv.width();
                querydiv.width(start);

                querydiv.animate({width: target},{
                    progress: function() { viewer.setLeft(querydiv.width());}
                }); 
                querydiv.resizable( "option", "disabled", false);

            } else { //close it 
                querydiv.animate({width: 0}, {
                    progress: function() { viewer.setLeft(querydiv.width());}
                }); 
                //viewer.setLeft(0);
                closer.addClass('pharmit_leftisclosed');
                closericon.addClass('ui-icon-carat-1-e');
                closericon.removeClass('ui-icon-carat-1-w');			
                querydiv.resizable( "option", "disabled", true );
            }
        });

        var header = $('<div>').appendTo(querydiv).addClass("pharmit_queryheader");

        //load features and load receptor
        var loaders = $('<div>').appendTo(header).addClass('pharmit_loaderdiv').addClass('pharmit_nowrap');
        var loadrec = $('<button>Load Receptor...</button>').button();
        var titletext = "Pharmacophore features can be provided as pharmacophore query formats (MOE, LigBuilder, LigandScout, PharmaGist, Pharmer) or automatically extracted from ligand structures (sdf, pdf, mol2, xyz). If a receptor is loaded, only interacting features will be enabled.";
        var loadfeatures = $('<button title="'+titletext+'">Load Features...</button>').button();

        //fileinput needs the file inputs in the dom
        element.append(querydiv);
        var loadrecfile = $('<input type="file">').appendTo(loaders).fileinput(loadrec).change(function(e) {readText(this, loadReceptor);});
        var loadfeaturesfile = $('<input type="file">').appendTo(loaders).fileinput(loadfeatures).change(function(e) {readText(this,loadFeatures);});		

        querydiv.detach();

        //query features
        body = $('<div>').appendTo(querydiv).addClass("pharmit_querybody");
        var featureheading = $('<div>Pharmacophore</div>').appendTo(body).addClass('pharmit_heading');
        featuregroup = $('<div>').appendTo(body);
        features = $('<div>').appendTo(featuregroup);
        features.accordion({header: "> div > h3", 
            animate: true, 
            active: false,
            collapsible: true,
            heightStyle:'content',
            beforeActivate: function( event, ui ) { 
                var fdiv = null;

                //deslect all features
                var fdivs = features.children();
                $.each(fdivs, function(key,fdiv) {
                    fdiv.feature.deselectFeature();
                });
                if(ui.newHeader.length > 0) { //being activated
                    fdiv = ui.newHeader.parent();
                    fdiv.get(0).feature.selectFeature();
                }

            }})
            .sortable({ //from jquery ui example
                axis: "y",
                handle: "h3",
                stop: function( event, ui ) {
                    // IE doesn't register the blur when sorting
                    // so trigger focusout handlers to remove .ui-state-focus
                    ui.item.children( "h3" ).triggerHandler( "focusout" );
                    // Refresh accordion to handle new order
                    $( this ).accordion( "refresh" );
                }
            });

        var buttondiv = $('<div>').appendTo(featuregroup).addClass('pharmit_featurebuttons');
        var addbutton = $('<button>Add</button>').appendTo(buttondiv)
        .button({text: true, icons: {secondary: "ui-icon-circle-plus"}})
        .click(function() {new Feature(viewer, features, defaultFeature);}); //feature adds a reference to itself in its container
        var sortbutton = $('<button>Sort</button>').appendTo(buttondiv).button({text: true, icons: {secondary: "ui-icon ui-icon-carat-2-n-s"}}).click(sortFeatures);

        createShapeQuery(body);

        createFilterQuery(body);

        prependModeSelect(header); //setup mode select after shape/pharmacophore elements are created

        //setup search button
        var buttons = $('<div>').addClass('pharmit_searchdiv');
        var run = $('<button id="pharmitsearchbutton" name="subset">Search</button>').appendTo(buttons).button();
        var select = $('<button>Select subset to search</button>').appendTo(buttons).button({text: false, icons: {primary: "ui-icon-triangle-1-s"}});
        select.tooltip().tooltip("disable"); //conflicts with menu
        run.val("");
        buttons.buttonset();

        //create a split button from a list of vendors and prepend it to header
        var createSearchButton = function(header,dbinfo) {
            var subsetinfo = dbinfo.standard;

            //default to what is already set (e.g., by loaded query) or molport
            var defaultdir = run.val();
            if(defaultdir == "") defaultdir = "molport";
            var founddefault = false;
            for(i = 0, n = subsetinfo.length; i < n; i++) {
                if(defaultdir == subsetinfo[i].subdir) {
                    founddefault = true;
                    break;
                }
            }
            //fallback to first dir if default not found
            if(!founddefault) {
                defaultdir = subsetinfo[0].subdir;
            }

            var ul = $('<ul>').appendTo($('body')).addClass('pharmit_floatmenu'); //can't be in query div because of overflow truncation
            var lis = [];
            var i, info, display;
            for(i = 0, n = subsetinfo.length; i < n; i++) {
                info = subsetinfo[i];
                display = escHTML(info.name);
                if(info.html) display = info.html; //optionally can provide html, but not from users
                lis[i] = '<li value='+i+' class="pharmit_subsetmenu">'+display+'<br>';
                lis[i] += '<span class="pharmit_subsetcnts">' + numeral(info.numConfs).format('0,0') + ' conformers of ' + numeral(info.numMols).format('0,0') + ' molecules</span>';
                lis[i] += '<span class="pharmit_subsettime">Updated: '+info.updated+'</span>';
                lis[i] += '</li>';

                //default to molport for now
                if(info.subdir == defaultdir) {
                    run.button("option",'label',"Search "+info.name);
                    run.val(info.subdir);
                }
            }

            ul.append(lis);
            ul.append($('<li> </li>'));
            var publicli = $('<li class="pharmit_contributed">Contributed Libraries</li>');
            var publicul =  $('<ul>').appendTo(publicli).addClass('pharmit_contributed_menu');
            var publiclis = [];
            var tosearch = [];
            var publicinfo = dbinfo.public;
            for(i = 0, n = publicinfo.length; i < n; i++) {
                info = publicinfo[i];
                if(info.numConfs === 0) continue;
                display = escHTML(info.name);
                var titlestr = "";
                tosearch[i] = display + info.updated; //can search by date
                if(info.description) {
                    titlestr = " title='"+escHTML(info.description)+"' ";
                    tosearch[i] += " " + escHTML(info.description);
                }
                tosearch[i] = tosearch[i].toLowerCase();
                publiclis[i] = '<li value='+subsetinfo.length+titlestr+' class="pharmit_subsetmenu">'+display+'<br>';
                subsetinfo.push(info);
                publiclis[i] += '<span class="pharmit_contrib_subsettime">Created: '+info.updated+'</span>';
                publiclis[i] += '<span class="pharmit_contrib_subsetcnts">' + numeral(info.numConfs).format('0,0') + ' conformers of ' + numeral(info.numMols).format('0,0') + ' molecules</span>';
                publiclis[i] += '</li>';
                publiclis[i] = $(publiclis[i]);
            }
            publicul.append(publiclis);
            publicdialog = $('<div class="pharmit_public_dialog" title="Contributed Libraries"> </div>');
            publicsearch = $('<input type="text" id="input" class="pharmit_contrib_search">');
            $('<div class="pharmit_search_div"><b>Search: </b></div>').append(publicsearch).appendTo(publicdialog);
            publicsearch.keyup(function (){
                var a, txtValue, filter;
                filter = publicsearch.val().toLowerCase();
                for(i in publiclis) {
                    if(tosearch[i].indexOf(filter) > -1) {
                        publiclis[i].css('display', "");
                    } else {
                        publiclis[i].css('display', "none");
                    }
                }
            });
            publicdialog.append(publicul);
            publicdialog.appendTo(body);
            publicdialog.dialog({
                autoOpen: false,
                height: 525,
                width: 600,
                modal: true,                
            });         
            $('<li class="pharmit_private">Contributed Libraries...</li>').appendTo(ul).click(
                    function() {
                        publicdialog.dialog("open");
                    });


            ul.append($('<li> </li>'));
            $('<li class="pharmit_private">Access Private Library...</li>').appendTo(ul).click(
                    function() {
                        privatedialog.dialog("open");
                    });
            publicul.menu({position:{
                my: "left top",
                at: "right top",
                collision: 'fit'
            }}).on('menuselect', function(event, ui) {
                var info = subsetinfo[ui.item.val()];
                run.button("option",'label',"Search "+info.name);
                run.val(info.subdir);
                publicdialog.dialog("close");
            });

            ul.hide().menu({position:{
                my: "left top",
                at: "right top",
                collision: 'fit'
            }}).on('menuselect', function(event, ui) {
                var info = subsetinfo[ui.item.val()];
                run.button("option",'label',"Search "+info.name);
                run.val(info.subdir);
            });

            //handlers
            run.click(doSearch);

            run.change(function() {
                //update button name
                var subdir = run.val();
                var i, n;
                for(i = 0, n = subsetinfo.length; i < n; i++) {
                    info = subsetinfo[i];
                    if(info.subdir == subdir) {
                        run.button("option",'label',"Search "+info.name);
                        break;
                    }
                }
                if(i == n) {
                    run.button("option",'label',"Unknown Subset");
                }				
            });
            select.click(
                    function() {
                        var menu = ul.show().position({
                            my: "left top",
                            at: "left bottom",
                            of: this,
                            collision: 'flipfit'
                        });
                        $(document).one('click', function() { menu.hide(); });
                        return false;
                    });

        };


        $.post(Pharmit.server, {cmd: "getsubsets"}, null, 'json').done(function(ret) {
            createSearchButton(header, ret);
        }).fail(function() {
            createSearchButton(header, {standard: [{name: "Error"}], public:[]});
            alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
        });

        header.prepend(buttons);


        //viewer settings
        var vizgroup = $('<div>').appendTo(body);
        $('<div>Visualization</div>').appendTo(vizgroup).addClass('pharmit_heading');
        var vizbody = $('<div>').appendTo(vizgroup).addClass('pharmit_vizdiv');
        viewer.appendViewerControls(vizbody);


        //load/save session
        var footer = $('<div>').appendTo(querydiv).addClass("pharmit_queryfooter");
        var bottomloaders = $('<div>').appendTo(footer).addClass("pharmit_bottomloaders").addClass('pharmit_nowrap');
        element.append(querydiv);

        var loadsession = $('<button>Load Session...</button>').button();

        var loadsessionfile = $('<input type="file">').appendTo(bottomloaders).fileinput(loadsession).change(function() {readText(this,loadSession);});	
        var savesession = $('<button>Save Session...</button>').appendTo(bottomloaders).button().click(saveSession);		

        viewer.setLeft(querydiv.width());

        //setup private dialog
        privatedialog = $('<div class="pharmit_private_dialog" title="Access Private Library">Enter your access code: </div>').appendTo(body);
        $('<input type="text" id="privateid" class="pharmit_privatetext">').appendTo(privatedialog);

        privatedialog.dialog({
            autoOpen: false,
            height: 125,
            width: 400,
            modal: true,
            buttons: {
                "Submit": function() {
                    //get info from server
                    $.post(Pharmit.server, {cmd: "getsubsets", subset: "Private/"+$('#privateid').val()}, null, 'json').done(function(ret) {
                        if(ret.error) {
                            alert(ret.error);
                        }
                        else {
                            var run = $('#pharmitsearchbutton');
                            run.button("option",'label',"Search "+ret.name);
                            run.val(ret.subdir);
                        }
                    }).fail(function() {
                        alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
                    });
                    privatedialog.dialog( "close" );
                },
                Cancel: function() {
                    privatedialog.dialog( "close" );
                }
            },
        });

    }

    return Query;
})();
