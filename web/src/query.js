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
	viewer.js
	Object responsible for maintaining molecular viewer.
	Query and Results both call have references to the viewer to set data.
	Query provides a callback for when pharmacophores are selected
*/

var Pharmit = Pharmit || {};

Pharmit.Query = (function() {
	
	var defaultFeature = {name:"Hydrophobic",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};

	
	function Query(element, viewer) {
		//private variables and functions
		var querydiv = $('<div>').addClass('pharmit_query');
		var vizgroup = null;
		var features = null;
		var featureheading = null;
		var receptorData = null;
		var receptorName = null; //filename (need ext)
		var ligandData = null;
		var ligandName = null;
		
		var doSearch = function() {
			
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
		
		var loadFeatures = function() {
			
		};
		
		var loadReceptor = function(data, fname) {
			receptorData = data;
			receptorName = fname;
			viewer.setReceptor(data, fname);
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
		
		var loadSession = function(data) {
			var query = $.parseJSON(data);
			
			features.detach();
			//replace features
			features.empty();
			if(query.points) {
				$.each(query.points, function(i, pt) {
					new Feature(viewer, features, pt);
				});
			}
			features.accordion("option","active",false);
			features.accordion("refresh");

			featureheading.after(features); 
			
			//get named settings, including visualization
			$.each(query, function(key,value) {
				var i = $('input[name='+key+']');
				if(i.length) {
					i.val(value);
				}
			});
			receptorData = query.receptor;
			receptorName = query.recname;			
			
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
			
			if(query.backgroundcolor) //this is the id of the correct radio button
				$('#'+query.backgroundcolor).prop('checked',true).change();
			viewer.setView(query.view);			

		};
		
		//return the query object
		var getQueryObj = function() {
			
			//get features
			var ret = {};
			ret.points = [];
			
			$.each(features.children(), function(key, fdiv) {
				ret.points.push(fdiv.feature.obj);
			});
			//everything with a name is something we want to save
			
			$.each($('[name]',querydiv), function(index, elem) {
				if(elem.name) ret[elem.name] = elem.value;
			});
			
			//radio buttons have to have the same name, so have to hack around background color
			ret.backgroundcolor = 'whiteBackground';
			if($('#blackBackground').prop('checked')) ret.backgroundcolor = 'blackBackground';
			
			//structures
			ret.ligand = ligandData;
			ret.ligandFormat = ligandName;
			ret.receptor = receptorData;
			ret.recname = receptorName;
			return ret;
		};
		
		var saveSession = function() {
			
			//IE doesn't support arbitrary data url's so much echo through a server
			//to download a file that is already on the client's machine
			// echo data back as a file to save
			var cmd = Pharmit.server+'?cmd=savedata&type="text%2Fphjson"&fname="pharmit.json"';
			var form = $('<form>', { 'action': cmd, 'method': 'post'});
			var qobj = getQueryObj();
			form.append($('<input>', {'name':"data",'type':"hidden",value:JSON.stringify(qobj,null,4)}));
			form.appendTo(document.body);
			form.submit();
			$(form).remove();
			

		};
		
		
		//create a split button from a list of vendors and prepend it to header
		var createSearchButton = function(header,vendors) {
			var buttons = $('<div>').addClass('searchdiv');
			var run = $('<button>Search '+vendors[0]+'</button>').appendTo(buttons).button();
			var select = $('<button>Select subset to search</button>').appendTo(buttons).button({text: false, icons: {primary: "ui-icon-triangle-1-s"}});
			
			buttons.buttonset();
			var ul = $('<ul>').appendTo($('body')).addClass('floatmenu'); //can't be in query div because of overflow truncation
			var lis = [];
			for(var i = 0, n = vendors.length; i < n; i++) {
				lis[i] = '<li>'+vendors[i]+'</li>';
			}
			ul.append(lis);
			ul.hide().menu().on('menuselect', function(event, ui) {
				run.button("option",'label',"Search "+ui.item.text());
			});
			
			//handlers
			run.click(doSearch);
			select.click(
					function() {
						var menu = ul.show().position({
							my: "left top",
							at: "left buttom",
							of: this
						});
						$(document).one('click', function() { menu.hide(); });
						return false;
					});
			
			header.prepend(buttons);
		};
		
		
		//public variables and functions
		
		
		//initialization code
		querydiv.resizable({handles: "e"});
		querydiv.disableSelection();
		
		var header = $('<div>').appendTo(querydiv).addClass("queryheader");
		createSearchButton(header,['MolPort','ZINC']);
		
		//load features and load receptor
		var loaders = $('<div>').appendTo(header);
		var loadfeatures = $('<button>Load Features...</button>').button();
		var loadrec = $('<button>Load Receptor...</button>').button();
		
		//fileinput needs the file inputs in the dom
		element.append(querydiv);
		var loadfeaturesfile = $('<input type="file">').appendTo(loaders).fileinput(loadfeatures).change(loadFeatures);		
		var loadrecfile = $('<input type="file">').appendTo(loaders).fileinput(loadrec).change(function(e) {readText(this, loadReceptor);});
		
		querydiv.detach();
		
		//query features
		var body = $('<div>').appendTo(querydiv).addClass("querybody");
		var featuregroup = $('<div>').appendTo(body);
		featureheading = $('<div>Pharmacophore</div>').appendTo(featuregroup).addClass('queryheading');
		features = $('<div>').appendTo(featuregroup);
		features.accordion({header: "> div > h3", animate: true, collapsible: true,heightStyle:'content'})
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
		
		var addbutton = $('<button>Add</button>').appendTo(featuregroup)
			.button({text: true, icons: {secondary: "ui-icon-circle-plus"}})
			.click(function() {new Feature(viewer, features, defaultFeature);}); //feature adds a reference to itself in its container
		var sortbutton = $('<button>Sort</button>').appendTo(featuregroup).button({text: true, icons: {secondary: "ui-icon ui-icon-carat-2-n-s"}}).click(sortFeatures);

		//filters
		var filtergroup = $('<div>').appendTo(body);
		$('<div>Filters</div>').appendTo(filtergroup).addClass('queryheading');
		var filters = $('<div>').appendTo(filtergroup);		
		
		var heading = $('<h3>Hit Reduction<br></h3>').appendTo(filters);
		var hitreductionsummary = $('<span class="headingsummary"></span>').appendTo(heading);

		var hitreductions = $('<div class="hitreduction"></div>').appendTo(filters);
		
		var row = $('<div>').addClass('filterrow').appendTo(hitreductions);
		row.append('<label title="Maximum number of orientations returned for each conformation" value="1" for="reduceorienttext">Max Hits per Conf:</label>');
		var maxorient = $('<input id="reduceorienttext" name="max-orient">').appendTo(row).spinner();
		
		row = $('<div>').addClass('filterrow').appendTo(hitreductions);
		row.append('<label title="Maximum number of conformations returned for each compound" value="1" for="reduceconfstext">Max Hits per Mol:</label>');
		var maxconfs = $('<input id="reduceconfstext" name="reduceConfs">').appendTo(row).spinner();
		
		row = $('<div>').addClass('filterrow').appendTo(hitreductions);
		row.append('<label title="Maximum number of hits returned" value="1" for="reducehitstext">Max Total Hits:</label>');
		var maxhits = $('<input id="reducehitstext" name="max-hits">').appendTo(row).spinner();
		
		
		heading = $('<h3>Hit Screening<br></h3>').appendTo(filters);
		var hitscreeningsummary = $('<span class="headingsummary"></span>').appendTo(heading);
		var hitscreening = $('<div class="hitscreening"></div>').appendTo(filters);
		row = $('<div>').addClass('filterrow').appendTo(hitscreening);
		var minweight = $('<input id="minmolweight" name="minMolWeight">').appendTo(row).spinner();
		row.append($('<label title="Minimum/maximum molecular weight (weights are approximate)" value="1" for="maxmolweight">&le;  Molecular Weight &le;</label>'));
		var maxweight = $('<input id="maxmolweight" name=maxMolWeight>').appendTo(row).spinner();

		row = $('<div>').addClass('filterrow').appendTo(hitscreening);
		var minnrot = $('<input id="minnrot" name="minrotbonds">').appendTo(row).spinner();
		row.append($('<label title="Minimum/maximum number of rotatable bonds" value="1" for="maxnrot"> &le;  Rotatable Bonds &le;</label>'));
		var maxnrot = $('<input id="maxnrot" name="maxrotbonds">').appendTo(row).spinner();

		filters.accordion({animate: true, collapsible: true,heightStyle:'content'});
		
		
		//viewer settings
		vizgroup = $('<div>').appendTo(body);
		$('<div>Visualization</div>').appendTo(vizgroup).addClass('queryheading');
		
		viewer.appendViewerControls(vizgroup);

		
		//load/save session
		var footer = $('<div>').appendTo(querydiv).addClass("queryfooter");
		element.append(querydiv);

		var loadsession = $('<button>Load Session...</button>').button();
		
		var loadsessionfile = $('<input type="file">').appendTo(footer).fileinput(loadsession).change(function() {readText(this,loadSession);});	
		var savesession = $('<button>Save Session...</button>').appendTo(footer).button().click(saveSession);		
				
	}

	return Query;
})();