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
	
	var defaultFeature = {name:"Hydrophobic",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null};

	/* object representing a single pharmacophore feature*/
	function Feature(features, fobj) {
		
		//setup html
		this.container = $('<div>').appendTo(features).addClass('featurediv');
		var heading = $('<h3>'+fobj.name+'<br></h3>').appendTo(this.container);
		this.summary = $('<span>').appendTo(heading).addClasS('featuresummary');
		var editdiv = $('<div>').appendTo(this.container);
		var locationdiv = $('<div>').appendTo(this.editdiv).addClass('locationdiv');
		
		$('<label>x:</label>').appendTo(locationdiv);
		this.x = $('<input>').appendTo(locationDiv).addClass('coordinput');
		this.x.spinner({step: 0.1, numberFormat: 'n'});
		
		$('<label>y:</label>').appendTo(locationdiv);
		this.y = $('<input>').appendTo(locationDiv).addClass('coordinput');
		this.y.spinner({step: 0.1, numberFormat: 'n'});

		$('<label>z:</label>').appendTo(locationdiv);
		this.z = $('<input>').appendTo(locationDiv).addClass('coordinput');
		this.z.spinner({step: 0.1, numberFormat: 'n'});

		$('<label>Radius:</label>').appendTo(locationdiv);
		this.radius = $('<input>').appendTo(locationDiv).addClass('radiusinput');
		this.radius.spinner({step: 0.1, numberFormat: 'n'});

		var orientdiv = this.orientdiv = $('<div>').appendTo(this.editdiv).addClass('orientdiv');
		
		$('<input type="checkbox"').appendTo(orientdiv);
		
		$('<label>&theta;:</label>').appendTo(orientdiv);
		this.theta = $('<input>').appendTo(orientdiv).addClass('orientinput');
		this.theta.spinner({step: 1, numberFormat: 'n'});
		
		$('<label>&phi;:</label>').appendTo(orientdiv);
		this.phi = $('<input>').appendTo(orientdiv).addClass('orientinput');
		this.phi.spinner({step: 1, numberFormat: 'n'});
		
		this.orientdiv.hide();
		
		this.setFeature(fobj);
	}

	//return feature values
	Feature.prototype.getFeature = function() {
		return this.obj;
	};
	
	//set the feature to fobj
	Feature.prototype.setFeature = function(fobj) {
		
	} ;
	//disable feature (remove from viewer, grey out)
	Feature.prototype.disable = function() {
		
	};
	
	//enable - show in viewer, style appropriately
	Feature.prototype.enable = function() {
		
	};
	
	//display in selected style
	Feature.prototype.select = function() {
		
	};
	
	//remove selection style
	Feature.prototype.deselect = function() {
		
	};
	
	//remove completely
	Feature.prototype.remove = function() {
		
	};
	
	function Query(element, viewer) {
		//private variables and functions
		var querydiv = $('<div>').addClass('pharmit_query');
		var features = null;

		
		var doSearch = function() {
			
		};
		
		var loadFeatures = function() {
			
		};
		
		var loadReceptor = function() {
			if(this.files.length > 0) {
				var file = this.files[0];
				var reader = new FileReader();
			    reader.onload = function(evt) {
			        viewer.setReceptor(evt.target.result,file.name);
			    };
			    reader.readAsText(file);
			}
		};
		
		var addFeature = function() {
			
		};
		
		var sortFeatures = function() {
			
		};
		
		var loadSession = function() {
			
		};
		
		var saveSession = function() {
			
		};
		
		
		//create a split button from a list of vendors and prepend it to header
		var createSearchButton = function(header,vendors) {
			var buttons = $('<div>').addClass('searchdiv');
			var run = $('<button>Search '+vendors[0]+'</button>').appendTo(buttons).button();
			var select = $('<button>Select subset to search</button>').appendTo(buttons).button({text: false, icons: {primary: "ui-icon-triangle-1-s"}});
			
			buttons.buttonset();
			var ul = $('<ul>').appendTo(buttons);
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
		
		//called when user picks a pharmacophore
		this.pickCallback = function(pharm) {
			
		};
		
		
		//initialization code
		querydiv.resizable({handles: "e"});
		var header = $('<div>').appendTo(querydiv).addClass("queryheader");
		createSearchButton(header,['MolPort','ZINC']);
		
		//load features and load receptor
		var loaders = $('<div>').appendTo(header);
		var loadfeatures = $('<button>Load Features...</button>').button();
		var loadrec = $('<button>Load Receptor...</button>').button();
		
		//fileinput needs the file inputs in the dom
		element.append(querydiv);
		var loadfeaturesfile = $('<input type="file">').appendTo(loaders).fileinput(loadfeatures).change(loadFeatures);		
		var loadrecfile = $('<input type="file">').appendTo(loaders).fileinput(loadrec).change(loadReceptor);
		
		querydiv.detach();
		
		//query features
		var body = $('<div>').appendTo(querydiv).addClass("querybody");
		var featuregroup = $('<div>').appendTo(body);
		$('<div>Pharmacophore</div>').appendTo(featuregroup).addClass('queryheading');
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
		
		var addbutton = $('<button>Add</button>').appendTo(featuregroup).button({text: true, icons: {secondary: "ui-icon-circle-plus"}}).click(addFeature);
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
		var vizgroup = $('<div>').appendTo(body);
		$('<div>Visualization</div>').appendTo(vizgroup).addClass('queryheading');
		
		viewer.appendViewerControls(vizgroup);

		
		//load/save session
		var footer = $('<div>').appendTo(querydiv).addClass("queryfooter");
		element.append(querydiv);

		var loadsession = $('<button>Load Session...</button>').button();
		
		var loadsessionfile = $('<input type="file">').appendTo(footer).fileinput(loadsession).change(loadSession);	
		var savesession = $('<button>Save Session...</button>').appendTo(footer).button().click(saveSession);		
				
	}

	return Query;
})();