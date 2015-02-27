
/* object representing a single pharmacophore feature*/
function Feature(viewer, features, fobj) {
	
	//setup html
	var F = this;
	this.obj = {}; //set in setFeature
	
	this.shape = null; //for viewer
	this.selected = false;
	//updateViewer will update the viewer's visualization of this feature,
	//mask the prototype so we don't make unnecessary calls
	//while constructing the object - after setFeature we will create the shape
	//and make updateViewer functional
	this.updateViewer = function() {};
	
	this.container = $('<div>').appendTo(features).addClass('featurediv');
	this.container.get(0).feature = this;
	this.container.disableSelection();
	
	//header has checkbox for enable, name, and a close icon
	var heading = $('<h3></h3>').appendTo(this.container);
	this.enabled = $('<div>').addClass('toggle toggle-light togglediv').appendTo(heading).toggles().on('toggle',
			function(e, active) {
				if(active) {
					F.obj.enabled = true;
					F.container.addClass("enabledfeature");
				}
				else {
					F.obj.enabled = false;	
					F.container.removeClass("enabledfeature");
				}
				F.updateViewer();
			}
			).click(function(e) {return false;}); //needed to stop clickthrue
	var namediv = $('<div>').addClass('featurenamediv').appendTo(heading);
	var namespan = $('<span>').addClass('featurenameheading').appendTo(namediv);
	var closediv = $('<div>').addClass('featureclose').appendTo(heading).click(function() {
		//remove from viewer
		if(F.shape) viewer.removeFeature(F.shape);
		//remove from dom
		F.container.feature = null;
		F.container.remove();
	});
	var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
	
	//summary has (x,y,z) Radius r
	var summary = $('<span>').appendTo(namediv).addClass('featuresummary');		
	summary.append($(document.createTextNode('(')));
	this.xsummary = $('<span>').appendTo(summary);
	summary.append($(document.createTextNode(',')));
	this.ysummary = $('<span>').appendTo(summary);
	summary.append($(document.createTextNode(',')));
	this.zsummary = $('<span>').appendTo(summary);
	summary.append($(document.createTextNode(') Radius ')));
	this.rsummary = $('<span>').appendTo(summary);

	var editdiv = $('<div>').appendTo(this.container);
	
	//feature kind selection (name)
	var select = this.select = $('<select name="featurename">').addClass('featureselect').appendTo(editdiv);
	$.each(this.featureNames, function(key,val) {
		$('<option value="'+val+'">'+val+'</option>').appendTo(select);
	});
	select.change(function() {
		var n = this.value;
		F.obj.name = n;
		namespan.text(n);
		if(n == "Aromatic" || n == "HydrogenDonor" || n == "HydrogenAcceptor") {
			F.orientdiv.show();
			F.sizediv.hide();
			F.obj.hasvec = true;  //indicate vector is actually relevant to feature
		}
		else if(n == "Hydrophobic") {
			F.orientdiv.hide();
			F.sizediv.show();
			F.obj.hasvec = false; //but in case we change kind of feature leave vector_on at current setting
		}
		else { //charges
			F.orientdiv.hide();
			F.sizediv.hide();
			F.obj.hasvec = false;
		}
		F.updateViewer();
		});
	select.selectmenu({width: "12em", change: function() {select.trigger('change');}});
	
	//position (x,y,z)
	var locationdiv = $('<div>').appendTo(editdiv).addClass('locationdiv');
	
	//because spinners are insane and don't trigger a change event on the
	//underlying element and split up spinning and stopping and changing...
	var spinObject = function(elem, init) {
		var change = function() {elem.change();};
		init.spin = init.stop = init.change = change;
		return init;
	};
	
	$('<label>x:</label>').appendTo(locationdiv);
	this.x = $('<input>').appendTo(locationdiv).addClass('coordinput').change(function() {
		//change x value
		var x = parseFloat(this.value);
		F.xsummary.text(numeral(x).format('0.[00]'));
		F.obj.x = x;
		F.updateViewer();
		});
	this.x.spinner(spinObject(F.x,{step: 0.1, numberFormat: 'n'}));
	
	
	$('<label>y:</label>').appendTo(locationdiv);
	this.y = $('<input>').appendTo(locationdiv).addClass('coordinput');
	this.y.spinner(spinObject(F.y,{step: 0.1, numberFormat: 'n'})).change(function() {
		//change y value
		var y = parseFloat(this.value);
		F.ysummary.text(numeral(y).format('0.[00]'));
		F.obj.y = y;
		F.updateViewer();
		});

	$('<label>z:</label>').appendTo(locationdiv);
	this.z = $('<input>').appendTo(locationdiv).addClass('coordinput');
	this.z.spinner(spinObject(F.z,{step: 0.1, numberFormat: 'n'})).change(function() {
		var z = parseFloat(this.value);
		F.zsummary.text(numeral(z).format('0.[00]'));
		F.obj.z = z;
		F.updateViewer();
		});

	//radius
	$('<label>Radius:</label>').appendTo(locationdiv);
	this.radius = $('<input>').appendTo(locationdiv).addClass('radiusinput');
	this.radius.spinner(spinObject(F.radius,{step: 0.1, numberFormat: 'n'})).change(function() {
		F.rsummary.text(numeral(this.value).format('0.[00]'));
		F.obj.radius = this.value;
		F.updateViewer();
		});

	//orientation (for hbonds and aromatic)
	var orientdiv = this.orientdiv = $('<div>').appendTo(editdiv).addClass('orientdiv');
	var theta = null, phi = null;
	
	
	var setVector = function(t, p) {
		//t and p are in degrees, create svector i nobj
		var theta = t*Math.PI/180;
		var psi = p*Math.PI/180;
		F.obj.svector = {
			x: Math.sin(theta)*Math.cos(psi),
			y: Math.sin(theta)*Math.sin(psi),
			z: Math.cos(theta)
		};
	};
	
	this.orientenabled = $('<input type="checkbox">').appendTo(orientdiv).change(
			function() {
				if($(this).prop( "checked")) {
					theta.spinner('option','disabled',false);
					phi.spinner('option','disabled',false);
					F.obj.vector_on = 0;
				}
				else {
					theta.spinner('option','disabled',true);
					phi.spinner('option','disabled',true);
					F.obj.vector_on = 1;
					setVector(theta.val(), phi.val());						
				}
				F.updateViewer();
			}
			);
	
	$('<label>&theta;:</label>').appendTo(orientdiv);
	theta = this.theta = $('<input>').appendTo(orientdiv).addClass('orientinput');
	this.theta.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(function() {
		setVector(theta.val(), phi.val());
		F.updateViewer();
		});
	
	$('<label>&phi;:</label>').appendTo(orientdiv);
	phi = this.phi = $('<input>').appendTo(orientdiv).addClass('orientinput');
	this.phi.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(function() {			
		setVector(theta.val(), phi.val());
		F.updateViewer();
		});
	
	this.orientdiv.hide();
	
	//size for hydrophobic
	var sizediv = this.sizediv = $('<div>').appendTo(editdiv).addClass('sizediv');		
	
	this.minsize = $('<input>').appendTo(sizediv).addClass('sizeinput');
	this.minsize.spinner(spinObject(F.minsize,{step: 1, numberFormat: 'n'})).change(function() {
		F.obj.minsize = this.value;
		F.updateViewer();
	});
	
	$('<label> &le; #Atoms &le;</label>').appendTo(sizediv);
	this.maxsize = $('<input>').appendTo(sizediv).addClass('sizeinput');
	this.maxsize.spinner(spinObject(F.maxsize,{step: 1, numberFormat: 'n'})).change(function() {
		F.obj.maxsize = this.value;
		F.updateViewer();
	});
	this.sizediv.hide();
	
	this.setFeature(fobj);
	
	delete this.updateViewer;
	this.viewer = viewer;
	this.updateViewer(); //call prototype to actually draw feature
	
	features.accordion( "refresh" ); //update accordian
	features.accordion("option","active",features.children().length-1);

}

//set the feature to fobj, fill in ui
Feature.prototype.setFeature = function(fobj) {
	//this.obj will be set by the change handlers
	this.select.val(fobj.name).trigger('change');
	this.select.selectmenu("refresh");
	this.x.val(fobj.x).trigger('change');
	this.y.val(fobj.y).trigger('change');
	this.z.val(fobj.z).trigger('change');
	this.radius.val(fobj.radius).trigger('change');	
	this.enabled.toggles(fobj.enabled);
	this.obj.enabled = !!fobj.enabled;
	this.orientenabled.prop('checked', fobj.vector_on == 1).trigger('change');
	
	if(!$.isNumeric(fobj.minsize))
		this.obj.minsize = "";
	else
		this.minsize.val(fobj.minsize).trigger('change');
	
	if(!$.isNumeric(fobj.maxsize))
		this.obj.maxsize = "";
	else
		this.maxsize.val(fobj.maxsize).trigger('change');
	
	if(!fobj.svector) {
		this.obj.svector = {x:1,y:0,z:0};
	}
	else {
		var x = fobj.svector.x;
		var y = fobj.svector.y;
		var z = fobj.svector.z;
		var theta = 180*Math.acos(z/(Math.sqrt(x*x+y*y+z*z)))/Math.PI;
		var phi = 180*Math.atan2(y,x)/Math.PI;
		this.theta.val(theta).trigger('change');
		this.phi.val(phi).trigger('change');
	}
};

Feature.prototype.featureNames = ['Aromatic','HydrogenDonor', 'HydrogenAcceptor', 
		'Hydrophobic', 'NegativeIon', 'PositiveIon'];

Feature.prototype.updateViewer = function() {
	//anything that changes the geometry requires a new shape 
	//(position, arrow orientation, radius)
	
	if(this.shape !== null) {
		this.viewer.removeFeature(this.shape);
	}
	if(this.obj.enabled) {
		var F = this;
		this.shape = this.viewer.addFeature(this.obj, function() {
			if(F.selected) {
				F.deselectFeature();
			} else {
				F.selectFeature();
			}
		});
	}
};


//display in selected style
Feature.prototype.selectFeature = function() {
	this.viewer.selectFeature(this.shape);
	this.container.addClass("selectedFeature");
	this.selected = true;
};

//remove selection style
Feature.prototype.deselectFeature = function() {
	this.viewer.unselectFeature(this.shape);
	this.container.removeClass("selectedFeature");
	this.selected = false;
};


/**
 * jQuery Fileinput Plugin v3.2.0
 *
 * Copyright 2013, Hannu Leinonen <hleinone@gmail.com>
 * Dual licensed under the MIT and GPL licenses:
 *   http://www.opensource.org/licenses/mit-license.php
 *   http://www.gnu.org/licenses/gpl.html
 */
(function($) {
    $.support.cssPseudoClasses = (function() {
        try {
            var input = $("<input type='checkbox' checked/>").appendTo('body');
            var style = $('<style type="text/css" />').appendTo('head');
            style.text("input:checked {width: 200px; display: none}");
            var support = input.css('width') == "200px";
            input.remove();
            style.remove();
            return support;
        } catch (e) {
            return false;
        }
    })();

    $.fn.fileinput = function(replacement) {
        var selector = this;
        if (!replacement) {
        	replacement = "<button class=\"fileinput\">Browse...</button>";
        }
        selector.each(function() {
            var element = $(this);
            element.wrap("<div class=\"fileinput-wrapper\" style=\" position: relative; display: inline-block;\" />");

            element.parent().mousemove(function(e) {
                var offL, offT, el = $(this);

                offL = el.offset().left;
                offT = el.offset().top;

                el.find("input").css({
                    "left":e.pageX - offL - el.find("input[type=file]").width() + 30,
                    "top":e.pageY - offT - 10
                });
            });

            element.attr("tabindex", "-1").css({filter: "alpha(opacity=0)", "-moz-opacity": 0, opacity: 0, position: "absolute", "z-index": -1});
            element.before(replacement);
            element.prev().addClass("fileinput");
            if (!$.support.cssPseudoClasses) {
                element.css({"z-index":"auto", "cursor":"pointer"});
                element.prev(".fileinput").css("z-index", -1);
                element.removeAttr("tabindex");
                element.prev(".fileinput").attr("tabindex", "-1");
                element.hover(function() {
                    $(this).prev(".fileinput").addClass("hover");
                }, function() {
                    $(this).prev(".fileinput").removeClass("hover");
                }).focusin(function() {
                    $(this).prev(".fileinput").addClass("focus");
                }).focusout(function() {
                    $(this).prev(".fileinput").removeClass("focus");
                }).mousedown(function() {
                    $(this).prev(".fileinput").addClass("active");
                }).mouseup(function() {
                    $(this).prev(".fileinput").removeClass("active");
                });
            } else {
                element.prev(".fileinput").click(function() {
                    element.click();
                });
                element.prev(":submit.fileinput").click(function(event) {
                    event.preventDefault();
                });
            }
        });
        return selector;
    };
})(jQuery);

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

var Pharmit = {};
$(document).ready(function() {
	var gup = function( name )
	{
	  name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
	  var regexS = "[\\?&]"+name+"=([^&#]*)";
	  var regex = new RegExp( regexS );
	  var results = regex.exec( window.location.href );
	  if( results === null )
	    return null;
	  else
	    return results[1];
	};
	
	
	Pharmit.server = '/fcgi-bin/pharmitserv.fcgi';
	var element = $('#pharmit');
	var viewer = new Pharmit.Viewer(element);
	var phresults = new Pharmit.PharmaResults(element, viewer);
	var query = new Pharmit.Query(element, viewer);
	
	//look for session in url
	if(gup('SESSION'))
	{		
		$.get(decodeURI(gup('SESSION')), function(ret) {
			query.loadSession(ret);
		});
	}

		
	
	//work around jquery bug
	$("button, input[type='button'], input[type='submit']").button()
    .bind('mouseup', function() {
        $(this).blur();     // prevent jquery ui button from remaining in the active state
    });	
});
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

Pharmit.PharmaResults = (function() {
	// private class variables and functions

	function PharmaResults(element) {
		//private variables and functions
		
		
		
		//public variables and functions
	}

	return PharmaResults;
})();
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
		
		var loadSession = this.loadSession = function(data) {
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
		querydiv.resizable({handles: "e",
			resize: function(event, ui) {
				viewer.setLeft(ui.size.width);
				}
			});
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
		features.accordion({header: "> div > h3", 
			animate: true, 
			collapsible: true,
			heightStyle:'content',
			beforeActivate: function( event, ui ) { 
				var fdiv = null;
				if(ui.newHeader.length > 0) { //being activated
					fdiv = ui.newHeader.parent();
					fdiv.get(0).feature.selectFeature();
				}
				if(ui.oldHeader.length > 0) {
					 fdiv = ui.oldHeader.parent();					
					fdiv.get(0).feature.deselectFeature();
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

		filters.accordion({animate: true, collapsible: true, heightStyle:'content'});
		
		
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
		
		viewer.setLeft(querydiv.width());

	}

	return Query;
})();
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

Pharmit.Viewer = (function() {
	"use strict";
	// private class variables and functions

	function Viewer(element) {
		//private variables and functions
		var ec = $3Dmol.elementColors;
		var colorStyles =  {
                none:  "None",
                rasmol:  "RasMol", 
                whiteCarbon: "White", 
                greenCarbon: "Green",
                cyanCarbon:  "Cyan", 
                magentaCarbon:"Magenta", 
                yellowCarbon:"Yellow", 
                orangeCarbon:"Orange", 
                purpleCarbon:"Purple", 
                blueCarbon:  "Blue"
		};
		
		var margins = {left: 0, right: 0}; 
		
		var featureColors = {'Aromatic': 'purple', 'HydrogenDonor': '0xeeeeee', 'HydrogenAcceptor': 'orange', 
							'Hydrophobic': 'green', 'NegativeIon': 'red', 'PositiveIon': 'blue'};

		
		var modelsAndStyles = {
				'Ligand': {model: null,
					style: {stick:{radius: 0.15}}
					},
				'Receptor': {model: null,
					style: {cartoon: {},
							line: {linewidth:3.0}}
					},
				'Results': {model: null,
						style: {stick: {}}
						}
		};
		var nonBondStyle = {sphere: {radius: 0.2}};
		var surface = null;
		var surfaceStyle = {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-0.8,0.8)}, opacity:0.8};
		var viewer = null;
		var shapes = [];
		
		var getExt = function(fname) {
			if(!fname) return "";
			var a = fname.split(".");
			if( a.length <= 1 ) {
			    return "";
			}
			return a.pop(); 
		};
		
		
		//create jquery selection object for picking molecule style
		var createStyleSelector = function(name, defaultval, vizgroup, callback) {
			var ret = $('<div>').appendTo(vizgroup);
			var id = name+"MolStyleSelect";
			$('<label for="'+id+'">'+name+' Style</label>').appendTo(ret).addClass('stylelabel');
			
			var select = $('<select name="'+id+'" id="'+id+'">').appendTo(ret).addClass('styleselector');
			$.each(colorStyles, function(key, value) {
				$('<option value="'+key+'">'+value+'</option>').appendTo(select);
			});
			
			select.val(defaultval);
			select.selectmenu({
				width: '8em', 
				appendTo: vizgroup, 
				change: function() {select.change();},
				position: {my: "left top", at: "left bottom", collision: "flip"}
			});
			select.change(function() {
				var scheme = this.value;
				var style = modelsAndStyles[name].style;
				$.each(style, function(key, substyle) {
					substyle.colorscheme = scheme;
				});
				var model = modelsAndStyles[name].model;
				if(model) model.setStyle({}, style);				
				select.selectmenu("refresh");
				viewer.render();
			});
			select.change();
		};
		
		//public variables and functions
		
		//add controls for change the styles of the div to the specified element
		this.appendViewerControls = function(vizgroup) {
			createStyleSelector("Ligand",'rasmol', vizgroup, null);
			createStyleSelector("Results",'rasmol', vizgroup, null);
			createStyleSelector("Receptor",'rasmol', vizgroup, null);
			
			//surface transparency
			$('<label for="surfaceopacity">Receptor Surface Opacity</label>').appendTo(vizgroup);
			var sliderdiv = $('<div>').addClass('surfacetransparencydiv').appendTo(vizgroup);
			$('<div id="surfaceopacity" name="surfaceopacity">').appendTo(sliderdiv)
				.slider({animate:'fast',step:0.05,'min':0,'max':1,'value':0.8,
					change: function(event, ui) { 
						surfaceStyle.opacity = ui.value;
						if(surface !== null) viewer.setSurfaceMaterialStyle(surface, surfaceStyle);
						viewer.render();
						}
				});
				
			
			//background color
			$('<label for="backgroundcolor">Background Color</label>').appendTo(vizgroup);
			var radiodiv = $('<div id="backgroundcolor">').appendTo(vizgroup);
			$('<input type="radio" id="whiteBackground" name="backgroundcolor"><label for="whiteBackground">White</label>').appendTo(radiodiv)
				.change(function() {
					if($(this).prop("checked")) {
						viewer.setBackgroundColor("white");
					}
					radiodiv.buttonset("refresh");
				}).prop("checked",true);
			$('<input type="radio" id="blackBackground" name="backgroundcolor"><label for="blackBackground">Black</label>').appendTo(radiodiv)
				.change(function() {
						if($(this).prop("checked")) {
							viewer.setBackgroundColor("black");
						}
						radiodiv.buttonset("refresh");
					});
			radiodiv.buttonset();
		};
		
		this.setReceptor = function(recstr, recname) {
			
			var receptor = modelsAndStyles.Receptor.model;

			//clear receptor
			if(receptor) viewer.removeModel(receptor);
			receptor = null;
			if(surface !== null) viewer.removeSurface(surface);
			surface = null;
			
			if(recstr) {
				var ext = getExt(recname);
				receptor = viewer.addModel(recstr, ext);
				receptor.setStyle({}, modelsAndStyles.Receptor.style);
				receptor.setStyle({bonds: 0}, nonBondStyle);
				modelsAndStyles.Receptor.model = receptor;
				viewer.render();
				//surface
				viewer.mapAtomProperties($3Dmol.applyPartialCharges,{model:receptor});
				surface = viewer.addSurface($3Dmol.SurfaceType.VDW, 
						surfaceStyle, 
						{model:receptor}, {bonds:0, invert:true});
				viewer.zoomTo();
			}
			viewer.render();
		};

		this.setLigand = function(ligstr, name) {
			
			var ligand = modelsAndStyles.Ligand.model;
			//lig receptor
			if(ligand) viewer.removeModel(ligand);
			ligand = null;
			
			if(ligstr) { 
				var ext = getExt(name);
				ligand = viewer.addModel(ligstr, ext);
				ligand.setStyle({}, modelsAndStyles.Ligand.style);
				modelsAndStyles.Ligand.model = ligand;
				viewer.zoomTo({model: ligand});
			}
			viewer.render();
		};

		
		this.setView = function(view) {
			if(view) viewer.setView(view);
		};
		
		this.getView = function() {
			return viewer.getView();
		};
		
		//add a feature as specified by fobj
		//returns an identifier for referencing the feature later (e.g., removeFeature)
		this.addFeature = function(fobj, clickHandler) {
			var sphere = {
				center: {x: fobj.x,
				y: fobj.y,
				z: fobj.z},
				radius: fobj.radius,
				color: featureColors[fobj.name],
				wireframe: true,
				linewidth: 1.5,
				clickable: true,
				callback: clickHandler
			};
			
			var shape = {sphere: null, arrows: [], label: null};
			shape.sphere = viewer.addSphere(sphere);
			viewer.render();
			shapes.push(shape);
			return shapes.length-1;
		};
		
		//change style of feature
		this.selectFeature = function(s) {
			var shape = shapes[s];
			if(shape && shape.sphere) {
				shape.sphere.updateStyle({wireframe: false});
				viewer.render();
			}
		};
		
		this.unselectFeature = function(s) {
			var shape = shapes[s];
			if(shape && shape.sphere) {
				shape.sphere.updateStyle({wireframe: true});
				viewer.render();
			}
		};
		
		this.removeFeature = function(s) {
			var shape = shapes[s];
			if(shape) {
				if(shape.sphere) viewer.removeShape(shape.sphere);
				$.each(shape.arrows, function(i, arrow) {
					viewer.removeShape(arrow);
				});
				if(shape.label) viewer.removeLabel(shape.label);
				viewer.render();
			}
			delete shapes[s];
			//clear back of array 
			while (shapes.length > 0 && typeof (shapes[shapes.length - 1]) === "undefined")
				shapes.pop();
		};
		
		//specify size of left div so we can move the center point of the viewer
		this.setLeft = function(x) {
			var dx = x-margins.left;
			margins.left = x;
			viewer.translate(dx, 0);
		};
		
		//specify size of right div so we can move the center point of the viewer
		this.setRight = function(x) {
			var dx = margins.right-x;
			margins.right = x;
			viewer.translate(dx, 0);
		};
		
		//initialization code
		viewer = new $3Dmol.GLViewer(element);
		viewer.setBackgroundColor('white');
		
	}

	return Viewer;
})();