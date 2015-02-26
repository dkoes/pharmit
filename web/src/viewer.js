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
			select.selectmenu({width: '8em', appendTo: vizgroup, change: function() {select.change();}});
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
			createStyleSelector("Receptor",'whiteCarbon', vizgroup, null);
			
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
				viewer.mapAtomProperties($3Dmol.partialCharges);
				surface = viewer.addSurface($3Dmol.SurfaceType.VDW, 
						surfaceStyle, 
						{model:receptor}, {bonds:0, invert:true});
			}
			viewer.zoomTo();
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
			}
			viewer.zoomTo({model: ligand});
			viewer.render();
		};

		
		this.setView = function(view) {
			if(view) viewer.setView(view);
		};
		
		this.getView = function() {
			return viewer.getView();
		};
		
		//initialization code
		viewer = new $3Dmol.GLViewer(element);
		viewer.setBackgroundColor('white');
		
	}

	return Viewer;
})();