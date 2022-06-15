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
		
		var margins = {left: 0, right: 0}; 
		
		var featureColors = {'Aromatic': 'purple', 'HydrogenDonor': '0xf0f0f0', 'HydrogenAcceptor': 'orange', 
							'Hydrophobic': 'green', 'NegativeIon': 'red', 'PositiveIon': 'blue', 'ExclusionSphere': 'grey',
							'InclusionSphere': 'yellow'};
        var featureQuality = {'Aromatic': 1.5, 'HydrogenDonor': 2, 'HydrogenAcceptor': 1, 
                            'Hydrophobic': 1, 'NegativeIon': 1.5, 'PositiveIon': 1.5, 'ExclusionSphere': 1.5,
                            'InclusionSphere': 1.5};                            

        var surface = null;
        var savedsurface = null;
        var surfaceStyle = {map:{prop:'partialCharge',scheme:new $3Dmol.Gradient.RWB(-0.8,0.8)}, opacity:0.8};
        var viewer = null;
        var shapes = [];
                		
		var modelsAndStyles = {
				'Ligand': {model: null,
					defaultColor: '#C8C8C8',
					colorscheme: $.extend({},$3Dmol.elementColors.defaultColors),
					selectedstyle: 'stick',
					styles: {
						stick: {name: "Stick",
							style: {stick:{radius: 0.15}},
							nonbond: {sphere: {radius: 0.15}}
						},
						wire: {name: "Wire",
							style: {line:{linewidth: 2.0}},
							nonbond: {sphere: {radius: 0.05}}
						},
						sphere: {name: "Sphere",
							style: {sphere:{}},
							nonbond: {sphere: {}}
						},
						none: {name: "None",
							style: {},
							nonbond: {}
						}
					}
					},
				'Receptor':{model: null,
					defaultColor: '#C8C8C8',
					colorscheme: $.extend({},$3Dmol.elementColors.defaultColors),
					backbone: null,
					selectedstyle: 'cartoonwire',
					styles: {
						stick: {name: "Stick",
							style: {stick:{radius: 0.1}},
							nonbond: {sphere: {radius: 0.1}}
						},
						wire: {name: "Wire",
							style: {line:{linewidth: 3.0}},
							nonbond: {sphere: {radius: 0.05}}
						},
						sphere: {name: "Sphere",
							style: {sphere:{}},
							nonbond: {sphere: {}}
						},
						cartoon: {name: "Cartoon",
							style: {cartoon:{}},
							nonbond: {sphere: {radius: 0.1}}
							},
						cartoonwire: {name: "Cartoon+Wire",
							style: {cartoon: {}, line: {linewidth: 3.0}},
							nonbond: {sphere: {radius: 0.2}}
							},
                        bindingsite: {name: "Binding Site",
                            style: {cartoon: {}, line: {linewidth: 3.0}},
                            nonbond: {sphere: {radius: 0.2}},
                            dynamic: function() {
                                var lig = modelsAndStyles.Ligand.model;
                                var rec = modelsAndStyles.Receptor.model;
                                if(lig && rec) {
                                    //obnoxiously, can't use within selector across models
                                    var ligatoms = lig.selectedAtoms(); 
                                    var close_to_lig = function(atom) {
                                      for(var i = 0; i < ligatoms.length; i++) {
                                        var a = ligatoms[i];
                                        var sqdist = (atom.x-a.x)*(atom.x-a.x) + (atom.y-a.y)*(atom.y-a.y) + (atom.z-a.z)*(atom.z-a.z);
                                        if(sqdist < 36.0) { //6A cutoff
                                            return true;
                                        }                                        
                                      }  
                                      return false;
                                    };
                                    rec.setStyle({model: rec, byres: true, predicate: close_to_lig},
                                                   {stick:{radius:0.1,colorscheme:modelsAndStyles.Receptor.colorscheme}},true);
                                                   
                                    if(surface !== null) {
                                        viewer.setSurfaceMaterialStyle(surface, {opacity:0});
                                        savedsurface = surface;   
                                    }
                                    var surfacepromise = viewer.addSurface($3Dmol.SurfaceType.VDW, 
                                            surfaceStyle, {model: rec, byres: true, predicate: close_to_lig},
                                            {model:rec, bonds: 0, invert:true});
                                    surface = surfacepromise.surfid;                                                   
                                } 
                              } 
                            },
						none: {name: "None",
							style: {},
							nonbond: {}
						}
					}
					}, 
				'Results': {model: null,
					defaultColor: 'gray',
					colorscheme: $.extend({},$3Dmol.elementColors.defaultColors),
					selectedstyle: 'stick',
					styles: {
						stick: {name: "Stick",
							style: {stick:{radius: 0.25}},
							nonbond: {sphere: {radius: 0.25}}
						},
						wire: {name: "Wire",
							style: {line:{linewidth: 4.0}},
							nonbond: {sphere: {radius: 0.1}}
						},
						sphere: {name: "Sphere",
							style: {sphere:{}},
							nonbond: {sphere: {}}
						},
						none: {name: "None",
							style: {},
							nonbond: {}
						}
					}
					}
		};
		
		var getExt = function(fname) { 
			if(!fname) return "";
			var a = fname.split(".");
			if( a.length <= 1 ) {
			    return "";
			}
			return a.pop(); 
		};
		
		//applies the current selectedstyle and colorscheme for name
		var updateStyle = function(name) {
			var rec = modelsAndStyles[name];
			var stylename = rec.selectedStyle;
			var s = rec.styles[stylename];
			var style = s.style;
			var nbond = s.nonbond;

			var model = rec.model;
            if(savedsurface) {
                //slight hack for binding site - if this exists reset surface
                if(surface) viewer.removeSurface(surface);
                surface = savedsurface;
                viewer.setSurfaceMaterialStyle(surface, surfaceStyle);
                savedsurface = null;
            }
			if(model) {
				model.setStyle({}, style);
				model.setStyle({bonds: 0}, nbond);
                if(s.dynamic) { //can run arbitrary commands
                    s.dynamic();
                }
			}				
			
			viewer.render();
		};
		
		//create jquery selection object for picking molecule style
		//adds a table row
		var createStyleSelector = function(name, table, callback) {
			var rec = modelsAndStyles[name];
			var ret = $('<tr>').appendTo(table).addClass('pharmit_styleselectrow');
			var id = name+"MolStyleSelect";
			$('<label for="'+id+'">'+name+':</label>').appendTo(ret).addClass('pharmit_stylelabel').appendTo($('<td>').appendTo(ret));
			
			var cell = $('<td nowrap>').appendTo(ret);
			var select = $('<select name="'+id+'" id="'+id+'">').appendTo(cell).addClass('pharmit_styleselector');
			$.each(rec.styles, function(key, value) {
				$('<option value="'+key+'">'+value.name+'</option>').appendTo(select);
			});
			
			select.val(rec.selectedstyle);
			select.selectmenu({
				width: '11em', 
				appendTo: table, 
				change: function() {select.change();},
				position: {my: "left top", at: "left bottom", collision: "flip"}
			});
			
			//workaround firefox bug - remove element style so css stylesheet takes effect
			select.selectmenu( "widget" ).css('width','');
			
			var colorscheme = rec.colorscheme;
			//give color scheme to all substyles, this is reference so change the original colorscheme should change the styles
			$.each(rec.styles, function(key, subrec) {
				$.each(subrec.style, function(key,value) {
					value.colorscheme = colorscheme;
				});
			});					
			
			select.change(function() {
				rec.selectedStyle = this.value;
				select.selectmenu("refresh");	       
				updateStyle(name);
			});
			
			var colorpicker = $('<input name="'+id+'color">').appendTo($('<td>').appendTo(ret));
			colorpicker.val(rec.defaultColor);
			colorpicker.change(function() {
				var c = this.value;
				colorpicker.spectrum("set",c);
				var color = parseInt(colorpicker.spectrum("get").toHex(),16); //number
				rec.colorscheme.C = color;
				updateStyle(name);
			});
			
			colorpicker.spectrum({
			    showPalette: true,
			    preferredFormat: "hex",
			    replacerClassName: 'ui-state-default ui-corner-all',
			    showPaletteOnly: true,
			    clickoutFiresChange: true,
			    palette: ['#C8C8C8', 'gray', 'white','green','cyan','magenta','yellow','orange','purple','blue'],
			    change: function(color) { 
			    	colorpicker.change();
			    }
			});		
			
			select.change();
			colorpicker.change();
			select.selectmenu("refresh");	       
	
		};
		
		//public variables and functions
		
		//add controls for change the styles of the div to the specified element
		this.appendViewerControls = function(vizgroup) {
			
			var table = $('<table>').appendTo(vizgroup);
			createStyleSelector("Ligand",  table, null);
			createStyleSelector("Results", table, null);
			createStyleSelector("Receptor",  table, null);
			
			//backbone color scheme
			var backdiv = $('<div>').addClass('pharmit_backbonediv').appendTo(vizgroup);
			$('<label for="receptorbackbone">Receptor Backbone:</label>').appendTo(backdiv);
			var bradiodiv = $('<div id="receptorbackbone">').appendTo(backdiv);
			$('<input type="radio" id="plainBackbone" name="receptorbackbone"><label for="plainBackbone">Plain</label>').appendTo(bradiodiv)
				.change(function() {
					if($(this).prop("checked")) {
						var rec = modelsAndStyles.Receptor;
						delete rec.styles.cartoon.style.cartoon.color;
						delete rec.styles.cartoonwire.style.cartoon.color;						
						updateStyle('Receptor');
					}
					bradiodiv.buttonset("refresh");
				}).prop("checked",true);
			$('<input type="radio" id="rainbowBackbone" name="receptorbackbone"><label for="rainbowBackbone">Rainbow</label>').appendTo(bradiodiv)
				.change(function() {
						if($(this).prop("checked")) {
							var rec = modelsAndStyles.Receptor;
							rec.styles.cartoon.style.cartoon.color = 'spectrum';
							rec.styles.cartoonwire.style.cartoon.color = 'spectrum';
							updateStyle('Receptor');
						}
						bradiodiv.buttonset("refresh");
					});
			bradiodiv.buttonset();
			
			//surface transparency
			var stdiv = $('<div>').addClass('pharmit_surfacetransparencydiv').appendTo(vizgroup);
			$('<label for="surfaceopacity">Receptor Surface Opacity:</label>').appendTo(stdiv);
			var sliderdiv = $('<div>').addClass('pharmit_surfacetransparencyslider').appendTo(stdiv);
			var surfaceinput = $('<input type="hidden" name="surfaceopacity">').appendTo(stdiv); //use named inputs for getting/setting values
			var defaultopacity = 0.8;
			surfaceinput.val(defaultopacity);
			var surfslider = $('<div id="surfaceopacity">').appendTo(sliderdiv)
				.slider({animate:'fast',step:0.05,'min':0,'max':1,'value':defaultopacity,
					change: function(event, ui) { 
						surfaceinput.val(ui.value).change();
						}
				});
			surfaceinput.change(function() {
				var val = surfaceinput.val();
                                surfaceStyle.opacity = val;
                                if(surface !== null) viewer.setSurfaceMaterialStyle(surface, surfaceStyle);
                                if(surfslider.slider("value") != val) surfslider.slider("value",val);
                                viewer.render();
			});	
			//background color
			var bcdiv = $('<div>').addClass('pharmit_backgroundcolordiv').appendTo(vizgroup);
			$('<label for="backgroundcolor">Background Color:</label>').appendTo(bcdiv);
			var radiodiv = $('<div id="backgroundcolor">').appendTo(bcdiv);
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
		
		//amount to offset viewer position by based on morgins
		var xoffset = function() {
			return (margins.left-margins.right)/2;
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
				modelsAndStyles.Receptor.model = receptor;
				updateStyle("Receptor");			
				
				//surface
				viewer.mapAtomProperties($3Dmol.applyPartialCharges,{model:receptor});
				var surfacepromise = viewer.addSurface($3Dmol.SurfaceType.VDW, 
						surfaceStyle, {model:receptor, bonds: 0, invert:true});
				surface = surfacepromise.surfid;
				viewer.zoomTo({});
                
                //hover labels
                var hover_label = null;
                receptor.setHoverable({},true, 
                    function(atom){  //hover
                      if(atom.resn && atom.resi) {
                        hover_label = viewer.addLabel(atom.resn+atom.resi+":"+atom.atom,
                                {position: atom, fontSize: 12, backgroundColor: "#2c1943", backgroundOpacity: 0.5, alignment: "bottomCenter"});
                        viewer.render();
                      }
                    },
                    function(){ viewer.removeLabel(hover_label);viewer.render();} //unhover
                    );
			}
			else
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
				modelsAndStyles.Ligand.model = ligand;
				updateStyle("Ligand");
                updateStyle('Receptor'); //in case binding site is set
				viewer.zoomTo({model: ligand});
			}
			else
				viewer.render();
		};

		this.setResult = function(molstr) { //assumed sdf
			//remove current result
			var mol = modelsAndStyles.Results.model;
			if(mol) viewer.removeModel(mol);
			
			if(molstr) {
				mol = viewer.addModel(molstr, "sdf");
				modelsAndStyles.Results.model = mol;
				updateStyle("Results");
				viewer.zoomTo({model: mol});
			}
			else
				viewer.render();
		};
		
		this.setView = function(view) {
			if(view) viewer.setView(view);
			else viewer.zoomTo(); //at least center on objects
		};
		
		this.getView = function() {
			return viewer.getView();
		};
		
		//add a mesh object, returns identifier
		this.addMesh = function(mesh, style) {
			if(style.kind == "exclusive") {
				mesh.color = "grey";
			} else if(style.kind == "inclusive") {
				mesh.color = "yellow";
			}
			$.extend(mesh, style);
			var ret = viewer.addCustom(mesh);
			viewer.render();
			return ret;
		};
		
		
		this.updateMesh = function(m, newstyle) {
			if(m) {
				if(m.updateStyle) m.updateStyle(newstyle);
				viewer.render();
			}
		};
		
		//removes previously created mesh object m
		this.removeMesh = function(m) {
			viewer.removeShape(m);
			viewer.render();
		};
		
		//hides, but does not remove, mesh object m
		this.hideMesh = function(m) {
			
		};
		
		//restores visibility of mesh object m
		this.unhideMesh = function(m) {
			
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
                quality: featureQuality[fobj.name],
				wireframe: true,
				linewidth: 1.5,
				clickable: true,
				callback: clickHandler
			};
			
			var shape = {sphere: null, arrows: [], label: null};
			shape.sphere = viewer.addSphere(sphere);
			if(fobj.selected)
				shape.sphere.updateStyle({wireframe: false});

			if(fobj.hasvec && fobj.vector_on && fobj.svector) {
				//draw arrow
				var vec = new $3Dmol.Vector3(fobj.svector.x, fobj.svector.y, fobj.svector.z);
				var len = fobj.radius+1.0;
				var mid = (len-0.5)/len; //where arrowhead starts as a ratio
				vec = vec.normalize();
				var start = vec.clone().multiplyScalar(fobj.radius).add(fobj);
				var end = vec.clone().multiplyScalar(len).add(fobj);
				var arrow = {
					start: start,
					end: end,
					radius: 0.075,
					radiusRatio: 2.0,
					mid: mid,
					wireframe: !fobj.selected,
					color: featureColors[fobj.name]
				};
				shape.arrows.push(viewer.addArrow(arrow));
				
				if(fobj.name == "Aromatic") { //double arrow
					start = vec.clone().multiplyScalar(-fobj.radius).add(fobj);
					end = vec.clone().multiplyScalar(-len).add(fobj);
					arrow.start = start;
					arrow.end = end;
					shape.arrows.push(viewer.addArrow(arrow));
				}
			}
			
			if(fobj.name == "Hydrophobic") {
				//may have size
				var label = fobj.minsize + ":" + fobj.maxsize;
				if(label != ":") {
					var lab = {
							position: {x: fobj.x, y: fobj.y, z: fobj.z},
							showBackground: true,
							fontColor: 'black',
							backgroundColor: featureColors[fobj.name],
							backgroundOpacity: 0.5,
							alignment: $3Dmol.SpriteAlignment.center
					};
					shape.label = viewer.addLabel(label, lab);
				}
			}
			viewer.render();
			shapes.push(shape);
			return shapes.length-1;
		};
		
		//change style of feature
		this.selectFeature = function(s) {
			var shape = shapes[s];
			if(shape && shape.sphere) {
				shape.sphere.updateStyle({wireframe: false});
				
				$.each(shape.arrows, function(i, arrow) {
					arrow.updateStyle({wireframe:false});
				});
				viewer.render();
			}
		};
		
		this.unselectFeature = function(s) {
			var shape = shapes[s];
			if(shape && shape.sphere) {
				shape.sphere.updateStyle({wireframe: true});
				$.each(shape.arrows, function(i, arrow) {
					arrow.updateStyle({wireframe:true});
				});
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
			viewer.translate(dx/2, 0);
		};
		
		//specify size of right div so we can move the center point of the viewer
		this.setRight = function(x) {
			var dx = margins.right-x;
			margins.right = x;
			viewer.translate(dx/2, 0);
		};
		
		var savedRender = null;
		this.disableRendering = function() {
			savedRender = viewer.render;
			viewer.render = function() {};
		};
		
		this.enableRendering = function() {
			if(savedRender) viewer.render = savedRender;
			viewer.render();
		};
		//initialization code
		viewer = new $3Dmol.GLViewer(element);
		viewer.setBackgroundColor('white');
        element.bind('mousemove touchmove', viewer._handleMouseMove); // need to grab even when things on top
		
	}

	return Viewer;
})();
