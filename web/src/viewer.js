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
	// private class variables and functions

	function Viewer(element) {
		//private variables and functions
		var ec = $3Dmol.elementColors;
		var colorStyles =  [//each element has name and color (3dmol)
		                    {name: "None", color: null},
		                    {name: "RasMol", color: ec.rasmol},
		                    {name: "White", color: ec.whiteCarbon},
		                    {name: "Green", color: ec.greenCarbon},
		                    {name: "Cyan", color: ec.cyanCarbon},
		                    {name: "Magenta", color: ec.magentaCarbon},
		                    {name: "Yellow", color: ec.yellowCarbon},
		                    {name: "Orange", color: ec.orangeCarbon},
		                    {name: "Purple", color: ec.purpleCarbon},
		                    {name: "Blue", color: ec.blueCarbon}
		                    ];
		
		var pickCallback = null;
		var viewer = null;
		var receptor = null;
		var ligand = null;
		var results = [];
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
		var createStyleSelector = function(name, defaultval, callback) {
			var ret = $('<div>');
			var id = name+"MolStyleSelect";
			$('<label for="'+id+'">'+name+' Style</label>').appendTo(ret).addClass('stylelabel');
			
			var select = $('<select name="'+id+'" id="'+id+'">').appendTo(ret).addClass('styleselector');
			for(var i = 0, n = colorStyles.length; i < n; i++) {
				$('<option value="'+colorStyles[i].name.toLowerCase()+'">'+colorStyles[i].name+'</option>').appendTo(select);
			}
			
			select.val(defaultval);
			select.selectmenu({width: 120});

			return ret;
		};
		
		//public variables and functions
		
		//add controls for change the styles of the div to the specified element
		this.appendViewerControls = function(vizgroup) {
			createStyleSelector("Ligand",1, null).appendTo(vizgroup);
			createStyleSelector("Results",1, null).appendTo(vizgroup);
			createStyleSelector("Receptor",2, null).appendTo(vizgroup);
			
			//surface transparency
			$('<label for="surfacetransparency">Receptor Surface Opacity</label>').appendTo(vizgroup);
			var sliderdiv = $('<div>').addClass('surfacetransparencydiv').appendTo(vizgroup);
			$('<div id="surfacetransparency">').appendTo(sliderdiv).slider({animate:'fast',step:0.05,'min':0,'max':1,'value':0.8});
		};
		
		this.setReceptor = function(recstr, recname) {
			
			if(!recstr) {
				//clear receptor
				if(receptor) viewer.removeModel(receptor);
				receptor = null;
			}
			else {
				var ext = getExt(recname);
				receptor = viewer.addModel(recstr, ext);
			}
			viewer.zoomTo();
			viewer.render();
		};

		this.setLigand = function(ligstr) {
			
			if(!recstr) {
				//clear receptor
				if(receptor) viewer.removeModel(receptor);
				receptor = null;
			}
			else {
				var ext = getExt(recname);
				receptor = viewer.addModel(recstr, ext);
			}
			viewer.zoomTo();
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