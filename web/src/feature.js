
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
	
	//for prettier display round values to 3 decimal places
	var round = function(x) { return Math.round(x*1000)/1000;};
	
	var c = $('<div>').appendTo(locationdiv).addClass('coorddiv');
	$('<label>x:</label>').appendTo(c);
	this.x = $('<input>').appendTo(c).addClass('coordinput').change(function() {
		//change x value
		var x = this.value = round(parseFloat(this.value));
		F.xsummary.text(numeral(x).format('0.[00]'));
		F.obj.x = round(x);
		F.updateViewer();
		});
	this.x.spinner(spinObject(F.x,{step: 0.1, numberFormat: 'n'}));
	
	c = $('<div>').appendTo(locationdiv).addClass('coorddiv');
	$('<label>y:</label>').appendTo(c);
	this.y = $('<input>').appendTo(c).addClass('coordinput');
	this.y.spinner(spinObject(F.y,{step: 0.1, numberFormat: 'n'})).change(function() {
		//change y value
		var y = this.value = round(parseFloat(this.value));
		F.ysummary.text(numeral(y).format('0.[00]'));
		F.obj.y = round(y);
		F.updateViewer();
		});

	c = $('<div>').appendTo(locationdiv).addClass('coorddiv');
	$('<label>z:</label>').appendTo(c);
	this.z = $('<input>').appendTo(c).addClass('coordinput');
	this.z.spinner(spinObject(F.z,{step: 0.1, numberFormat: 'n'})).change(function() {
		var z = this.value = round(parseFloat(this.value));
		F.zsummary.text(numeral(z).format('0.[00]'));
		F.obj.z = z;
		F.updateViewer();
		});

	//radius
	c = $('<div>').appendTo(locationdiv).addClass('coorddiv');
	$('<label>Radius:</label>').appendTo(c);
	this.radius = $('<input>').appendTo(c).addClass('radiusinput');
	this.radius.spinner(spinObject(F.radius,{step: 0.1, numberFormat: 'n'})).change(function() {
		F.rsummary.text(numeral(this.value).format('0.[00]'));
		F.obj.radius = this.value = round(this.value);
		F.updateViewer();
		});

	//orientation (for hbonds and aromatic)
	var orientdiv = this.orientdiv = $('<div>').appendTo(editdiv).addClass('orientdiv');
	var theta = null, phi = null;
	
	
	var updateVector = function () {
		//parse text, round to integers, and calculate vector for orientation
		var p = Math.round(parseFloat(phi.val())) % 360;
		var t = Math.round(parseFloat(theta.val())) % 360;
		if(isNaN(p)) p = 0;
		if(isNaN(t)) t = 0;
		phi.val(p);
		theta.val(t);
		//convert to radians
		t = t*Math.PI/180;
		p = p*Math.PI/180;
		F.obj.svector = {
			x: Math.sin(t)*Math.cos(p),
			y: Math.sin(t)*Math.sin(p),
			z: Math.cos(t)
		};
		F.updateViewer();

	};
	
	this.orientenabled = $('<input type="checkbox">').appendTo(orientdiv).change(
			function() {
				if($(this).prop( "checked")) {
					theta.spinner('option','disabled',false);
					phi.spinner('option','disabled',false);
					F.obj.vector_on = 1;
				}
				else {
					theta.spinner('option','disabled',true);
					phi.spinner('option','disabled',true);
					F.obj.vector_on = 0;
					updateVector();
				}
				F.updateViewer();
			}
			);
	
	var nowrap = $('<span>').addClass('nowrap').appendTo(orientdiv);
	$('<label>&theta;:</label>').appendTo(nowrap);
	theta = this.theta = $('<input>').appendTo(nowrap).addClass('orientinput');
	this.theta.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(updateVector);

	nowrap = $('<span>').addClass('nowrap').appendTo(orientdiv);
	$('<label>&phi;:</label>').appendTo(nowrap);
	phi = this.phi = $('<input>').appendTo(nowrap).addClass('orientinput');
	this.phi.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(updateVector);
	
	this.orientdiv.hide();
	
	//size for hydrophobic
	var sizediv = this.sizediv = $('<div>').appendTo(editdiv).addClass('sizediv');		
	
	this.minsize = $('<input>').appendTo(sizediv).addClass('sizeinput');
	this.minsize.spinner(spinObject(F.minsize,{step: 1, numberFormat: 'n'})).change(function() {
		F.obj.minsize = this.value;
		F.updateViewer();
	});
	
	$('<label> &le; #Atoms &le;</label>').appendTo(sizediv).addClass('nowrap');
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
	this.obj.selected = true;
};

//remove selection style
Feature.prototype.deselectFeature = function() {
	this.viewer.unselectFeature(this.shape);
	this.container.removeClass("selectedFeature");
	this.selected = false;
	this.obj.selected = false;
};

