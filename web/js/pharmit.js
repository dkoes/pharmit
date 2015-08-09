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
	feature.js
	Represents a single pharmacophore feature.
*/

/* object representing a single pharmacophore feature*/
function Feature(viewer, features, fobj, type) {
	
	var featureNames = null;
	if(type == Feature.INCLUSIVESHAPE) {
		featureNames = ['InclusionSphere'];
	} else if(type == Feature.EXCLUSIVESHAPE) {
		featureNames = ['ExclusionSphere'];
	} else { //pharmacophores
		featureNames = ['Aromatic','HydrogenDonor', 'HydrogenAcceptor', 
	                          		'Hydrophobic', 'NegativeIon', 'PositiveIon'];
	}
	
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
	
	this.container = $('<div>').appendTo(features).addClass('pharmit_featurediv');
	this.container.get(0).feature = this;
	this.container.disableSelection();
	
	//header has checkbox for enable, name, and a close icon
	var heading = $('<h3></h3>').appendTo(this.container);
	this.enabled = $('<div>').addClass('toggle toggle-light pharmit_togglediv').appendTo(heading).toggles().on('toggle',
			function(e, active) {
				if(active) {
					F.obj.enabled = true;
					F.container.addClass("pharmit_enabledfeature");
				}
				else {
					F.obj.enabled = false;	
					F.container.removeClass("pharmit_enabledfeature");
				}
				F.updateViewer();
			}
			).click(function(e) {return false;}); //needed to stop clickthrue
	
	var namediv = $('<div>').addClass('pharmit_featurenamediv').appendTo(heading);
	var namespan = $('<span>').addClass('pharmit_featurenameheading').appendTo(namediv);
	var closediv = $('<div>').addClass('pharmit_featureclose').appendTo(heading).click(function() {
		F.deleteFeature();
	});
	var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
	
	//summary has (x,y,z) Radius r
	var summary = $('<span>').appendTo(namediv).addClass('pharmit_featuresummary');		
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
	var select = this.select = $('<select>').addClass('pharmit_featureselect').appendTo(editdiv);
	$.each(featureNames, function(key,val) {
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
	select.selectmenu({width: "15em", change: function() {select.trigger('change');}});
	
	if(featureNames.length <= 1) {
		select.next('.ui-selectmenu-button').hide(); //don't bother showing if no choices
	}
	
	//position (x,y,z)
	var locationdiv = $('<div>').appendTo(editdiv).addClass('pharmit_locationdiv');
	
	//because spinners are insane and don't trigger a change event on the
	//underlying element and split up spinning and stopping and changing...
	var spinObject = function(elem, init) {
		var change = function() {elem.change();};
		init.spin = init.stop = init.change = change;
		return init;
	};
	
	//for prettier display round values to 3 decimal places
	var round = function(x) { return Math.round(x*1000)/1000;};
	
	//intelligently apply rounding and return final result
	var updateNumber = function(elem, text) {
		if($.isNumeric(elem.value)) {
			var x = parseFloat(elem.value);
			if(text) text.text(numeral(x).format('0.0[0]'));
			F.updateViewer();
			return round(x);
		}
	};
	
	var c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>x:</label>').appendTo(c);
	this.x = $('<input>').appendTo(c).addClass('pharmit_coordinput').change(function() {
			//change x value
			F.obj.x = updateNumber(this, F.xsummary);
		});
	this.x.spinner(spinObject(F.x,{step: 0.1, numberFormat: 'n'}));
	
	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>y:</label>').appendTo(c);
	this.y = $('<input>').appendTo(c).addClass('pharmit_coordinput');
	this.y.spinner(spinObject(F.y,{step: 0.1, numberFormat: 'n'})).change(function() {
			//change y value
			F.obj.y = updateNumber(this, F.ysummary);
		});

	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>z:</label>').appendTo(c);
	this.z = $('<input>').appendTo(c).addClass('pharmit_coordinput');
	this.z.spinner(spinObject(F.z,{step: 0.1, numberFormat: 'n'})).change(function() {
			F.obj.z = updateNumber(this, F.zsummary);
		});

	//radius
	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>Radius:</label>').appendTo(c);
	this.radius = $('<input>').appendTo(c).addClass('pharmit_radiusinput');
	this.radius.spinner(spinObject(F.radius,{step: 0.1, numberFormat: 'n'})).change(function() {
		F.obj.radius = updateNumber(this, F.rsummary);
		});

	//orientation (for hbonds and aromatic)
	var orientdiv = this.orientdiv = $('<div>').appendTo(editdiv).addClass('pharmit_orientdiv');
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
	
	var nowrap = $('<span>').addClass('pharmit_nowrap').appendTo(orientdiv);
	$('<label>&theta;:</label>').appendTo(nowrap);
	theta = this.theta = $('<input>').appendTo(nowrap).addClass('pharmit_orientinput');
	this.theta.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(updateVector);

	nowrap = $('<span>').addClass('pharmit_nowrap').appendTo(orientdiv);
	$('<label>&phi;:</label>').appendTo(nowrap);
	phi = this.phi = $('<input>').appendTo(nowrap).addClass('pharmit_orientinput');
	this.phi.spinner(spinObject(F.theta,{step: 1, numberFormat: 'n'})).change(updateVector);
	
	this.orientdiv.hide();
	
	//size for hydrophobic
	var sizediv = this.sizediv = $('<div>').appendTo(editdiv).addClass('pharmit_sizediv');		
	
	this.minsize = $('<input>').appendTo(sizediv).addClass('pharmit_sizeinput');
	this.minsize.spinner(spinObject(F.minsize,{step: 1, numberFormat: 'n'})).change(function() {
		F.obj.minsize = this.value;
		F.updateViewer();
	});
	
	$('<label> &le; #Atoms &le;</label>').appendTo(sizediv).addClass('pharmit_nowrap');
	this.maxsize = $('<input>').appendTo(sizediv).addClass('pharmit_sizeinput');
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

Feature.INCLUSIVESHAPE = 1;
Feature.EXCLUSIVESHAPE = 2;

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

//don't show in viewer, but don't change enabled state
Feature.prototype.hideFeature = function() {
	this.obj.hidden = true;
	if(this.shape !== null) this.viewer.removeFeature(this.shape);
	this.updateViewer();
};

//undow hid
Feature.prototype.unhideFeature = function() {
	this.obj.hidden = false;
	this.updateViewer();
};


Feature.prototype.updateViewer = function() {
	//anything that changes the geometry requires a new shape 
	//(position, arrow orientation, radius)
	
	if(this.shape !== null) {
		this.viewer.removeFeature(this.shape);
		this.shape = null;
	}
	if(this.obj.enabled && !this.obj.hidden) {
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
	this.container.addClass("pharmit_selectedFeature");
	this.selected = true;
	this.obj.selected = true;
};

//remove selection style
Feature.prototype.deselectFeature = function() {
	this.viewer.unselectFeature(this.shape);
	this.container.removeClass("pharmit_selectedFeature");
	this.selected = false;
	this.obj.selected = false;
};

//remove completely
Feature.prototype.deleteFeature = function() {
	
	//remove from viewer
	if(this.shape !== null) this.viewer.removeFeature(this.shape);
	//remove from dom
	this.container.feature = null;
	this.container.remove();
};


/*
 * A JavaScript implementation of the RSA Data Security, Inc. MD5 Message
 * Digest Algorithm, as defined in RFC 1321.
 * Version 2.2 Copyright (C) Paul Johnston 1999 - 2009
 * Other contributors: Greg Holt, Andrew Kepert, Ydnar, Lostinet
 * Distributed under the BSD License
 * See http://pajhome.org.uk/crypt/md5 for more info.
 */

var Pharmit = Pharmit || {};
//dkoes - wrap
(function() {
/*
 * Configurable variables. You may need to tweak these to be compatible with
 * the server-side, but the defaults work in most cases.
 */

var hexcase = 0;   /* hex output format. 0 - lowercase; 1 - uppercase        */
var b64pad  = "";  /* base-64 pad character. "=" for strict RFC compliance   */

/*
 * These are the functions you'll usually want to call
 * They take string arguments and return either hex or base-64 encoded strings
 */
function hex_md5(s)    { return rstr2hex(rstr_md5(str2rstr_utf8(s))); }
function b64_md5(s)    { return rstr2b64(rstr_md5(str2rstr_utf8(s))); }
function any_md5(s, e) { return rstr2any(rstr_md5(str2rstr_utf8(s)), e); }
function hex_hmac_md5(k, d)
  { return rstr2hex(rstr_hmac_md5(str2rstr_utf8(k), str2rstr_utf8(d))); }
function b64_hmac_md5(k, d)
  { return rstr2b64(rstr_hmac_md5(str2rstr_utf8(k), str2rstr_utf8(d))); }
function any_hmac_md5(k, d, e)
  { return rstr2any(rstr_hmac_md5(str2rstr_utf8(k), str2rstr_utf8(d)), e); }

/*
 * Perform a simple self-test to see if the VM is working
 */
function md5_vm_test()
{
  return hex_md5("abc").toLowerCase() == "900150983cd24fb0d6963f7d28e17f72";
}

/*
 * Calculate the MD5 of a raw string
 */
function rstr_md5(s)
{
  return binl2rstr(binl_md5(rstr2binl(s), s.length * 8));
}

/*
 * Calculate the HMAC-MD5, of a key and some data (raw strings)
 */
function rstr_hmac_md5(key, data)
{
  var bkey = rstr2binl(key);
  if(bkey.length > 16) bkey = binl_md5(bkey, key.length * 8);

  var ipad = Array(16), opad = Array(16);
  for(var i = 0; i < 16; i++)
  {
    ipad[i] = bkey[i] ^ 0x36363636;
    opad[i] = bkey[i] ^ 0x5C5C5C5C;
  }

  var hash = binl_md5(ipad.concat(rstr2binl(data)), 512 + data.length * 8);
  return binl2rstr(binl_md5(opad.concat(hash), 512 + 128));
}

/*
 * Convert a raw string to a hex string
 */
function rstr2hex(input)
{
  var hex_tab = hexcase ? "0123456789ABCDEF" : "0123456789abcdef";
  var output = "";
  var x;
  for(var i = 0; i < input.length; i++)
  {
    x = input.charCodeAt(i);
    output += hex_tab.charAt((x >>> 4) & 0x0F)  +  hex_tab.charAt( x & 0x0F);
  }
  return output;
}

/*
 * Convert a raw string to a base-64 string
 */
function rstr2b64(input)
{
  var tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var output = "";
  var len = input.length;
  for(var i = 0; i < len; i += 3)
  {
    var triplet = (input.charCodeAt(i) << 16)
                | (i + 1 < len ? input.charCodeAt(i+1) << 8 : 0)
                | (i + 2 < len ? input.charCodeAt(i+2)      : 0);
    for(var j = 0; j < 4; j++)
    {
      if(i * 8 + j * 6 > input.length * 8) output += b64pad;
      else output += tab.charAt((triplet >>> 6*(3-j)) & 0x3F);
    }
  }
  return output;
}

/*
 * Convert a raw string to an arbitrary string encoding
 */
function rstr2any(input, encoding)
{
  var divisor = encoding.length;
  var i, j, q, x, quotient;

  /* Convert to an array of 16-bit big-endian values, forming the dividend */
  var dividend = Array(Math.ceil(input.length / 2));
  for(i = 0; i < dividend.length; i++)
  {
    dividend[i] = (input.charCodeAt(i * 2) << 8) | input.charCodeAt(i * 2 + 1);
  }

  /*
   * Repeatedly perform a long division. The binary array forms the dividend,
   * the length of the encoding is the divisor. Once computed, the quotient
   * forms the dividend for the next step. All remainders are stored for later
   * use.
   */
  var full_length = Math.ceil(input.length * 8 /
                                    (Math.log(encoding.length) / Math.log(2)));
  var remainders = Array(full_length);
  for(j = 0; j < full_length; j++)
  {
    quotient = Array();
    x = 0;
    for(i = 0; i < dividend.length; i++)
    {
      x = (x << 16) + dividend[i];
      q = Math.floor(x / divisor);
      x -= q * divisor;
      if(quotient.length > 0 || q > 0)
        quotient[quotient.length] = q;
    }
    remainders[j] = x;
    dividend = quotient;
  }

  /* Convert the remainders to the output string */
  var output = "";
  for(i = remainders.length - 1; i >= 0; i--)
    output += encoding.charAt(remainders[i]);

  return output;
}

/*
 * Encode a string as utf-8.
 * For efficiency, this assumes the input is valid utf-16.
 */
function str2rstr_utf8(input)
{
  var output = "";
  var i = -1;
  var x, y;

  while(++i < input.length)
  {
    /* Decode utf-16 surrogate pairs */
    x = input.charCodeAt(i);
    y = i + 1 < input.length ? input.charCodeAt(i + 1) : 0;
    if(0xD800 <= x && x <= 0xDBFF && 0xDC00 <= y && y <= 0xDFFF)
    {
      x = 0x10000 + ((x & 0x03FF) << 10) + (y & 0x03FF);
      i++;
    }

    /* Encode output as utf-8 */
    if(x <= 0x7F)
      output += String.fromCharCode(x);
    else if(x <= 0x7FF)
      output += String.fromCharCode(0xC0 | ((x >>> 6 ) & 0x1F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0xFFFF)
      output += String.fromCharCode(0xE0 | ((x >>> 12) & 0x0F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0x1FFFFF)
      output += String.fromCharCode(0xF0 | ((x >>> 18) & 0x07),
                                    0x80 | ((x >>> 12) & 0x3F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
  }
  return output;
}

/*
 * Encode a string as utf-16
 */
function str2rstr_utf16le(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode( input.charCodeAt(i)        & 0xFF,
                                  (input.charCodeAt(i) >>> 8) & 0xFF);
  return output;
}

function str2rstr_utf16be(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode((input.charCodeAt(i) >>> 8) & 0xFF,
                                   input.charCodeAt(i)        & 0xFF);
  return output;
}

/*
 * Convert a raw string to an array of little-endian words
 * Characters >255 have their high-byte silently ignored.
 */
function rstr2binl(input)
{
  var output = Array(input.length >> 2);
  for(var i = 0; i < output.length; i++)
    output[i] = 0;
  for( i = 0; i < input.length * 8; i += 8)
    output[i>>5] |= (input.charCodeAt(i / 8) & 0xFF) << (i%32);
  return output;
}

/*
 * Convert an array of little-endian words to a string
 */
function binl2rstr(input)
{
  var output = "";
  for(var i = 0; i < input.length * 32; i += 8)
    output += String.fromCharCode((input[i>>5] >>> (i % 32)) & 0xFF);
  return output;
}

/*
 * Calculate the MD5 of an array of little-endian words, and a bit length.
 */
function binl_md5(x, len)
{
  /* append padding */
  x[len >> 5] |= 0x80 << ((len) % 32);
  x[(((len + 64) >>> 9) << 4) + 14] = len;

  var a =  1732584193;
  var b = -271733879;
  var c = -1732584194;
  var d =  271733878;

  for(var i = 0; i < x.length; i += 16)
  {
    var olda = a;
    var oldb = b;
    var oldc = c;
    var oldd = d;

    a = md5_ff(a, b, c, d, x[i+ 0], 7 , -680876936);
    d = md5_ff(d, a, b, c, x[i+ 1], 12, -389564586);
    c = md5_ff(c, d, a, b, x[i+ 2], 17,  606105819);
    b = md5_ff(b, c, d, a, x[i+ 3], 22, -1044525330);
    a = md5_ff(a, b, c, d, x[i+ 4], 7 , -176418897);
    d = md5_ff(d, a, b, c, x[i+ 5], 12,  1200080426);
    c = md5_ff(c, d, a, b, x[i+ 6], 17, -1473231341);
    b = md5_ff(b, c, d, a, x[i+ 7], 22, -45705983);
    a = md5_ff(a, b, c, d, x[i+ 8], 7 ,  1770035416);
    d = md5_ff(d, a, b, c, x[i+ 9], 12, -1958414417);
    c = md5_ff(c, d, a, b, x[i+10], 17, -42063);
    b = md5_ff(b, c, d, a, x[i+11], 22, -1990404162);
    a = md5_ff(a, b, c, d, x[i+12], 7 ,  1804603682);
    d = md5_ff(d, a, b, c, x[i+13], 12, -40341101);
    c = md5_ff(c, d, a, b, x[i+14], 17, -1502002290);
    b = md5_ff(b, c, d, a, x[i+15], 22,  1236535329);

    a = md5_gg(a, b, c, d, x[i+ 1], 5 , -165796510);
    d = md5_gg(d, a, b, c, x[i+ 6], 9 , -1069501632);
    c = md5_gg(c, d, a, b, x[i+11], 14,  643717713);
    b = md5_gg(b, c, d, a, x[i+ 0], 20, -373897302);
    a = md5_gg(a, b, c, d, x[i+ 5], 5 , -701558691);
    d = md5_gg(d, a, b, c, x[i+10], 9 ,  38016083);
    c = md5_gg(c, d, a, b, x[i+15], 14, -660478335);
    b = md5_gg(b, c, d, a, x[i+ 4], 20, -405537848);
    a = md5_gg(a, b, c, d, x[i+ 9], 5 ,  568446438);
    d = md5_gg(d, a, b, c, x[i+14], 9 , -1019803690);
    c = md5_gg(c, d, a, b, x[i+ 3], 14, -187363961);
    b = md5_gg(b, c, d, a, x[i+ 8], 20,  1163531501);
    a = md5_gg(a, b, c, d, x[i+13], 5 , -1444681467);
    d = md5_gg(d, a, b, c, x[i+ 2], 9 , -51403784);
    c = md5_gg(c, d, a, b, x[i+ 7], 14,  1735328473);
    b = md5_gg(b, c, d, a, x[i+12], 20, -1926607734);

    a = md5_hh(a, b, c, d, x[i+ 5], 4 , -378558);
    d = md5_hh(d, a, b, c, x[i+ 8], 11, -2022574463);
    c = md5_hh(c, d, a, b, x[i+11], 16,  1839030562);
    b = md5_hh(b, c, d, a, x[i+14], 23, -35309556);
    a = md5_hh(a, b, c, d, x[i+ 1], 4 , -1530992060);
    d = md5_hh(d, a, b, c, x[i+ 4], 11,  1272893353);
    c = md5_hh(c, d, a, b, x[i+ 7], 16, -155497632);
    b = md5_hh(b, c, d, a, x[i+10], 23, -1094730640);
    a = md5_hh(a, b, c, d, x[i+13], 4 ,  681279174);
    d = md5_hh(d, a, b, c, x[i+ 0], 11, -358537222);
    c = md5_hh(c, d, a, b, x[i+ 3], 16, -722521979);
    b = md5_hh(b, c, d, a, x[i+ 6], 23,  76029189);
    a = md5_hh(a, b, c, d, x[i+ 9], 4 , -640364487);
    d = md5_hh(d, a, b, c, x[i+12], 11, -421815835);
    c = md5_hh(c, d, a, b, x[i+15], 16,  530742520);
    b = md5_hh(b, c, d, a, x[i+ 2], 23, -995338651);

    a = md5_ii(a, b, c, d, x[i+ 0], 6 , -198630844);
    d = md5_ii(d, a, b, c, x[i+ 7], 10,  1126891415);
    c = md5_ii(c, d, a, b, x[i+14], 15, -1416354905);
    b = md5_ii(b, c, d, a, x[i+ 5], 21, -57434055);
    a = md5_ii(a, b, c, d, x[i+12], 6 ,  1700485571);
    d = md5_ii(d, a, b, c, x[i+ 3], 10, -1894986606);
    c = md5_ii(c, d, a, b, x[i+10], 15, -1051523);
    b = md5_ii(b, c, d, a, x[i+ 1], 21, -2054922799);
    a = md5_ii(a, b, c, d, x[i+ 8], 6 ,  1873313359);
    d = md5_ii(d, a, b, c, x[i+15], 10, -30611744);
    c = md5_ii(c, d, a, b, x[i+ 6], 15, -1560198380);
    b = md5_ii(b, c, d, a, x[i+13], 21,  1309151649);
    a = md5_ii(a, b, c, d, x[i+ 4], 6 , -145523070);
    d = md5_ii(d, a, b, c, x[i+11], 10, -1120210379);
    c = md5_ii(c, d, a, b, x[i+ 2], 15,  718787259);
    b = md5_ii(b, c, d, a, x[i+ 9], 21, -343485551);

    a = safe_add(a, olda);
    b = safe_add(b, oldb);
    c = safe_add(c, oldc);
    d = safe_add(d, oldd);
  }
  return Array(a, b, c, d);
}

/*
 * These functions implement the four basic operations the algorithm uses.
 */
function md5_cmn(q, a, b, x, s, t)
{
  return safe_add(bit_rol(safe_add(safe_add(a, q), safe_add(x, t)), s),b);
}
function md5_ff(a, b, c, d, x, s, t)
{
  return md5_cmn((b & c) | ((~b) & d), a, b, x, s, t);
}
function md5_gg(a, b, c, d, x, s, t)
{
  return md5_cmn((b & d) | (c & (~d)), a, b, x, s, t);
}
function md5_hh(a, b, c, d, x, s, t)
{
  return md5_cmn(b ^ c ^ d, a, b, x, s, t);
}
function md5_ii(a, b, c, d, x, s, t)
{
  return md5_cmn(c ^ (b | (~d)), a, b, x, s, t);
}

/*
 * Add integers, wrapping at 2^32. This uses 16-bit operations internally
 * to work around bugs in some JS interpreters.
 */
function safe_add(x, y)
{
  var lsw = (x & 0xFFFF) + (y & 0xFFFF);
  var msw = (x >> 16) + (y >> 16) + (lsw >> 16);
  return (msw << 16) | (lsw & 0xFFFF);
}

/*
 * Bitwise rotate a 32-bit number to the left.
 */
function bit_rol(num, cnt)
{
  return (num << cnt) | (num >>> (32 - cnt));
}

Pharmit.hex_md5 = hex_md5;
})();
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
/*
 * Pretty much all html is dynamically created into the provided element.
 * Assumes jquery, 3Dmol, DataTables, jquery-toggles, and numeral.js are available..
 */
var Pharmit = Pharmit ||  {};
$(document).ready(function() {
	
	Pharmit.server = '/fcgi-bin/pharmitserv.fcgi';
	Pharmit.email = 'dkoes+pharmit@pitt.edu';
	
	//global variable checking - we should add nothing but Pharmit
	var globalsBefore = {};
    for (var key in window)
         globalsBefore[key] = true;

    $( document ).tooltip();
    Pharmit.inFormSubmit = false; // unfortunately there doesn't seem to be a better way to distinguish beforeunload events due to forms
	Pharmit.checkGlobals = function() {
	    var leaked = [];
        for (var key in window)
            if (!(key in globalsBefore))
                leaked.push(key);
        if (leaked.length > 0)
            console.log('Leaked global variables: [' + leaked.join(', ') + ']');
    };
	
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

	
	var element = $('#pharmit').addClass('pharmit_main');
	var viewer = new Pharmit.Viewer(element);
	var results = new Pharmit.Results(element, viewer);
	var query = new Pharmit.Query(element, viewer, results);
	
	//look for session in url
	if(gup('SESSION'))
	{		
		$.get(decodeURI(gup('SESSION')), function(ret) {
			query.loadSession(ret);
		});
	}

		
	Pharmit.checkGlobals();
	//work around jquery bug
	$("button, input[type='button'], input[type='submit']").button()
    .bind('mouseup', function() {
        $(this).blur();     // prevent jquery ui button from remaining in the active state
    });
	
	//message passing, can send a ligand and/or receptor to pharmit
	//to have it loaded, but 
	
	var receiveMessage = function(event) {
		//console.log("receivemsg "+event.data);
		if(event.data == "ack") { //acks let us verify that we are listening
			event.source.postMessage("ack2","*");
		}
		else if(event.data == "ack2") {
			//ignore, message should handle these
		}
		else {
			try {
				var obj = $.parseJSON(event.data);
				query.setLigandAndReceptor(obj.ligand, obj.ligandFormat, obj.receptor, obj.recname);
			}
			catch(e) {
				alert("Communication error: "+e);
			}
		}
	};
	window.addEventListener("message", receiveMessage);

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
	minresults.js
	This is a div for managing minimization results.
	
*/


var Pharmit = Pharmit || {};

Pharmit.MinResults = (function() {
	// private class variables and functions
	

	function MinResults(results, viewer) {
		//private variables and functions
		var mindiv = null;
		var qid = null;
		var save = null;
		var onclose = null; //what to call when close button is clicked
		var table =null; //results table
		var maxscore = null; //filters
		var maxRMSD = null; 
		var singleConfs = null;
		var startTotal = 0;
		var timeout = null;
		var query = null;
		
		var processData = function(data) {
			
			if(data.status === 0)
				return {error: data.msg};
			
			var ret = data.data;
			
			for(var i = 0; i < ret.length; i++) {
				//round floats
				ret[i][3] = numeral(ret[i][3]).format('0.00');
				ret[i][4] = numeral(ret[i][4]).format('0.000');
				ret[i][2] = results.mangleName(ret[i][2]);

			}
			return ret;
		};
		
		//construct objects for fetching smina data
		var ajaxData = function(postData, settings) {
			
			postData.cmd = "getsminadata";
			postData.qid = qid;
	    	
            //add filters if present
            if(maxRMSD.val() !== '') {
            	postData.maxRMSD = maxRMSD.val();
            }
            if(maxscore.val() !== '') {
            	postData.maxScore = maxscore.val();
            }
            if(singleConfs.prop('checked')) {
            	postData.unique = 1;
            }
            return postData;
		};
		
		
		//public variables and functions
		var cancel = this.cancel = function() {
			
			clearTimeout(timeout);
			if($.fn.DataTable.isDataTable(table) && qid > 0) {
				table.DataTable().clear();
				$.post(Pharmit.server, {
					cmd: 'cancelsmina',
					qid: qid
				});				
			}			
			qid = null;
			viewer.setResult();
			save.button( "option", "disabled", true );
			mindiv.hide();
		};
		
		//perform the query
		this.minimize = function(q, qobj, nummols, closer) {
			startTotal = nummols;
			onclose = closer;
			qid = q;
			query = qobj;
			var postData = {cmd: 'startsmina',
					qid: qid,
					receptorid: qobj.receptorid,
					recname: qobj.recname
			};
			
			$('.pharmit_mincontainer .pharmit_resbody').css({opacity: 0}); //don't show table and associated crust until it is setup
			//start provided query
			$.post(Pharmit.server, postData, null, 'json').done(function(ret) {
				if(ret.status) { //success
					//setup table
					sminaid = ret.sminaid;
					var numrows = Math.floor((body.height()-100)/29); //magic numbers!
					
					table.dataTable({
						searching: false,
						pageLength: numrows,
						destroy: true, //replace any existing table
						lengthChange: false,
						order: [[ 3, "asc" ]],
						orderMulti: false,
						columnDefs: [
						                {
						                    targets: [ 0 ], //position
						                    visible: false
						                },
						                {
						                    targets: [ 1 ], //orig position
						                    visible: false
						                },
						                {
						                    targets: [ 2 ], //name
						                    className: "pharmit_minname",
						                    searchable: false,
						                    sortable: false
						                },
						                {
						                    targets: [ 3 ], //score
						                    className: "pharmit_minscore",
						                    searchable: false
						                },
						                {
						                	 targets: [ 4 ], //rmsd
							                 className: "pharmit_minrmsd",
							                 searchable: false
						                }
						            ],
						 language: {
							 	emptyTable: "Minimizing...",
							 	infoEmpty: " ",
							 	infoFiltered: "<br>(filtered from _MAX_ total hits)",
							 	info: "Minimizing..."
						 },
						 serverSide: true,
						 processing: false,
						 ajax: {
						    	url: Pharmit.server,
						    	data: ajaxData,
						    	dataSrc: processData
						 }

					});
					table.show();
					$('.pharmit_mincontainer .pharmit_resbody').css({opacity: 100}); 
				} else {
					cancel();
					if(onclose) onclose();
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				cancel();
				if(onclose) onclose();
				alert("Error contacting minimization server.  Please inform "+Pharmit.email+ " if this problem persists.");
			});	
			
			//show div
			mindiv.show();
		};
		
		//download results
		var saveResults = function() {
			
			Pharmit.inFormSubmit = true;			
			//have to use stupid form trick - mostly because of IE and safari
			var cmd = Pharmit.server+'?cmd=savesmina&qid='+qid;
			
            //add filters if present
            if(maxRMSD.val() !== '') {
            	cmd += "&maxRMSD=" + maxRMSD.val();
            }
            if(maxscore.val() !== '') {
            	cmd += "&maxScore=" + maxscore.val();
            }
            if(singleConfs.prop('checked')) {
            	cmd += "&unique=1";
            }
            
            var order = table.DataTable().order();
            cmd += "&order[0][column]="+order[0][0];
            cmd += "&order[0][dir]="+order[0][1];    		
			
			var form = $('<form>', { 'action': cmd, 'method': 'post'});
			form.appendTo(document.body);
			form.submit();
			$(form).remove();		
		};
		
		
		//initialization code
		mindiv = $('<div>').appendTo(results.div).addClass('pharmit_rescontainer pharmit_mincontainer');
		//header
		var header = $('<div>').appendTo(mindiv).addClass("pharmit_resheader");
		var title = $('<div>Minimization Results</div>').appendTo(header).addClass('pharmit_heading').addClass("pharmit_rightheading");
		var closediv = $('<div>').addClass('pharmit_resclose').appendTo(title).click(function() {
			//cancel the current query 
			cancel();
			//do what our caller told us to do on close
			if(onclose) onclose();
		});
		var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
		
		
		//body, should stretch to fill
		var body = $('<div>').appendTo(mindiv).addClass("pharmit_resbody");

		//skeleton of datatable
		table = $('<table width="100%" class="display compact" cellspacing="0">').addClass('pharmit_mintable').appendTo(body);
		var headrow = $('<tr>').appendTo($('<thead>').appendTo(table));
		$('<th>mid</th>').appendTo(headrow);
		$('<th>origpos</th>').appendTo(headrow);
		$('<th>Name</th>').appendTo(headrow);
		$('<th title="Computed using the AutoDock Vina scoring function">Score</th>').appendTo(headrow);
		$('<th title="The minimized RMSD, the heavy atom difference between the pharmacophore aligned pose and the minimized pose.  Small values imply a closer agreement to the query pharmacophore.">mRMSD</th>').appendTo(headrow);
		$('<tbody>').appendTo(table);
		
		//event handler for loading data, only install this once
		table.on('xhr.dt', function(e, settings, json) {
			var lang = table.DataTable().settings()[0].oLanguage;

			if(json.finished) {
				save.button( "option", "disabled", false );			
				lang.sInfo = "Showing _START_ to _END_ of _TOTAL_ entries";
			} 
	        else if(json.status === 0) {
	        	alert(json.msg);
	        }
			else {
	            viewer.setResult(); //clear in case clicked on

	            lang.sInfo = numeral(json.recordsTotal).format('0,0') + "/" + numeral(startTotal).format('0,0');
	            
	            clearTimeout(timeout); //push next poll farther into the future
	            timeout = setTimeout(function() {
					if($.fn.DataTable.isDataTable(table) && qid > 0) {
						table.DataTable().ajax.reload();
					}
				}, 1000);
			}					 
		});	
		
		table.on('draw.dt', function() {
			$('.pharmit_namecol span').powerTip({mouseOnToPopup:true,placement:'s',smartPlacement:true});
		});
		
		$('tbody',table).on( 'click', 'tr', function () {
			var mid = table.DataTable().row(this).data()[0];
			var r = this;
			$(".pharmit_iterate_button").remove();

	        if ( $(this).hasClass('selected') ) {
	            $(this).removeClass('selected');
	            viewer.setResult(); //clear
	        }
	        else {
	            table.DataTable().$('tr.selected').removeClass('selected');
	            $(this).addClass('selected');
	            
	            $.post(Pharmit.server, {cmd: 'getsminamol',
            		 qid: qid,
            		 molid: mid
            		}).done(function(ret) {
	            			if( $(r).hasClass('selected')) //still selected
	            				viewer.setResult(ret);
	            				var ibutton = $('<div class="pharmit_iterate_button" title="Start new pharmit session around selected ligand">').appendTo($('td',r).last());
	            				ibutton.button({ icons: {primary: "ui-icon-arrowthickstop-1-e"}, text: false});
								ibutton.tooltip({show: {delay: 500}});
	            				ibutton.click(function(event) {
	            					event.stopPropagation();
	            					//create new window around this molecule
	            					var win = window.open("search.html");
	            					var data = {ligand: ret, ligandFormat: mid+".sdf", receptor: query.receptor, recname: query.recname};
	            					var msg = new Message(JSON.stringify(data), win, '*');
	            				});	            			
	            		});
	        }
	    });
		
		//footer
		var footer = $('<div>').appendTo(mindiv).addClass("pharmit_resfooter");
		var paramdiv = $('<div>').appendTo(footer).addClass("pharmit_minparams").disableSelection();
		//filters for minimization
		var filtertable = $('<table>').appendTo(paramdiv);
		var row = $('<tr>').appendTo(filtertable).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max Score</label>'));
		var cell = $('<td>').appendTo(row);
		maxscore = $('<input name="sminamaxscore">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		singleConfs = $('<input type="checkbox" name="sminaunique">');
		$('<td>').appendTo(row).addClass('pharmit_checkcell').append(singleConfs);
		$('<td>').appendTo(row).append($('<label>Single<br>conformer</label>'));
		
		
		row = $('<tr>').appendTo(filtertable).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max mRMSD</label>'));
		cell = $('<td>').appendTo(row);
		maxRMSD = $('<input name="sminamaxrmsd">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		cell = $('<td colspan=2>').appendTo(row).addClass('pharmit_applycell');
		$('<button>Apply</button>').appendTo(cell).button().click(function() {table.DataTable().ajax.reload();});
		
		//save button
		var bottomrow = $('<div>').appendTo(footer).addClass('pharmit_minbottom');
		save = $('<button>Save...</button>').appendTo(bottomrow).button().click(saveResults);				

		mindiv.hide();
		
		var lastheight = 0;
		var lastnum = 0;
		$(window).resize(function() {
			if(qid > 0) {
				var total = body.height();
				if(total != lastheight) {
					lastheight = total;
					total -= $('thead', table).height();
					total -= $('.dataTables_info', body).height();
					total -= $('.dataTables_paginate', body).height();
					total -= 20; //padding
					var single = $('tr.odd',table).first().height()+1; //include border
					var num = Math.floor(total/single);
					if(num != lastnum) { //really only do draw calls when needed
						table.DataTable().page.len(num);
						table.DataTable().draw();
						lastnum = num;
					}
				}
			}
		});
	}

	return MinResults;
})();

//object for sending messages to a window, but only after we receive an ack
function Message(data, w, dest) {
	var curWindow = w;
	var curDest = dest;
	var curMsg = data;
	var isAcked = 0;
	
	function receiveMessage(event)
	{
		if(event.data == "ack2")
		{
			isAcked = 1;
		}
	}
	
	function check() {
		if(isAcked) {
			curWindow.postMessage(curMsg,curDest);
			curDest ="";
			curMsg = "";
			curWindow = null;
			isAcked = 0;
			window.removeEventListener("message", receiveMessage);
		}
		else if(curWindow) {
			curWindow.postMessage("ack", curDest);
			setTimeout(check, 250);
		}
	}
	
	window.addEventListener("message", receiveMessage);
	w.postMessage("ack",dest);		
	setTimeout(check, 250);
}


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
	phresults.js
	This is a div for managing pharmacophore results.
	
*/


var Pharmit = Pharmit || {};

Pharmit.PhResults = (function() {
	// private class variables and functions

	function PhResults(r, viewer, minresults) {
		//private variables and functions
		var phdiv = null;
		var qid = null;
		var query = null;
		var table = null;
		var body = null;
		var minimize = null;
		var save = null;
		var timeout = null;
		var results = r;
		var receptor = null;
		
		//format that provided data (mangle the names appropriately)		
		var processData = function(data) {
			
			if(data.status === 0)
				return {error: data.msg};
			
			var ret = data.data;
			
			for(var i = 0; i < ret.length; i++) {
				//round rmsd
				ret[i][1] = numeral(ret[i][1]).format('0.000');
				ret[i][0] = results.mangleName(ret[i][0]);
			}
			
			return ret;
		};
		
		
		var lastheight = 0;
		var lastnum = 0;
		var resize = function() {
			if(qid > 0) {
				var total = body.height();
				if(total != lastheight) {
					lastheight = total;
					total -= $('thead', table).height();
					total -= $('.dataTables_info', body).height();
					total -= $('.dataTables_paginate', body).height();
					total -= 24; //padding
					var single = $('tr.odd',table).first().height()+1; //include border
					var num = Math.floor(total/single);
					if(num != lastnum) { //really only do draw calls when needed
						table.DataTable().page.len(num);
						table.DataTable().draw();
						lastnum = num;
					}
				}
			}
		};
		
		//public variables and functions				
		$.fn.DataTable.ext.pager.numbers_length = 5;
		//perform the query
		this.query = function(qobj) {
			query = $.extend({}, qobj);
			//don't need receptor or ligand structures and they are big
			receptor = query.receptor; //save for iteration
			delete query.receptor;
			delete query.ligand;
			//start provided query
			var postData = {cmd: 'startquery',
					json: JSON.stringify(query)
			};
			
			if(qid !== null) postData.oldqid = qid;
			
			$.post(Pharmit.server, postData, null, 'json').done(function(ret) {
				if(ret.status) { //success
					
					//setup table
					qid = ret.qid;
					var numrows = Math.floor((body.height()-120)/28); //magic numbers!
					table.dataTable({
						searching: false,
						pageLength: numrows,
						destroy: true, //replace any existing table
						lengthChange: false,
						order: [[ 1, "asc" ]],
						orderMulti: false,
						columnDefs: [
						                {
						                    targets: [ 0 ], //name
						                    className: "pharmit_namecol",
						                    searchable: false,
						                    sortable: false
						                },
						                {
						                    targets: [ 1 ], //rmsd
						                    className: "pharmit_rmsdcol",
						                    searchable: false
						                },
						                {
						                    targets: [ 2 ], //mass
						                    className: "pharmit_masscol",
						                    searchable: false
						                },
						                {
						                    targets: [ 3 ], //bonds
						                    className: "pharmit_bondscol",
						                    searchable: false
						                },
						                {
						                    targets: [ 4 ], //id
						                    visible: false
						                }
						            ],
						 language: {
							 emptyTable: "Searching "+numeral(ret.numMols).format('0,0')+
							 	" molecules and "+numeral(ret.numConfs).format('0,0')+" conformers...",
							 	infoFiltered: '',
							 	infoEmpty: "",
							 	info: "<span class='pharmit_pulse'>Searching...</span>"
						 },
						 serverSide: true,
						 processing: false,
						 ajax: {
						    	url: Pharmit.server,
						    	data: {
						    		cmd: "getdata",
						    		qid: qid
						    	},
						    	dataSrc: processData
						 }

					});												
					
				} else {
					cancel();
					results.close();
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				cancel();
				results.close();
				alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
			});

			phdiv.show(); //make sure we're showing
		};
		
		//cancel any query. clear out the table, and hide the div
		var cancel = this.cancel = function() {
			
			minresults.cancel();
			clearTimeout(timeout); //stop polling
			if(qid !== null) {
				$.post(Pharmit.server, 
						{cmd: 'cancelquery',
						oldqid: qid});
			}
			
			if($.fn.DataTable.isDataTable(table)) {
				table.DataTable().clear();
			}
			qid = null;
			minimize.button( "option", "disabled", true );
			save.button( "option", "disabled", true );
			
		};
		
		//download and save results
		var saveResults = function() {
			//have to use stupid form trick - mostly because of IE and safari
			var cmd = Pharmit.server+'?cmd=saveres&qid='+qid;
			var form = $('<form>', { 'action': cmd, 'method': 'post'});
			form.appendTo(document.body);
			Pharmit.inFormSubmit = true;			
			form.submit();
			$(form).remove();		
		};
		
		//initiate minimization
		var minimizeResults = function() {
			//hide us, show minresults
			phdiv.hide();
			var cnt = table.DataTable().ajax.json().recordsTotal;
			
			var qobj = $.extend({}, query);
			qobj.receptor = receptor;
			minresults.minimize(qid, qobj, cnt, function() {
				phdiv.show();				
			});
		};
		
		//initialization code
		phdiv = $('<div>').appendTo(results.div).addClass('pharmit_rescontainer');
		//header
		var header = $('<div>').appendTo(phdiv).addClass("pharmit_resheader");
		var heading = $('<div>Pharmacophore Results</div>').appendTo(header).addClass('pharmit_heading').addClass("pharmit_rightheading");
		var closediv = $('<div>').addClass('pharmit_resclose').appendTo(heading).click(function() {
			//cancel the current query 
			cancel();
			//close our parent
			results.close();
		});
		var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
		
		//body, should stretch to fill
		body = $('<div>').appendTo(phdiv).addClass("pharmit_resbody");
		
		//skeleton of datatable
		table = $('<table width="100%" class="display compact" cellspacing="0">').addClass('pharmit_phtable').appendTo(body);
		var headrow = $('<tr>').appendTo($('<thead>').appendTo(table));
		$('<th>Name</th>').appendTo(headrow);
		$('<th>RMSD</th>').appendTo(headrow);
		$('<th>Mass</th>').appendTo(headrow);
		$('<th>RBnds</th>').appendTo(headrow);
		$('<th>mid</th>').appendTo(headrow);
		$('<tbody>').appendTo(table);
		

		//setup event handlers for table - this should be done but once
		table.on('xhr.dt', function(e, settings, json) {
			if(json.finished) {
				var lang = table.DataTable().settings()[0].oLanguage;
				if(json.recordsTotal === 0) {
					lang.emptyTable = lang.sEmptyTable = "No results found";
				} else {
					lang.sInfo = "Showing _START_ to _END_ of _TOTAL_ hits";								
					minimize.button( "option", "disabled", false );
					save.button( "option", "disabled", false );
				}
				
			} 
			else if(json.status === 0) {
				//alert(json.msg);
			}
			else {
	            viewer.setResult(); //clear in case clicked on
				//poll server
	            clearTimeout(timeout); //no more than one at once
				timeout = setTimeout(function() {
					if(qid > 0) {
						table.DataTable().ajax.reload();
					}
				}, 1000);	
			}					 
		});	
		
		table.on('draw.dt', function() {
			$('.pharmit_namecol span').powerTip({mouseOnToPopup:true,placement:'s',smartPlacement:true});
		});
		
		$('tbody',table).on( 'click', 'tr', function () {
			var r = this;
			var mid = table.DataTable().row(r).data()[4];
			$(".pharmit_iterate_button").remove();
        	        if ( $(r).hasClass('selected') ) {
        	            $(r).removeClass('selected');
        	            viewer.setResult(); //clear
        	        }
        	        else {
        	            table.DataTable().$('tr.selected').removeClass('selected');
        	            $(r).addClass('selected');
        	            
        	            $.post(Pharmit.server,
        	            		{cmd: 'getmol',
        	            		 qid: qid,
        	            		 loc: mid
        	            		}).done(function(ret) {
        	            			if( $(r).hasClass('selected')) { //still selected
        	            				viewer.setResult(ret);
        	            				var ibutton = $('<div class="pharmit_iterate_button" title="Start new pharmit session around selected ligand">').appendTo($('td',r).last());
        	            				ibutton.button({ icons: {primary: "ui-icon-arrowthickstop-1-e"}, text: false});					
        						ibutton.tooltip({show: {delay: 500}});
        	            				ibutton.click(function(event) {
        	            					event.stopPropagation();
        	            					//create new window around this molecule
        	            					var win = window.open("search.html");
        	            					var data = {ligand: ret, ligandFormat: mid+".sdf", receptor: receptor, recname: query.recname};
        	            					var msg = new Message(JSON.stringify(data), win, '*');
        	            				});
        	            			}
        	            		});
        	        }
	        });		
		
		//footer
		var footer = $('<div>').appendTo(phdiv).addClass("pharmit_resfooter");
		//minimize and save buttons
		var bottomloaders = $('<div>').appendTo(footer).addClass("pharmit_bottomloaders").addClass('pharmit_nowrap');

		minimize = $('<button>Minimize</button>').appendTo(bottomloaders).button({disabled: true}).click(minimizeResults);
		save = $('<button>Save...</button>').appendTo(bottomloaders).button({disabled: true}).click(saveResults);		
		
		//resize event - set number of rows
		$(window).resize(resize);		

	}

	return PhResults;
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
	query.js
	Left div that manages the query information.
*/

var Pharmit = Pharmit || {};

Pharmit.Query = (function() {
	
	var defaultFeature = {name:"Hydrophobic",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
	var defaultInShapeFeature = {name:"InclusionSphere",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
	var defaultExShapeFeature = {name:"ExclusionSphere",x:0,y:0,z:0,radius:1.0,enabled:true,vector_on:0,minsize:"",maxsize:"",svector:null,hasvec:false};
	var pharmaGistRegEx = /@<TRIPOS>MOLECULE[^@]*?@<TRIPOS>ATOM\n(\s*\d+\s*(ACC|DON|CAT|ANI|HYD|AR).*)*\n@<TRIPOS>BOND\n/g;
	var privatedialog = null;
	var endsWith = function(str, suffix) {
	    return str.indexOf(suffix, str.length - suffix.length) !== -1;
	};
	
	function Query(element, viewer, results) {
		//private variables and functions
		var querydiv = $('<div>').addClass('pharmit_query pharmit_overlay');
		var features = null;
		var inshapefeatures = null;
		var exshapefeatures = null;
		
		var featuregroup = null;
		var featurenone = null;
		var receptorData = null;
		var receptorName = null; //filename (need ext)
		var receptorKey = null; //md5 key to avoid transfering full structure
		var ligandData = null;
		var ligandName = null;
		
		var doSearch = function() {
			var qobj = getQueryObj();
			//remove receptor and ligand data to reduce size of query
			if(qobj.receptor && qobj.reckey) {
				delete qobj.receptor; //server can get data w/reckey if need be
			}
			
			//results manages queries
			results.phquery(qobj);
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
				}
			});
			
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
			
			//radio buttons have to have the same name, so have to hack around background color
			ret.backgroundcolor = 'whiteBackground';
			if($('#blackBackground').prop('checked')) ret.backgroundcolor = 'blackBackground';
			
			//structures
			ret.ligand = ligandData;
			ret.ligandFormat = ligandName;
			ret.receptor = receptorData;
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
			var qobj = getQueryObj();
			form.append($('<input>', {'name':"data",'type':"hidden",value:JSON.stringify(qobj,null,4)}));
			form.appendTo(document.body);
			form.submit();
			$(form).remove();			

		};
	
		//escape special characters to avoid injections
		var escHTML = function(str) { return $('<div/>').text(str).html(); };	
		
		var updateShapeMesh = function(sel) {
			//send query to server to get mesh
			var qobj = getQueryObj();
			var kind = sel.meshstyle.kind;
			if(qobj.receptor && qobj.receptorid) delete qobj.receptor; //don't send receptor if it is cached
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
			
			var instylediv = $('<div>').appendTo(liganddiv).addClass('pharmit_meshstylediv');
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
			
			var exstylediv = $('<div>').appendTo(recdiv).addClass('pharmit_meshstylediv');
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
			var shapemodediv = $('<div>').prependTo(header).addClass('pharmit_shapemodediv');	
			
			var shapemodeid = "ShapeModeSelect";
			
			var shapeselect = $('<select name="'+shapemodeid+'" id="'+shapemodeid+'">').addClass('pharmit_styleselector').appendTo(shapemodediv);
			
			$('<option value="filter">Pharmacophore Search / Shape Filter</option>').appendTo(shapeselect);
			$('<option value="search">Shape Search / Pharmacophore Filter</option>').appendTo(shapeselect);	
			
			shapeselect.val("filter");
			shapeselect.selectmenu({
				width: '11em', 
				appendTo: header, 
				change: function() {
					shapeselect.change();
				},
				position: {my: "left top", at: "left bottom", collision: "flip"}
			});
			
			//workaround firefox bug - remove element style so css stylesheet takes effect
			shapeselect.selectmenu( "widget" ).css('width','');
			
			shapeselect.change(function() {
				if(this.value == 'search') {
					$('.pharmit_shapefiltertext').hide();
					$('.pharmit_shapesearchtext').show();
					featuregroup.hide();
					featurenone.show();
				} else { //filter
					$('.pharmit_shapefiltertext').show();
					$('.pharmit_shapesearchtext').hide();
					featuregroup.show();
					featurenone.hide();
				}
				shapeselect.selectmenu("refresh");	       

			});
			
			shapeselect.change(); //set proper visibility of text after it is created

		};
		
		//create a split button from a list of vendors and prepend it to header
		var createSearchButton = function(header,dbinfo) {
			var subsetinfo = dbinfo.standard;
			var buttons = $('<div>').addClass('pharmit_searchdiv');
			var run = $('<button id="pharmitsearchbutton" name="subset">Search '+subsetinfo[0].name+'</button>').appendTo(buttons).button();
			var select = $('<button>Select subset to search</button>').appendTo(buttons).button({text: false, icons: {primary: "ui-icon-triangle-1-s"}});
			select.tooltip().tooltip("disable"); //conflicts with menu
			run.val(subsetinfo[0].subdir);
			
			buttons.buttonset();
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
				if(info.subdir == 'molport') {
					run.button("option",'label',"Search "+info.name);
					run.val(info.subdir);
				}
			}
			ul.append(lis);
			ul.append($('<li> </li>'));
			var publicli = $('<li class="pharmit_contributed">Contributed Libraries</li>');
			var publicul =  $('<ul>').appendTo(publicli);
			var publiclis = [];
			var publicinfo = dbinfo.public;
			for(i = 0, n = publicinfo.length; i < n; i++) {
				info = publicinfo[i];
				display = escHTML(info.name);
				var titlestr = "";
				if(info.description) titlestr = " title='"+escHTML(info.description)+"' ";
				publiclis[i] = '<li value='+subsetinfo.length+titlestr+' class="pharmit_subsetmenu">'+display+'<br>';
				subsetinfo.push(info);
				publiclis[i] += '<span class="pharmit_subsetcnts">' + numeral(info.numConfs).format('0,0') + ' conformers of ' + numeral(info.numMols).format('0,0') + ' molecules</span>';
				publiclis[i] += '<span class="pharmit_subsettime">Created: '+info.updated+'</span>';
				publiclis[i] += '</li>';
			}
			publicul.append(publiclis);
			ul.append(publicli);
			ul.append($('<li> </li>'));

			$('<li class="pharmit_private">Access Private Library</li>').appendTo(ul).click(
					function() {
						privatedialog.dialog("open");
					});
			
			ul.hide().menu().on('menuselect', function(event, ui) {
				var info = subsetinfo[ui.item.val()];
				run.button("option",'label',"Search "+info.name);
				run.val(info.subdir);
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
		
		$.post(Pharmit.server, {cmd: "getsubsets"}, null, 'json').done(function(ret) {
			createSearchButton(header, ret);
		}).fail(function() {
			createSearchButton(header, {standard: [{name: "Error"}], public:[]});
			alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
		});
		
		
		//load features and load receptor
		var loaders = $('<div>').appendTo(header).addClass('pharmit_loaderdiv').addClass('pharmit_nowrap');
		var loadrec = $('<button>Load Receptor...</button>').button();
		var titletext = "Pharmacophore features can be provided as pharmacophore query formats (MOE, LigBuilder, LigandScout, Pharmer) or automatically extracted from ligand structures (sdf, pdf, mol2, xyz). If a receptor is loaded, only interacting features will be enabled.";
		var loadfeatures = $('<button title="'+titletext+'">Load Features...</button>').button();
		
		//fileinput needs the file inputs in the dom
		element.append(querydiv);
		var loadrecfile = $('<input type="file">').appendTo(loaders).fileinput(loadrec).change(function(e) {readText(this, loadReceptor);});
		var loadfeaturesfile = $('<input type="file">').appendTo(loaders).fileinput(loadfeatures).change(function(e) {readText(this,loadFeatures);});		
		
		querydiv.detach();
		
		//query features
		var body = $('<div>').appendTo(querydiv).addClass("pharmit_querybody");
		var featureheading = $('<div>Pharmacophore</div>').appendTo(body).addClass('pharmit_heading');
		featuregroup = $('<div>').appendTo(body);
		featurenone = $('<div>Pharmacophore filtering of aligned shape hits currently not implemented.</div>').addClass('pharmit_disabledmessage').appendTo(body);
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
	results.js
	Object responsible for maintaining results from pharmacophore search
	and energy minimization.  Placed on right side of screen.
	Completely hidden if no active query.  Can be collapsed.
	Closing the div is equivalent to canceling the query (make this event is fired beforeunload)
	Pharmacophore results and minimization results are seperate divs.
	
*/

var Pharmit = Pharmit || {};

Pharmit.Results = (function() {
	// private class variables and functions
	

	function Results(element, viewer) {
		//private variables and functions
		var resultsdiv = this.div = $('<div>').addClass('pharmit_results pharmit_overlay').appendTo(element);
		var phresults = null;
		var minresults = null;
		
		
		//public variables and functions
		
		//perform the query
		this.phquery = function(qobj) {
			// cancel current query first
			phresults.cancel();
			//start provided query
			phresults.query(qobj);						
			//show div
			this.show();
		};
		
		
		//show panel, updating viwer
		this.show = function() {
			resultsdiv.show();
			viewer.setRight(resultsdiv.width());
		};
		
		//completely hide panel
		this.close = function() {
			resultsdiv.hide();
			viewer.setRight(0);
		};
		
		//if a compound name can be mapped to a url, return it
		var getNameURL = function(name) {
		    var m = null;
		    if((m = name.match(/PubChem-?(\d+)/))) {
		        return "https://pubchem.ncbi.nlm.nih.gov/compound/"+m[1];
		    }
		    else if((m = name.match(/MolPort-?(\S+)/))) {
		        return "https://www.molport.com/shop/moleculelink/about-this-molecule/"+m[1].replace(/-/g,"");
		    }
		    else if((m = name.match(/CHEMBL/))) {
		        return "https://www.ebi.ac.uk/chembldb/index.php/compound/inspect/"+name;
		    }
		    else if((m = name.match(/ChemDiv-?(\S+)/))) {
		        return "http://chemistryondemand.com:8080/eShop/search_results.jsp?idnumber="+m[1];
		    }
		    
		    return null;
		};
		
		//convert a mol name into something more presentable
		this.mangleName = function(name) {
			var names = name.split(" ");
			var tip = "";
			for(var i = 0; i < names.length; i++) {
			    var n = names[i];
			    var url = getNameURL(n);
			    if(url) {
			        tip += "<a target='_blank' href='"+url+"'>"+n+"</a><br>";
			    }
			    else {
			        tip += n+"<br>";
			    }
			}
			var ret = '<span data-powertip="' + tip +
						'">'+names[0]+'</span>';
			
			return ret;
		};
		
		//initialization code
		var closer = $('<div>').appendTo(resultsdiv).addClass('pharmit_rightclose');
		var closericon = $('<span>').addClass("ui-icon ui-icon-carat-1-e").appendTo(closer);
		
		//initialization code
		resultsdiv.resizable({handles: "w",
			resize: function(event, ui) {
				viewer.setRight(ui.size.width);
			    $(this).css("left", ''); //workaround for chrome/jquery bug
			}
		});
		resultsdiv.disableSelection();
		

		closer.click(function() { //todo, refactor w/query
			if(closer.hasClass('pharmit_rightisclosed')) {
				closer.removeClass('pharmit_rightisclosed');
				closericon.removeClass('ui-icon-carat-1-w');
				closericon.addClass('ui-icon-carat-1-e');
				var start = resultsdiv.width();
				resultsdiv.css('width', ''); //restore stylesheet width	
				var target = resultsdiv.width();
				resultsdiv.width(start);
				
				resultsdiv.animate({width: target},{
					progress: function() { viewer.setRight(resultsdiv.width());}
				}); 
				resultsdiv.resizable( "option", "disabled", false);

			} else { //close it 
				resultsdiv.animate({width: 0}, {
					progress: function() { viewer.setRight(resultsdiv.width());}
					}); 
				//viewer.setLeft(0);
				closer.addClass('pharmit_rightisclosed');
				closericon.addClass('ui-icon-carat-1-w');
				closericon.removeClass('ui-icon-carat-1-e');			
				resultsdiv.resizable( "option", "disabled", true );
			}
		});
		
		
		
		//minimization results
		minresults = new Pharmit.MinResults(this, viewer);
		
		//pharmacophore results
		phresults = new Pharmit.PhResults(this, viewer, minresults);		

		resultsdiv.hide(); //wait for query
		if(resultsdiv.is(":visible")) {
			viewer.setRight(resultsdiv.width());
		}
		
		//be nice and cancel queries when finishing
		$(window).on('beforeunload', function(){
			if(!Pharmit.inFormSubmit) {
				phresults.cancel();
			}
			else {
				Pharmit.inFormSubmit = false; //no longer
			}
		});
	}

	return Results;
})();
// Spectrum Colorpicker v1.6.0
// https://github.com/bgrins/spectrum
// Author: Brian Grinstead
// License: MIT

(function (factory) {
    "use strict";

    if (typeof define === 'function' && define.amd) { // AMD
        define(['jquery'], factory);
    }
    else if (typeof exports == "object" && typeof module == "object") { // CommonJS
        module.exports = factory;
    }
    else { // Browser
        factory(jQuery);
    }
})(function($, undefined) {
    "use strict";

    var defaultOpts = {

        // Callbacks
        beforeShow: noop,
        move: noop,
        change: noop,
        show: noop,
        hide: noop,

        // Options
        color: false,
        flat: false,
        showInput: false,
        allowEmpty: false,
        showButtons: true,
        clickoutFiresChange: false,
        showInitial: false,
        showPalette: false,
        showPaletteOnly: false,
        hideAfterPaletteSelect: false,
        togglePaletteOnly: false,
        showSelectionPalette: true,
        localStorageKey: false,
        appendTo: "body",
        maxSelectionSize: 7,
        cancelText: "cancel",
        chooseText: "choose",
        togglePaletteMoreText: "more",
        togglePaletteLessText: "less",
        clearText: "Clear Color Selection",
        noColorSelectedText: "No Color Selected",
        preferredFormat: false,
        className: "", // Deprecated - use containerClassName and replacerClassName instead.
        containerClassName: "",
        replacerClassName: "",
        showAlpha: false,
        theme: "sp-light",
        palette: [["#ffffff", "#000000", "#ff0000", "#ff8000", "#ffff00", "#008000", "#0000ff", "#4b0082", "#9400d3"]],
        selectionPalette: [],
        disabled: false,
        offset: null
    },
    spectrums = [],
    IE = !!/msie/i.exec( window.navigator.userAgent ),
    rgbaSupport = (function() {
        function contains( str, substr ) {
            return !!~('' + str).indexOf(substr);
        }

        var elem = document.createElement('div');
        var style = elem.style;
        style.cssText = 'background-color:rgba(0,0,0,.5)';
        return contains(style.backgroundColor, 'rgba') || contains(style.backgroundColor, 'hsla');
    })(),
    inputTypeColorSupport = (function() {
        var colorInput = $("<input type='color' value='!' />")[0];
        return colorInput.type === "color" && colorInput.value !== "!";
    })(),
    replaceInput = [
        "<div class='sp-replacer'>",
            "<div class='sp-preview'><div class='sp-preview-inner'></div></div>",
            "<div class='sp-dd'>&#9660;</div>",
        "</div>"
    ].join(''),
    markup = (function () {

        // IE does not support gradients with multiple stops, so we need to simulate
        //  that for the rainbow slider with 8 divs that each have a single gradient
        var gradientFix = "";
        if (IE) {
            for (var i = 1; i <= 6; i++) {
                gradientFix += "<div class='sp-" + i + "'></div>";
            }
        }

        return [
            "<div class='sp-container sp-hidden'>",
                "<div class='sp-palette-container'>",
                    "<div class='sp-palette sp-thumb sp-cf'></div>",
                    "<div class='sp-palette-button-container sp-cf'>",
                        "<button type='button' class='sp-palette-toggle'></button>",
                    "</div>",
                "</div>",
                "<div class='sp-picker-container'>",
                    "<div class='sp-top sp-cf'>",
                        "<div class='sp-fill'></div>",
                        "<div class='sp-top-inner'>",
                            "<div class='sp-color'>",
                                "<div class='sp-sat'>",
                                    "<div class='sp-val'>",
                                        "<div class='sp-dragger'></div>",
                                    "</div>",
                                "</div>",
                            "</div>",
                            "<div class='sp-clear sp-clear-display'>",
                            "</div>",
                            "<div class='sp-hue'>",
                                "<div class='sp-slider'></div>",
                                gradientFix,
                            "</div>",
                        "</div>",
                        "<div class='sp-alpha'><div class='sp-alpha-inner'><div class='sp-alpha-handle'></div></div></div>",
                    "</div>",
                    "<div class='sp-input-container sp-cf'>",
                        "<input class='sp-input' type='text' spellcheck='false'  />",
                    "</div>",
                    "<div class='sp-initial sp-thumb sp-cf'></div>",
                    "<div class='sp-button-container sp-cf'>",
                        "<a class='sp-cancel' href='#'></a>",
                        "<button type='button' class='sp-choose'></button>",
                    "</div>",
                "</div>",
            "</div>"
        ].join("");
    })();

    function paletteTemplate (p, color, className, opts) {
        var html = [];
        for (var i = 0; i < p.length; i++) {
            var current = p[i];
            if(current) {
                var tiny = tinycolor(current);
                var c = tiny.toHsl().l < 0.5 ? "sp-thumb-el sp-thumb-dark" : "sp-thumb-el sp-thumb-light";
                c += (tinycolor.equals(color, current)) ? " sp-thumb-active" : "";
                var formattedString = tiny.toString(opts.preferredFormat || "rgb");
                var swatchStyle = rgbaSupport ? ("background-color:" + tiny.toRgbString()) : "filter:" + tiny.toFilter();
                html.push('<span title="' + formattedString + '" data-color="' + tiny.toRgbString() + '" class="' + c + '"><span class="sp-thumb-inner" style="' + swatchStyle + ';" /></span>');
            } else {
                var cls = 'sp-clear-display';
                html.push($('<div />')
                    .append($('<span data-color="" style="background-color:transparent;" class="' + cls + '"></span>')
                        .attr('title', opts.noColorSelectedText)
                    )
                    .html()
                );
            }
        }
        return "<div class='sp-cf " + className + "'>" + html.join('') + "</div>";
    }

    function hideAll() {
        for (var i = 0; i < spectrums.length; i++) {
            if (spectrums[i]) {
                spectrums[i].hide();
            }
        }
    }

    function instanceOptions(o, callbackContext) {
        var opts = $.extend({}, defaultOpts, o);
        opts.callbacks = {
            'move': bind(opts.move, callbackContext),
            'change': bind(opts.change, callbackContext),
            'show': bind(opts.show, callbackContext),
            'hide': bind(opts.hide, callbackContext),
            'beforeShow': bind(opts.beforeShow, callbackContext)
        };

        return opts;
    }

    function spectrum(element, o) {

        var opts = instanceOptions(o, element),
            flat = opts.flat,
            showSelectionPalette = opts.showSelectionPalette,
            localStorageKey = opts.localStorageKey,
            theme = opts.theme,
            callbacks = opts.callbacks,
            resize = throttle(reflow, 10),
            visible = false,
            dragWidth = 0,
            dragHeight = 0,
            dragHelperHeight = 0,
            slideHeight = 0,
            slideWidth = 0,
            alphaWidth = 0,
            alphaSlideHelperWidth = 0,
            slideHelperHeight = 0,
            currentHue = 0,
            currentSaturation = 0,
            currentValue = 0,
            currentAlpha = 1,
            palette = [],
            paletteArray = [],
            paletteLookup = {},
            selectionPalette = opts.selectionPalette.slice(0),
            maxSelectionSize = opts.maxSelectionSize,
            draggingClass = "sp-dragging",
            shiftMovementDirection = null;

        var doc = element.ownerDocument,
            body = doc.body,
            boundElement = $(element),
            disabled = false,
            container = $(markup, doc).addClass(theme),
            pickerContainer = container.find(".sp-picker-container"),
            dragger = container.find(".sp-color"),
            dragHelper = container.find(".sp-dragger"),
            slider = container.find(".sp-hue"),
            slideHelper = container.find(".sp-slider"),
            alphaSliderInner = container.find(".sp-alpha-inner"),
            alphaSlider = container.find(".sp-alpha"),
            alphaSlideHelper = container.find(".sp-alpha-handle"),
            textInput = container.find(".sp-input"),
            paletteContainer = container.find(".sp-palette"),
            initialColorContainer = container.find(".sp-initial"),
            cancelButton = container.find(".sp-cancel"),
            clearButton = container.find(".sp-clear"),
            chooseButton = container.find(".sp-choose"),
            toggleButton = container.find(".sp-palette-toggle"),
            isInput = boundElement.is("input"),
            isInputTypeColor = isInput && inputTypeColorSupport && boundElement.attr("type") === "color",
            shouldReplace = isInput && !flat,
            replacer = (shouldReplace) ? $(replaceInput).addClass(theme).addClass(opts.className).addClass(opts.replacerClassName) : $([]),
            offsetElement = (shouldReplace) ? replacer : boundElement,
            previewElement = replacer.find(".sp-preview-inner"),
            initialColor = opts.color || (isInput && boundElement.val()),
            colorOnShow = false,
            preferredFormat = opts.preferredFormat,
            currentPreferredFormat = preferredFormat,
            clickoutFiresChange = !opts.showButtons || opts.clickoutFiresChange,
            isEmpty = !initialColor,
            allowEmpty = opts.allowEmpty && !isInputTypeColor;

        function applyOptions() {

            if (opts.showPaletteOnly) {
                opts.showPalette = true;
            }

            toggleButton.text(opts.showPaletteOnly ? opts.togglePaletteMoreText : opts.togglePaletteLessText);

            if (opts.palette) {
                palette = opts.palette.slice(0);
                paletteArray = $.isArray(palette[0]) ? palette : [palette];
                paletteLookup = {};
                for (var i = 0; i < paletteArray.length; i++) {
                    for (var j = 0; j < paletteArray[i].length; j++) {
                        var rgb = tinycolor(paletteArray[i][j]).toRgbString();
                        paletteLookup[rgb] = true;
                    }
                }
            }

            container.toggleClass("sp-flat", flat);
            container.toggleClass("sp-input-disabled", !opts.showInput);
            container.toggleClass("sp-alpha-enabled", opts.showAlpha);
            container.toggleClass("sp-clear-enabled", allowEmpty);
            container.toggleClass("sp-buttons-disabled", !opts.showButtons);
            container.toggleClass("sp-palette-buttons-disabled", !opts.togglePaletteOnly);
            container.toggleClass("sp-palette-disabled", !opts.showPalette);
            container.toggleClass("sp-palette-only", opts.showPaletteOnly);
            container.toggleClass("sp-initial-disabled", !opts.showInitial);
            container.addClass(opts.className).addClass(opts.containerClassName);

            reflow();
        }

        function initialize() {

            if (IE) {
                container.find("*:not(input)").attr("unselectable", "on");
            }

            applyOptions();

            if (shouldReplace) {
                boundElement.after(replacer).hide();
            }

            if (!allowEmpty) {
                clearButton.hide();
            }

            if (flat) {
                boundElement.after(container).hide();
            }
            else {

                var appendTo = opts.appendTo === "parent" ? boundElement.parent() : $(opts.appendTo);
                if (appendTo.length !== 1) {
                    appendTo = $("body");
                }

                appendTo.append(container);
            }

            updateSelectionPaletteFromStorage();

            offsetElement.bind("click.spectrum touchstart.spectrum", function (e) {
                if (!disabled) {
                    toggle();
                }

                e.stopPropagation();

                if (!$(e.target).is("input")) {
                    e.preventDefault();
                }
            });

            if(boundElement.is(":disabled") || (opts.disabled === true)) {
                disable();
            }

            // Prevent clicks from bubbling up to document.  This would cause it to be hidden.
            container.click(stopPropagation);

            // Handle user typed input
            textInput.change(setFromTextInput);
            textInput.bind("paste", function () {
                setTimeout(setFromTextInput, 1);
            });
            textInput.keydown(function (e) { if (e.keyCode == 13) { setFromTextInput(); } });

            cancelButton.text(opts.cancelText);
            cancelButton.bind("click.spectrum", function (e) {
                e.stopPropagation();
                e.preventDefault();
                revert();
                hide();
            });

            clearButton.attr("title", opts.clearText);
            clearButton.bind("click.spectrum", function (e) {
                e.stopPropagation();
                e.preventDefault();
                isEmpty = true;
                move();

                if(flat) {
                    //for the flat style, this is a change event
                    updateOriginalInput(true);
                }
            });

            chooseButton.text(opts.chooseText);
            chooseButton.bind("click.spectrum", function (e) {
                e.stopPropagation();
                e.preventDefault();

                if (isValid()) {
                    updateOriginalInput(true);
                    hide();
                }
            });

            toggleButton.text(opts.showPaletteOnly ? opts.togglePaletteMoreText : opts.togglePaletteLessText);
            toggleButton.bind("click.spectrum", function (e) {
                e.stopPropagation();
                e.preventDefault();

                opts.showPaletteOnly = !opts.showPaletteOnly;

                // To make sure the Picker area is drawn on the right, next to the
                // Palette area (and not below the palette), first move the Palette
                // to the left to make space for the picker, plus 5px extra.
                // The 'applyOptions' function puts the whole container back into place
                // and takes care of the button-text and the sp-palette-only CSS class.
                if (!opts.showPaletteOnly && !flat) {
                    container.css('left', '-=' + (pickerContainer.outerWidth(true) + 5));
                }
                applyOptions();
            });

            draggable(alphaSlider, function (dragX, dragY, e) {
                currentAlpha = (dragX / alphaWidth);
                isEmpty = false;
                if (e.shiftKey) {
                    currentAlpha = Math.round(currentAlpha * 10) / 10;
                }

                move();
            }, dragStart, dragStop);

            draggable(slider, function (dragX, dragY) {
                currentHue = parseFloat(dragY / slideHeight);
                isEmpty = false;
                if (!opts.showAlpha) {
                    currentAlpha = 1;
                }
                move();
            }, dragStart, dragStop);

            draggable(dragger, function (dragX, dragY, e) {

                // shift+drag should snap the movement to either the x or y axis.
                if (!e.shiftKey) {
                    shiftMovementDirection = null;
                }
                else if (!shiftMovementDirection) {
                    var oldDragX = currentSaturation * dragWidth;
                    var oldDragY = dragHeight - (currentValue * dragHeight);
                    var furtherFromX = Math.abs(dragX - oldDragX) > Math.abs(dragY - oldDragY);

                    shiftMovementDirection = furtherFromX ? "x" : "y";
                }

                var setSaturation = !shiftMovementDirection || shiftMovementDirection === "x";
                var setValue = !shiftMovementDirection || shiftMovementDirection === "y";

                if (setSaturation) {
                    currentSaturation = parseFloat(dragX / dragWidth);
                }
                if (setValue) {
                    currentValue = parseFloat((dragHeight - dragY) / dragHeight);
                }

                isEmpty = false;
                if (!opts.showAlpha) {
                    currentAlpha = 1;
                }

                move();

            }, dragStart, dragStop);

            if (!!initialColor) {
                set(initialColor);

                // In case color was black - update the preview UI and set the format
                // since the set function will not run (default color is black).
                updateUI();
                currentPreferredFormat = preferredFormat || tinycolor(initialColor).format;

                addColorToSelectionPalette(initialColor);
            }
            else {
                updateUI();
            }

            if (flat) {
                show();
            }

            function paletteElementClick(e) {
                if (e.data && e.data.ignore) {
                    set($(e.target).closest(".sp-thumb-el").data("color"));
                    move();
                }
                else {
                    set($(e.target).closest(".sp-thumb-el").data("color"));
                    move();
                    updateOriginalInput(true);
                    if (opts.hideAfterPaletteSelect) {
                      hide();
                    }
                }

                return false;
            }

            var paletteEvent = IE ? "mousedown.spectrum" : "click.spectrum touchstart.spectrum";
            paletteContainer.delegate(".sp-thumb-el", paletteEvent, paletteElementClick);
            initialColorContainer.delegate(".sp-thumb-el:nth-child(1)", paletteEvent, { ignore: true }, paletteElementClick);
        }

        function updateSelectionPaletteFromStorage() {

            if (localStorageKey && window.localStorage) {

                // Migrate old palettes over to new format.  May want to remove this eventually.
                try {
                    var oldPalette = window.localStorage[localStorageKey].split(",#");
                    if (oldPalette.length > 1) {
                        delete window.localStorage[localStorageKey];
                        $.each(oldPalette, function(i, c) {
                             addColorToSelectionPalette(c);
                        });
                    }
                }
                catch(e) { }

                try {
                    selectionPalette = window.localStorage[localStorageKey].split(";");
                }
                catch (e) { }
            }
        }

        function addColorToSelectionPalette(color) {
            if (showSelectionPalette) {
                var rgb = tinycolor(color).toRgbString();
                if (!paletteLookup[rgb] && $.inArray(rgb, selectionPalette) === -1) {
                    selectionPalette.push(rgb);
                    while(selectionPalette.length > maxSelectionSize) {
                        selectionPalette.shift();
                    }
                }

                if (localStorageKey && window.localStorage) {
                    try {
                        window.localStorage[localStorageKey] = selectionPalette.join(";");
                    }
                    catch(e) { }
                }
            }
        }

        function getUniqueSelectionPalette() {
            var unique = [];
            if (opts.showPalette) {
                for (var i = 0; i < selectionPalette.length; i++) {
                    var rgb = tinycolor(selectionPalette[i]).toRgbString();

                    if (!paletteLookup[rgb]) {
                        unique.push(selectionPalette[i]);
                    }
                }
            }

            return unique.reverse().slice(0, opts.maxSelectionSize);
        }

        function drawPalette() {

            var currentColor = get();

            var html = $.map(paletteArray, function (palette, i) {
                return paletteTemplate(palette, currentColor, "sp-palette-row sp-palette-row-" + i, opts);
            });

            updateSelectionPaletteFromStorage();

            if (selectionPalette) {
                html.push(paletteTemplate(getUniqueSelectionPalette(), currentColor, "sp-palette-row sp-palette-row-selection", opts));
            }

            paletteContainer.html(html.join(""));
        }

        function drawInitial() {
            if (opts.showInitial) {
                var initial = colorOnShow;
                var current = get();
                initialColorContainer.html(paletteTemplate([initial, current], current, "sp-palette-row-initial", opts));
            }
        }

        function dragStart() {
            if (dragHeight <= 0 || dragWidth <= 0 || slideHeight <= 0) {
                reflow();
            }
            container.addClass(draggingClass);
            shiftMovementDirection = null;
            boundElement.trigger('dragstart.spectrum', [ get() ]);
        }

        function dragStop() {
            container.removeClass(draggingClass);
            boundElement.trigger('dragstop.spectrum', [ get() ]);
        }

        function setFromTextInput() {

            var value = textInput.val();

            if ((value === null || value === "") && allowEmpty) {
                set(null);
                updateOriginalInput(true);
            }
            else {
                var tiny = tinycolor(value);
                if (tiny.isValid()) {
                    set(tiny);
                    updateOriginalInput(true);
                }
                else {
                    textInput.addClass("sp-validation-error");
                }
            }
        }

        function toggle() {
            if (visible) {
                hide();
            }
            else {
                show();
            }
        }

        function show() {
            var event = $.Event('beforeShow.spectrum');

            if (visible) {
                reflow();
                return;
            }

            boundElement.trigger(event, [ get() ]);

            if (callbacks.beforeShow(get()) === false || event.isDefaultPrevented()) {
                return;
            }

            hideAll();
            visible = true;

            $(doc).bind("click.spectrum", clickout);
            $(window).bind("resize.spectrum", resize);
            replacer.addClass("sp-active");
            container.removeClass("sp-hidden");

            reflow();
            updateUI();

            colorOnShow = get();

            drawInitial();
            callbacks.show(colorOnShow);
            boundElement.trigger('show.spectrum', [ colorOnShow ]);
        }

        function clickout(e) {
            // Return on right click.
            if (e.button == 2) { return; }

            if (clickoutFiresChange) {
                updateOriginalInput(true);
            }
            else {
                revert();
            }
            hide();
        }

        function hide() {
            // Return if hiding is unnecessary
            if (!visible || flat) { return; }
            visible = false;

            $(doc).unbind("click.spectrum", clickout);
            $(window).unbind("resize.spectrum", resize);

            replacer.removeClass("sp-active");
            container.addClass("sp-hidden");

            callbacks.hide(get());
            boundElement.trigger('hide.spectrum', [ get() ]);
        }

        function revert() {
            set(colorOnShow, true);
        }

        function set(color, ignoreFormatChange) {
            if (tinycolor.equals(color, get())) {
                // Update UI just in case a validation error needs
                // to be cleared.
                updateUI();
                return;
            }

            var newColor, newHsv;
            if (!color && allowEmpty) {
                isEmpty = true;
            } else {
                isEmpty = false;
                newColor = tinycolor(color);
                newHsv = newColor.toHsv();

                currentHue = (newHsv.h % 360) / 360;
                currentSaturation = newHsv.s;
                currentValue = newHsv.v;
                currentAlpha = newHsv.a;
            }
            updateUI();

            if (newColor && newColor.isValid() && !ignoreFormatChange) {
                currentPreferredFormat = preferredFormat || newColor.getFormat();
            }
        }

        function get(opts) {
            opts = opts || { };

            if (allowEmpty && isEmpty) {
                return null;
            }

            return tinycolor.fromRatio({
                h: currentHue,
                s: currentSaturation,
                v: currentValue,
                a: Math.round(currentAlpha * 100) / 100
            }, { format: opts.format || currentPreferredFormat });
        }

        function isValid() {
            return !textInput.hasClass("sp-validation-error");
        }

        function move() {
            updateUI();

            callbacks.move(get());
            boundElement.trigger('move.spectrum', [ get() ]);
        }

        function updateUI() {

            textInput.removeClass("sp-validation-error");

            updateHelperLocations();

            // Update dragger background color (gradients take care of saturation and value).
            var flatColor = tinycolor.fromRatio({ h: currentHue, s: 1, v: 1 });
            dragger.css("background-color", flatColor.toHexString());

            // Get a format that alpha will be included in (hex and names ignore alpha)
            var format = currentPreferredFormat;
            if (currentAlpha < 1 && !(currentAlpha === 0 && format === "name")) {
                if (format === "hex" || format === "hex3" || format === "hex6" || format === "name") {
                    format = "rgb";
                }
            }

            var realColor = get({ format: format }),
                displayColor = '';

             //reset background info for preview element
            previewElement.removeClass("sp-clear-display");
            previewElement.css('background-color', 'transparent');

            if (!realColor && allowEmpty) {
                // Update the replaced elements background with icon indicating no color selection
                previewElement.addClass("sp-clear-display");
            }
            else {
                var realHex = realColor.toHexString(),
                    realRgb = realColor.toRgbString();

                // Update the replaced elements background color (with actual selected color)
                if (rgbaSupport || realColor.alpha === 1) {
                    previewElement.css("background-color", realRgb);
                }
                else {
                    previewElement.css("background-color", "transparent");
                    previewElement.css("filter", realColor.toFilter());
                }

                if (opts.showAlpha) {
                    var rgb = realColor.toRgb();
                    rgb.a = 0;
                    var realAlpha = tinycolor(rgb).toRgbString();
                    var gradient = "linear-gradient(left, " + realAlpha + ", " + realHex + ")";

                    if (IE) {
                        alphaSliderInner.css("filter", tinycolor(realAlpha).toFilter({ gradientType: 1 }, realHex));
                    }
                    else {
                        alphaSliderInner.css("background", "-webkit-" + gradient);
                        alphaSliderInner.css("background", "-moz-" + gradient);
                        alphaSliderInner.css("background", "-ms-" + gradient);
                        // Use current syntax gradient on unprefixed property.
                        alphaSliderInner.css("background",
                            "linear-gradient(to right, " + realAlpha + ", " + realHex + ")");
                    }
                }

                displayColor = realColor.toString(format);
            }

            // Update the text entry input as it changes happen
            if (opts.showInput) {
                textInput.val(displayColor);
            }

            if (opts.showPalette) {
                drawPalette();
            }

            drawInitial();
        }

        function updateHelperLocations() {
            var s = currentSaturation;
            var v = currentValue;

            if(allowEmpty && isEmpty) {
                //if selected color is empty, hide the helpers
                alphaSlideHelper.hide();
                slideHelper.hide();
                dragHelper.hide();
            }
            else {
                //make sure helpers are visible
                alphaSlideHelper.show();
                slideHelper.show();
                dragHelper.show();

                // Where to show the little circle in that displays your current selected color
                var dragX = s * dragWidth;
                var dragY = dragHeight - (v * dragHeight);
                dragX = Math.max(
                    -dragHelperHeight,
                    Math.min(dragWidth - dragHelperHeight, dragX - dragHelperHeight)
                );
                dragY = Math.max(
                    -dragHelperHeight,
                    Math.min(dragHeight - dragHelperHeight, dragY - dragHelperHeight)
                );
                dragHelper.css({
                    "top": dragY + "px",
                    "left": dragX + "px"
                });

                var alphaX = currentAlpha * alphaWidth;
                alphaSlideHelper.css({
                    "left": (alphaX - (alphaSlideHelperWidth / 2)) + "px"
                });

                // Where to show the bar that displays your current selected hue
                var slideY = (currentHue) * slideHeight;
                slideHelper.css({
                    "top": (slideY - slideHelperHeight) + "px"
                });
            }
        }

        function updateOriginalInput(fireCallback) {
            var color = get(),
                displayColor = '',
                hasChanged = true; //dkoes - we don't wait for close to use the color so always update

            if (color) {
                displayColor = color.toString(currentPreferredFormat);
                // Update the selection palette with the current color
                addColorToSelectionPalette(color);
            }

            if (isInput) {
                boundElement.val(displayColor);
            }

            if (fireCallback && hasChanged) {
                callbacks.change(color);
                boundElement.trigger('change', [ color ]);
            }
        }

        function reflow() {
            dragWidth = dragger.width();
            dragHeight = dragger.height();
            dragHelperHeight = dragHelper.height();
            slideWidth = slider.width();
            slideHeight = slider.height();
            slideHelperHeight = slideHelper.height();
            alphaWidth = alphaSlider.width();
            alphaSlideHelperWidth = alphaSlideHelper.width();

            if (!flat) {
                container.css("position", "absolute");
                if (opts.offset) {
                    container.offset(opts.offset);
                } else {
                    container.offset(getOffset(container, offsetElement));
                }
            }

            updateHelperLocations();

            if (opts.showPalette) {
                drawPalette();
            }

            boundElement.trigger('reflow.spectrum');
        }

        function destroy() {
            boundElement.show();
            offsetElement.unbind("click.spectrum touchstart.spectrum");
            container.remove();
            replacer.remove();
            spectrums[spect.id] = null;
        }

        function option(optionName, optionValue) {
            if (optionName === undefined) {
                return $.extend({}, opts);
            }
            if (optionValue === undefined) {
                return opts[optionName];
            }

            opts[optionName] = optionValue;
            applyOptions();
        }

        function enable() {
            disabled = false;
            boundElement.attr("disabled", false);
            offsetElement.removeClass("sp-disabled");
        }

        function disable() {
            hide();
            disabled = true;
            boundElement.attr("disabled", true);
            offsetElement.addClass("sp-disabled");
        }

        function setOffset(coord) {
            opts.offset = coord;
            reflow();
        }

        initialize();

        var spect = {
            show: show,
            hide: hide,
            toggle: toggle,
            reflow: reflow,
            option: option,
            enable: enable,
            disable: disable,
            offset: setOffset,
            set: function (c) {
                set(c);
                updateOriginalInput();
            },
            get: get,
            destroy: destroy,
            container: container
        };

        spect.id = spectrums.push(spect) - 1;

        return spect;
    }

    /**
    * checkOffset - get the offset below/above and left/right element depending on screen position
    * Thanks https://github.com/jquery/jquery-ui/blob/master/ui/jquery.ui.datepicker.js
    */
    function getOffset(picker, input) {
        var extraY = 0;
        var dpWidth = picker.outerWidth();
        var dpHeight = picker.outerHeight();
        var inputHeight = input.outerHeight();
        var doc = picker[0].ownerDocument;
        var docElem = doc.documentElement;
        var viewWidth = docElem.clientWidth + $(doc).scrollLeft();
        var viewHeight = docElem.clientHeight + $(doc).scrollTop();
        var offset = input.offset();
        offset.top += inputHeight;

        offset.left -=
            Math.min(offset.left, (offset.left + dpWidth > viewWidth && viewWidth > dpWidth) ?
            Math.abs(offset.left + dpWidth - viewWidth) : 0);

        offset.top -=
            Math.min(offset.top, ((offset.top + dpHeight > viewHeight && viewHeight > dpHeight) ?
            Math.abs(dpHeight + inputHeight - extraY) : extraY));

        return offset;
    }

    /**
    * noop - do nothing
    */
    function noop() {

    }

    /**
    * stopPropagation - makes the code only doing this a little easier to read in line
    */
    function stopPropagation(e) {
        e.stopPropagation();
    }

    /**
    * Create a function bound to a given object
    * Thanks to underscore.js
    */
    function bind(func, obj) {
        var slice = Array.prototype.slice;
        var args = slice.call(arguments, 2);
        return function () {
            return func.apply(obj, args.concat(slice.call(arguments)));
        };
    }

    /**
    * Lightweight drag helper.  Handles containment within the element, so that
    * when dragging, the x is within [0,element.width] and y is within [0,element.height]
    */
    function draggable(element, onmove, onstart, onstop) {
        onmove = onmove || function () { };
        onstart = onstart || function () { };
        onstop = onstop || function () { };
        var doc = document;
        var dragging = false;
        var offset = {};
        var maxHeight = 0;
        var maxWidth = 0;
        var hasTouch = ('ontouchstart' in window);

        var duringDragEvents = {};
        duringDragEvents.selectstart = prevent;
        duringDragEvents.ragstart = prevent;
        duringDragEvents["touchmove mousemove"] = move;
        duringDragEvents["touchend mouseup"] = stop;

        function prevent(e) {
            if (e.stopPropagation) {
                e.stopPropagation();
            }
            if (e.preventDefault) {
                e.preventDefault();
            }
            e.returnValue = false;
        }

        function move(e) {
            if (dragging) {
                // Mouseup happened outside of window
                if (IE && doc.documentMode < 9 && !e.button) {
                    return stop();
                }

                var touches = e.originalEvent && e.originalEvent.touches;
                var pageX = touches ? touches[0].pageX : e.pageX;
                var pageY = touches ? touches[0].pageY : e.pageY;

                var dragX = Math.max(0, Math.min(pageX - offset.left, maxWidth));
                var dragY = Math.max(0, Math.min(pageY - offset.top, maxHeight));

                if (hasTouch) {
                    // Stop scrolling in iOS
                    prevent(e);
                }

                onmove.apply(element, [dragX, dragY, e]);
            }
        }

        function start(e) {
            var rightclick = (e.which) ? (e.which == 3) : (e.button == 2);

            if (!rightclick && !dragging) {
                if (onstart.apply(element, arguments) !== false) {
                    dragging = true;
                    maxHeight = $(element).height();
                    maxWidth = $(element).width();
                    offset = $(element).offset();

                    $(doc).bind(duringDragEvents);
                    $(doc.body).addClass("sp-dragging");

                    if (!hasTouch) {
                        move(e);
                    }

                    prevent(e);
                }
            }
        }

        function stop() {
            if (dragging) {
                $(doc).unbind(duringDragEvents);
                $(doc.body).removeClass("sp-dragging");
                onstop.apply(element, arguments);
            }
            dragging = false;
        }

        $(element).bind("touchstart mousedown", start);
    }

    function throttle(func, wait, debounce) {
        var timeout;
        return function () {
            var context = this, args = arguments;
            var throttler = function () {
                timeout = null;
                func.apply(context, args);
            };
            if (debounce) clearTimeout(timeout);
            if (debounce || !timeout) timeout = setTimeout(throttler, wait);
        };
    }

    /**
    * Define a jQuery plugin
    */
    var dataID = "spectrum.id";
    $.fn.spectrum = function (opts, extra) {

        if (typeof opts == "string") {

            var returnValue = this;
            var args = Array.prototype.slice.call( arguments, 1 );

            this.each(function () {
                var spect = spectrums[$(this).data(dataID)];
                if (spect) {
                    var method = spect[opts];
                    if (!method) {
                        throw new Error( "Spectrum: no such method: '" + opts + "'" );
                    }

                    if (opts == "get") {
                        returnValue = spect.get();
                    }
                    else if (opts == "container") {
                        returnValue = spect.container;
                    }
                    else if (opts == "option") {
                        returnValue = spect.option.apply(spect, args);
                    }
                    else if (opts == "destroy") {
                        spect.destroy();
                        $(this).removeData(dataID);
                    }
                    else {
                        method.apply(spect, args);
                    }
                }
            });

            return returnValue;
        }

        // Initializing a new instance of spectrum
        return this.spectrum("destroy").each(function () {
            var options = $.extend({}, opts, $(this).data());
            var spect = spectrum(this, options);
            $(this).data(dataID, spect.id);
        });
    };

    $.fn.spectrum.load = true;
    $.fn.spectrum.loadOpts = {};
    $.fn.spectrum.draggable = draggable;
    $.fn.spectrum.defaults = defaultOpts;

    $.spectrum = { };
    $.spectrum.localization = { };
    $.spectrum.palettes = { };

    $.fn.spectrum.processNativeColorInputs = function () {
        if (!inputTypeColorSupport) {
            $("input[type=color]").spectrum({
                preferredFormat: "hex6"
            });
        }
    };

    // TinyColor v1.1.1
    // https://github.com/bgrins/TinyColor
    // Brian Grinstead, MIT License

    (function() {

    var trimLeft = /^[\s,#]+/,
        trimRight = /\s+$/,
        tinyCounter = 0,
        math = Math,
        mathRound = math.round,
        mathMin = math.min,
        mathMax = math.max,
        mathRandom = math.random;

    var tinycolor = function tinycolor (color, opts) {

        color = (color) ? color : '';
        opts = opts || { };

        // If input is already a tinycolor, return itself
        if (color instanceof tinycolor) {
           return color;
        }
        // If we are called as a function, call using new instead
        if (!(this instanceof tinycolor)) {
        	/* jshint newcap: false */
            return new tinycolor(color, opts);
        }

        var rgb = inputToRGB(color);
        this._originalInput = color;
        this._r = rgb.r;
        this._g = rgb.g;
        this._b = rgb.b;
        this._a = rgb.a;
        this._roundA = mathRound(100*this._a)/100;
        this._format = opts.format || rgb.format;
        this._gradientType = opts.gradientType;

        // Don't let the range of [0,255] come back in [0,1].
        // Potentially lose a little bit of precision here, but will fix issues where
        // .5 gets interpreted as half of the total, instead of half of 1
        // If it was supposed to be 128, this was already taken care of by `inputToRgb`
        if (this._r < 1) { this._r = mathRound(this._r); }
        if (this._g < 1) { this._g = mathRound(this._g); }
        if (this._b < 1) { this._b = mathRound(this._b); }

        this._ok = rgb.ok;
        this._tc_id = tinyCounter++;
    };

    tinycolor.prototype = {
        isDark: function() {
            return this.getBrightness() < 128;
        },
        isLight: function() {
            return !this.isDark();
        },
        isValid: function() {
            return this._ok;
        },
        getOriginalInput: function() {
          return this._originalInput;
        },
        getFormat: function() {
            return this._format;
        },
        getAlpha: function() {
            return this._a;
        },
        getBrightness: function() {
            var rgb = this.toRgb();
            return (rgb.r * 299 + rgb.g * 587 + rgb.b * 114) / 1000;
        },
        setAlpha: function(value) {
            this._a = boundAlpha(value);
            this._roundA = mathRound(100*this._a) / 100;
            return this;
        },
        toHsv: function() {
            var hsv = rgbToHsv(this._r, this._g, this._b);
            return { h: hsv.h * 360, s: hsv.s, v: hsv.v, a: this._a };
        },
        toHsvString: function() {
            var hsv = rgbToHsv(this._r, this._g, this._b);
            var h = mathRound(hsv.h * 360), s = mathRound(hsv.s * 100), v = mathRound(hsv.v * 100);
            return (this._a == 1) ?
              "hsv("  + h + ", " + s + "%, " + v + "%)" :
              "hsva(" + h + ", " + s + "%, " + v + "%, "+ this._roundA + ")";
        },
        toHsl: function() {
            var hsl = rgbToHsl(this._r, this._g, this._b);
            return { h: hsl.h * 360, s: hsl.s, l: hsl.l, a: this._a };
        },
        toHslString: function() {
            var hsl = rgbToHsl(this._r, this._g, this._b);
            var h = mathRound(hsl.h * 360), s = mathRound(hsl.s * 100), l = mathRound(hsl.l * 100);
            return (this._a == 1) ?
              "hsl("  + h + ", " + s + "%, " + l + "%)" :
              "hsla(" + h + ", " + s + "%, " + l + "%, "+ this._roundA + ")";
        },
        toHex: function(allow3Char) {
            return rgbToHex(this._r, this._g, this._b, allow3Char);
        },
        toHexString: function(allow3Char) {
            return '#' + this.toHex(allow3Char);
        },
        toHex8: function() {
            return rgbaToHex(this._r, this._g, this._b, this._a);
        },
        toHex8String: function() {
            return '#' + this.toHex8();
        },
        toRgb: function() {
            return { r: mathRound(this._r), g: mathRound(this._g), b: mathRound(this._b), a: this._a };
        },
        toRgbString: function() {
            return (this._a == 1) ?
              "rgb("  + mathRound(this._r) + ", " + mathRound(this._g) + ", " + mathRound(this._b) + ")" :
              "rgba(" + mathRound(this._r) + ", " + mathRound(this._g) + ", " + mathRound(this._b) + ", " + this._roundA + ")";
        },
        toPercentageRgb: function() {
            return { r: mathRound(bound01(this._r, 255) * 100) + "%", g: mathRound(bound01(this._g, 255) * 100) + "%", b: mathRound(bound01(this._b, 255) * 100) + "%", a: this._a };
        },
        toPercentageRgbString: function() {
            return (this._a == 1) ?
              "rgb("  + mathRound(bound01(this._r, 255) * 100) + "%, " + mathRound(bound01(this._g, 255) * 100) + "%, " + mathRound(bound01(this._b, 255) * 100) + "%)" :
              "rgba(" + mathRound(bound01(this._r, 255) * 100) + "%, " + mathRound(bound01(this._g, 255) * 100) + "%, " + mathRound(bound01(this._b, 255) * 100) + "%, " + this._roundA + ")";
        },
        toName: function() {
            if (this._a === 0) {
                return "transparent";
            }

            if (this._a < 1) {
                return false;
            }

            return hexNames[rgbToHex(this._r, this._g, this._b, true)] || false;
        },
        toFilter: function(secondColor) {
            var hex8String = '#' + rgbaToHex(this._r, this._g, this._b, this._a);
            var secondHex8String = hex8String;
            var gradientType = this._gradientType ? "GradientType = 1, " : "";

            if (secondColor) {
                var s = tinycolor(secondColor);
                secondHex8String = s.toHex8String();
            }

            return "progid:DXImageTransform.Microsoft.gradient("+gradientType+"startColorstr="+hex8String+",endColorstr="+secondHex8String+")";
        },
        toString: function(format) {
            var formatSet = !!format;
            format = format || this._format;

            var formattedString = false;
            var hasAlpha = this._a < 1 && this._a >= 0;
            var needsAlphaFormat = !formatSet && hasAlpha && (format === "hex" || format === "hex6" || format === "hex3" || format === "name");

            if (needsAlphaFormat) {
                // Special case for "transparent", all other non-alpha formats
                // will return rgba when there is transparency.
                if (format === "name" && this._a === 0) {
                    return this.toName();
                }
                return this.toRgbString();
            }
            if (format === "rgb") {
                formattedString = this.toRgbString();
            }
            if (format === "prgb") {
                formattedString = this.toPercentageRgbString();
            }
            if (format === "hex" || format === "hex6") {
                formattedString = this.toHexString();
            }
            if (format === "hex3") {
                formattedString = this.toHexString(true);
            }
            if (format === "hex8") {
                formattedString = this.toHex8String();
            }
            if (format === "name") {
                formattedString = this.toName();
            }
            if (format === "hsl") {
                formattedString = this.toHslString();
            }
            if (format === "hsv") {
                formattedString = this.toHsvString();
            }

            return formattedString || this.toHexString();
        },

        _applyModification: function(fn, args) {
            var color = fn.apply(null, [this].concat([].slice.call(args)));
            this._r = color._r;
            this._g = color._g;
            this._b = color._b;
            this.setAlpha(color._a);
            return this;
        },
        lighten: function() {
            return this._applyModification(lighten, arguments);
        },
        brighten: function() {
            return this._applyModification(brighten, arguments);
        },
        darken: function() {
            return this._applyModification(darken, arguments);
        },
        desaturate: function() {
            return this._applyModification(desaturate, arguments);
        },
        saturate: function() {
            return this._applyModification(saturate, arguments);
        },
        greyscale: function() {
            return this._applyModification(greyscale, arguments);
        },
        spin: function() {
            return this._applyModification(spin, arguments);
        },

        _applyCombination: function(fn, args) {
            return fn.apply(null, [this].concat([].slice.call(args)));
        },
        analogous: function() {
            return this._applyCombination(analogous, arguments);
        },
        complement: function() {
            return this._applyCombination(complement, arguments);
        },
        monochromatic: function() {
            return this._applyCombination(monochromatic, arguments);
        },
        splitcomplement: function() {
            return this._applyCombination(splitcomplement, arguments);
        },
        triad: function() {
            return this._applyCombination(triad, arguments);
        },
        tetrad: function() {
            return this._applyCombination(tetrad, arguments);
        }
    };

    // If input is an object, force 1 into "1.0" to handle ratios properly
    // String input requires "1.0" as input, so 1 will be treated as 1
    tinycolor.fromRatio = function(color, opts) {
        if (typeof color == "object") {
            var newColor = {};
            for (var i in color) {
                if (color.hasOwnProperty(i)) {
                    if (i === "a") {
                        newColor[i] = color[i];
                    }
                    else {
                        newColor[i] = convertToPercentage(color[i]);
                    }
                }
            }
            color = newColor;
        }

        return tinycolor(color, opts);
    };

    // Given a string or object, convert that input to RGB
    // Possible string inputs:
    //
    //     "red"
    //     "#f00" or "f00"
    //     "#ff0000" or "ff0000"
    //     "#ff000000" or "ff000000"
    //     "rgb 255 0 0" or "rgb (255, 0, 0)"
    //     "rgb 1.0 0 0" or "rgb (1, 0, 0)"
    //     "rgba (255, 0, 0, 1)" or "rgba 255, 0, 0, 1"
    //     "rgba (1.0, 0, 0, 1)" or "rgba 1.0, 0, 0, 1"
    //     "hsl(0, 100%, 50%)" or "hsl 0 100% 50%"
    //     "hsla(0, 100%, 50%, 1)" or "hsla 0 100% 50%, 1"
    //     "hsv(0, 100%, 100%)" or "hsv 0 100% 100%"
    //
    function inputToRGB(color) {

        var rgb = { r: 0, g: 0, b: 0 };
        var a = 1;
        var ok = false;
        var format = false;

        if (typeof color == "string") {
            color = stringInputToObject(color);
        }

        if (typeof color == "object") {
            if (color.hasOwnProperty("r") && color.hasOwnProperty("g") && color.hasOwnProperty("b")) {
                rgb = rgbToRgb(color.r, color.g, color.b);
                ok = true;
                format = String(color.r).substr(-1) === "%" ? "prgb" : "rgb";
            }
            else if (color.hasOwnProperty("h") && color.hasOwnProperty("s") && color.hasOwnProperty("v")) {
                color.s = convertToPercentage(color.s);
                color.v = convertToPercentage(color.v);
                rgb = hsvToRgb(color.h, color.s, color.v);
                ok = true;
                format = "hsv";
            }
            else if (color.hasOwnProperty("h") && color.hasOwnProperty("s") && color.hasOwnProperty("l")) {
                color.s = convertToPercentage(color.s);
                color.l = convertToPercentage(color.l);
                rgb = hslToRgb(color.h, color.s, color.l);
                ok = true;
                format = "hsl";
            }

            if (color.hasOwnProperty("a")) {
                a = color.a;
            }
        }

        a = boundAlpha(a);

        return {
            ok: ok,
            format: color.format || format,
            r: mathMin(255, mathMax(rgb.r, 0)),
            g: mathMin(255, mathMax(rgb.g, 0)),
            b: mathMin(255, mathMax(rgb.b, 0)),
            a: a
        };
    }


    // Conversion Functions
    // --------------------

    // `rgbToHsl`, `rgbToHsv`, `hslToRgb`, `hsvToRgb` modified from:
    // <http://mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript>

    // `rgbToRgb`
    // Handle bounds / percentage checking to conform to CSS color spec
    // <http://www.w3.org/TR/css3-color/>
    // *Assumes:* r, g, b in [0, 255] or [0, 1]
    // *Returns:* { r, g, b } in [0, 255]
    function rgbToRgb(r, g, b){
        return {
            r: bound01(r, 255) * 255,
            g: bound01(g, 255) * 255,
            b: bound01(b, 255) * 255
        };
    }

    // `rgbToHsl`
    // Converts an RGB color value to HSL.
    // *Assumes:* r, g, and b are contained in [0, 255] or [0, 1]
    // *Returns:* { h, s, l } in [0,1]
    function rgbToHsl(r, g, b) {

        r = bound01(r, 255);
        g = bound01(g, 255);
        b = bound01(b, 255);

        var max = mathMax(r, g, b), min = mathMin(r, g, b);
        var h, s, l = (max + min) / 2;

        if(max == min) {
            h = s = 0; // achromatic
        }
        else {
            var d = max - min;
            s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
            switch(max) {
                case r: h = (g - b) / d + (g < b ? 6 : 0); break;
                case g: h = (b - r) / d + 2; break;
                case b: h = (r - g) / d + 4; break;
            }

            h /= 6;
        }

        return { h: h, s: s, l: l };
    }

    // `hslToRgb`
    // Converts an HSL color value to RGB.
    // *Assumes:* h is contained in [0, 1] or [0, 360] and s and l are contained [0, 1] or [0, 100]
    // *Returns:* { r, g, b } in the set [0, 255]
    function hslToRgb(h, s, l) {
        var r, g, b;

        h = bound01(h, 360);
        s = bound01(s, 100);
        l = bound01(l, 100);

        function hue2rgb(p, q, t) {
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        }

        if(s === 0) {
            r = g = b = l; // achromatic
        }
        else {
            var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
            var p = 2 * l - q;
            r = hue2rgb(p, q, h + 1/3);
            g = hue2rgb(p, q, h);
            b = hue2rgb(p, q, h - 1/3);
        }

        return { r: r * 255, g: g * 255, b: b * 255 };
    }

    // `rgbToHsv`
    // Converts an RGB color value to HSV
    // *Assumes:* r, g, and b are contained in the set [0, 255] or [0, 1]
    // *Returns:* { h, s, v } in [0,1]
    function rgbToHsv(r, g, b) {

        r = bound01(r, 255);
        g = bound01(g, 255);
        b = bound01(b, 255);

        var max = mathMax(r, g, b), min = mathMin(r, g, b);
        var h, s, v = max;

        var d = max - min;
        s = max === 0 ? 0 : d / max;

        if(max == min) {
            h = 0; // achromatic
        }
        else {
            switch(max) {
                case r: h = (g - b) / d + (g < b ? 6 : 0); break;
                case g: h = (b - r) / d + 2; break;
                case b: h = (r - g) / d + 4; break;
            }
            h /= 6;
        }
        return { h: h, s: s, v: v };
    }

    // `hsvToRgb`
    // Converts an HSV color value to RGB.
    // *Assumes:* h is contained in [0, 1] or [0, 360] and s and v are contained in [0, 1] or [0, 100]
    // *Returns:* { r, g, b } in the set [0, 255]
     function hsvToRgb(h, s, v) {

        h = bound01(h, 360) * 6;
        s = bound01(s, 100);
        v = bound01(v, 100);

        var i = math.floor(h),
            f = h - i,
            p = v * (1 - s),
            q = v * (1 - f * s),
            t = v * (1 - (1 - f) * s),
            mod = i % 6,
            r = [v, q, p, p, t, v][mod],
            g = [t, v, v, q, p, p][mod],
            b = [p, p, t, v, v, q][mod];

        return { r: r * 255, g: g * 255, b: b * 255 };
    }

    // `rgbToHex`
    // Converts an RGB color to hex
    // Assumes r, g, and b are contained in the set [0, 255]
    // Returns a 3 or 6 character hex
    function rgbToHex(r, g, b, allow3Char) {

        var hex = [
            pad2(mathRound(r).toString(16)),
            pad2(mathRound(g).toString(16)),
            pad2(mathRound(b).toString(16))
        ];

        // Return a 3 character hex if possible
        if (allow3Char && hex[0].charAt(0) == hex[0].charAt(1) && hex[1].charAt(0) == hex[1].charAt(1) && hex[2].charAt(0) == hex[2].charAt(1)) {
            return hex[0].charAt(0) + hex[1].charAt(0) + hex[2].charAt(0);
        }

        return hex.join("");
    }
        // `rgbaToHex`
        // Converts an RGBA color plus alpha transparency to hex
        // Assumes r, g, b and a are contained in the set [0, 255]
        // Returns an 8 character hex
        function rgbaToHex(r, g, b, a) {

            var hex = [
                pad2(convertDecimalToHex(a)),
                pad2(mathRound(r).toString(16)),
                pad2(mathRound(g).toString(16)),
                pad2(mathRound(b).toString(16))
            ];

            return hex.join("");
        }

    // `equals`
    // Can be called with any tinycolor input
    tinycolor.equals = function (color1, color2) {
        if (!color1 || !color2) { return false; }
        return tinycolor(color1).toRgbString() == tinycolor(color2).toRgbString();
    };
    tinycolor.random = function() {
        return tinycolor.fromRatio({
            r: mathRandom(),
            g: mathRandom(),
            b: mathRandom()
        });
    };


    // Modification Functions
    // ----------------------
    // Thanks to less.js for some of the basics here
    // <https://github.com/cloudhead/less.js/blob/master/lib/less/functions.js>

    function desaturate(color, amount) {
        amount = (amount === 0) ? 0 : (amount || 10);
        var hsl = tinycolor(color).toHsl();
        hsl.s -= amount / 100;
        hsl.s = clamp01(hsl.s);
        return tinycolor(hsl);
    }

    function saturate(color, amount) {
        amount = (amount === 0) ? 0 : (amount || 10);
        var hsl = tinycolor(color).toHsl();
        hsl.s += amount / 100;
        hsl.s = clamp01(hsl.s);
        return tinycolor(hsl);
    }

    function greyscale(color) {
        return tinycolor(color).desaturate(100);
    }

    function lighten (color, amount) {
        amount = (amount === 0) ? 0 : (amount || 10);
        var hsl = tinycolor(color).toHsl();
        hsl.l += amount / 100;
        hsl.l = clamp01(hsl.l);
        return tinycolor(hsl);
    }

    function brighten(color, amount) {
        amount = (amount === 0) ? 0 : (amount || 10);
        var rgb = tinycolor(color).toRgb();
        rgb.r = mathMax(0, mathMin(255, rgb.r - mathRound(255 * - (amount / 100))));
        rgb.g = mathMax(0, mathMin(255, rgb.g - mathRound(255 * - (amount / 100))));
        rgb.b = mathMax(0, mathMin(255, rgb.b - mathRound(255 * - (amount / 100))));
        return tinycolor(rgb);
    }

    function darken (color, amount) {
        amount = (amount === 0) ? 0 : (amount || 10);
        var hsl = tinycolor(color).toHsl();
        hsl.l -= amount / 100;
        hsl.l = clamp01(hsl.l);
        return tinycolor(hsl);
    }

    // Spin takes a positive or negative amount within [-360, 360] indicating the change of hue.
    // Values outside of this range will be wrapped into this range.
    function spin(color, amount) {
        var hsl = tinycolor(color).toHsl();
        var hue = (mathRound(hsl.h) + amount) % 360;
        hsl.h = hue < 0 ? 360 + hue : hue;
        return tinycolor(hsl);
    }

    // Combination Functions
    // ---------------------
    // Thanks to jQuery xColor for some of the ideas behind these
    // <https://github.com/infusion/jQuery-xcolor/blob/master/jquery.xcolor.js>

    function complement(color) {
        var hsl = tinycolor(color).toHsl();
        hsl.h = (hsl.h + 180) % 360;
        return tinycolor(hsl);
    }

    function triad(color) {
        var hsl = tinycolor(color).toHsl();
        var h = hsl.h;
        return [
            tinycolor(color),
            tinycolor({ h: (h + 120) % 360, s: hsl.s, l: hsl.l }),
            tinycolor({ h: (h + 240) % 360, s: hsl.s, l: hsl.l })
        ];
    }

    function tetrad(color) {
        var hsl = tinycolor(color).toHsl();
        var h = hsl.h;
        return [
            tinycolor(color),
            tinycolor({ h: (h + 90) % 360, s: hsl.s, l: hsl.l }),
            tinycolor({ h: (h + 180) % 360, s: hsl.s, l: hsl.l }),
            tinycolor({ h: (h + 270) % 360, s: hsl.s, l: hsl.l })
        ];
    }

    function splitcomplement(color) {
        var hsl = tinycolor(color).toHsl();
        var h = hsl.h;
        return [
            tinycolor(color),
            tinycolor({ h: (h + 72) % 360, s: hsl.s, l: hsl.l}),
            tinycolor({ h: (h + 216) % 360, s: hsl.s, l: hsl.l})
        ];
    }

    function analogous(color, results, slices) {
        results = results || 6;
        slices = slices || 30;

        var hsl = tinycolor(color).toHsl();
        var part = 360 / slices;
        var ret = [tinycolor(color)];

        for (hsl.h = ((hsl.h - (part * results >> 1)) + 720) % 360; --results; ) {
            hsl.h = (hsl.h + part) % 360;
            ret.push(tinycolor(hsl));
        }
        return ret;
    }

    function monochromatic(color, results) {
        results = results || 6;
        var hsv = tinycolor(color).toHsv();
        var h = hsv.h, s = hsv.s, v = hsv.v;
        var ret = [];
        var modification = 1 / results;

        while (results--) {
            ret.push(tinycolor({ h: h, s: s, v: v}));
            v = (v + modification) % 1;
        }

        return ret;
    }

    // Utility Functions
    // ---------------------

    tinycolor.mix = function(color1, color2, amount) {
        amount = (amount === 0) ? 0 : (amount || 50);

        var rgb1 = tinycolor(color1).toRgb();
        var rgb2 = tinycolor(color2).toRgb();

        var p = amount / 100;
        var w = p * 2 - 1;
        var a = rgb2.a - rgb1.a;

        var w1;

        if (w * a == -1) {
            w1 = w;
        } else {
            w1 = (w + a) / (1 + w * a);
        }

        w1 = (w1 + 1) / 2;

        var w2 = 1 - w1;

        var rgba = {
            r: rgb2.r * w1 + rgb1.r * w2,
            g: rgb2.g * w1 + rgb1.g * w2,
            b: rgb2.b * w1 + rgb1.b * w2,
            a: rgb2.a * p  + rgb1.a * (1 - p)
        };

        return tinycolor(rgba);
    };


    // Readability Functions
    // ---------------------
    // <http://www.w3.org/TR/AERT#color-contrast>

    // `readability`
    // Analyze the 2 colors and returns an object with the following properties:
    //    `brightness`: difference in brightness between the two colors
    //    `color`: difference in color/hue between the two colors
    tinycolor.readability = function(color1, color2) {
        var c1 = tinycolor(color1);
        var c2 = tinycolor(color2);
        var rgb1 = c1.toRgb();
        var rgb2 = c2.toRgb();
        var brightnessA = c1.getBrightness();
        var brightnessB = c2.getBrightness();
        var colorDiff = (
            Math.max(rgb1.r, rgb2.r) - Math.min(rgb1.r, rgb2.r) +
            Math.max(rgb1.g, rgb2.g) - Math.min(rgb1.g, rgb2.g) +
            Math.max(rgb1.b, rgb2.b) - Math.min(rgb1.b, rgb2.b)
        );

        return {
            brightness: Math.abs(brightnessA - brightnessB),
            color: colorDiff
        };
    };

    // `readable`
    // http://www.w3.org/TR/AERT#color-contrast
    // Ensure that foreground and background color combinations provide sufficient contrast.
    // *Example*
    //    tinycolor.isReadable("#000", "#111") => false
    tinycolor.isReadable = function(color1, color2) {
        var readability = tinycolor.readability(color1, color2);
        return readability.brightness > 125 && readability.color > 500;
    };

    // `mostReadable`
    // Given a base color and a list of possible foreground or background
    // colors for that base, returns the most readable color.
    // *Example*
    //    tinycolor.mostReadable("#123", ["#fff", "#000"]) => "#000"
    tinycolor.mostReadable = function(baseColor, colorList) {
        var bestColor = null;
        var bestScore = 0;
        var bestIsReadable = false;
        for (var i=0; i < colorList.length; i++) {

            // We normalize both around the "acceptable" breaking point,
            // but rank brightness constrast higher than hue.

            var readability = tinycolor.readability(baseColor, colorList[i]);
            var readable = readability.brightness > 125 && readability.color > 500;
            var score = 3 * (readability.brightness / 125) + (readability.color / 500);

            if ((readable && ! bestIsReadable) ||
                (readable && bestIsReadable && score > bestScore) ||
                ((! readable) && (! bestIsReadable) && score > bestScore)) {
                bestIsReadable = readable;
                bestScore = score;
                bestColor = tinycolor(colorList[i]);
            }
        }
        return bestColor;
    };


    // Big List of Colors
    // ------------------
    // <http://www.w3.org/TR/css3-color/#svg-color>
    var names = tinycolor.names = {
        aliceblue: "f0f8ff",
        antiquewhite: "faebd7",
        aqua: "0ff",
        aquamarine: "7fffd4",
        azure: "f0ffff",
        beige: "f5f5dc",
        bisque: "ffe4c4",
        black: "000",
        blanchedalmond: "ffebcd",
        blue: "00f",
        blueviolet: "8a2be2",
        brown: "a52a2a",
        burlywood: "deb887",
        burntsienna: "ea7e5d",
        cadetblue: "5f9ea0",
        chartreuse: "7fff00",
        chocolate: "d2691e",
        coral: "ff7f50",
        cornflowerblue: "6495ed",
        cornsilk: "fff8dc",
        crimson: "dc143c",
        cyan: "0ff",
        darkblue: "00008b",
        darkcyan: "008b8b",
        darkgoldenrod: "b8860b",
        darkgray: "a9a9a9",
        darkgreen: "006400",
        darkgrey: "a9a9a9",
        darkkhaki: "bdb76b",
        darkmagenta: "8b008b",
        darkolivegreen: "556b2f",
        darkorange: "ff8c00",
        darkorchid: "9932cc",
        darkred: "8b0000",
        darksalmon: "e9967a",
        darkseagreen: "8fbc8f",
        darkslateblue: "483d8b",
        darkslategray: "2f4f4f",
        darkslategrey: "2f4f4f",
        darkturquoise: "00ced1",
        darkviolet: "9400d3",
        deeppink: "ff1493",
        deepskyblue: "00bfff",
        dimgray: "696969",
        dimgrey: "696969",
        dodgerblue: "1e90ff",
        firebrick: "b22222",
        floralwhite: "fffaf0",
        forestgreen: "228b22",
        fuchsia: "f0f",
        gainsboro: "dcdcdc",
        ghostwhite: "f8f8ff",
        gold: "ffd700",
        goldenrod: "daa520",
        gray: "808080",
        green: "008000",
        greenyellow: "adff2f",
        grey: "808080",
        honeydew: "f0fff0",
        hotpink: "ff69b4",
        indianred: "cd5c5c",
        indigo: "4b0082",
        ivory: "fffff0",
        khaki: "f0e68c",
        lavender: "e6e6fa",
        lavenderblush: "fff0f5",
        lawngreen: "7cfc00",
        lemonchiffon: "fffacd",
        lightblue: "add8e6",
        lightcoral: "f08080",
        lightcyan: "e0ffff",
        lightgoldenrodyellow: "fafad2",
        lightgray: "d3d3d3",
        lightgreen: "90ee90",
        lightgrey: "d3d3d3",
        lightpink: "ffb6c1",
        lightsalmon: "ffa07a",
        lightseagreen: "20b2aa",
        lightskyblue: "87cefa",
        lightslategray: "789",
        lightslategrey: "789",
        lightsteelblue: "b0c4de",
        lightyellow: "ffffe0",
        lime: "0f0",
        limegreen: "32cd32",
        linen: "faf0e6",
        magenta: "f0f",
        maroon: "800000",
        mediumaquamarine: "66cdaa",
        mediumblue: "0000cd",
        mediumorchid: "ba55d3",
        mediumpurple: "9370db",
        mediumseagreen: "3cb371",
        mediumslateblue: "7b68ee",
        mediumspringgreen: "00fa9a",
        mediumturquoise: "48d1cc",
        mediumvioletred: "c71585",
        midnightblue: "191970",
        mintcream: "f5fffa",
        mistyrose: "ffe4e1",
        moccasin: "ffe4b5",
        navajowhite: "ffdead",
        navy: "000080",
        oldlace: "fdf5e6",
        olive: "808000",
        olivedrab: "6b8e23",
        orange: "ffa500",
        orangered: "ff4500",
        orchid: "da70d6",
        palegoldenrod: "eee8aa",
        palegreen: "98fb98",
        paleturquoise: "afeeee",
        palevioletred: "db7093",
        papayawhip: "ffefd5",
        peachpuff: "ffdab9",
        peru: "cd853f",
        pink: "ffc0cb",
        plum: "dda0dd",
        powderblue: "b0e0e6",
        purple: "800080",
        rebeccapurple: "663399",
        red: "f00",
        rosybrown: "bc8f8f",
        royalblue: "4169e1",
        saddlebrown: "8b4513",
        salmon: "fa8072",
        sandybrown: "f4a460",
        seagreen: "2e8b57",
        seashell: "fff5ee",
        sienna: "a0522d",
        silver: "c0c0c0",
        skyblue: "87ceeb",
        slateblue: "6a5acd",
        slategray: "708090",
        slategrey: "708090",
        snow: "fffafa",
        springgreen: "00ff7f",
        steelblue: "4682b4",
        tan: "d2b48c",
        teal: "008080",
        thistle: "d8bfd8",
        tomato: "ff6347",
        turquoise: "40e0d0",
        violet: "ee82ee",
        wheat: "f5deb3",
        white: "fff",
        whitesmoke: "f5f5f5",
        yellow: "ff0",
        yellowgreen: "9acd32"
    };

    // Make it easy to access colors via `hexNames[hex]`
    var hexNames = tinycolor.hexNames = flip(names);


    // Utilities
    // ---------

    // `{ 'name1': 'val1' }` becomes `{ 'val1': 'name1' }`
    function flip(o) {
        var flipped = { };
        for (var i in o) {
            if (o.hasOwnProperty(i)) {
                flipped[o[i]] = i;
            }
        }
        return flipped;
    }

    // Return a valid alpha value [0,1] with all invalid values being set to 1
    function boundAlpha(a) {
        a = parseFloat(a);

        if (isNaN(a) || a < 0 || a > 1) {
            a = 1;
        }

        return a;
    }

    // Take input from [0, n] and return it as [0, 1]
    function bound01(n, max) {
        if (isOnePointZero(n)) { n = "100%"; }

        var processPercent = isPercentage(n);
        n = mathMin(max, mathMax(0, parseFloat(n)));

        // Automatically convert percentage into number
        if (processPercent) {
            n = parseInt(n * max, 10) / 100;
        }

        // Handle floating point rounding errors
        if ((math.abs(n - max) < 0.000001)) {
            return 1;
        }

        // Convert into [0, 1] range if it isn't already
        return (n % max) / parseFloat(max);
    }

    // Force a number between 0 and 1
    function clamp01(val) {
        return mathMin(1, mathMax(0, val));
    }

    // Parse a base-16 hex value into a base-10 integer
    function parseIntFromHex(val) {
        return parseInt(val, 16);
    }

    // Need to handle 1.0 as 100%, since once it is a number, there is no difference between it and 1
    // <http://stackoverflow.com/questions/7422072/javascript-how-to-detect-number-as-a-decimal-including-1-0>
    function isOnePointZero(n) {
        return typeof n == "string" && n.indexOf('.') != -1 && parseFloat(n) === 1;
    }

    // Check to see if string passed in is a percentage
    function isPercentage(n) {
        return typeof n === "string" && n.indexOf('%') != -1;
    }

    // Force a hex value to have 2 characters
    function pad2(c) {
        return c.length == 1 ? '0' + c : '' + c;
    }

    // Replace a decimal with it's percentage value
    function convertToPercentage(n) {
        if (n <= 1) {
            n = (n * 100) + "%";
        }

        return n;
    }

    // Converts a decimal to a hex value
    function convertDecimalToHex(d) {
        return Math.round(parseFloat(d) * 255).toString(16);
    }
    // Converts a hex value to a decimal
    function convertHexToDecimal(h) {
        return (parseIntFromHex(h) / 255);
    }

    var matchers = (function() {

        // <http://www.w3.org/TR/css3-values/#integers>
        var CSS_INTEGER = "[-\\+]?\\d+%?";

        // <http://www.w3.org/TR/css3-values/#number-value>
        var CSS_NUMBER = "[-\\+]?\\d*\\.\\d+%?";

        // Allow positive/negative integer/number.  Don't capture the either/or, just the entire outcome.
        var CSS_UNIT = "(?:" + CSS_NUMBER + ")|(?:" + CSS_INTEGER + ")";

        // Actual matching.
        // Parentheses and commas are optional, but not required.
        // Whitespace can take the place of commas or opening paren
        var PERMISSIVE_MATCH3 = "[\\s|\\(]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")\\s*\\)?";
        var PERMISSIVE_MATCH4 = "[\\s|\\(]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")\\s*\\)?";

        return {
            rgb: new RegExp("rgb" + PERMISSIVE_MATCH3),
            rgba: new RegExp("rgba" + PERMISSIVE_MATCH4),
            hsl: new RegExp("hsl" + PERMISSIVE_MATCH3),
            hsla: new RegExp("hsla" + PERMISSIVE_MATCH4),
            hsv: new RegExp("hsv" + PERMISSIVE_MATCH3),
            hsva: new RegExp("hsva" + PERMISSIVE_MATCH4),
            hex3: /^([0-9a-fA-F]{1})([0-9a-fA-F]{1})([0-9a-fA-F]{1})$/,
            hex6: /^([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})$/,
            hex8: /^([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})$/
        };
    })();

    // `stringInputToObject`
    // Permissive string parsing.  Take in a number of formats, and output an object
    // based on detected format.  Returns `{ r, g, b }` or `{ h, s, l }` or `{ h, s, v}`
    function stringInputToObject(color) {

        color = color.replace(trimLeft,'').replace(trimRight, '').toLowerCase();
        var named = false;
        if (names[color]) {
            color = names[color];
            named = true;
        }
        else if (color == 'transparent') {
            return { r: 0, g: 0, b: 0, a: 0, format: "name" };
        }

        // Try to match string input using regular expressions.
        // Keep most of the number bounding out of this function - don't worry about [0,1] or [0,100] or [0,360]
        // Just return an object and let the conversion functions handle that.
        // This way the result will be the same whether the tinycolor is initialized with string or object.
        var match;
        if ((match = matchers.rgb.exec(color))) {
            return { r: match[1], g: match[2], b: match[3] };
        }
        if ((match = matchers.rgba.exec(color))) {
            return { r: match[1], g: match[2], b: match[3], a: match[4] };
        }
        if ((match = matchers.hsl.exec(color))) {
            return { h: match[1], s: match[2], l: match[3] };
        }
        if ((match = matchers.hsla.exec(color))) {
            return { h: match[1], s: match[2], l: match[3], a: match[4] };
        }
        if ((match = matchers.hsv.exec(color))) {
            return { h: match[1], s: match[2], v: match[3] };
        }
        if ((match = matchers.hsva.exec(color))) {
            return { h: match[1], s: match[2], v: match[3], a: match[4] };
        }
        if ((match = matchers.hex8.exec(color))) {
            return {
                a: convertHexToDecimal(match[1]),
                r: parseIntFromHex(match[2]),
                g: parseIntFromHex(match[3]),
                b: parseIntFromHex(match[4]),
                format: named ? "name" : "hex8"
            };
        }
        if ((match = matchers.hex6.exec(color))) {
            return {
                r: parseIntFromHex(match[1]),
                g: parseIntFromHex(match[2]),
                b: parseIntFromHex(match[3]),
                format: named ? "name" : "hex"
            };
        }
        if ((match = matchers.hex3.exec(color))) {
            return {
                r: parseIntFromHex(match[1] + '' + match[1]),
                g: parseIntFromHex(match[2] + '' + match[2]),
                b: parseIntFromHex(match[3] + '' + match[3]),
                format: named ? "name" : "hex"
            };
        }

        return false;
    }

    window.tinycolor = tinycolor;
    })();


    $(function () {
        if ($.fn.spectrum.load) {
            $.fn.spectrum.processNativeColorInputs();
        }
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
		
		//applies the current selectedstyle and colorscheme for name
		var updateStyle = function(name) {
			var rec = modelsAndStyles[name];
			var stylename = rec.selectedStyle;
			var s = rec.styles[stylename];
			var style = s.style;
			var nbond = s.nonbond;

			var model = rec.model;
			if(model) {
				model.setStyle({}, style);
				model.setStyle({bonds: 0}, nbond);
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
				surface = viewer.addSurface($3Dmol.SurfaceType.VDW, 
						surfaceStyle, {model:receptor, bonds: 0, invert:true});
				viewer.zoomTo({});
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
			if(m.updateStyle) m.updateStyle(newstyle);
			viewer.render();
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
		
	}

	return Viewer;
})();
