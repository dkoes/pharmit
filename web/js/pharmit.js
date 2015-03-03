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
	Represents a single pharmacophore feature.
*/

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
	var select = this.select = $('<select name="featurename">').addClass('pharmit_featureselect').appendTo(editdiv);
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
	
	var c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>x:</label>').appendTo(c);
	this.x = $('<input>').appendTo(c).addClass('pharmit_coordinput').change(function() {
		//change x value
		var x = this.value = round(parseFloat(this.value));
		F.xsummary.text(numeral(x).format('0.[00]'));
		F.obj.x = round(x);
		F.updateViewer();
		});
	this.x.spinner(spinObject(F.x,{step: 0.1, numberFormat: 'n'}));
	
	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>y:</label>').appendTo(c);
	this.y = $('<input>').appendTo(c).addClass('pharmit_coordinput');
	this.y.spinner(spinObject(F.y,{step: 0.1, numberFormat: 'n'})).change(function() {
		//change y value
		var y = this.value = round(parseFloat(this.value));
		F.ysummary.text(numeral(y).format('0.[00]'));
		F.obj.y = round(y);
		F.updateViewer();
		});

	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>z:</label>').appendTo(c);
	this.z = $('<input>').appendTo(c).addClass('pharmit_coordinput');
	this.z.spinner(spinObject(F.z,{step: 0.1, numberFormat: 'n'})).change(function() {
		var z = this.value = round(parseFloat(this.value));
		F.zsummary.text(numeral(z).format('0.[00]'));
		F.obj.z = z;
		F.updateViewer();
		});

	//radius
	c = $('<div>').appendTo(locationdiv).addClass('pharmit_coorddiv');
	$('<label>Radius:</label>').appendTo(c);
	this.radius = $('<input>').appendTo(c).addClass('pharmit_radiusinput');
	this.radius.spinner(spinObject(F.radius,{step: 0.1, numberFormat: 'n'})).change(function() {
		F.rsummary.text(numeral(this.value).format('0.[00]'));
		F.obj.radius = this.value = round(this.value);
		F.updateViewer();
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
		'Hydrophobic', 'NegativeIon', 'PositiveIon','ExclusionSphere'];

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
	
	//global variable checking - we should add nothing but Pharmit
	var globalsBefore = {};
    for (var key in window)
         globalsBefore[key] = true;

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
	
	
	Pharmit.server = '/fcgi-bin/pharmitserv.fcgi';
	Pharmit.email = 'dkoes+pharmit@pitt.edu';
	
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
		var onclose = null; //what to call when close button is clicked
		
		
		//public variables and functions
		
		//perform the query
		this.minimize = function(qid, closer) {
			onclose = closer;
			//if we aren't hidden, need to cancel current query first
			
			//start provided query
			
			//show div
			mindiv.show();
		};
		
		//download results
		var saveResults = function() {
			
		};
		
		//get results filtered
		var applyFilters = function() {
			
		};
		
		//initialization code
		mindiv = $('<div>').appendTo(results.div).addClass('pharmit_rescontainer');
		//header
		var header = $('<div>').appendTo(mindiv).addClass("pharmit_resheader");
		var title = $('<div>Minimization Results</div>').appendTo(header).addClass('pharmit_heading').addClass("pharmit_rightheading");
		var closediv = $('<div>').addClass('pharmit_resclose').appendTo(title).click(function() {
			//cancel the current query 
			
			//hide
			mindiv.hide();
			//do what our caller told us to do on close
			if(onclose) onclose();
		});
		var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
		
		
		//body, should stretch to fill
		var body = $('<div>').appendTo(mindiv).addClass("pharmit_resbody");

		//footer
		var footer = $('<div>').appendTo(mindiv).addClass("pharmit_resfooter");
		var paramdiv = $('<div>').appendTo(footer).addClass("pharmit_minparams").disableSelection();
		//filters for minimization
		var table = $('<table>').appendTo(paramdiv);
		var row = $('<tr>').appendTo(table).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max Score</label>'));
		var cell = $('<td>').appendTo(row);
		$('<input name="sminamaxscore">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		$('<td>').appendTo(row).addClass('pharmit_checkcell').append($('<input type="checkbox" name="sminaunique">'));
		$('<td>').appendTo(row).append($('<label>Single<br>conformer</label>'));
		
		
		row = $('<tr>').appendTo(table).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max mRMSD</label>'));
		cell = $('<td>').appendTo(row);
		$('<input name="sminamaxrmsd">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		cell = $('<td colspan=2>').appendTo(row).addClass('pharmit_applycell');
		$('<button>Apply</button>').appendTo(cell).button().click(applyFilters);
		
		//save button
		var bottomrow = $('<div>').appendTo(footer).addClass('pharmit_minbottom');
		var save = $('<button>Save...</button>').appendTo(bottomrow).button().click(saveResults);				

		mindiv.hide();
	}

	return MinResults;
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
	phresults.js
	This is a div for managing pharmacophore results.
	
*/


var Pharmit = Pharmit || {};

Pharmit.PhResults = (function() {
	// private class variables and functions
	var phdiv = null;
	var table = null;
	var body = null;
	
	function PhResults(results, viewer, minresults) {
		//private variables and functions
		var phdiv = null;
		var qid = null;
		
		//format that provided data (mangle the names appropriately)		
		var processData = function(data) {
			var ret = data.data;
			
			for(var i = 0; i < ret.length; i++) {
				//round rmsd
				ret[i][1] = numeral(ret[i][1]).format('0.000');
			}
			return ret;
		};
		
		
		var lastheight = 0;
		var lastnum = 0;
		var resize = function() {
			if($.fn.DataTable.isDataTable(table)) {
				var total = body.height();
				if(total != lastheight) {
					lastheight = total;
					total -= $('thead', table).height();
					total -= $('.dataTables_info', body).height();
					total -= $('.dataTables_paginate', body).height();
					total -= 10; //padding
					var single = $('tr.odd',table).first().height()+1; //include border
					var num = Math.floor(total/single);
					if(num != lastnum) { //really only do draw calls when needed
						table.DataTable().page.len(num);
						table.DataTable().draw();
					}
				}
			}
		};
		
		//public variables and functions				
		$.fn.DataTable.ext.pager.numbers_length = 5;
		//perform the query
		this.query = function(qobj) {
			//don't need receptor or ligand structures and they are big
			delete qobj.receptor;
			delete qobj.ligand;
			//start provided query
			var postData = {cmd: 'startquery',
					json: JSON.stringify(qobj)
			};
			
			if(qid !== null) postData.oldqid = qid;
			
			$.post(Pharmit.server, postData, null, 'json').done(function(ret) {
				if(ret.status) { //success
					
					//setup table
					qid = ret.qid;
					var numrows = Math.floor((body.height()-75)/28); //magic numbers!
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
							 	info: "Searching..."
						 },
						 serverSide: true,
						 processing: true,
						 ajax: {
						    	url: Pharmit.server,
						    	data: {
						    		cmd: "getdata",
						    		qid: qid
						    	},
						    	dataSrc: processData
						 }

					});
										
					//event handler for loading data, keep loading until done
					table.on('xhr.dt', function(e, settings, json) {
						if(json.finished) {
							var lang = table.DataTable().settings()[0].oLanguage;
							if(json.recordsTotal === 0) {
								lang.emptyTable = lang.sEmptyTable = "No results found";
							} else {
								lang.sInfo = "Showing _START_ to _END_ of _TOTAL_ entries";								
							}
							
						} 
						else {
							//keep polling server
							setTimeout(function() {
								table.DataTable().ajax.reload();
							}, 1000);
						}					 
					});	
					
					
					$('tbody',table).on( 'click', 'tr', function () {
				        if ( $(this).hasClass('selected') ) {
				            $(this).removeClass('selected');
				        }
				        else {
				            table.DataTable().$('tr.selected').removeClass('selected');
				            $(this).addClass('selected');
				        }
				    });
					
				} else {
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
			});

		};
		
		//cancel any query. clear out the table, and hide the div
		//note that quiting is always controlled by Results
		var quit = this.quit = function() {
			
			if(qid !== null) {
				$.post(Pharmit.server, 
						{cmd: 'cancelquery',
						oldqid: qid});
			}
			
			if($.fn.DataTable.isDataTable(table)) {
				table.DataTable().destroy();
			}
		};
		
		//download and save results
		var saveResults = function() {
			
		};
		
		//initiate minimization
		var minimizeResults = function() {
			//hide us, show minresults
			phdiv.hide();
			minresults.minimize(0, function() {
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
			quit();
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
		

		//footer
		var footer = $('<div>').appendTo(phdiv).addClass("pharmit_resfooter");
		//minimize and save buttons
		var bottomloaders = $('<div>').appendTo(footer).addClass("pharmit_bottomloaders").addClass('pharmit_nowrap');

		var minimize = $('<button>Minimize</button>').appendTo(bottomloaders).button().click(minimizeResults);
		var save = $('<button>Save...</button>').appendTo(bottomloaders).button().click(saveResults);		
		
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

	
	function Query(element, viewer, results) {
		//private variables and functions
		var querydiv = $('<div>').addClass('pharmit_query pharmit_overlay');
		var features = null;
		var featureheading = null;
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
			delete qobj.ligand;
			
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
		
		//take an array of pharmacophore features (query.points) and
		//put them in the query view
		var setFeatures = function(featurearray) {
			var start = new Date().getTime();
			
			viewer.disableRendering();
			//while we're removing/adding bunches of features, don't bother rendering until the end
			
			features.detach();
			//replace features
			var old = features.children();
			$.each(old, function(i, fdiv) {
				fdiv.feature.deleteFeature();
			});
			
			features.empty();
			if(featurearray) {
				$.each(featurearray, function(i, pt) {
					new Feature(viewer, features, pt);
				});
			}
			features.accordion("option","active",false);
			features.accordion("refresh");

			featureheading.after(features); 
			
			viewer.enableRendering();
			var end = new Date().getTime();
			var time = end - start;
			console.log('setFeatures time: ' + time);
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
					if(ret.mol) {
						//this was molecular data, save it
						ligandName = lname;
						ligandData = data;
						viewer.setLigand(data, lname);						
					}
					setFeatures(ret.points);					
					
				} else {
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				alert("Error contacting server.  Please inform "+Pharmit.email+ " if this problem persists.");
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
			ret.receptorid = receptorKey;
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
			var buttons = $('<div>').addClass('pharmit_searchdiv');
			var run = $('<button>Search '+vendors[0]+'</button>').appendTo(buttons).button();
			var select = $('<button>Select subset to search</button>').appendTo(buttons).button({text: false, icons: {primary: "ui-icon-triangle-1-s"}});
			
			buttons.buttonset();
			var ul = $('<ul>').appendTo($('body')).addClass('pharmit_floatmenu'); //can't be in query div because of overflow truncation
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
		createSearchButton(header,['MolPort','ZINC']);
		
		//load features and load receptor
		var loaders = $('<div>').appendTo(header).addClass('pharmit_loaderdiv').addClass('pharmit_nowrap');
		var loadrec = $('<button>Load Receptor...</button>').button();
		var loadfeatures = $('<button>Load Features...</button>').button();
		
		//fileinput needs the file inputs in the dom
		element.append(querydiv);
		var loadrecfile = $('<input type="file">').appendTo(loaders).fileinput(loadrec).change(function(e) {readText(this, loadReceptor);});
		var loadfeaturesfile = $('<input type="file">').appendTo(loaders).fileinput(loadfeatures).change(function(e) {readText(this,loadFeatures);});		
		
		querydiv.detach();
		
		//query features
		var body = $('<div>').appendTo(querydiv).addClass("pharmit_querybody");
		var featuregroup = $('<div>').appendTo(body);
		featureheading = $('<div>Pharmacophore</div>').appendTo(featuregroup).addClass('pharmit_heading');
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

		//filters
		var filtergroup = $('<div>').appendTo(body);
		$('<div>Filters</div>').appendTo(filtergroup).addClass('pharmit_heading');
		var filters = $('<div>').appendTo(filtergroup);		
		
		var heading = $('<h3>Hit Reduction<br></h3>').appendTo(filters);
		var hitreductionsummary = $('<span>').appendTo(heading).addClass('pharmit_headingsummary');

		var hitreductions = $('<div>').addClass("pharmit_hitreduction").appendTo(filters);
		var reducetable = $('<table>').appendTo(hitreductions);
		var row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
		$('<td>').append('<label title="Maximum number of orientations returned for each conformation" value="1" for="reduceorienttext">Max Hits per Conf:</label>').appendTo(row);
		var cell = $('<td>').appendTo(row);
		$('<input id="reduceorienttext" name="max-orient">').appendTo(cell).spinner();
		
		row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
		$('<td>').append('<label title="Maximum number of conformations returned for each compound" value="1" for="reduceconfstext">Max Hits per Mol:</label>').appendTo(row);
		cell = $('<td>').appendTo(row);
		$('<input id="reduceconfstext" name="reduceConfs">').appendTo(cell).spinner();
		
		row = $('<tr>').addClass('pharmit_filterrow').appendTo(reducetable);
		$('<td>').append('<label title="Maximum number of hits returned" value="1" for="reducehitstext">Max Total Hits:</label>').appendTo(row);
		cell = $('<td>').appendTo(row);
		$('<input id="reducehitstext" name="max-hits">').appendTo(cell).spinner();
		
		
		heading = $('<h3>Hit Screening<br></h3>').appendTo(filters);
		var hitscreeningsummary = $('<span>').appendTo(heading).addClass('pharmit_headingsummary');
		var hitscreening = $('<div>').appendTo(filters).addClass('pharmit_hitscreening');
		var screentable = $('<table>').appendTo(hitscreening);
		
		row = $('<tr>').addClass('pharmit_filterrow').appendTo(screentable);
		cell = $('<td>').appendTo(row);
		$('<input id="minmolweight" name="minMolWeight">').appendTo(cell).spinner();
		$('<td>').appendTo(row).append($('<label title="Minimum/maximum molecular weight (weights are approximate)" value="1" for="maxmolweight">&le;  MolWeight &le;</label>'));
		cell = $('<td>').appendTo(row);
		$('<input id="maxmolweight" name=maxMolWeight>').appendTo(cell).spinner();

		row = $('<tr>').addClass('pharmit_filterrow').appendTo(screentable);
		cell = $('<td>').appendTo(row);
		$('<input id="minnrot" name="minrotbonds">').appendTo(cell).spinner();
		$('<td>').appendTo(row).append($('<label title="Minimum/maximum number of rotatable bonds" value="1" for="maxnrot"> &le;  RotBonds &le;</label>'));
		cell = $('<td>').appendTo(row);
		$('<input id="maxnrot" name="maxrotbonds">').appendTo(cell).spinner();

		filters.accordion({animate: true, active: false, collapsible: true, heightStyle:'content'});
		
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
			phresults.quit();
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
	}

	return Results;
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
		
		var featureColors = {'Aromatic': 'purple', 'HydrogenDonor': '0xf0f0f0', 'HydrogenAcceptor': 'orange', 
							'Hydrophobic': 'green', 'NegativeIon': 'red', 'PositiveIon': 'blue', 'ExclusionSphere': 'grey'};

		
		var modelsAndStyles = {
				'Ligand': {model: null,
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
				updateStyle(name);
			});
			
			var colorpicker = $('<input name="'+id+'color">').appendTo($('<td>').appendTo(ret));
			colorpicker.val('#C8C8C8');
			colorpicker.change(function() {
				var c = this.value;
				colorpicker.spectrum("set",c);
				var color = parseInt(colorpicker.spectrum("get").toHex(),16); //number
				rec.colorscheme.C = color;
				updateStyle(name);
			});
			
			colorpicker.spectrum({
				color: '#C8C8C8',
			    showPalette: true,
			    replacerClassName: 'ui-state-default ui-corner-all',
			    showPaletteOnly: true,
			    clickoutFiresChange: true,
			    palette: ['#C8C8C8', 'white','green','cyan','magenta','yellow','orange','purple','blue'],
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
			$('<div id="surfaceopacity" name="surfaceopacity">').appendTo(sliderdiv)
				.slider({animate:'fast',step:0.05,'min':0,'max':1,'value':0.8,
					change: function(event, ui) { 
						surfaceStyle.opacity = ui.value;
						if(surface !== null) viewer.setSurfaceMaterialStyle(surface, surfaceStyle);
						viewer.render();
						}
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
						surfaceStyle, 
						{model:receptor}, {bonds:0, invert:true});
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