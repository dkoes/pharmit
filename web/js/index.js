(function(d){var h=[];d.loadImages=function(a,e){"string"==typeof a&&(a=[a]);for(var f=a.length,g=0,b=0;b<f;b++){var c=document.createElement("img");c.onload=function(){g++;g==f&&d.isFunction(e)&&e()};c.src=a[b];h.push(c)}}})(window.jQuery||window.Zepto);
 $.fn.hasAttr = function(name) { var attr = $(this).attr(name); return typeof attr !== typeof undefined && attr !== false; };

$(document).ready(function() {
r = function() {
$('.img').attr('src', (window.devicePixelRatio > 1) ? ((window.devicePixelRatio > 2) ? 'images/pasted-image-1344.png' : 'images/pasted-image-896.png') : 'images/pasted-image-448.png');};
$(window).resize(r);
r();

//get counts from server
	$.post('/fcgi-bin/pharmitserv.fcgi', {cmd: "getsubsets"}, null, 'json').done(function(ret) {
		$('#numberstandard').text(numeral(ret.standard.length).format('0,0'));
		var numconfs = 0;
		var nummols = 0;
		for(var i = 0; i < ret.standard.length; i++) {
			numconfs += ret.standard[i].numConfs;
			nummols += ret.standard[i].numMols;
		}
		$('#numconfsstandard').text(numeral(numconfs).format('0,0'));
		$('#nummolsstandard').text(numeral(nummols).format('0,0'));
		
		$('#numberpublic').text(numeral(ret.public.length).format('0,0'));
		numconfs = 0;
		nummols = 0;
		for(var i = 0; i < ret.public.length; i++) {
			numconfs += ret.public[i].numConfs;
			nummols += ret.public[i].numMols;
		}
		$('#numconfspublic').text(numeral(numconfs).format('0,0'));
		$('#nummolspublic').text(numeral(nummols).format('0,0'));
		
	}).fail(function() {
		$('#numberstandard').text('X');
		$('#numconfsstandard').text('X');
		$('#nummolsstandard').text('X');
		$('#numberpublic').text('X');
		$('#numconfspublic').text('X');
		$('#nummolspublic').text('X');
		
	});
});