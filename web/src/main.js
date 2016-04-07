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
	Pharmit.email = 'dkoes@pitt.edu';
	
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
		else if(event.data.type && event.data.type == "init") {
			//actually have no idea where this event comes from, but it pops up in chrome
		}
		else if(event.contexts) {
			//same here
		}
		else {
			try {
				var obj = $.parseJSON(event.data);
				query.setLigandAndReceptor(obj.ligand, obj.ligandFormat, obj.receptor, obj.recname);
			}
			catch(e) {
				console.log("Communication error: "+e);
			}
		}
	};
	window.addEventListener("message", receiveMessage);

});
