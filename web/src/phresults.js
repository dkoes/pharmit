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

	function PhResults(results, viewer, minresults) {
		//private variables and functions
		var phdiv = null;
		
		//public variables and functions
		
		//perform the query
		this.query = function(qobj) {
			//if we aren't hidden, need to cancel current query first
			
			//start provided query
			
		};
		
		//cancel any query. clear out the table, and hide the div
		//note that quiting is always controlled by Results
		this.quit = function() {
			
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
			//hide our parent
			results.hide();
		});
		var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
		
		//body, should stretch to fill
		var body = $('<div>').appendTo(phdiv).addClass("pharmit_resbody");

		//footer
		var footer = $('<div>').appendTo(phdiv).addClass("pharmit_resfooter");
		//minimize and save buttons
		var bottomloaders = $('<div>').appendTo(footer).addClass("pharmit_bottomloaders").addClass('pharmit_nowrap');

		var minimize = $('<button>Minimize</button>').appendTo(bottomloaders).button().click(minimizeResults);
		var save = $('<button>Save...</button>').appendTo(bottomloaders).button().click(saveResults);		
		
	}

	return PhResults;
})();