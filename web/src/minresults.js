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