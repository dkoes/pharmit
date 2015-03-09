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
		
		//convert a mol name into something more presentable
		this.mangleName = function(name) {
			var names = name.split(" ");
			var ret = '<span title="' + names.join("\n") +
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