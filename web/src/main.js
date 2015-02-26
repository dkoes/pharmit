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

var Pharmit = {};
$(document).ready(function() {
	
	Pharmit.server = '/fcgi-bin/pharmitserv.fcgi';
	var element = $('#pharmit');
	var viewer = new Pharmit.Viewer(element);
	var phresults = new Pharmit.PharmaResults(element, viewer);
	var query = new Pharmit.Query(element, viewer);
		
	//work around jquery bug
	$("button, input[type='button'], input[type='submit']").button()
    .bind('mouseup', function() {
        $(this).blur();     // prevent jquery ui button from remaining in the active state
    });	
});