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
		var save = null;
		var onclose = null; //what to call when close button is clicked
		
		var processData = function(data) {
			
			if(data.status === 0)
				return {error: data.msg};
			
			var ret = data.data;
			
			for(var i = 0; i < ret.length; i++) {
				//round floats
				ret[i][3] = numeral(ret[i][3]).format('0.00');
				ret[i][4] = numeral(ret[i][4]).format('0.000');
			}
			return ret;
		};
		
		//public variables and functions
		var cancel = this.cancel = function() {
			
			if($.fn.DataTable.isDataTable(table)) {
				table.DataTable().destroy();
			}
			
			save.button( "option", "disabled", true );
			mindiv.hide();
		};
		
		//perform the query
		this.minimize = function(qid, qobj, closer) {
			onclose = closer;
			
			var postData = {cmd: 'startsmina',
					qid: qid,
					receptorid: qobj.receptorid,
					recname: qobj.recname
			};
			
			//start provided query
			$.post(Pharmit.server, postData, null, 'json').done(function(ret) {
				if(ret.status) { //success
					
					//setup table
					sminaid = ret.sminaid;
					var numrows = Math.floor((body.height()-85)/28); //magic numbers!
					table.dataTable({
						searching: false,
						pageLength: numrows,
						destroy: true, //replace any existing table
						lengthChange: false,
						order: [[ 1, "asc" ]],
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
							 emptyTable: "Nada, zip, zilcho",
							 	infoFiltered: '',
							 	infoEmpty: "",
							 	info: "Minimizing..."
						 },
						 serverSide: true,
						 processing: true,
						 ajax: {
						    	url: Pharmit.server,
						    	data: {
						    		cmd: "getsminadata",
						    		qid: qid
						    	},
						    	dataSrc: processData
						 }

					});
										
					//event handler for loading data, keep loading until done
					table.on('xhr.dt', function(e, settings, json) {
						if(json.finished) {
							save.button( "option", "disabled", false );														
						} 
				        else if(json.status === 0) {
				        	alert(json.msg);
				        }
						else {
				            viewer.setResult(); //clear in case clicked on

							//keep polling server
							setTimeout(function() {
								table.DataTable().ajax.reload();
							}, 1000);
						}					 
					});	
					
					
					$('tbody',table).on( 'click', 'tr', function () {
						var mid = table.DataTable().row(this).data()[0];
						var r = this;
				        if ( $(this).hasClass('selected') ) {
				            $(this).removeClass('selected');
				            viewer.setResult(); //clear
				        }
				        else {
				            table.DataTable().$('tr.selected').removeClass('selected');
				            $(this).addClass('selected');
				            
				            $.post(Pharmit.server,
				            		{cmd: 'getsminamol',
				            		 qid: qid,
				            		 molid: mid
				            		}).done(function(ret) {
				            			if( $(r).hasClass('selected')) //still selected
				            				viewer.setResult(ret);
				            		});
				        }
				    });
				} else {
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				alert("Error contacting minimization server.  Please inform "+Pharmit.email+ " if this problem persists.");
			});	
			
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
		save = $('<button>Save...</button>').appendTo(bottomrow).button().click(saveResults);				

		mindiv.hide();
	}

	return MinResults;
})();