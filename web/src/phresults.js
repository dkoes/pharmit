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