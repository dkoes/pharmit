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
	SearchResults.js
	This is a div for managing pharmacophore/shape results.
	
*/


var Pharmit = Pharmit || {};

Pharmit.SearchResults = (function() {
	// private class variables and functions

	function SearchResults(r, viewer, minresults, type) {
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
			receptor = query.receptor; //save for iteration
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
						order: type == 'shape' ? [[ 1, "desc" ]] : [[ 1, "asc" ]],
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
		
		this.hide = function() {
			phdiv.hide();
		};
		
		this.show = function() {
			phdiv.show();
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
		var headingname = "Pharmacophore Results";
		if(type == "shape") {
			headingname = "Aligned Shape Results";
		}
		var heading = $('<div>'+headingname+'</div>').appendTo(header).addClass('pharmit_heading').addClass("pharmit_rightheading");
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
		
		if(type == 'shape') {
			$('<th>Sim</th>').appendTo(headrow);			
		} else {
			$('<th>RMSD</th>').appendTo(headrow);
		}
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

	return SearchResults;
})();
