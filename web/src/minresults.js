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
		var qid = null;
		var save = null;
		var onclose = null; //what to call when close button is clicked
		var table =null; //results table
		var maxscore = null; //filters
		var maxRMSD = null; 
		var singleConfs = null;
		var startTotal = 0;
		var timeout = null;
		var query = null;
		
		var processData = function(data) {
			
			if(data.status === 0)
				return {error: data.msg};
			
			var ret = data.data;
			
			for(var i = 0; i < ret.length; i++) {
				//round floats
				ret[i][3] = numeral(ret[i][3]).format('0.00');
				ret[i][4] = numeral(ret[i][4]).format('0.000');
				ret[i][2] = results.mangleName(ret[i][2]);

			}
			return ret;
		};
		
		//construct objects for fetching smina data
		var ajaxData = function(postData, settings) {
			
			postData.cmd = "getsminadata";
			postData.qid = qid;
	    	
            //add filters if present
            if(maxRMSD.val() !== '') {
            	postData.maxRMSD = maxRMSD.val();
            }
            if(maxscore.val() !== '') {
            	postData.maxScore = maxscore.val();
            }
            if(singleConfs.prop('checked')) {
            	postData.unique = 1;
            }
            return postData;
		};
		
		
		//public variables and functions
		var cancel = this.cancel = function() {
			
			clearTimeout(timeout);
			if($.fn.DataTable.isDataTable(table) && qid > 0) {
				table.DataTable().clear();
				$.post(Pharmit.server, {
					cmd: 'cancelsmina',
					qid: qid
				});				
			}			
			qid = null;
			viewer.setResult();
			save.button( "option", "disabled", true );
			mindiv.hide();
		};
		
		//perform the query
		this.minimize = function(q, qobj, nummols, closer) {
			startTotal = nummols;
			onclose = closer;
			qid = q;
			query = qobj;
			var postData = {cmd: 'startsmina',
					qid: qid,
					receptorid: qobj.receptorid,
					recname: qobj.recname
			};
			
			table.hide(); //don't show it until it is setup
			//start provided query
			$.post(Pharmit.server, postData, null, 'json').done(function(ret) {
				if(ret.status) { //success
					table.show();
					//setup table
					sminaid = ret.sminaid;
					var numrows = Math.floor((body.height()-95)/28); //magic numbers!
					
					table.dataTable({
						searching: false,
						pageLength: numrows,
						destroy: true, //replace any existing table
						lengthChange: false,
						order: [[ 3, "asc" ]],
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
							 	emptyTable: "Minimizing...",
							 	infoEmpty: " ",
							 	infoFiltered: "<br>(filtered from _MAX_ total hits)",
							 	info: "Minimizing..."
						 },
						 serverSide: true,
						 processing: false,
						 ajax: {
						    	url: Pharmit.server,
						    	data: ajaxData,
						    	dataSrc: processData
						 }

					});


				} else {
					cancel();
					if(onclose) onclose();
					alert("Error: "+ret.msg);
				}
			}).fail(function() {
				cancel();
				if(onclose) onclose();
				alert("Error contacting minimization server.  Please inform "+Pharmit.email+ " if this problem persists.");
			});	
			
			//show div
			mindiv.show();
		};
		
		//download results
		var saveResults = function() {
			
			Pharmit.inFormSubmit = true;			
			//have to use stupid form trick - mostly because of IE and safari
			var cmd = Pharmit.server+'?cmd=savesmina&qid='+qid;
			var form = $('<form>', { 'action': cmd, 'method': 'post'});
			form.appendTo(document.body);
			form.submit();
			$(form).remove();		
		};
		
		
		//initialization code
		mindiv = $('<div>').appendTo(results.div).addClass('pharmit_rescontainer pharmit_mincontainer');
		//header
		var header = $('<div>').appendTo(mindiv).addClass("pharmit_resheader");
		var title = $('<div>Minimization Results</div>').appendTo(header).addClass('pharmit_heading').addClass("pharmit_rightheading");
		var closediv = $('<div>').addClass('pharmit_resclose').appendTo(title).click(function() {
			//cancel the current query 
			cancel();
			//do what our caller told us to do on close
			if(onclose) onclose();
		});
		var close = $('<span>').addClass('ui-icon-circle-close ui-icon').appendTo(closediv);
		
		
		//body, should stretch to fill
		var body = $('<div>').appendTo(mindiv).addClass("pharmit_resbody");

		//skeleton of datatable
		table = $('<table width="100%" class="display compact" cellspacing="0">').addClass('pharmit_mintable').appendTo(body);
		var headrow = $('<tr>').appendTo($('<thead>').appendTo(table));
		$('<th>mid</th>').appendTo(headrow);
		$('<th>origpos</th>').appendTo(headrow);
		$('<th>Name</th>').appendTo(headrow);
		$('<th title="Computed using the AutoDock Vina scoring function">Score</th>').appendTo(headrow);
		$('<th title="The minimized RMSD, the heavy atom difference between the pharmacophore aligned pose and the minimized pose.  Small values imply a closer agreement to the query pharmacophore.">mRMSD</th>').appendTo(headrow);
		$('<tbody>').appendTo(table);
		
		//event handler for loading data, only install this once
		table.on('xhr.dt', function(e, settings, json) {
			var lang = table.DataTable().settings()[0].oLanguage;

			if(json.finished) {
				save.button( "option", "disabled", false );			
				lang.sInfo = "Showing _START_ to _END_ of _TOTAL_ entries";
			} 
	        else if(json.status === 0) {
	        	alert(json.msg);
	        }
			else {
	            viewer.setResult(); //clear in case clicked on

	            lang.sInfo = numeral(json.recordsTotal).format('0,0') + "/" + numeral(startTotal).format('0,0');
	            
	            clearTimeout(timeout); //push next poll farther into the future
	            timeout = setTimeout(function() {
					if($.fn.DataTable.isDataTable(table) && qid > 0) {
						table.DataTable().ajax.reload();
					}
				}, 1000);
			}					 
		});	
		
		table.on('draw.dt', function() {
			$('.pharmit_namecol span').tooltip({position:{my: 'left-5 top+5', at: 'left bottom', collision: 'flipfit'}});
		});
		
		$('tbody',table).on( 'click', 'tr', function () {
			var mid = table.DataTable().row(this).data()[0];
			var r = this;
			$(".pharmit_iterate_button").remove();

	        if ( $(this).hasClass('selected') ) {
	            $(this).removeClass('selected');
	            viewer.setResult(); //clear
	        }
	        else {
	            table.DataTable().$('tr.selected').removeClass('selected');
	            $(this).addClass('selected');
	            
	            $.post(Pharmit.server, {cmd: 'getsminamol',
            		 qid: qid,
            		 molid: mid
            		}).done(function(ret) {
	            			if( $(r).hasClass('selected')) //still selected
	            				viewer.setResult(ret);
	            				var ibutton = $('<div class="pharmit_iterate_button" title="Start new pharmit session around selected ligand">').appendTo($('td',r).last());
	            				ibutton.button({ icons: {primary: "ui-icon-arrowthickstop-1-e"}, text: false});
								ibutton.tooltip({show: {delay: 500}});
	            				ibutton.click(function(event) {
	            					event.stopPropagation();
	            					//create new window around this molecule
	            					var win = window.open("search.html");
	            					var data = {ligand: ret, ligandFormat: mid+".sdf", receptor: query.receptor, recname: query.recname};
	            					var msg = new Message(JSON.stringify(data), win, '*');
	            				});	            			
	            		});
	        }
	    });
		
		//footer
		var footer = $('<div>').appendTo(mindiv).addClass("pharmit_resfooter");
		var paramdiv = $('<div>').appendTo(footer).addClass("pharmit_minparams").disableSelection();
		//filters for minimization
		var filtertable = $('<table>').appendTo(paramdiv);
		var row = $('<tr>').appendTo(filtertable).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max Score</label>'));
		var cell = $('<td>').appendTo(row);
		maxscore = $('<input name="sminamaxscore">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		singleConfs = $('<input type="checkbox" name="sminaunique">');
		$('<td>').appendTo(row).addClass('pharmit_checkcell').append(singleConfs);
		$('<td>').appendTo(row).append($('<label>Single<br>conformer</label>'));
		
		
		row = $('<tr>').appendTo(filtertable).addClass('pharmit_paramrow');
		$('<td>').appendTo(row).append($('<label>Max mRMSD</label>'));
		cell = $('<td>').appendTo(row);
		maxRMSD = $('<input name="sminamaxrmsd">').appendTo(cell).addClass('pharmit_sminainput').spinner();
		cell = $('<td colspan=2>').appendTo(row).addClass('pharmit_applycell');
		$('<button>Apply</button>').appendTo(cell).button().click(function() {table.DataTable().ajax.reload();});
		
		//save button
		var bottomrow = $('<div>').appendTo(footer).addClass('pharmit_minbottom');
		save = $('<button>Save...</button>').appendTo(bottomrow).button().click(saveResults);				

		mindiv.hide();
		
		var lastheight = 0;
		var lastnum = 0;
		$(window).resize(function() {
			if(qid > 0) {
				var total = body.height();
				if(total != lastheight) {
					lastheight = total;
					total -= $('thead', table).height();
					total -= $('.dataTables_info', body).height();
					total -= $('.dataTables_paginate', body).height();
					total -= 20; //padding
					var single = $('tr.odd',table).first().height()+1; //include border
					var num = Math.floor(total/single);
					if(num != lastnum) { //really only do draw calls when needed
						table.DataTable().page.len(num);
						table.DataTable().draw();
						lastnum = num;
					}
				}
			}
		});
	}

	return MinResults;
})();
