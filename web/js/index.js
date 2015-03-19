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
	
	//setup handlers for pdb search
	$('#pdbtext').keyup(function(event) {
		var val = event.target.value;
		if(val.length == 4 || (val.length <= 6 && val[4] === ':')) {
			//get ligand info
			val = val.substring(0,4);
			var sel = $('#pdbligand').empty();
			$('<option disabled selected>...</option>').appendTo(sel);
			$.get('http://www.rcsb.org/pdb/rest/ligandInfo?structureId='+val).done(function(ret) {
				var ligands = $(ret).find('ligand');
				sel.empty();
				$.each(ligands, function(k,v) {
					var lname = $(v).attr('chemicalID');
					$('<option value="'+lname+'">'+lname+"</option>").appendTo(sel);
				});
				$('#pdbligand').prop('disabled',false);
				$('#pdbsubmit').prop('disabled',false);
			});
		} else {
			//assume anything else invalid
			$('#pdbligand').prop('disabled',true);
			$('#pdbsubmit').prop('disabled',true);
		}
		
	});
	
	//return true if pos [x,y,z] has a squared distance less than sqDist
	//to at least one position in coords
	var isclose = function(pos, coords, sqDist) {
		for(var i = 0, n = coords.length; i < n; i++) {
			var a = coords[i];
			var sum = 0;
			for(var j = 0; j < 3; j++) {
				var diff = pos[j]-a[j];
				sum += diff*diff;
			}
			if(sum < sqDist) return true;
		}
		return false;
	}
	
	$('#pdbform').submit(function(event) {
		var pdb = $('#pdbtext').val();
		var ligname = $('#pdbligand').val();
		var water = $('#pdbwater').val();
		var win = window.open("search.html");
		
		var requestedChain = null;
		if(pdb.length > 4) { //presumably has chain
			requestedChain = pdb[5];
			pdb = pdb.substring(0,4);
		}

		//download pdb
		$.get('http://www.pdb.org/pdb/files/'+pdb+".pdb").done(function(mol) {
			var lines = mol.split('\n');
			var reclines = []; //anthing that starts with ATOM
			var waterlines = []; //HOH or WAT resname
			var ligandlines = []; //first occurance only
			var ligchain = null;
			var maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
			var minx = Infinity, miny = Infinity, minz = Infinity;
			var x,y,z, i, n;
			var line;
			var ligcoords = [];
			for(i = 0, n = lines.length; i < n; i++) {
				line = lines[i];
				var resname = line.substring(17,20).trim();
				var chain = line.substring(21,22);
				
				if(requestedChain != null && chain != requestedChain) {
					continue;
				}

				if(resname == ligname) {
					//only retain the ligand(s) on the first chain
					if(ligchain == null) ligchain = chain;
					if(chain == ligchain) {
						ligandlines.push(line);
						//get coordinates
						x = parseFloat(line.substring(30,38));
						y = parseFloat(line.substring(38,46));
						z = parseFloat(line.substring(46,54));
						ligcoords.push([x,y,z]);
						minx = Math.min(x,minx);
						miny = Math.min(y,miny);
						minz = Math.min(z,minz);
						maxx = Math.max(x,maxx);
						maxy = Math.max(y,maxy);
						maxz = Math.max(z,maxz);
					}
				}
				else if(line.lastIndexOf('ATOM ',0) === 0) {
					reclines.push(line);
				} else if(resname == "HOH" || resname == "WAT") {						
					waterlines.push(line);
				}
			}				
			
			if(water != 'ignore') {
				//only consider waters in the binding site
				maxx += 4; maxy +=4; maxz += 4;
				minx -= 4; miny -= 4; minz -= 4;
				var bindingwaters = [];
				for(i = 0, n = waterlines.length; i < n; i++) {
					line = waterlines[i];
					x = parseFloat(line.substring(30,38));
					y = parseFloat(line.substring(38,46));
					z = parseFloat(line.substring(46,54));
					if(x > minx && x < maxx && y > miny && y < maxy && z > minz && z < maxz) {
						//compare to ligcoords
						if(isclose([x,y,z], ligcoords, 16)) {
							bindingwaters.push(line);
						}
					}
				}
				
				if(water == 'rec') {
					$.merge(reclines, bindingwaters);
				} else {
					$.merge(ligandlines, bindingwaters);
				}
			}
			
			var receptor = reclines.join('\n');
			var ligand = ligandlines.join('\n');
			
			if(ligand.length == 0) {
				win.close();
				alert("Unabled to extract ligand.  Check chain identifier.");
			} else {
				var data = {ligand: ligand, ligandFormat: ligname+".pdb", receptor: receptor, recname: pdb+".pdb"};
				var msg = new Message(JSON.stringify(data), win, '*');
			}
			
		}).fail(function() {
			alert("Could not fetch pdb "+pdb);
		});
	});
});
