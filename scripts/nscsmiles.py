#!/usr/local/bin/python

# extract smiles for DTP/NCI compounds from PubChem with NSC identifiers
# if you are extracting all of pubchem anyway, it might be more efficient
# to use nscids.py and pubchemsmiles.py

import sys, urllib.request, urllib.error, urllib.parse, json
import nscavail

sids = urllib.request.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/DTP.NCI/sids/TXT').read().split('\n')

# I couldn't figure out a way to get the NSC identifier and the isomeric smiels at the same time in one download

for sidline in sids:
    sid = sidline.strip()
    try:
        sidjson = json.load(urllib.request.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%s/JSON' % sid))
        if len(sidjson) > 0:
            record = sidjson['PC_Substances'][0]
            nscid = 'NSC' + record['source']['db']['source_id']['str']
            if nscavail.nscavail(nscid):
                clist = record['compound']
                for cl in clist:
                    if cl['id']['type'] == 1:
                        cid = cl['id']['id']['cid']
                        smile = urllib.request.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/IsomericSMILES/TXT' % cid).read().strip()
                        print(smile, nscid)
    except Exception as e:
        sys.stderr.write('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%s/JSON\n' % sid)        
