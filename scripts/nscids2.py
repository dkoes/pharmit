#!/usr/local/bin/python

#get cids and names of nsc compounds
#
#since this is a large set and I can't do the necessary query
#with pubchem rest, it is faster to scan all of pubchem

import json,sys,tempfile,gzip,re,urllib2
import nscavail
from Bio import Entrez

Entrez.email = "dkoes@pitt.edu"
records = Entrez.read(Entrez.esearch(db='pcsubstance',term="NSC"))

total = int(records['Count'])
batchsize = 100000

sids = []
for start in xrange(0,total,batchsize):
	records = Entrez.read(Entrez.esearch(db='pcsubstance',term="NSC",retmax=batchsize,retstart=start))
	for sid in records['IdList']:
		sids.append(int(sid))

for sid in sids:
    try:
        #query each one individually! hopefully more robust..
        record = urllib2.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%d/json' % sid).read()
        record = json.loads(record)
        record = record['PC_Substances'][0]
        name = ''
        for syn in record['synonyms']:
            if syn.startswith('NSC'):
                name = syn
                break
        #standardize on NSC-[num]
        if re.match(r'NSC\d+',name):
            name = re.sub(r'NSC','NSC-',name)  
        if re.match(r'NSC\s+\d+',name):
            name = re.sub(r'NSC\s+','NSC-',name)  
        if name != '' and 'compound' in record:
            #get cmpd id
            cid = 0
            for rec in record['compound']:
                if 'id' in rec and rec['id']['type'] == 1 and 'id' in rec['id'] and 'cid' in rec['id']['id']:
                    cid = rec['id']['id']['cid']
                    break
        if cid != 0 and nscavail.nscavail(name):
            print cid,name
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        sys.stderr.write("Error with sid %d\n"%sid)
        sys.stderr.write(str(e))
                
