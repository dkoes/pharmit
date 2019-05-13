#!/usr/local/bin/python

#get names a smi of putative NSC compounds
#this doesn't query nsc, since that kept getting stuck

import json,sys,tempfile,gzip,re,urllib.request,urllib.error,urllib.parse
import nscavail
from Bio import Entrez
from rdkit.Chem import AllChem as Chem
import io

Entrez.email = "dkoes@pitt.edu"
records = Entrez.read(Entrez.esearch(db='pcsubstance',term="NSC"))

total = int(records['Count'])
batchsize = 100000

sids = []
for start in range(0,total,batchsize):
	records = Entrez.read(Entrez.esearch(db='pcsubstance',term="NSC",retmax=batchsize,retstart=start))
	for sid in records['IdList']:
		sids.append(int(sid))

sys.stderr.write("Found %d sids\n"%len(sids))
for sid in sids:
    try:
        #query each one individually! hopefully more robust..
        sys.stderr.write("sid %d\n"%sid)
        sdf = urllib.request.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%d/sdf' % sid,timeout=30).read()        
        mol = next(Chem.ForwardSDMolSupplier(io.StringIO(sdf)))
        if mol.GetNumAtoms() == 0:
            continue
        name = ''
        for syn in mol.GetProp('PUBCHEM_SUBSTANCE_SYNONYM').split('\n'):
            m = re.match(r'^NSC.*?(\d+)',syn,re.IGNORECASE)
            if m:
                name = 'NSC%s' % m.group(1)
                #standardize on NSC[num]
                print(Chem.MolToSmiles(mol),name)
                sys.stderr.write(name+"\n")
                sys.stdout.flush()
                break

    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        sys.stderr.write("Error with sid %d\n"%sid)
        sys.stderr.write(str(e))
                
