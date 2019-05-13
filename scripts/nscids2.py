#!/usr/local/bin/python

#get cids and names of nsc compounds
#
#since this is a large set and I can't do the necessary query
#with pubchem rest, it is faster to scan all of pubchem

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

for sid in sids:
    try:
        #query each one individually! hopefully more robust..
        sdf = urllib.request.urlopen('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/%d/sdf' % sid).read()        
        mol = next(Chem.ForwardSDMolSupplier(io.StringIO(sdf)))
        if mol.GetNumAtoms() == 0:
            continue
        name = ''
        for syn in mol.GetProp('PUBCHEM_SUBSTANCE_SYNONYM').split():
            if syn.startswith('NSC'):
                name = syn
                break
        #standardize on NSC[num]
        if re.match(r'NSC-\d+',name):
            name = re.sub(r'NSC-','NSC',name)  
        if re.match(r'NSC\s+\d+',name):
            name = re.sub(r'NSC\s+','NSC',name)  
        if nscavail.nscavail(name):
            print(Chem.MolToSmiles(mol),name)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        sys.stderr.write("Error with sid %d\n"%sid)
        sys.stderr.write(str(e))
                
