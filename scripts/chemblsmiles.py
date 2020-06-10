#!/usr/bin/env python3

#why does everyone make it so hard to get smiles
#defaults to getting everything, will eventually support getting
#latest snapshots

#optionally take two arguments: a cid/names file and the name of an output file
#if provided, the cid/names file will be used to output a subset to the specified output

import ftplib,sys,tempfile,gzip,re
from rdkit.Chem import AllChem as Chem

        
ftp = ftplib.FTP('ftp.ebi.ac.uk')
ftp.login()
ftp.cwd('pub/databases/chembl/ChEMBLdb/latest/')
files = ftp.nlst()

for f in files:
    if not f.endswith('.sdf.gz'):
        continue
    #todo - get version and do something with it
        
    #it would be nice to be fancy and stream download and parsing,
    #but for simplicity we will download each file whole and then parse
    #which requires sufficient disk space in /tmp
    temp = tempfile.TemporaryFile(mode='r+b')
    ftp.retrbinary('RETR %s' % f, temp.write)
    temp.seek(0)
    data = gzip.GzipFile(fileobj=temp)
   
    suppl = Chem.ForwardSDMolSupplier(data);
    for mol in suppl:
        if mol is None: 
            continue
        smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
        name = mol.GetProp("_Name")
        print(smiles,name)
    temp.close()
    
ftp.close()
        
