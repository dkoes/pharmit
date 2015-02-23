#!/usr/local/bin/python

#extract smiles from PubChem (not sure why this isn't easier)
#defaults to getting everything, will eventually supprot getting
#latest snapshots

import ftplib,sys,tempfile,gzip,re

ftp = ftplib.FTP('ftp.ncbi.nih.gov')
ftp.login()
ftp.cwd('pubchem/Compound/CURRENT-Full/SDF')
files = ftp.nlst()

for f in files:
    #it would be nice to be fancy and stream download and parsing,
    #but for simplicity we will download each file whole and then parse
    #which requires sufficient disk space in /tmp
    temp = tempfile.TemporaryFile(mode='r+b')
    ftp.retrbinary('RETR %s' % f, temp.write)
    temp.seek(0)
    data = gzip.GzipFile(fileobj=temp)
    cid = None
    smile = None
    line = data.readline()
    while line:
        #look for cid or smiles data tag and then grab next line
        if re.search(r'PUBCHEM_COMPOUND_CID',line):
            cid = data.readline().strip()
        elif re.search(r'PUBCHEM_OPENEYE_ISO_SMILES',line):
            smile = data.readline().strip()
        elif line.startswith('$$$$'):
            print '%s\tPubChem-%s'% (smile,cid) 
            smile = None
            cid = None
        line = data.readline()
        