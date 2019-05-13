#!/usr/local/bin/python

#extract smiles from PubChem (not sure why this isn't easier)
#defaults to getting everything, will eventually support getting
#latest snapshots

#optionally take two arguments: a cid/names file and the name of an output file
#if provided, the cid/names file will be used to output a subset to the specified output

import ftplib,sys,tempfile,gzip,re

nscids = {}
nscout = None
if len(sys.argv) == 3:
    #providing an nsc ids files (cid,nscname) and output 
    nscf = open(sys.argv[1])
    nscout = open(sys.argv[2],'w')
    for line in nscf:
        (cid,name) = line.strip().split()
        nscids[cid] = name
        
ftp = ftplib.FTP('ftp.ncbi.nih.gov')
ftp.login()
ftp.cwd('pubchem/Compound/CURRENT-Full/SDF')
files = ftp.nlst()

for f in files:
    if not f.endswith('.sdf.gz'):
        continue
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
            print('%s\tPubChem-%s'% (smile,cid)) 
            if cid in nscids:
                nscout.write('%s\t%s' % (smile,nscids[cid]))
                sys.stderr.write('%s\t%s' % (smile,nscids[cid]))
            smile = None
            cid = None
        line = data.readline()
    temp.close()
    
ftp.close()
        
