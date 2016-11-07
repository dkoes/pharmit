#!/usr/local/bin/python

#get cids and names of nsc compounds
#since this is a large set and I can't do the necessary query
#with pubchem rest, it is faster to scan all of pubchem

import ftplib,sys,tempfile,gzip,re
import nscavail

ftp = ftplib.FTP('ftp.ncbi.nih.gov')
ftp.login()
ftp.cwd('pubchem/Substance/CURRENT-Full/SDF')
files = ftp.nlst()


for f in files:
    if not f.endswith('.sdf.gz'):
        continue
    sys.stderr.write("%s\n"%f)
    #it would be nice to be fancy and stream download and parsing,
    #but for simplicity we will download each file whole and then parse
    #which requires sufficient disk space in /tmp
    temp = tempfile.TemporaryFile(mode='r+b')
    try:
       ftp = ftplib.FTP('ftp.ncbi.nih.gov')
       ftp.login()
       ftp.cwd('pubchem/Substance/CURRENT-Full/SDF')
       ftp.retrbinary('RETR %s' % f, temp.write)
    except:
      sys.stderr.write("Failure with %s\n" % f)
      continue
    temp.seek(0)
    data = gzip.GzipFile(fileobj=temp)
    cid = None
    name = None
    line = data.readline()
    while line:
        #look for cid or smiles data tag and then grab next line
        if re.search(r'PUBCHEM_CID_ASSOCIATIONS',line):
            cids = data.readline().strip().split()
            if len(cids) > 1:
                cid = cids[0]
        elif re.search(r'PUBCHEM_SUBSTANCE_SYNONYM',line):
            while line.strip() != '':
                line = data.readline().strip()
                if line.startswith('NSC'):
                    name = re.sub(r'\s+','',line) #name can't have whitespace
                    #standardize on NSC-[num]
                    if re.match(r'NSC\d+',name):
                        name = re.sub(r'NSC','NSC-',name)
                    break
        elif line.startswith('$$$$'):
            if cid != None and name != None:
                try:
                    if nscavail.nscavail(name):
                        print cid,name
                except:
                    sys.stderr.write("avail fail %s,%s\n" % (cid,name))
            name = None
            cid = None
        line = data.readline()
    temp.close()
    
ftp.close()
        
