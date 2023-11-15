#!/usr/bin/env python3
# A continually running loop that builds libraries from user submitted information.
#It is assumed only one instance of buildlibs is running (if not, need to fix race condition getting pending submission).

#Needs to be able to call createconfs.py and pharmitserver 

import time, sys, os, MySQLdb, subprocess, json, psutil, signal, gzip, traceback, io
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser

#set an error message for the library specified by row
def setError(conn, which, msg):
    conn.ping(True)
    c = conn.cursor()
    c.execute("UPDATE `databases` SET status=%s, message=%s WHERE id = %s",("Error", msg, which))

def setStatus(conn, which, status, msg):
    conn.ping(True)
    c = conn.cursor()
    c.execute("UPDATE `databases` SET status=%s, message=%s WHERE id = %s",(status, msg, which))

#send a signal to any running pharmitserver processes indicating that they should reload database info
def reset_server():
    for proc in psutil.process_iter():
        if proc.name() == 'pharmitserver' or proc.name() == 'pharmit':
            proc.send_signal(signal.SIGUSR1)
           
class BytesIOWrapper(io.RawIOBase):
    '''Because some people like to upload non-ASCII files'''
    def __init__(self, file, encoding='utf-8', errors='strict'):
        self.file, self.encoding, self.errors = file, encoding, errors
        self.buf = b''
    def readinto(self, buf):
        if not self.buf:
            self.buf = self.file.read(4096).encode(self.encoding, self.errors)
            if not self.buf:
                return 0
        length = min(len(buf), len(self.buf))
        buf[:length] = self.buf[:length]
        self.buf = self.buf[length:]
        return length
    def readable():
        return True

def create_sdf_ligs(conn, libraryid, cprefixes):
    #break up molecules into individual files, assume molecules with the same name
    #are the same conformer; return true if successfull
    #outputs an ligs.in file in the current directory
    whichprefix = -1
    lastsmi = ''
    fname = False
    confconn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")

    try:
        ligout = open("ligs.in", 'wt')        
        #throw out non-ASCII characters
        infile = BytesIOWrapper(gzip.open('input.sdf.gz','rt',encoding='utf-8',errors='ignore'))
        mols = Chem.ForwardSDMolSupplier(infile)
        molcnt = 0
        for mol in mols:
            try:
                if mol is None: continue
                Chem.SanitizeMol(mol)
                can = Chem.MolToSmiles(mol,isomericSmiles=True)
                
                if can != lastsmi: #a different molecule
                    if len(can) > 250: #way too big
                        continue
                    molcnt += 1
                    #get/assign a unique id
                    cursor = confconn.cursor()
                    cursor.execute('SELECT id FROM structures WHERE smile = %s', (can,))
                    row = cursor.fetchone()
                    if row == None:
                        #insert without sdfs to get unique id 
                        cursor.execute('INSERT INTO structures (smile,weight) VALUES(%s,%s) ', (can, Chem.CalcExactMolWt(mol)))
                        cursor.execute('SELECT id FROM structures WHERE smile = %s', (can,))
                        row = cursor.fetchone()
    
                    uniqueid = row[0]
                    #we do not store the user supplied conformers, but leave sdfloc blank
                    if fname and not out.closed:
                        writer.close() #we have a file we have previously opened
                        out.close()
                        fname = False
                    whichprefix = (whichprefix+1)%len(cprefixes)
                    subdir = "%s/user/%s/" % (cprefixes[whichprefix],libraryid)
                    if not os.path.isdir(subdir):
                        os.makedirs(subdir)
                    fname = "%s/%d.sdf.gz" % (subdir,uniqueid)
                    out = gzip.open(fname, 'wt')
                    writer = Chem.SDWriter(out)
                    
                    if mol.HasProp('_Name'):
                        name = mol.GetProp('_Name')
                    else:
                        name = str(molcnt)
                    
                    ligout.write('%s %d %s\n' % (fname,uniqueid, name))
                    
                #have file setup for this molecule, may be conformer
                writer.write(mol)
                
                        
            except: #catch rdkit issues
                traceback.print_exc()
                if fname and not out.closed:
                    fname = False
                    writer.close()
                    out.close()
                continue
        
        if fname and not out.closed:
            writer.close()
            out.close()
        return True
    except:
        traceback.print_exc()
        return False;
            
            
def make_libraries(conn, dbprefixfile,row,numactives=0):
      #create database info - assumes current directory has ligs.in in it
    which = row['id']
    setStatus(conn, which, "MakeLib","Creating search index")

    #make dbinfo from database row
    dbinfo = dict()
    for (k,v) in row.items():
        dbinfo[k] = str(v)
    dbinfo['fromuser'] = True #indicate this is a contributed library
    if int(row['isprivate']):
        dbinfo['subdir'] = 'Private/%s' % which
    else:
        dbinfo['subdir'] = 'Public/%s' % which
    if numactives > 0:
        dbinfo['numactives'] = numactives
    jsonfile = "dbinfo.json"
    f = open(jsonfile,'wt')
    f.write(json.dumps(dbinfo,indent=4))
    f.close()
    #build libraries
    print(os.getcwd())
    cmd = "%s dbcreateserverdir -ligs %s -prefixes %s -dbinfo %s > pharmit.out" %(options.pharmit, "ligs.in", dbprefixfile, jsonfile)
    print(cmd)
    ret = subprocess.call(cmd,shell=True)
    if ret != 0:
        print("Problem generating databases",cmd)
        setError(conn,which, "Problem generating databases")
        return
    
    reset_server()
    #get number of molecules
    numConfs = 0
    numMols = 0
    for line in open(dbprefixfile):        
        dbinfofile = '%s/%s/dbinfo.json' % (line.strip(), dbinfo['subdir'])
        info = json.loads(open(dbinfofile).read())
        numConfs += info["numConfs"]
        numMols += info["numMols"]
    conn.ping(True)
    c = conn.cursor()
    c.execute("UPDATE `databases` SET status=%s, message=%s, completed=NOW(), nummols=%s, numconfs=%s WHERE id = %s",("Completed", "Library successfully created", numMols, numConfs ,which))        
    
    
    
if __name__ == '__main__':
    
    parser = OptionParser(usage="Usage: %prog [options]")
    parser.add_option("--maxconfs", dest="maxconfs",action="store",
                      help="maximum number of conformers to generate per a molecule (default 20)", default="20", type="int", metavar="CNT")
    parser.add_option('--confprefixes', dest="confprefixfile", action="store", 
            help="file containing path prefixes for storing conformers", default="",metavar="FILE")
    parser.add_option('--dbprefixes', dest="dbprefixfile", action="store", 
            help="file containing path prefixes for storing databases", default="",metavar="FILE")    
    parser.add_option('--pharmit', dest="pharmit", action="store", 
            help="pharmit executable for creating libraries", default="pharmitserver",metavar="FILE")
    parser.add_option('--createconfs', dest="createconfs", action="store", 
            help="createconfs.py script for generating conformers", default="createconfs.py",metavar="FILE")
        
    (options, args) = parser.parse_args()

    cprefixfile = os.path.abspath(options.confprefixfile)
    if not cprefixfile or not os.path.isfile(cprefixfile):
        print("Require prefix file for storing structures")
        sys.exit(-1)
    cprefixes = open(cprefixfile).read().splitlines()
        
    dbprefixfile = os.path.abspath(options.dbprefixfile)
    if not dbprefixfile or not os.path.isfile(dbprefixfile):
        print("Require prefix file for storing databases")
        sys.exit(-1)        
    
    #look for pending libraries
    while True:
        try:
            conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="pharmit")
            conn.autocommit(True)
            conn.query("SELECT * FROM `databases` WHERE status = 'Pending' ORDER BY submitted ASC")
            rows = conn.store_result().fetch_row(how=1,maxrows=0)
            for row in rows:
                #validate values
                dir = row["directory"]
                which = row["id"]
                if not os.path.isdir(dir):
                    setError(conn, which, "Problem reading input directory")
                if os.path.exists(dir+"/input.sdf.gz"): #sdf
                    os.chdir(dir)
                    setStatus(conn, which, "GenConf","Processing conformers")
                    if not create_sdf_ligs(conn, which, cprefixes):
                        setError(conn, which, "Problem reading sdf file")
                        continue
                    make_libraries(conn, dbprefixfile, row)
                elif os.path.exists(dir+"/input.smi"):
                    os.chdir(dir)
                    setStatus(conn, which, "GenConf","Generating conformers")
                    #check for a "benchmark" set with actives marked
                    f = open(dir+"/input.smi")
                    numactives = 0
                    for line in f:
                        if "active" in line:
                            numactives += 1
                    cmd = "%s --nonames --prefixes %s input.smi > ligs.in 2> conf.err" % (options.createconfs,cprefixfile)
                    print(os.getcwd())
                    print(cmd)
                    ret = subprocess.call(cmd,shell=True)
                    if ret != 0:
                        setError(conn,which,"Problem generating conformers")
                        continue
                    
                    make_libraries(conn, dbprefixfile, row, numactives)                    

                else:
                    setError(conn, row, "Problem reading input files")
            
            conn.close()
            time.sleep(1)
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()
            time.sleep(1)
