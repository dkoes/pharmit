#!/usr/bin/env python
# A continually running loop that builds libraries from user submitted information.
#It is assumed only one instance of buildlibs is running (if not, need to fix race condition getting pending submission).

#Needs to be able to call createconfs.py and pharmitserver 

import time, sys, os, MySQLdb, subprocess, json, psutil, signal
from rdkit import Chem
from optparse import OptionParser

#set an error message for the library specified by row
def setError(conn, which, msg):
    conn.ping(True)
    conn.query("UPDATE `databases` SET status=%s, message=%s WHERE id = %s",("Error", msg, which))

def setStatus(conn, which, status, msg):
    conn.ping(True)
    conn.query("UPDATE `databases` SET status=%s, message=%s WHERE id = %s",(status, msg, which))

#send a signal to any running pharmitserver processes indicating that they should reload database info
def reset_server():
    for proc in psutil.process_iter():
        if proc.name == 'pharmitserver' or proc.name == 'pharmit':
            proc.send_signal(signal.SIGUSR1)
            
def create_sdf_ligs(conn, libraryid, cprefixes):
    #break up molecules into individual files, assume molecules with the same name
    #are the same conformer; return true if successfull
    #outputs an in.ligs file in the current directory
    whichprefix = -1
    lastsmi = ''
    fname = False
    try:
        ligout = open("ligs.in", 'w')        
        infile = gzip.open('input.sdf.gz')
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
                    cursor = conn.cursor()
                    cursor.execute('SELECT id FROM structures WHERE smile = %s', (can,))
                    row = cursor.fetchone()
                    if row == None:
                        #insert without sdfs to get unique id 
                        cursor.execute('INSERT INTO structures (smile,weight) VALUES(%s,%s) ', (can, Chem.CalcExactMolWt(mol)))
                        cursor.execute('SELECT id FROM structures WHERE smile = %s', (can,))
                        row = cursor.fetchone()
    
                    uniqueid = row[0]
                    #we do not store the user supplied conformers, but leave sdfloc blank
                    if fname:
                        writer.close() #we have a file we have previously opened
                        out.close()
                    whichprefix = (whichprefix+1)%len(cprefixes)
                    subdir = "%s/user/%s/" % (cprefixes[whichprefix],libraryid)
                    if not os.path.isdir(subdir):
                        os.makedirs(subdir)
                    fname = "%s/%d.sdf.gz" % (subdir,uniqueid)
                    out = gzip.open(fname, 'w')
                    writer = Chem.SDWriter(out)
                    
                    if mol.HasProp('_Name'):
                        name = mol.GetProp('_Name')
                    else:
                        name = str(molcnt)
                    
                    ligout.write('%s %d %s\n' % (fname,uniqueid, name))
                    
                #have file setup for this molecule, may be conformer
                writer.write(mol)
                
                        
            except: #catch rdkit issues
                continue
        
        if fname:
            writer.close()
            out.close()
        return True
    except:
        print sys.exc_info()
        return False;
            
            
def make_libraries(conn, dbprefixes,row):
      #create database info - assumes current directory has in.ligs in it
    which = row['id']
    row['fromuser'] = True #indicate this is a contributed library
    if int(row['isprivate']):
        row['subdir'] = 'Private/%s' % which
    else:
        row['subdir'] = 'Public/%s' % which
    jsonfile = "dbinfo.json"
    f = open(jsonfile,'w')
    f.write(json.dumps(row,indent=4))
    f.close()
    #build libraries
    ret = subprocess.call("%s dbcreateserverdir -ligs %s -prefixes %s -dbinfo %s > pharmit.out" %(options.pharmit, "in.ligs", dbprefixes, jsonfile))
    if ret != 0:
        setError(conn,which, "Problem generating databases")
        continue
    
    reset_server()
    conn.ping(True)
    conn.query("UPDATE `databases` SET status=%s, message=%s, completed=NOW() WHERE id = %s",("Completed", "Library successfully created", which))        
    
    
    
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
        
    cprefixes = os.path.abspath(options.confprefixfile)
    if not cprefixes or not os.path.isfile(cprefixes):
        print "Require prefix file for storing structures"
        sys.exit(-1)
        
    dbprefixes = os.path.abspath(options.dbprefixfile)
    if not dbprefixes or not os.path.isfile(dbprefixes):
        print "Require prefix file for storing databases"
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
                    setStatus(conn, which, "GenConf","Processing conformers")
                    if not create_sdf_ligs(conn, which, cprefixes):
                        setError(conn, which, "Problem reading sdf file")
                        continue
                    make_libraries(conn, dbprefixes, row)
                elif os.path.exists(dir+"/input.smi"):
                    os.chdir(dir)
                    setStatus(conn, which, "GenConf","Generating conformers")
                    ret = subprocess.call("%s --nonames --prefixes %s input.smi > in.ligs" % (options.createconfs,cprefixes),shell=True)
                    if ret != 0:
                        setError(conn,which,"Problem generating conformers")
                        continue
                    
                    make_libraries(conn, dbprefixes, row)                    

                else:
                    setError(conn, row, "Problem reading input files")
            
            conn.close()
            time.sleep(1)
        except KeyboardInterrupt:
            raise
        except:
            print sys.exc_info()
            time.sleep(1)