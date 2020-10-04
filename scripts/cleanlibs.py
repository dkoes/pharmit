#!/usr/bin/env python3
# Actually remove user libraries marked for removal.  This should probably be
# run as a cronjob, but the pharmit server needs to be restarted after it is
# done, so am running manually for now.


import time, sys, os, MySQLdb, json, shutil, glob
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser


if __name__ == '__main__':
    
    parser = OptionParser(usage="Usage: %prog [options]")
    parser.add_option('--dbprefixes', dest="dbprefixfile", action="store", 
            help="file containing path prefixes for storing databases", default="",metavar="FILE")    

        
    (options, args) = parser.parse_args()

    dbprefixfile = os.path.abspath(options.dbprefixfile)
    if not dbprefixfile or not os.path.isfile(dbprefixfile):
        print("Require prefix file for storing databases")
        sys.exit(-1)        
    prefixes = [line.rstrip() for line in open(dbprefixfile)]
    #look for removed libraries
    conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="pharmit")
    conn.autocommit(True)
    conn.query("SELECT * FROM `databases` WHERE status = 'REMOVE' ORDER BY submitted ASC")
    rows = conn.store_result().fetch_row(how=1,maxrows=0)
    for row in rows:
        #validate values
        dir = row["directory"]
        which = row["id"]
        if not os.path.isdir(dir):
            print("Problem reading input directory",dir)
        else:
            info = json.load(open(dir+'/dbinfo.json'))
            for p in prefixes:
                path = p+'/'+info['subdir']
                if os.path.exists(path) and ('Public' in path or 'Private' in path):
                    realdir = glob.glob(path+'-*')
                    if len(realdir) == 1:
                        realdir = realdir[0]
                        print(path,realdir)
                        os.remove(path) # this is a symlink
                        shutil.rmtree(realdir)
        
    
    conn.close()

