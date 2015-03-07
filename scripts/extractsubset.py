#!/usr/local/bin/python

#Given a prefix, extract a ligand file from the database of all the compounds
#that have a name that starts with the prefix.

import sys,subprocess, re, MySQLdb, os, collections

def sortNames(prefix, names):
    #sort alphabetically, but with prefixed names first
    prefixed = []
    unprefixed = []
    for n in names:
        if n.startswith(prefix):
            prefixed.append(n)
        else:
            unprefixed.append(n)
            
    prefixed.sort()
    unprefixed.sort()
    
    return prefixed+unprefixed
    

if len(sys.argv) < 2:
    print "Need prefix"
    sys.exit(-1)
    
prefix = sys.argv[1]
conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")
    
cursor = conn.cursor()
cursor.execute("SELECT smile,name FROM names WHERE name LIKE %s", (prefix+'%',))

compounds = collections.defaultdict(list) #smile -> [names]
rows = cursor.fetchall()
for row in rows:
    smile = row[0]
    name = row[1].strip()
    if name.startswith(prefix): #should be redundant
        compounds[smile].append(name)
        
#for each compound, get sdf location and id, output with sorted names
for smile in compounds.iterkeys():
    cursor.execute("SELECT id, sdfloc FROM structures WHERE smile = %s", (smile,))
    (i, sdfloc) = cursor.fetchall()[0]
    print sdfloc,i,' '.join(sortNames(prefix,compounds[smile]))
