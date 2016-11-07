#!/usr/local/bin/python

#Given a prefix, extract a ligand file from the database of all the compounds
#that have a name that starts with the prefix.

import sys,subprocess, re, MySQLdb, os, collections
import itertools

def sortNames(prefix, names):
    #sort alphabetically, but with prefixed names first
    prefixed = []
    unprefixed = []
    for n in names:
        if n.startswith('MolPort') and '_' in n:
            n = n.split('_')[0] #remove meta data accidentally included in name
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

smiles = set()
rows = cursor.fetchall()
for row in rows:
    smile = row[0]
    name = row[1].strip()
    if name.startswith(prefix): #should be redundant
        smiles.add(smile)
        
#for each compound, get sdf location and id, output with sorted names
for smile in smiles:
    cursor.execute("SELECT id, sdfloc FROM structures WHERE smile = %s", (smile,))
    rows = cursor.fetchall() #should be one
    if len(rows) >= 1:
        (i, sdfloc) = rows[0]
        #get all names
        cursor.execute("SELECT name FROM names WHERE smile = %s", (smile,))
        names = cursor.fetchall()
        names = list(itertools.chain.from_iterable(names)) 
        bigname =' '.join(sortNames(prefix,names))
        bigname = bigname.replace('\n','')
        print sdfloc,i,bigname
