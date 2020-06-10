#!/usr/bin/env python3

#Given a smiles file and prefix, extract a ligand file from the database of all the compounds
#with those smiles that have a name that starts with the prefix.

import sys,subprocess, re, MySQLdb, os, collections
import itertools
from rdkit.Chem import AllChem as Chem

def sortNames(prefix, names):
    #sort alphabetically, but with prefixed names first, make sure there are no duplicates
    prefixed = set()
    unprefixed = set()
    for n in names:
        n = n.strip();
        if n.startswith('MolPort') and '_' in n: #workaround for bad molport names
            n = n.split('_')[0]
        if n.startswith(prefix):
            prefixed.add(n)
        else:
            unprefixed.add(n)
            
    prefixed = list(set(prefixed))
    unprefixed = list(set(unprefixed))
    prefixed.sort()
    unprefixed.sort()
    
    return prefixed+unprefixed
    

if len(sys.argv) < 3:
    print("Need smiles and prefix")
    sys.exit(-1)

smilesf = sys.argv[1]
prefix = sys.argv[2]

#read in smiles
f = open(smilesf)
smiles = set()
for line in f:
    #read in the smiles
    vals = line.split(None,1)
    if len(vals) == 0:
        continue

    #remove salts from compound
    cmpds = vals[0].split('.')
    smile = max(cmpds, key=len) #take largest component by smiles length

    try: #catch any rdkit problems
        if len(smile) > 275:
            continue #conveniently skips molecule that hangs rdkit
        mol = Chem.MolFromSmiles(smile)
        Chem.SanitizeMol(mol)
        #to be sure, canonicalize smile (with iso)
        can = Chem.MolToSmiles(mol,isomericSmiles=True)
        smiles.add(can)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        sys.stderr.write("%s %s\n" %(e,smile))
        
            
            
            
conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")
    
cursor = conn.cursor()
        
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
        if len(names) < 1000: #avoid what are likely solvents
            print(sdfloc,i,bigname)
