#!/usr/local/bin/python

#Turns smiles into conformers and deposits them in a local database if needed
#Each line of the smiles file is of the form 
#[smiles] [name]
#and name is assumed to be of a form that uniquely indicates the vendor

import sys,subprocess, re, MySQLdb, cStringIO
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser
import multiprocessing

def getRMS(mol, c1,c2):
    (rms,trans) = Chem.GetAlignmentTransform(mol,mol,c1,c2)
    return rms

def createconfs(uniqueid, smile, mol, options):
    maxconfs = options.maxconfs
    sample = int(options.sample*options.maxconfs)
    rmsdcut = options.rms
    energycut = options.energy
    cids = Chem.EmbedMultipleConfs(mol, sample,randomSeed=1301979)
    cenergy = []           
    output = cStringIO.StringIO()
    sdwriter = Chem.SDWriter(output) 
    mol.SetProp("_Name",str(uniqueid)) #the id should be the name 

    for conf in cids:
        #not passing confID only minimizes the first conformer
        converged = not Chem.UFFOptimizeMolecule(mol,confId=conf)
        cenergy.append(Chem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
    
    sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
    if len(sortedcids) > 0:
        mine = cenergy[sortedcids[0]]
    else:
        mine = 0
    #filter out via rmsd cutoff
    written = {}
    for conf in sortedcids:
        if len(written) >= maxconfs:
            break
        #check rmsd and energy
        passed = True
        if cenergy[conf]-mine > energycut:
            break
        for seenconf in written.iterkeys():
            rms = getRMS(mol,seenconf,conf) 
            if rms < rmsdcut:
                passed = False
                break
        if(passed):
            written[conf] = True
            sdwriter.write(mol,conf)

    sdwriter.close()
    #add conformers to database
    conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")
    cursor = conn.cursor()
    cursor.execute('UPDATE structures SET sdfs=COMPRESS(%s) WHERE `id`=%s AND smile=%s',(output.getvalue(),uniqueid, smile))
    conn.commit()
    cursor.close()
    conn.close()

#end creatconfs

#reads inputs from multiprocessing queue
def dowork(queue):
    while True:
        ins = queue.get()
        if not ins:
            return #al ldone
        try:
            createconfs(*ins)
        except Exception as e:
            print e, ins

if __name__ == '__main__':
	
    parser = OptionParser(usage="Usage: %prog [options] <input>.smi")
    parser.add_option("--maxconfs", dest="maxconfs",action="store",
                      help="maximum number of conformers to generate per a molecule (default 20)", default="20", type="int", metavar="CNT")
    parser.add_option("--sample_multiplier", dest="sample",action="store",
                      help="sample N*maxconfs conformers and choose the maxconformers with lowest energy (default 1)", default="1", type="float", metavar="N")
    parser.add_option("--seed", dest="seed",action="store",
                      help="random seed (default 9162006)", default="9162006", type="int", metavar="s")
    parser.add_option("--rms_threshold", dest="rms",action="store",
                      help="filter based on rms (default 0.7)", default="0.7", type="float", metavar="R")
    parser.add_option("--energy_window", dest="energy",action="store",
                      help="filter based on energy difference with lowest energy conformer", default="10", type="float", metavar="E")
    parser.add_option("-r","--replace", dest="replace",action="store_true",default=False,
                      help="replace already computed conformers")
    parser.add_option("--threads", dest="threads",action="store",
          help="number of threads to use", default="0", type="int", metavar="N")
    parser.add_option("-v","--verbose", dest="verbose",action="store_true",default=False,
                  help="verbose output")
    
    (options, args) = parser.parse_args()
    if(len(args) != 1):
        parser.error("Need input smiles")
        sys.exit(-1)
    
    conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")
    
    try:
        f = open(args[0])
    except IOError:
        print "Could not read file",sys.argv[1]
        
    #setup multiprocessing queues
    numt = multiprocessing.cpu_count()
    if options.threads > 0:
        numt = options.threads
        
    queue = multiprocessing.Queue(numt)
    for _ in xrange(numt):
        multiprocessing.Process(target=dowork, args=(queue,)).start()

    if options.verbose:
        print "Running with",numt,"threads"
    for line in f:
    #read in the smiles
        vals = line.split(None,1)
        if len(vals) != 2:
            #must have name 
            print "Missing name of",vals[0]
            continue
        name = vals[1]    
        #remove salts from compound
        cmpds = vals[0].split('.')
        smile = max(cmpds, key=len) #take largest component by smiles length
    
        try: #catch any rdkit problems
            mol = Chem.MolFromSmiles(smile)
            Chem.SanitizeMol(mol)
            #to be sure, canonicalize smile
            can = Chem.MolToSmiles(mol)
            cursor = conn.cursor()
            cursor.execute('SELECT * FROM structures WHERE smile = %s', (can,))
            #if smile is not in structures
            isnew = cursor.fetchone() == None
            if isnew:
                #insert without sdfs to get unique id 
                cursor.execute('INSERT INTO structures (smile,weight,sdfs) VALUES(%s,%s,\'\') ', (can, Chem.CalcExactMolWt(mol)))
                
            #get unique id
            cursor.execute('SELECT id FROM structures WHERE smile = %s', (can,))
            result = cursor.fetchone();
            uniqueid = result[0]
            
            #we always update the name
            cursor.execute('INSERT IGNORE INTO names (smile,name) VALUES(%s,%s)', (can,name))
            conn.commit()
            if options.verbose:
                print uniqueid,can
                            
            if isnew or options.replace:
                #create conformers and insert them
                queue.put((uniqueid, can, mol, options))
                
            
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            print e
    
    #clear out queues
    for _ in xrange(numt):
        queue.put(None)





#regardless, add name to names
