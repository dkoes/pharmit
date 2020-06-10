#!/usr/bin/env python3

#I forgot to strip the newline from compound names, so now have to go through
#the whole database and fix it

import sys,subprocess, re, MySQLdb, os
    
conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="conformers")
    
cursor = conn.cursor()
cursor.execute("SELECT smile,name FROM names WHERE name LIKE '%%\n'")

rows = cursor.fetchall()
 
for row in rows:
    smile = row[0]
    name = row[1]
    newname = name.strip()
    print(smile,newname, name) 
    if newname != name: #should be redundant
        c = conn.cursor()
        c.execute("UPDATE names SET name = %s WHERE smile = %s AND name = %s", (newname,smile,name))
        conn.commit()
