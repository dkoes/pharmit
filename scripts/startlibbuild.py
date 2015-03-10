#!/usr/local/bin/python

#takes an fcgi request to create a compound database, sets up the diretory 
#with the input files and inserts the request into the Database

#another script will poll the database looking for databases that need to be
#built (this way only one is done at a time); this script just sets up the input directory

import sys, flup, MySQLdb, gzip, cgi, os,random,string

def setError(conn, id, err):
    #set the error state with a message in the database
    c = conn.cursor()
    c.execute("REPLACE INTO `databases` (status, message, id) VALUES('Error', %s, %s)", (err, id))
    
def application(environ, start_response):
    form = cgi.FieldStorage(fp=environ['wsgi.input'], environ=environ, keep_blank_values=True)
    conn = MySQLdb.connect (host = "localhost",user = "pharmit",db="pharmit")
    conn.autocommit(True)
    output = "Okay"
    try:
        name = form['dbname'].value
        description = form['description'].value
        isprivate = form['access'].value == 'private'
        email = form['email'].value
        file = form['compounds']
        
        #generate a random id
        id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(20))

        #setup input directory 
        os.mkdir(id)
        dir = os.path.abspath(id)
        
        #copy file into directory as input.[smi|sdf.gz]
        if file.filename.endswith('.sdf.gz'):
            open(id+'/input.sdf.gz', 'wb').write(file.file.read())
        elif file.filename.endswith('.smi.gz'):
            open(id+'/input.smi','w').write(gzip.GzipFile(mode='r',fileobj=file.file))
        elif file.filename.endswith('.smi'): #store smis uncompressed
            open(id + '/input.smi','w').write(file.file.read())
        elif file.filename.endswith('.sdf'): #store sdfs compressed
            gzip.open(id+'/input.sdf.gz','wb').write(file.file.read())
        else:
            output = "Error\nUnsupported file format in file %s"%file.filename
        #insert row in databases table
        c = conn.cursor()
        c.execute("REPLACE INTO `databases` (email, name, description, id, isprivate, status, message, directory) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)",
                   (email, name, description, id, isprivate, "Pending", "Your submission is pending in the queue.",dir))
    except:
        output = "Error\n"+str(sys.exc_info())
    status = '200 OK'
    response_headers = [('Content-type', 'text/plain'),
                        ('Content-Length', str(len(output)))]
    start_response(status, response_headers)
    return [output]



    
if __name__ == '__main__':    
    from flup.server.fcgi import WSGIServer
    from optparse import OptionParser

    parser = OptionParser(usage="Usage: %prog --userdir <directory>")
    parser.add_option("--userdir", dest="userdir",action="store",
                      help="directory to deposit submited data into", default=".", type="string", metavar="DIR")
    parser.add_option("--port", dest="port",action="store",
                  help="port to listen on", default=11111, type="int", metavar="P")
    
    (options, args) = parser.parse_args()
    os.chdir(options.userdir)
    WSGIServer(application, bindAddress=('localhost',options.port)).run()