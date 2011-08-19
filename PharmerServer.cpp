/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "PharmerServer.h"
#include <unistd.h>
#include <set>

#include "cgicc/Cgicc.h"
#include "cgicc/HTTPHTMLHeader.h"
#include "cgicc/HTTPPlainHeader.h"
#include "cgicc/HTTPResponseHeader.h"
#include "cgicc/HTMLClasses.h"
#include <cgicc/CgiEnvironment.h>

#include "cgi.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <setjmp.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <json/json.h>
#include <boost/date_time.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/assign/list_of.hpp>
#include <curl/curl.h>

#include "pharmerdb.h"
#include "pharmarec.h"
#include "queryparsers.h"
#include "Timer.h"
#include <openbabel/obconversion.h>

#include <vector>
#include <cstdio>
#include "PharmerServerCommands.h"

using namespace std;
using namespace cgicc;
using namespace boost;
using namespace OpenBabel;

typedef struct sockaddr SA;
#define LISTENQ 1024
#define SERVERTHREADS 16


//define parsers both here and in PharmaQuery so we can be flexibly about
//dealing with file mappings
static TextQueryParser textParser;
static JSonQueryParser jsonParser;
static PH4Parser ph4Parser;
static PMLParser pmlParser;

static unordered_map<string, QueryParser*> parsers = assign::map_list_of
		("",(QueryParser*) &textParser)
		(".txt", (QueryParser*) &textParser)
		(".json", (QueryParser*) &jsonParser)
		(".query", (QueryParser*) &jsonParser)
		(".ph4", (QueryParser*) &ph4Parser)
		(".pml", (QueryParser*) &pmlParser);

/*
 * open_listenfd - open and return a listening socket on port
 *     Returns -1 and sets errno on Unix error.
 *     Shamelessly stolen from Introduction to Computer Systems, a wonderful
 *     undergraduate textbook
 */
int open_listenfd(int port)
{
	int listenfd, optval = 1;
	struct sockaddr_in serveraddr;

	/* Create a socket descriptor */
	if ((listenfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
		return -1;

	/* Eliminates "Address already in use" error from bind. */
	if (setsockopt(listenfd, SOL_SOCKET, SO_REUSEADDR, (const void *) &optval,
			sizeof(int)) < 0)
		return -1;

	/* Listenfd will be an endpoint for all requests to port
	 on any IP address for this host */
	bzero((char *) &serveraddr, sizeof(serveraddr));
	serveraddr.sin_family = AF_INET;
	serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
	serveraddr.sin_port = htons((unsigned short) port);
	if (::bind(listenfd, (SA *) &serveraddr, sizeof(serveraddr)) < 0)
		return -1;

	/* Make it a listening socket ready to accept connection requests */
	if (listen(listenfd, LISTENQ) < 0)
		return -1;
	return listenfd;
}


static void server_thread(unsigned listenfd, unordered_map<string, shared_ptr<Command> >& commands)
{
	FCGX_Request request;
	FCGX_InitRequest(&request, listenfd, 0);


	while (FCGX_Accept_r(&request) == 0)
	{
		try
		{
			Timer time;
			FastCgiIO IO(request);

			Cgicc CGI(&IO);
			//handle the cmd
			string cmdname = cgiGetString(CGI, "cmd");
			transform(cmdname.begin(), cmdname.end(), cmdname.begin(), ::tolower);
			if (commands.count(cmdname) == 0)
			{
				IO << HTTPResponseHeader("HTTP/1.1", 433,
						"Invalid query syntax.  Invalid command.");
			}
			else
			{
				commands[cmdname]->execute(CGI, IO);

				if (!Quiet && !commands[cmdname]->isFrequent())
				{
					cout << "Finished processing query: " << cgiGetString(CGI,
							"cmd") << " in " << time.elapsed() << "s ("
							<< time.elapsedUser() << "s) ";
					posix_time::ptime now(
							posix_time::second_clock::local_time());
					cout << " at " << posix_time::to_simple_string(now) << "\n";
				}
			}
		}
		catch (const std::exception& e)
		{
			// handle error condition
			cout << "Exception thrown " << e.what();
		}

		FCGX_Finish_r(&request);
	}
}


//start up a server, assumes databases are segregrated by molweight
void pharmer_server(unsigned port, vector< vector<MolWeightDatabase> >& databases,
		const string& logdir, unsigned totalConfs, unsigned totalMols)
{
	FILE *LOG;
	filesystem::path logdirpath;
	SpinMutex logmutex;
	SpinMutex recmutex;
	const Pharmas *pharmas = &databases.back().back().db->getPharmas();

	logdirpath = filesystem::path(logdir);
	OBConversion conv; //load plugins

	string logname = (logdirpath / "LOG").file_string();
	LOG = fopen(logname.c_str(), "a");
	if (LOG == NULL)
	{
		cerr << "Could not open log directory " << logname << "\n";
		exit(-1);
	}
	setlinebuf(LOG);

	WebQueryManager queries(databases);

	FCGX_Init();
	curl_global_init(CURL_GLOBAL_ALL);

	int listenfd = open_listenfd(port);
	if (listenfd < 0)
	{
		cerr << strerror(errno);
		exit(-1);
	}
	cout << "Listening on port " << port << "\n";

	cout << totalConfs << " conformations of " << totalMols << " compounds; "
			<< (double) totalConfs / (double) totalMols << " average per mol\n";

	unordered_map<string, shared_ptr<Command> > commands = assign::map_list_of
			("startquery",shared_ptr<Command>(new StartQuery(queries, LOG, logmutex, logdirpath,pharmas, totalConfs, totalMols)))
					("hasreceptor", shared_ptr<Command>(new HasReceptor(LOG, logmutex,logdirpath)))
					("setreceptor", shared_ptr<Command>(new SetReceptor(LOG, logmutex,logdirpath)))
					("cancelquery",shared_ptr<Command>(new CancelQuery(LOG, logmutex, queries)))
					("getdata",shared_ptr<Command>(new GetData(LOG, logmutex, queries)))
					("getpharma",shared_ptr<Command>(new GetPharma(LOG, logmutex, pharmas, parsers, logdirpath)))
					("getmol",shared_ptr<Command>(new GetMol(LOG, logmutex, queries)))
					("saveres",shared_ptr<Command>(new SaveRes(LOG, logmutex, queries)))
					("ping", shared_ptr<Command>(new Ping(LOG, logmutex)))
					("receptor",shared_ptr<Command>(new Receptor(LOG, logmutex)))
					("echo",shared_ptr<Command>(new Echo(LOG, logmutex)))
					("savedata",shared_ptr<Command>(new SaveData(LOG, logmutex)))
					("getstatus",shared_ptr<Command>(new GetStatus(LOG, logmutex, queries)));


	for (unsigned i = 0; i < SERVERTHREADS; i++)
	{
		thread server(server_thread, listenfd, ref(commands));
	}

	{
		//startup message
		posix_time::ptime t(posix_time::second_clock::local_time());
		SpinLock lock(logmutex);
		fprintf(LOG, "startup %s\n", posix_time::to_simple_string(t).c_str());
		fflush(LOG);
	}

	while(true)
	{
		this_thread::sleep(posix_time::time_duration(0,3,0,0));
		unsigned npurged = queries.purgeOldQueries();
		if(npurged > 0)
		{
			posix_time::ptime t(posix_time::second_clock::local_time());
			SpinLock lock(logmutex);
			fprintf(LOG, "purged %s %d\n", posix_time::to_simple_string(t).c_str(), npurged);

			//return mem to system
			MallocExtension::instance()->ReleaseFreeMemory();
			fflush(LOG);
		}
	}
}



//add a query
//first parse the text and return 0 if invalid
unsigned WebQueryManager::add(const Pharmas& pharmas, Json::Value& data,
		const QueryParameters& qp, unsigned oldqid)
{
	vector<PharmaPoint> queryPoints;

	readPharmaPointsJSON(pharmas, data, queryPoints);

	//check result - need at least 3 points to define a triangle
	if (queryPoints.size() < 3)
	{
		return 0;
	}

	unique_lock<mutex>(lock);

	//TODO: optimize successive but similar queries
	if (oldqid > 0 && queries.count(oldqid) > 0)
	{
		PharmerQuery *oldq = queries[oldqid];
		oldq->cancel();
		if(oldq->finished())
		{
			//immediately deallocate so we can reuse mem for this query
			delete oldq;
			queries.erase(oldqid);
		}
	}
	unsigned id = nextID++;
	queries[id] = new PharmerQuery(databases, queryPoints, qp, databases.size());
	queries[id]->execute(false); //don't wait for result
	return id;
}

WebQueryHandle WebQueryManager::get(unsigned qid)
{
	unique_lock<mutex>(lock);

	if (queries.count(qid) == 0)
		return WebQueryHandle(NULL);
	PharmerQuery *q = queries[qid];
	if (q == NULL)
		return WebQueryHandle(NULL);

	q->access();
	return WebQueryHandle(q);
}

#define TIMEOUT (60*30)
//delete any queries that are older than a timeout
unsigned WebQueryManager::purgeOldQueries()
{
	unique_lock<mutex>(lock);
	vector<unsigned> toErase;
	for (QueryMap::iterator itr = queries.begin(), end = queries.end(); itr
			!= end; itr++)
	{
		PharmerQuery *q = itr->second;
		if (q->idle() > TIMEOUT || q->cancelled())
		{
			if(q->finished())
			{
				delete q;

				queries[itr->first] = NULL;
				toErase.push_back(itr->first); //this way can bypass iterator invalidation issues
			}
			else
			{
				q->cancel();
			}
		}
	}

	BOOST_FOREACH(unsigned qid, toErase)
	{
		queries.erase(qid);
	}

	return toErase.size();
}

//count types of queries
void WebQueryManager::getCounts(unsigned& active, unsigned& inactive,
		unsigned& defunct)
{
	active = inactive = defunct = 0;
	shared_lock<shared_mutex> L(mutex);

	for (QueryMap::iterator itr = queries.begin(), end = queries.end(); itr
			!= end; itr++)
	{
		PharmerQuery *q = itr->second;
		if (!q->finished())
			active++;
		else if (q->idle() > TIMEOUT)
			defunct++;
		else
			inactive++;
	}
}
