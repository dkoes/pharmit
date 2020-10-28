/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

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
#include <csignal>
#include "pharmerdb.h"
#include "pharmarec.h"
#include "queryparsers.h"
#include "Timer.h"
#include <openbabel/obconversion.h>

#include <vector>
#include <cstdio>
#include "PharmerServerCommands.h"
#include "dbloader.h"

using namespace std;
using namespace cgicc;
using namespace OpenBabel;

typedef struct sockaddr SA;
#define LISTENQ 1024

//define parsers both here and in PharmaQuery so we can be flexibly about
//dealing with file mappings
static TextQueryParser textParser;
static JSonQueryParser jsonParser;
static PH4Parser ph4Parser;
static PMLParser pmlParser;

static boost::unordered_map<string, QueryParser*> parsers = assign::map_list_of
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


static void server_thread(unsigned listenfd, boost::unordered_map<string, std::shared_ptr<Command> >& commands)
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


static WebQueryManager *queriesptr = NULL;
static void signalhandler(int sig)
{
  if (sig == SIGUSR1 && queriesptr != NULL)
  {
	  //look for new directories
	  queriesptr->addUserDirectories();
  }
}





//start up a server
void pharmer_server(unsigned port, const vector<filesystem::path>& prefixpaths,
		boost::unordered_map<string, StripedSearchers >& databases,
		const string& logdir, const string& minServer, unsigned minPort)
{
	FILE *LOG;
	filesystem::path logdirpath;
	SpinMutex logmutex;
	SpinMutex recmutex;

	OpenBabel::OBPlugin::LoadAllPlugins();

	if (databases.size() == 0)
	{
		cerr << "No valid databases specified\n";
		exit(-1);
	}
	//these all better be the same...
	const Pharmas *pharmas = &databases.begin()->second.stripes.back()->getPharmas();

	logdirpath = filesystem::path(logdir);
	OBConversion conv; //load plugins

	string logname = (logdirpath / "LOG").string();
	LOG = fopen(logname.c_str(), "a");
	if (LOG == NULL)
	{
		cerr << "Could not open log directory " << logname << "\n";
		exit(-1);
	}
	setlinebuf(LOG);

	WebQueryManager queries(databases, prefixpaths);
	queriesptr = &queries; //for signal handler

	FCGX_Init();

	int listenfd = open_listenfd(port);
	if (listenfd < 0)
	{
		cerr << strerror(errno);
		exit(-1);
	}
	cout << "Listening on port " << port << "\n";

	boost::unordered_map<string, std::shared_ptr<Command> > commands = assign::map_list_of
			("startquery",std::shared_ptr<Command>(new StartQuery(queries, LOG, logmutex, logdirpath,pharmas)))
					("hasreceptor", std::shared_ptr<Command>(new HasReceptor(LOG, logmutex,logdirpath)))
					("setreceptor", std::shared_ptr<Command>(new SetReceptor(LOG, logmutex,logdirpath)))
					("getmesh", std::shared_ptr<Command>(new GetMesh(LOG, logmutex,logdirpath)))
					("cancelquery",std::shared_ptr<Command>(new CancelQuery(LOG, logmutex, queries)))
					("getdata",std::shared_ptr<Command>(new GetData(LOG, logmutex, queries)))
					("getpharma",std::shared_ptr<Command>(new GetPharma(LOG, logmutex, pharmas, parsers, logdirpath)))
					("getsubsets",std::shared_ptr<Command>(new GetSubsets(queries, LOG, logmutex)))
					("getmol",std::shared_ptr<Command>(new GetMol(LOG, logmutex, queries)))
					("saveres",std::shared_ptr<Command>(new SaveRes(LOG, logmutex, queries)))
					("ping", std::shared_ptr<Command>(new Ping(LOG, logmutex)))
					("receptor",std::shared_ptr<Command>(new Receptor(LOG, logmutex)))
					("echo",std::shared_ptr<Command>(new Echo(LOG, logmutex)))
					("savedata",std::shared_ptr<Command>(new SaveData(LOG, logmutex)))
					("getstatus",std::shared_ptr<Command>(new GetStatus(LOG, logmutex, queries)))
					("startsmina",std::shared_ptr<Command>(new StartSmina(LOG, logmutex, queries, logdirpath, minServer, minPort)))
					("cancelsmina",std::shared_ptr<Command>(new CancelSmina(LOG, logmutex, queries)))
					("getsminadata",std::shared_ptr<Command>(new GetSminaData(LOG, logmutex, queries, minServer, minPort)))
					("getsminamol", std::shared_ptr<Command>(new GetSminaMol(LOG, logmutex, queries, minServer, minPort)))
					("savesmina", std::shared_ptr<Command>(new SaveSmina(LOG, logmutex, queries, minServer, minPort)));


	for (unsigned i = 0; i < SERVERTHREADS; i++)
	{
		thread server(server_thread, listenfd, boost::ref(commands));
	}

	{
		//startup message
		posix_time::ptime t(posix_time::second_clock::local_time());
		SpinLock lock(logmutex);
		fprintf(LOG, "startup %s\n", posix_time::to_simple_string(t).c_str());
		fflush(LOG);
	}

	//load user libraries after startup
	queries.addUserDirectories();

	signal(SIGUSR1, signalhandler);

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


void WebQueryManager::addUserDirectories()
{
	DBMap newpublic;
	loadNewFromPrefixes(publicPrefixes, newpublic, publicDatabases);

	DBMap newprivate;
	loadNewFromPrefixes(privatePrefixes, newprivate, privateDatabases);

	if(newpublic.size() > 0 || newprivate.size() > 0)
	{
		//maps aren't thread-safe so protect
		boost::unique_lock<boost::mutex>(lock);

		publicDatabases.insert(newpublic.begin(), newpublic.end());
		privateDatabases.insert(newprivate.begin(), newprivate.end());
		setupJSONInfo();
	}

}

//add a query
//first parse the text and return 0 if invalid
unsigned WebQueryManager::add(const Pharmas& pharmas, Json::Value& data,
		const QueryParameters& qp, unsigned oldqid,
		unsigned& totalMols, unsigned& totalConfs, string& msg)
{
	vector<PharmaPoint> queryPoints;

	readPharmaPointsJSON(pharmas, data, queryPoints);

	ShapeConstraints excluder;
	excluder.readJSONExclusion(data);

	//check result - need at least 3 points to define a triangle
	if (queryPoints.size() < 3 && !qp.isshape)
	{
		msg = "Invalid query.  Three features are required.";
		return 0;
	}

	if(qp.isshape && queryPoints.size() > 0 && !excluder.isMeaningful())
	{ //is this constraint necessary? ...apparently since otherwise people will search all of PubChem  with no shape constraints set at all
		msg = "Please provide more expressive shape constraints to reduce the number of hits for pharmacophore filtering.";
		return 0;
	}

	//identify databases to search
	vector< std::shared_ptr<PharmerDatabaseSearcher> > dbs;
	unsigned numslices = 1;
	StripedSearchers *searchers = NULL;
	if(databases.count(qp.subset))
	{
		searchers = &databases[qp.subset];
	}
	else
	{
		boost::unique_lock<boost::mutex> L(lock); //public/private database may change underneath us
		if(publicDatabases.count(qp.subset))
		{
			searchers = &publicDatabases[qp.subset];
		}
		else if(privateDatabases.count(qp.subset))
		{
			searchers = &privateDatabases[qp.subset];
		}
		else
		{
			msg = "Unknown subset.";
			return 0;
		}
	}

	if(qp.isshape && !searchers->hasShape)
	{
		msg = "Database is missing shape information.";
		return 0;
	}

	dbs = searchers->stripes;
	numslices = min(thread::hardware_concurrency(),(unsigned)dbs.size()); //how many threads we should run, don't do more than available
	totalMols = searchers->totalMols;
	totalConfs = searchers->totalConfs;

	boost::unique_lock<boost::mutex> L(lock);

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
	queries[id] = new PharmerQuery(dbs, queryPoints, qp, excluder, numslices);
	if(!queries[id]->isValid(msg)) return 0;
	L.unlock();

	queries[id]->execute(false); //don't wait for result
	return id;
}

//order json's
static bool jsonInfoSorter(const Json::Value& a, const Json::Value& b)
{
	return a["name"] < b["name"];
}

//se json info to describe databases
void WebQueryManager::setupJSONInfo()
{
	vector<Json::Value> jsons;
	BOOST_FOREACH(DBMap::value_type i, databases )
	{
		jsons.push_back(i.second.getJSON());
	}

	sort(jsons.begin(), jsons.end(), jsonInfoSorter);

	Json::Value ret;
	ret["standard"].resize(jsons.size());
	for(unsigned i = 0, n = jsons.size(); i < n; i++)
	{
		ret["standard"][i] = jsons[i];
	}

	//public databases
	jsons.clear();
	BOOST_FOREACH(DBMap::value_type i, publicDatabases )
	{
		jsons.push_back(i.second.getJSON());
	}

	sort(jsons.begin(), jsons.end(), jsonInfoSorter);

	ret["public"].resize(jsons.size());
	for(unsigned i = 0, n = jsons.size(); i < n; i++)
	{
		ret["public"][i] = jsons[i];
	}

	json = ret;

	//private should not get sent back
	jsons.clear();
	BOOST_FOREACH(DBMap::value_type i, privateDatabases )
	{
		jsons.push_back(i.second.getJSON());
	}

	sort(jsons.begin(), jsons.end(), jsonInfoSorter);

	Json::Value priv;
	priv.resize(jsons.size());
	for(unsigned i = 0, n = jsons.size(); i < n; i++)
	{
		priv[i] = jsons[i];
	}
	privatejson = priv;
}

//return the json dbinfo for the specified db
//this is primarily used for private databases
Json::Value WebQueryManager::getSingleJSON(const string& id)
{
	if(privateDatabases.count(id))
		return privateDatabases[id].getJSON();
	if(databases.count(id))
		return databases[id].getJSON();
	if(publicDatabases.count(id))
		return publicDatabases[id].getJSON();

	//otherwise error
	Json::Value ret;
	ret["error"] = "Invalid access code";
	return ret;
}


WebQueryHandle WebQueryManager::get(unsigned qid)
{
	boost::unique_lock<boost::mutex>(lock);

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
	boost::unique_lock<boost::mutex>(lock);
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
	boost::unique_lock<boost::mutex> L(lock);

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
