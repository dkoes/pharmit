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

/*
 * PharmerServerCommands.h
 *
 *  Created on: Jan 25, 2011
 *      Author: dkoes
 *
 *  Interface for processing commands from a web client
 *
 */

#ifndef PHARMITSERVER_PHARMERSERVERCOMMANDS_H_
#define PHARMITSERVER_PHARMERSERVERCOMMANDS_H_

#include "cgi.h"
#include "PharmerQuery.h"
#include "PharmerServer.h"
#include "pharmarec.h"
#include <string>
#include <cstdio>
#include <gperftools/malloc_extension.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>

#ifndef SKIP_REGISTERZINC
#include <curl/curl.h>
#endif

using namespace std;
using namespace cgicc;

class Command
{
protected:

	FILE *LOG;
	SpinMutex& logmutex;

	//send a full JSON error message
	void sendError(FastCgiIO& IO, Cgicc& CGI, const char *msg)
	{
		using namespace boost;
		IO << HTTPPlainHeader();
		IO << "{\"status\" : 0, \"msg\": \"";
		IO << msg << "\"}\n";

		if (LOG)
		{
			posix_time::ptime t(posix_time::second_clock::local_time());
			fprintf(LOG, "error %s %s %s\n",
					posix_time::to_simple_string(t).c_str(),
					CGI.getEnvironment().getRemoteAddr().c_str(), msg);
			fprintf(LOG, "error cgi %s\n", cgiDump(CGI).c_str());
			fflush(LOG);
		}
	}

	//acknowledge a command that has no resonse to avoid http errors
	void sendAck(FastCgiIO& IO, const char *msg)
	{
		using namespace boost;
		IO << HTTPPlainHeader();
		IO << msg << "\n";
	}

public:
	Command(FILE * l, SpinMutex& logm) :
			LOG(l), logmutex(logm)
	{
	}
	virtual ~Command()
	{
	}

	virtual void execute(Cgicc& CGI, FastCgiIO& IO)
	{
	}

	//don't every dump
	virtual bool isFrequent()
	{
		return false;
	}
};

//any command that processes an existing query that needs to be validated
class QueryCommand: public Command
{
protected:
	WebQueryManager& queries;

	WebQueryHandle getQuery(Cgicc& CGI, FastCgiIO& IO)
	{
		if (!cgiTagExists(CGI, "qid"))
		{
			sendError(IO, CGI, "Bad data request. No query id.");
			return WebQueryHandle();
		}
		else
		{
			unsigned qid = cgiGetInt(CGI, "qid");
			if (qid == 0)
			{
				sendError(IO, CGI, "Bad data request. Invalid query id.");
				return WebQueryHandle();
			}
			else
			{
				WebQueryHandle query(queries.get(qid));
				if (!query)
				{
					//presumably garbage collected after a timeout
					sendError(IO, CGI, "Query data no longer available.");
				}
				return query;
			}
		}
	}
public:
	QueryCommand(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			Command(l, lm), queries(qs)
	{
	}
};


//set receptor key
class SetReceptor: public Command
{
	boost::filesystem::path logdirpath;
	SpinMutex recmutex;

public:
	SetReceptor(FILE * l, SpinMutex& lm, const boost::filesystem::path& ldp) :
			Command(l, lm), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		string key = cgiGetString(CGI, "key");
		string recstr = cgiGetString(CGI, "receptor");
		boost::filesystem::path rname = logdirpath / key;
		IO << HTTPPlainHeader();

		SpinLock lock(recmutex);
		if (boost::filesystem::exists(rname))
		{
			//exists already, but set it again anyway
			IO << "exists";
			ofstream rec(rname.string().c_str());
			rec << recstr;
		}
		else if (recstr.length() > 0)
		{
			//set
			ofstream rec(rname.string().c_str());
			rec << recstr;
			IO << "saved";
		}
	}
};

//return subset info, if a specific subset is provided, just provide its info (even if private)
class GetSubsets: public Command
{
	WebQueryManager& queries;

public:
	GetSubsets(WebQueryManager& qs, FILE * l, SpinMutex& lm):
		Command(l, lm), queries(qs) {}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		IO << HTTPPlainHeader();
		Json::FastWriter writer;

		if(!cgiTagExists(CGI, "subset"))
		{
			Json::Value info = queries.getJSONInfo();
			IO << writer.write(info);
		}
		else
		{
			string id = cgiGetString(CGI, "subset");
			Json::Value info = queries.getSingleJSON(id);
			IO << writer.write(info);
		}
	}
};


class StartQuery: public Command
{
	WebQueryManager& queries;
	boost::filesystem::path logdirpath;
	SpinMutex recmutex;

	const Pharmas *pharmas;

public:
	StartQuery(WebQueryManager& qs, FILE * l, SpinMutex& lm,
			const boost::filesystem::path& ldp, const Pharmas *ph) :
			Command(l, lm), queries(qs), logdirpath(ldp), pharmas(ph)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		namespace filesystem = boost::filesystem;
		if (!cgiTagExists(CGI, "json"))
		{
			//no query
			sendError(IO, CGI, "Invalid query syntax. No query data.");
		}
		else
		{
			Json::Value root; // will contains the root value after parsing.
			Json::Reader reader;
			bool parsingSuccessful = reader.parse(cgiGetString(CGI,
					"json"), root);
			if (!parsingSuccessful)
			{
				sendError(IO, CGI,
						"Invalid query. Could not parse query data.");
			}
			else
			{
				//check for memoized receptor
				if (root.isMember("receptorid"))
				{
					if (!root["receptor"].isString()
							|| root["receptor"].asString().length() == 0)
					{
						string recstr;
						filesystem::path rname = logdirpath
								/ root["receptorid"].asString();
						SpinLock lock(recmutex);
						if (filesystem::exists(rname))
						{
							ifstream rec(rname.string().c_str());
							stringstream str;
							str << rec.rdbuf();
							root["receptor"] = str.str();
						}
					}
				}

				unsigned totalMols = 0, totalConfs = 0;
				string msg;
				unsigned oldqid = cgiGetInt(CGI, "oldqid");

				//DEBUG code - output queries
				if(true)
				{
					SpinLock lock(logmutex);

					posix_time::ptime t(posix_time::second_clock::local_time());
					fprintf(LOG, "startq %s %s %s\n",
							posix_time::to_simple_string(t).c_str(),
							CGI.getEnvironment().getRemoteAddr().c_str(), cgiGetString(CGI,"json").c_str());
					fflush(LOG);
					lock.release();
				}
				if(true)
				{
					posix_time::ptime t(posix_time::second_clock::local_time());
					string lastname = (logdirpath/ "lastq").string();
					FILE *last = fopen(lastname.c_str(), "w");
					fprintf(last, "%s\n%s\n%s\n",
                                                        posix_time::to_simple_string(t).c_str(),
                                                        CGI.getEnvironment().getRemoteAddr().c_str(), cgiGetString(CGI,"json").c_str());
					fclose(last);
				}



				unsigned qid = queries.add(*pharmas, root,
						QueryParameters(root), oldqid, totalMols, totalConfs, msg);
				if (qid == 0) //invalid query
				{
					sendError(IO, CGI, msg.c_str());
				}
				else
				{
					//write log
					SpinLock lock(logmutex);

					posix_time::ptime t(posix_time::second_clock::local_time());
					fprintf(LOG, "query %s %s %d ",
							posix_time::to_simple_string(t).c_str(),
							CGI.getEnvironment().getRemoteAddr().c_str(), qid);

					if (root.isMember("receptorid"))
					{
						fprintf(LOG, "%s ",
								root["receptorid"].asString().c_str());
					}
					fprintf(LOG, "\n");
					fflush(LOG);
					lock.release();


					//output json with query id and size of database
					IO << HTTPPlainHeader();
					IO << "{\"status\": 1, \"qid\": " << qid
							<< ", \"numMols\": " << totalMols
							<< ", \"numConfs\": " << totalConfs
							<< "}";
				}
			}
		}
	}
};

class HasReceptor: public Command
{
	boost::filesystem::path logdirpath;

public:
	HasReceptor(FILE * l, SpinMutex& lm, const boost::filesystem::path& ldp) :
			Command(l, lm), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		string key = cgiGetString(CGI, "key");
		IO << HTTPPlainHeader();
		if (boost::filesystem::exists(logdirpath / key))
		{
			IO << "{\"recavail\": 1 }";
		}
		else
		{
			IO << "{\"recavail\": 0 }";
		}
	}
};

class GetMesh: public Command
{
	boost::filesystem::path logdirpath;

public:
	GetMesh(FILE * l, SpinMutex& lm, const boost::filesystem::path& ldp) :
			Command(l, lm), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		namespace filesystem = boost::filesystem;
		if (!cgiTagExists(CGI, "json"))
		{
			//no query
			sendError(IO, CGI, "Invalid syntax. No query data.");
		}
		else
		{
			Json::Value root; // will contains the root value after parsing.
			Json::Reader reader;
			bool parsingSuccessful = reader.parse(cgiGetString(CGI,
					"json"), root);
			if (!parsingSuccessful)
			{
				sendError(IO, CGI,
						"Invalid query. Could not parse query data.");
			}
			else
			{
				//check for memoized receptor
				if (root.isMember("receptorid"))
				{
					if (!root["receptor"].isString()
							|| root["receptor"].asString().length() == 0)
					{
						string recstr;
						filesystem::path rname = logdirpath
								/ root["receptorid"].asString();
						//not sure what good a mutex would do here..
						if (filesystem::exists(rname))
						{
							ifstream rec(rname.string().c_str());
							stringstream str;
							str << rec.rdbuf();
							root["receptor"] = str.str();
						}
					}

					ShapeConstraints excluder;
					excluder.readJSONExclusion(root);
					Json::Value mesh;

					if(cgiGetString(CGI, "type") == "exclusive")
					{
						mesh["tolerance"] = root["extolerance"];
						excluder.getExclusiveMesh(mesh);
					}
					else if(cgiGetString(CGI, "type") == "inclusive")
					{
						mesh["tolerance"] = root["intolerance"];
						excluder.getInclusiveMesh(mesh);
					}

					IO << HTTPPlainHeader();
					Json::FastWriter writer;
					IO << writer.write(mesh);
				}
			}
		}
	}
};


class CancelQuery: public QueryCommand
{
public:
	CancelQuery(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		{
			unsigned oldqid = cgiGetInt(CGI, "oldqid");
			WebQueryHandle query(queries.get(oldqid));
			if (query)
			{
				query->cancel();
			}
		}
		//make sure the handle is out of scope before purging
		queries.purgeOldQueries();
		sendAck(IO,"");
	}
};

class CancelSmina: public QueryCommand
{
public:
	CancelSmina(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		unsigned qid = cgiGetInt(CGI, "qid");
		WebQueryHandle query(queries.get(qid));
		if (query)
		{
			query->cancelSmina();
		}
		sendAck(IO,"");
	}
};

class Ping: public Command
{
public:
	Ping(FILE *l, SpinMutex& lm) :
			Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		IO << HTTPPlainHeader(); //just acknowledge existance
	}
};

class GetData: public QueryCommand
{

	//load params with paramaters from data
	//we are using the jquery datatables communication protocal
	void initDataParams(Cgicc& data, DataParameters& params)
	{
		params.drawCode = cgiGetInt(data, "draw");
		params.start = cgiGetInt(data, "start");
		params.num = cgiGetInt(data, "length");
		params.sort = (SortType::SortType) cgiGetInt(data, "order[0][column]");
		if (cgiGetString(data, "order[0][dir]") == "desc")
			params.reverseSort = true;
		;
		params.extraInfo = true;
	}

public:
	GetData(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			IO << HTTPPlainHeader();
			DataParameters params;
			initDataParams(CGI, params);
			Json::Value val;
			query->setDataJSON(params, val);
			IO << val;
		}
	}

	virtual bool isFrequent()
	{
		return true;
	}

};

class GetPharma: public Command
{
	const Pharmas *pharmas;
	boost::unordered_map<string, QueryParser*>& parsers;
	const boost::filesystem::path logdirpath;

public:
	GetPharma(FILE * l, SpinMutex& lm, const Pharmas *ph,
			boost::unordered_map<string, QueryParser*>& p,
			const boost::filesystem::path& ldp) :
			Command(l, lm), pharmas(ph), parsers(p), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		namespace filesystem = boost::filesystem;

		using namespace OpenBabel;
		vector<PharmaPoint> points;
		if (!cgiTagExists(CGI, "ligand") || !cgiTagExists(CGI, "ligandname"))
		{
			sendError(IO, CGI, "No ligand for pharmacophore identification.");
			return;
		}
		else
		{
			Json::Value val;
			string filedata = cgiGetString(CGI, "ligand");
			string filename = cgiGetString(CGI, "ligandname");
			OBFormat* format = OBConversion::FormatFromExt(filename.c_str());
			string ext = filesystem::extension(filename);
			if (parsers.count(ext))
			{
				//a pharmacophore query format
				stringstream str(filedata);
				vector<PharmaPoint> points;
				ShapeConstraints excluder;
				if (parsers[ext]->parse(*pharmas, str, points, excluder))
				{
					if (points.size() > 50)
					{
						//this many points is a pointless query and at some point
						//you just generate too many triplets
						sendError(IO, CGI,
								"I'm sorry, your query has too many points (more than 50!).");
					}
					else
					{
						convertPharmaJson(val, points);
						excluder.addToJSON(val);
						IO << HTTPPlainHeader();
						Json::FastWriter writer;
						val["status"] = 1;
						val["mol"] = false;
						//output the json as one line
						IO << writer.write(val);
					}
				}
				else
				{
					sendError(IO, CGI,
							"Error parsing query format file. Note that third party file formats (e.g. pml, ph4) are reverse engineered. Please submit any examples that do not work so we can improve the parser.");
				}
			}
			else //molecular data
			{
				if (format == NULL)
				{
					stringstream err;
					err << "Could not understand molecular data in " << filename
							<<
							" (OpenBabel compatible format required). Pharmacophore features must be in ph4, pml, or json format.";
					sendError(IO, CGI, err.str().c_str());
					return;
				}

				//check for receptor
				string receptor;
				OBFormat* rformat = NULL;
				if (cgiTagExists(CGI, "reckey") && cgiTagExists(CGI, "recname"))
				{
					string key = cgiGetString(CGI, "reckey");
					filesystem::path rname = logdirpath / key;
					if (filesystem::exists(rname))
					{
						ifstream rec(rname.string().c_str());
						stringstream str;
						str << rec.rdbuf();
						receptor = str.str();
						string fname = cgiGetString(CGI, "recname");
						rformat = OBConversion::FormatFromExt(fname.c_str());
					}
				}
				//if we didn't get a receptor from a reckey, see if we were sent the whole thing
				if (receptor.size() == 0 && cgiTagExists(CGI, "receptor")
						&& cgiTagExists(CGI, "recname"))
				{
					receptor = cgiGetString(CGI, "receptor");
					string fname = cgiGetString(CGI, "recname");
					rformat = OBConversion::FormatFromExt(fname.c_str());
				}

				if (!jsonPharmaQuery(*pharmas, val, filedata, format, receptor,
						rformat))
				{
					sendError(IO, CGI, "Could not parse molecular data.");
					return;
				}

				IO << HTTPPlainHeader();
				Json::FastWriter writer;
				//val filled in above
				val["status"] = 1;
				val["mol"] = true;
				//output the json as one line
				IO << writer.write(val);
			}
		}
	}
};

class GetMol: public QueryCommand
{
public:
	GetMol(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			IO << HTTPPlainHeader();
			unsigned index = cgiGetInt(CGI, "loc");
			query->outputMol(index, IO, false, cgiTagExists(CGI, "minimize"));
		}
	}
};

class SaveRes: public QueryCommand
{
public:
	SaveRes(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		using namespace boost;
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			unsigned nR = query->numResults(); //do before lock
			SpinLock lock(logmutex);
			posix_time::ptime t(
					posix_time::second_clock::local_time());
			fprintf(LOG, "save %s %s %lu %u\n", posix_time::to_simple_string(
					t).c_str(), CGI.getEnvironment().getRemoteHost().c_str(),
					cgiGetInt(CGI, "qid"), nR);
			fflush(LOG);
			lock.release();

			IO << "Content-Type: text/sdf\n";
			IO
			<< "Content-Disposition: attachment; filename=query_results.sdf\n";
			IO << endl;
			query->outputMols(IO);
		}
	}

};


class Receptor: public Command
{
public:
	Receptor(FILE * l, SpinMutex& lm) :
			Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		//TODO: cache receptor here (have to generate identical hash as client)
		IO << HTTPHTMLHeader();
		IO << html() << body();
		BOOST_FOREACH(const FormFile& f, CGI.getFiles())
		{
			IO << f.getData();
		}
		IO << body() << html();
	}
};

class Echo: public Command
{
public:
	Echo(FILE * l, SpinMutex& lm) :
			Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		//response goes into iframe and textarea, the body of which is parsed out
		//as long as </textarea> isn't in the file this works
		IO << HTTPHTMLHeader();
		IO << html() << body();
		IO << textarea();
		BOOST_FOREACH(const FormFile& f, CGI.getFiles())
		{
			IO << f.getData();
		}
		IO << textarea() << body() << html();
	}
};

//send back non-file upload stuff to save
class SaveData: public Command
{
public:
	SaveData(FILE * l, SpinMutex& lm) :
			Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		string type = "text/plain";
		if (cgiTagExists(CGI, "type"))
			type = cgiGetString(CGI, "type");
		string fname = "result.txt";
		if (cgiTagExists(CGI, "fname"))
			fname = cgiGetString(CGI, "fname");

		IO << "Content-Type: " << type << "\n";
		IO << "Content-Disposition: attachment; filename=" << fname << "\n";
		IO << endl;
		IO << cgiGetString(CGI, "data");
	}
};

class GetStatus: public QueryCommand
{
public:
	GetStatus(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
			QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		IO << HTTPPlainHeader();
		unsigned active, inactive, defunct;
		queries.getCounts(active, inactive, defunct);
		double load = 0;

		ifstream ldfile("/proc/loadavg");
		ldfile >> load;
		size_t mem = 0;
		MallocExtension::instance()->GetNumericProperty(
				"generic.current_allocated_bytes", &mem);
		double gb = round(100.0 * mem / (1024.0 * 1024 * 1024))
				/ 100.0;
		IO << "{\"msg\": \"Active: " << active << " Inactive: "
				<< inactive << " Defunct: " << defunct
				<< " Memory: " << gb << "GB"
						" Load: " << load << " TotalQ: "
				<< queries.processedQueries() << "\"}\n";
	}

	virtual bool isFrequent()
	{
		return true;
	}

};

//smina commands
class SminaCommand: public QueryCommand
{
protected:
	//do a buffered copy
	static void copy_stream(ostream& out, istream& in)
	{
		const int BUFFSZ = 4096;
		char buffer[BUFFSZ];
		unsigned read = 0;
		do
		{
			in.read(buffer, BUFFSZ);
			read = in.gcount();
			if (read)
				out.write(buffer, read);
		}
		while (read > 0 && in);
	}

	string server;
	string port;

public:
	SminaCommand(FILE * l, SpinMutex& lm, WebQueryManager& qs,
			const string& s, unsigned p) :
			QueryCommand(l, lm, qs),
					server(s), port(boost::lexical_cast<string>(p))
	{

	}
};

//start a smina minimization using the current receptor and query results
class StartSmina: public SminaCommand
{
	boost::asio::io_service my_io_service;
	boost::filesystem::path logdirpath;

public:
	StartSmina(FILE * l, SpinMutex& lm, WebQueryManager& qs,
			const boost::filesystem::path& ldp, const string& s, unsigned p) :
			SminaCommand(l, lm, qs, s, p), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		using namespace boost;
		using namespace boost::asio;
		using namespace boost::asio::ip;
		using namespace OpenBabel;
		namespace filesystem = boost::filesystem;

		if (server.length() == 0)
		{
			sendError(IO, CGI, "No minimization server configured");
			return;
		}
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			unsigned nR = query->numResults();
			SpinLock lock(logmutex);
			unsigned max = cgiGetInt(CGI, "num");
			posix_time::ptime t(posix_time::second_clock::local_time());
			fprintf(LOG, "startmin %s %s %u %u %lu\n",
					posix_time::to_simple_string(t).c_str(),
					CGI.getEnvironment().getRemoteHost().c_str(),
					nR, max, cgiGetInt(CGI, "qid"));
			fflush(LOG);
			lock.release();

			//need receptor
			string receptor;
			string rid = cgiGetString(CGI, "receptorid");
			if (rid.size() == 0)
			{
				sendError(IO, CGI, "Invalid receptor id.");
				return;
			}

			//use the passed receptor if it is there
			filesystem::path rname = logdirpath / rid;
			if (filesystem::exists(rname))
			{
				ifstream rec(rname.string().c_str());
				stringstream str;
				str << rec.rdbuf();
				receptor = str.str();
			}

			if (receptor.size() == 0)
			{
				sendError(IO, CGI, "Missing receptor.");
				return;
			}

			string fname = cgiGetString(CGI, "recname"); //need for file extension
			if (fname.length() == 0)
			{
				sendError(IO, CGI, "Missing receptor filename.");
				return;
			}

			OBFormat* rformat = OBConversion::FormatFromExt(fname.c_str());
			OBConversion conv;
			conv.SetInFormat(rformat);
			conv.SetOutFormat("PDBQT");
			conv.SetOptions("r", OBConversion::OUTOPTIONS);  //rigid
			OBMol rec;
			conv.ReadString(&rec, receptor);

			rec.AddHydrogens(true);
			FOR_ATOMS_OF_MOL(a, rec)
			{
			  a->GetPartialCharge();
			}

			string pdbqtrec = conv.WriteString(&rec);

			//connect to server
			stream_ptr minstrm(new tcp::iostream(server, port));
			if (!minstrm)
			{
				sendError(IO, CGI, "Could not connect to minimization server.");
				return;
			}

			*minstrm << "startmin\n";
			*minstrm << query->getSminaID() << "\n";
			*minstrm << "receptor " << pdbqtrec.length() << " 1\n"; //ispdbqt
			*minstrm << pdbqtrec;
			*minstrm << "1 0\n"; //reorient on, frag off

			unsigned sminaid = 0;
			*minstrm >> sminaid;

			//start asynchronous sending of ligand data
			query->sendSminaResults(server, port, minstrm, sminaid, max);
			IO << HTTPPlainHeader();
			IO << "{\"status\": 1, \"sminaid\": " << sminaid << "}";

		}

	}

};

class GetSminaMol: public SminaCommand
{
public:
	GetSminaMol(FILE * l, SpinMutex& lm, WebQueryManager& qs, const string& s,
			unsigned p) :
			SminaCommand(l, lm, qs, s, p)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			if (!cgiTagExists(CGI, "molid"))
			{
				sendError(IO, CGI, "Missing id.");
				return;
			}
			unsigned mid = cgiGetInt(CGI, "molid");
			IO << HTTPPlainHeader();
			unsigned sminaid = query->getSminaID();

			boost::asio::ip::tcp::iostream strm(server, port);
			if (!strm)
			{
				IO
						<< "{\"status\" : 0, \"msg\" : \"Could not connect to minimization server.\"}";
				return;
			}
			strm << "getmol\n" << sminaid << " " << mid << "\n";
			copy_stream(IO, strm); //returns mol
		}

	}
};

//parameters for retrieving minimization results
struct SminaParameters
{
	double maxScore;
	double maxRMSD;

	unsigned start; //were to start
	unsigned num; //how many to include, if zero then all

	unsigned sortType; //score=0,rmsd=1,origpos=2

	bool reverseSort;

	bool unique;

	//need to default to finie numbers for parsing
	SminaParameters() :
			maxScore(99999), maxRMSD(99999), start(0), num(0), sortType(0), reverseSort(
					false), unique(false)
	{
	}

	SminaParameters(Cgicc& data) :
			maxScore(99999), maxRMSD(99999), start(0), num(0), sortType(0), reverseSort(
					false), unique(false)
	{
		start = cgiGetInt(data, "start");
		num = cgiGetInt(data, "length");

		int sortcol = cgiGetInt(data, "order[0][column]");
		//column 3->score, 4->rmsd
		if (sortcol == 3)
			sortType = 0;
		else if (sortcol == 4)
			sortType = 1;
		else
			sortType = 2;

		string dir = cgiGetString(data, "order[0][dir]");
		if (dir == "desc" || dir == "dsc")
			reverseSort = true;
		else
			reverseSort = false;

		if (cgiTagExists(data, "maxscore"))
			maxScore = cgiGetDouble(data, "maxscore");
		if (cgiTagExists(data, "maxrmsd"))
			maxRMSD = cgiGetDouble(data, "maxrmsd");

		if (cgiTagExists(data, "unique"))
			unique = cgiGetInt(data, "unique");

	}

	void write(ostream& out) const
			{
		out << maxRMSD << " ";
		out << maxScore << " ";
		out << start << " ";
		out << num << " ";
		out << sortType << " ";
		out << reverseSort << " ";
		out << unique << "\n";
	}
};

class SaveSmina: public SminaCommand
{
public:
	SaveSmina(FILE * l, SpinMutex& lm, WebQueryManager& qs, const string& s,
			unsigned p) :
			SminaCommand(l, lm, qs, s, p)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		using namespace boost;
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			unsigned nR = query->numResults();
			SpinLock lock(logmutex);
			unsigned max = cgiGetInt(CGI, "num");
			posix_time::ptime t(posix_time::second_clock::local_time());
			fprintf(LOG, "save %s %s %u %u %lu\n",
					posix_time::to_simple_string(t).c_str(),
					CGI.getEnvironment().getRemoteHost().c_str(),
					nR, max, cgiGetInt(CGI, "qid"));
			fflush(LOG);
			lock.release();

			IO << "Content-Type: text/sdf\n";
			IO
					<< "Content-Disposition: attachment; filename=minimized_results.sdf.gz\n";
			IO << endl;

			boost::asio::ip::tcp::iostream strm(server, port);
			if (!strm)
			{
				IO << "ERROR contacting minimization server.\n";
				return;
			}
			strm << "getmols\n" << query->getSminaID() << "\n";
			SminaParameters param(CGI);
			param.write(strm);

			copy_stream(IO, strm);
		}

	}

};

//return smina minimization data (scores) for a query
class GetSminaData: public SminaCommand
{
public:
	GetSminaData(FILE * l, SpinMutex& lm, WebQueryManager& qs, const string& s,
			unsigned p) :
			SminaCommand(l, lm, qs, s, p)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			IO << HTTPPlainHeader();
			unsigned sminaid = query->getSminaID();
			SminaParameters param(CGI);
			boost::asio::ip::tcp::iostream strm(server, port);
			if (!strm)
			{
				IO
						<< "{\"status\" : 0, \"error\" : \"Could not connect to minimization server.\"}";
				return;
			}

			int drawcode = cgiGetInt(CGI, "draw");

			strm << "getjsonscores\n";
			strm << sminaid << " " << drawcode << "\n";
			param.write(strm);

			//start the json output and provide draw code
			copy_stream(IO, strm);
		}

	}

	virtual bool isFrequent()
	{
		return true;
	}

};

#endif /* PHARMITSERVER_PHARMERSERVERCOMMANDS_H_ */
