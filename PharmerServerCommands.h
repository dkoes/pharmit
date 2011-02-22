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

/*
 * PharmerServerCommands.h
 *
 *  Created on: Jan 25, 2011
 *      Author: dkoes
 *
 *  Interface for processing commands from a web client
 *
 */

#ifndef PHARMERSERVERCOMMANDS_H_
#define PHARMERSERVERCOMMANDS_H_

#include "cgi.h"
#include "PharmerQuery.h"
#include "PharmerServer.h"
#include "pharmarec.h"
#include <string>
#include <cstdio>
#include <google/malloc_extension.h>

using namespace std;
using namespace cgicc;
using namespace boost;
using namespace OpenBabel;

class Command
{
protected:

	FILE *LOG;
	SpinMutex& logmutex;

	//send a full JSON error message
	void sendError(FastCgiIO& IO, const char *msg)
	{
		IO << HTTPPlainHeader();
		IO << "{\"status\" : 0, \"msg\": \"";
		IO << msg << "\"}\n";

		if(LOG)
		{
			fprintf(LOG,"error %s\n", msg);
		}
	}

public:
	Command(FILE * l, SpinMutex& logm): LOG(l), logmutex(logm) {}
	virtual ~Command() {}

	virtual void execute(Cgicc& CGI, FastCgiIO& IO) { }

	//don't every dump
	virtual bool isFrequent() { return false; }
};

//any command that processes an existing query that needs to be validated
class QueryCommand : public Command
{
protected:
	WebQueryManager& queries;

	WebQueryHandle getQuery(Cgicc& CGI, FastCgiIO& IO)
	{
		form_iterator itr = CGI["qid"];
		if (!cgiTagExists(CGI, "qid"))
		{
			sendError(IO, "Bad data request. No query id.");
			return WebQueryHandle();
		}
		else
		{
			unsigned qid = cgiGetInt(CGI, "qid");
			if (qid == 0)
			{
				sendError(IO, "Bad data request. Invalid query id.");
				return WebQueryHandle();
			}
			else
			{
				WebQueryHandle query(queries.get(qid));
				if (!query)
				{
					//presumably garbage collected after a timeout
					sendError(IO, "Query data no longer available.");
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

class StartQuery : public Command
{
	WebQueryManager& queries;
	filesystem::path logdirpath;
	SpinMutex recmutex;

	const Pharmas *pharmas;
	unsigned totalConfs;
	unsigned totalMols;

public:
	StartQuery(WebQueryManager& qs, FILE * l, SpinMutex& lm,
			const filesystem::path& ldp, const Pharmas *ph, unsigned tc, unsigned tm) :
		Command(l, lm), queries(qs), logdirpath(ldp), pharmas(ph), totalConfs(tc), totalMols(tm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		if (!cgiTagExists(CGI, "json"))
		{
			//no query
			sendError(IO, "Invalid query syntax. No query data.");
		}
		else
		{
			Json::Value root; // will contains the root value after parsing.
			Json::Reader reader;
			bool parsingSuccessful = reader.parse(cgiGetString(CGI,
					"json"), root);
			if (!parsingSuccessful)
			{
				sendError(IO,
						"Invalid query. Could not parse query data.");
			}
			else
			{
				//check for memoized receptor
				if(root.isMember("receptorid"))
				{
					string recstr;
					filesystem::path rname = logdirpath / root["receptorid"].asString();
					SpinLock lock(recmutex);
					if(filesystem::exists(rname))
					{
						ifstream rec(rname.string().c_str());
						stringstream str;
						str << rec.rdbuf();
						root["receptor"] = str.str();
					}
					else if(root["receptor"].isString() && root["receptor"].asString().length() > 0)
					{
						//set
						ofstream rec(rname.string().c_str());
						rec <<  root["receptor"].asString();
					}
				}

				unsigned oldqid = cgiGetInt(CGI, "oldqid");
				unsigned qid = queries.add(*pharmas, root, QueryParameters(
						root), oldqid);
				if (qid == 0) //invalid query
				{
					sendError(IO, "Invalid query.");
				}
				else
				{
					//write log
					SpinLock lock(logmutex);

					posix_time::ptime t(posix_time::second_clock::local_time());
					fprintf(LOG, "query %s %s %d ", posix_time::to_simple_string(t).c_str(),
							CGI.getEnvironment().getRemoteAddr().c_str(), qid);

					if(root.isMember("receptorid"))
					{
						fprintf(LOG, "%s ", root["receptorid"].asString().c_str());
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

class HasReceptor : public Command
{
	filesystem::path logdirpath;

public:
	HasReceptor(FILE * l, SpinMutex& lm, const filesystem::path& ldp) :
		Command(l, lm), logdirpath(ldp)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		string key = cgiGetString(CGI,"key");
		IO << HTTPPlainHeader();
		if(filesystem::exists(logdirpath / key))
		{
			IO << "{\"recavail\": 1 }";
		}
		else
		{
			IO << "{\"recavail\": 0 }";
		}
	}
};

class CancelQuery : public QueryCommand
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
	}
};


class Ping: public Command
{
public:
	Ping(FILE *l, SpinMutex& lm): Command(l, lm) {}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		IO << HTTPPlainHeader(); //just acknowledge existance
	}
};

class GetData : public QueryCommand
{

	//load params with paramaters from data
	void initDataParams(Cgicc& data, DataParameters& params)
	{
		params.start = cgiGetInt(data, "startIndex");
		params.num = cgiGetInt(data, "results");
		string s = cgiGetString(data, "sort");
		if(s == "RMSD")
			params.sort = SortType::RMSD;
		else if(s == "MolWeight")
			params.sort = SortType::MolWeight;
		else if(s == "RBnds")
			params.sort = SortType::NRBnds;
		params.reverseSort = cgiGetInt(data, "reverse");
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
		if(query)
		{
			IO << HTTPPlainHeader();
			DataParameters params;
			initDataParams(CGI, params);
			query->outputData(params, IO, true);
		}
	}

	virtual bool isFrequent() { return true; }

};

class GetPharma : public Command
{
	const Pharmas *pharmas;
	unordered_map<string, QueryParser*>& parsers;
public:
	GetPharma(FILE * l, SpinMutex& lm, const Pharmas *ph,
			unordered_map<string, QueryParser*>& p) :
		Command(l, lm), pharmas(ph), parsers(p)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		vector<PharmaPoint> points;
		const std::vector<FormFile>& files = CGI.getFiles();
		if (files.size() == 0)
		{
			sendError(IO,
					"No file for pharmacophore identification.");
		}
		else //just process first file
		{
			Json::Value val;
			string filedata = files[0].getData();
			string filename = files[0].getFilename();

			OBFormat* format = OBConversion::FormatFromExt(filename.c_str());
			string ext = filesystem::extension(filename);

			if (parsers.count(ext))
			{
				//a pharmacophore query format
				stringstream str(filedata);
				vector<PharmaPoint> points;
				if (parsers[ext]->parse(*pharmas, str, points))
				{
					convertPharmaJson(val, points);
					IO << HTTPHTMLHeader();
					IO << html() << body();
					Json::FastWriter writer;
					val["status"] = 1;
					val["mol"] = false;
					//output the json as one line
					IO << writer.write(val);
					IO << body() << html();
				}
				else
				{
					sendError(IO, "Error parsing query format file.");
				}
			}
			else if (format == NULL || !jsonPharmaQuery(*pharmas, val, filedata, format))
			{
				sendError(IO, "Could not understand molecular data.");
			}
			else
			{
				IO << HTTPHTMLHeader();
				IO << html() << body();
				Json::FastWriter writer;
				val["status"] = 1;
				val["mol"] = true;
				//output the json as one line
				IO << writer.write(val);
				//followed by mol
				IO << filedata << "\n";
				IO << body() << html();
			}
		}
	}
};

class GetMol : public QueryCommand
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
			query->outputMol(index, IO, true, cgiTagExists(CGI, "minimize"));
		}
	}
};

class SaveRes : public QueryCommand
{
public:
	SaveRes(FILE * l, SpinMutex& lm, WebQueryManager& qs) :
		QueryCommand(l, lm, qs)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		WebQueryHandle query = getQuery(CGI, IO);
		if (query)
		{
			SpinLock lock(logmutex);
			posix_time::ptime t(
					posix_time::second_clock::local_time());
			fprintf(LOG, "save %s %s %lu\n", posix_time::to_simple_string(
					t).c_str(), CGI.getEnvironment().getRemoteHost().c_str(), cgiGetInt(CGI, "qid"));
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

class Receptor : public Command
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

class Echo : public Command
{
public:
	Echo(FILE * l, SpinMutex& lm) :
		Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		//response goes into iframe, the body of which is parsed out
		IO << HTTPHTMLHeader();
		IO << html() << body();
		BOOST_FOREACH(const FormFile& f, CGI.getFiles())
		{
			IO << f.getData();
		}
		IO << body() << html();
	}
};

//send back non-file upload stuff to save
class SaveData : public Command
{
public:
	SaveData(FILE * l, SpinMutex& lm) :
		Command(l, lm)
	{
	}

	void execute(Cgicc& CGI, FastCgiIO& IO)
	{
		string type = "text/plain";
		if(cgiTagExists(CGI, "type")) type = cgiGetString(CGI, "type");
		string fname = "result.txt";
		if(cgiTagExists(CGI,"fname")) fname = cgiGetString(CGI, "fname");

		IO << "Content-Type: " << type << "\n";
		IO << "Content-Disposition: attachment; filename=" << fname << "\n";
		IO << endl;
		IO << cgiGetString(CGI, "data");
	}
};


class GetStatus : public QueryCommand
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

	virtual bool isFrequent() { return true; }

};



#endif /* PHARMERSERVERCOMMANDS_H_ */
