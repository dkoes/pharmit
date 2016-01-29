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
 * dbloader.cpp
 *
 *  Created on: Mar 11, 2015
 *      Author: dkoes
 */

#include "dbloader.h"
#include <glob.h>
#include <boost/algorithm/string/predicate.hpp>
using namespace boost;
using namespace std;

//thread class for loading database info
struct LoadDatabase
{
	unsigned totalConf;
	unsigned totalMols;
	bool hasShape;

	LoadDatabase() :
			totalConf(0), totalMols(0), hasShape(true)
	{

	}

	void operator()( boost::shared_ptr<PharmerDatabaseSearcher>& database, unsigned i,
			filesystem::path dbpath)
	{
		shared_ptr<PharmerDatabaseSearcher> db(new PharmerDatabaseSearcher(dbpath));

		if (!db->isValid())
		{
			cerr << "Error reading database " << dbpath;
			exit(-1);
		}
		totalConf += db->numConformations();
		totalMols += db->numMolecules();
		hasShape &= db->hasShape();
		database = db;
	}
};

//load databases based on commandline arguments
void loadDatabases(vector<filesystem::path>& dbpaths, StripedSearchers& databases)
{
	databases.totalConfs = 0;
	databases.totalMols = 0;
	databases.stripes.resize(dbpaths.size());
	vector<LoadDatabase> loaders(dbpaths.size());
	thread_group loading_threads;
	for (unsigned i = 0, n = dbpaths.size(); i < n; i++)
	{
		if (!filesystem::is_directory(dbpaths[i]))
		{
			cerr << "Invalid database directory path: " << dbpaths[i] << "\n";
			exit(-1);
		}

		loading_threads.add_thread(
				new thread(ref(loaders[i]), ref(databases.stripes[i]), i, dbpaths[i]));
	}
	loading_threads.join_all();

	BOOST_FOREACH(const LoadDatabase& ld, loaders)
	{
		databases.totalConfs += ld.totalConf;
		databases.totalMols += ld.totalMols;
		databases.hasShape &= ld.hasShape;
	}
}

//from stack overflow
inline vector<string> glob(const std::string& pat){
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

//load striped databases from the specified prefixes
//get the keys for each database from the database json file
void loadFromPrefixes(vector<filesystem::path>& prefixes, unordered_map<string, StripedSearchers >& databases)
{
	unordered_map<string, StripedSearchers > blank;
	loadNewFromPrefixes(prefixes, databases, blank);
}

//only load databases that don't already have keys in olddatabases
void loadNewFromPrefixes(vector<filesystem::path>& prefixes,
		unordered_map<string, StripedSearchers >& databases,
		const unordered_map<string, StripedSearchers >& olddatabases)
{
	assert(prefixes.size() > 0);
	filesystem::path jsons = prefixes[0] / "*" / "dbinfo.json";
	vector<string> infos = glob(jsons.c_str());

	for(unsigned i = 0, n = infos.size(); i < n; i++)
	{
		filesystem::path subdir(infos[i]);
		subdir = subdir.remove_filename();
		filesystem::path name = subdir.filename();

		Json::Value json;
		Json::Reader reader;
		ifstream info(infos[i].c_str());
		if(!reader.parse(info, json)) {
			cerr << "Error reading database info " << infos[i] << "\n";
			continue;
		}
		if(!json.isMember("subdir")) {
			cerr << "Missing subdir int database info " << infos[i] << "\n";
			continue;
		}
		string specified = json["subdir"].asString();

		if(!algorithm::ends_with(specified, name.string())) //ignore prefixed subdirs like Public, the specified subdir must match the actual subdir
		{
			cerr << "Ignoring " << name << "\n";
			continue;
		}
		else if(olddatabases.count(specified) == 0)
		{
			vector<filesystem::path> dbpaths(prefixes.size());
			bool badsubdir = false;

			if(json.isMember("splitdirs") && json["splitdirs"].isArray()) {
				//we have a really large database (like pubchem) that needs to be
				//broken up further
				Json::Value splits = json["splitdirs"];
				unsigned nsplits = splits.size();
				dbpaths.resize(0);
				dbpaths.reserve(prefixes.size() * nsplits);
				for(unsigned p = 0, np = prefixes.size(); p < np; p++)
				{
					filesystem::path dir = prefixes[p] / name / splits[i].asString();
					dbpaths.push_back(dir);

					filesystem::path infofile = dir / "info";
					if(!filesystem::exists(infofile) || filesystem::file_size(infofile) == 0)
					{
						badsubdir = true;
						cerr << "Invalid subdir info: " << infofile << "\n";
					}
				}
			}
			else {
				for(unsigned p = 0, np = prefixes.size(); p < np; p++)
				{
					filesystem::path dir = prefixes[p] / name;
					dbpaths[p] = dir;
				}
			}

			if(!badsubdir)
				loadDatabases(dbpaths, databases[specified]);
		}
	}
}
