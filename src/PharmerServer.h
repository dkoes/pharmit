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
 * PharmerServer.h
 *
 *  Created on: Sep 10, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_PHARMERSERVER_H_
#define PHARMITSERVER_PHARMERSERVER_H_

#include <vector>
#include <string>
#include <boost/unordered_map.hpp>
#include "SpinLock.h"
#include "PharmerQuery.h"

using namespace std;

#define SERVERTHREADS 16
void pharmer_server(unsigned port, const vector<boost::filesystem::path>& prefixpaths,
		boost::unordered_map<string, StripedSearchers>& databases,
		const string& logdir, const string& minServer,unsigned minPort);


//wrapper for webquery that automatically takes care of usecnt
class WebQueryHandle
{
	PharmerQuery *ptr;
public:
	WebQueryHandle(): ptr(NULL) {}
	WebQueryHandle(PharmerQuery *q): ptr(q)
	{
		if(ptr) ptr->incrementUseCnt();
	}

	WebQueryHandle(const WebQueryHandle& rhs): ptr(rhs.ptr)
	{
		if(ptr) ptr->incrementUseCnt();
	}

	~WebQueryHandle()
	{
		if(ptr) ptr->decrementUseCnt();
	}

	PharmerQuery& operator*() { return *ptr; }
	PharmerQuery* operator->() { return ptr; }

	operator bool() { return ptr != NULL; }
};

//an instance of this classes manages a whole bunch of pharmer queries being
//served over the web (so have to deal with concurrency and timeout dellocation)
//each query is assigned a unique id for later reference
class WebQueryManager
{
	unsigned nextID; //counter to generate unique IDs
	typedef boost::unordered_map<unsigned, PharmerQuery*> QueryMap;
	QueryMap queries;


	typedef boost::unordered_map<string, StripedSearchers> DBMap;
	DBMap databases; //standard
	DBMap publicDatabases; //user supplied
	DBMap privateDatabases;


	vector<boost::filesystem::path> publicPrefixes;
	vector<boost::filesystem::path> privatePrefixes;

	Json::Value json;
	Json::Value privatejson;

	boost::mutex lock;
public:
	WebQueryManager(boost::unordered_map<string, StripedSearchers>& dbs,
			const vector<boost::filesystem::path>& prefixes): nextID(1), databases(dbs)
	{
		//setup private and public prefixes, which should be subdirs of the
		//standard library prefixes
		for(unsigned i = 0, n = prefixes.size(); i < n; i++)
		{
			publicPrefixes.push_back(prefixes[i] / "Public");
			privatePrefixes.push_back(prefixes[i] / "Private");
		}
		setupJSONInfo();
	}

	//add a query
	//first parse the text and return 0 if invalid
	//if oldqid is set, then deallocate/reuse it
	//also set the number of conformers/mols that will be searched
	//if there is an error, put the message in msg
	unsigned add(const Pharmas& pharma, Json::Value& data, const QueryParameters& qp,
			unsigned oldqid, unsigned& totalMols, unsigned& totalConfs, string& msg);

	//check public and private prefixes for new databases and load them
	void addUserDirectories();

	//return pointer to query qid (or null if not present)
	//increments inUse on query, caller must decrement
	WebQueryHandle get(unsigned qid);

	unsigned purgeOldQueries();

	void getCounts(unsigned& active, unsigned& inactive, unsigned& defunct);
	unsigned processedQueries() const { return nextID-1; }

	void setupJSONInfo();
	Json::Value getJSONInfo() { return json; }
	Json::Value getSingleJSON(const string& id);
};

#endif /* PHARMITSERVER_PHARMERSERVER_H_ */
