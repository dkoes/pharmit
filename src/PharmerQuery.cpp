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
 * PharmerQuery.cpp
 *
 *  Created on: Aug 5, 2010
 *      Author: dkoes
 */

#include <gperftools/malloc_extension.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/assign/list_of.hpp>
#include "PharmerQuery.h"
#include "Corresponder.h"
#include "queryparsers.h"
#include "Timer.h"
#include "ShapeResults.h"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace OpenBabel;

static TextQueryParser textParser;
static JSonQueryParser jsonParser;
static PH4Parser ph4Parser;
static PMLParser pmlParser;

static const int MAX_FEATURES = 64;

static unordered_map<string, QueryParser*> parsers = assign::map_list_of("",
		(QueryParser*) &textParser)(".txt", (QueryParser*) &textParser)(".json",
		(QueryParser*) &jsonParser)(".query", (QueryParser*) &jsonParser)(
		".ph4", (QueryParser*) &ph4Parser)(".pml", (QueryParser*) &pmlParser);

bool PharmerQuery::validFormat(const string& ext)
{
	return parsers.count(ext) > 0;
}

void PharmerQuery::initializeTriplets()
{
	//create all possible triplets
	const unsigned n = points.size();
	coralloc.setSize(n);
	triplets.reserve(math::binomial_coefficient<double>(n, 3));
	tripIndex.resize(extents[n][n][n]);
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = i + 1; j < n; j++)
		{
			for (unsigned k = j + 1; k < n; k++)
			{
				unsigned index = triplets.size();
				triplets.push_back(QueryTriplet(points, i, j, k));
				tripIndex[i][j][k] = index;
				tripIndex[i][k][j] = index;
				tripIndex[j][i][k] = index;
				tripIndex[j][k][i] = index;
				tripIndex[k][i][j] = index;
				tripIndex[k][j][i] = index;
			}
		}
	}

	valid = true;
}

PharmerQuery::PharmerQuery(
		const vector< std::shared_ptr<PharmerDatabaseSearcher> >& dbs,
		istream& in, const string& ext, const QueryParameters& qp, unsigned nth) :
		databases(dbs), params(qp), valid(false), stopQuery(false), tripletMatchThread(
		NULL), shapeMatchThread(NULL), lastAccessed(time(NULL)), corrsQs(dbs.size()), currsort(
				SortType::Undefined), currrev(false), nthreads(nth), dbcnt(0), inUseCnt(0), numactives(0),
				totalmols(0), sminaid(0)
{
	if (dbs.size() == 0)
	{
		errorStr = "No databases provided.";
		return;
	}

	//parse query file
	if (parsers.count(ext) == 0)
	{
		errorStr = "Invalid format extension.";
		return;
	}
	if (!parsers[ext]->parse(dbs.back()->getPharmas(), in, points,
			excluder))
	{
		errorStr = "Could not parse query.";
		return;
	}


	dbcnt = dbs.size();
	if (!params.isshape && points.size() < 3)
	{
		errorStr = "Require at least three pharmacophore points.";
		return;
	}
	if(!params.isshape && points.size() >  MAX_FEATURES) 
	{
		errorStr = "Too many pharmacophore features specified.";
		return;
	}
	initializeTriplets();
}

PharmerQuery::PharmerQuery(
		const vector< std::shared_ptr<PharmerDatabaseSearcher> >& dbs,
		const vector<PharmaPoint>& pts, const QueryParameters& qp,
		const ShapeConstraints& ex, unsigned nth) :
		databases(dbs), points(pts), params(qp), excluder(ex), valid(false), stopQuery(
				false), tripletMatchThread(NULL), shapeMatchThread(NULL), lastAccessed(time(NULL)), corrsQs(
				dbs.size()), currsort(SortType::Undefined), currrev(false), nthreads(nth), dbcnt(
				0), inUseCnt(0), numactives(0), totalmols(0), sminaid(0)
{
	if (dbs.size() == 0)
	{
		errorStr = "No databases provided.";
		return;
	}


	dbcnt = dbs.size();

	if(dbcnt > 0)
	{
		Json::FastWriter writer;
		const Json::Value& dbinfo = dbs[0]->getJSON();

		if(dbinfo["numactives"].isNumeric())
		{ //the presence of this field indicates this is a benchmark database with compounds marked as active (in the title)
			numactives = dbinfo["numactives"].asUInt();
			totalmols = 0;
			for(unsigned i = 0; i < dbcnt; i++)
			{
				totalmols += dbs[i]->numMolecules();
			}
		}
	}

	if (points.size() < 3)
	{
		errorStr = "Require at least three pharmacophore points.";
		return;
	}
	if(points.size() >  MAX_FEATURES)
        {
                errorStr = "Too many pharmacophore features specified.";
                return;
        }

	initializeTriplets();
}

static void incrementDegrees(const Triplet& t, unsigned *degrees)
{
	for (unsigned i = 0; i < 3; i++)
	{
		//assume always connecting from previous triangle, so a nonzero node
		//is going to be sharing an edge
		if (degrees[t.getPoints()[i].index] == 0)
			degrees[t.getPoints()[i].index] += 2;
		else
			degrees[t.getPoints()[i].index]++;
	}
}

//generate the query triplets
//need a "triangulation"; a set of triangles where every
//point is connected to at least 3 other points
//it is not clear (to me at least) what the "optimal" choice is
//we want a small number of uncommon triangles, but we would
//also like the triangles to build off of previous triangles so we
//can throw out bad correspondences early

//so here's the idea: triangles are sorted by decreasing "goodness"
//always build onto current collection by connecting a point of zero degree
//to that last chosen triangle, using at least one point of degree 2
//
//we expand triplets so we have all possible legal orderings based on overlap
void PharmerQuery::generateQueryTriplets(PharmerDatabaseSearcher& pharmdb,
		vector<vector<QueryTriplet> >& expandedtrips)
{
	vector<QueryTriplet> trips;
	trips.reserve(points.size() - 1);
	//have the database d rank the triplets (lower is better, since this
	//is presumably correlated to frequency)
	vector<double> ranking;
	unsigned mintrip = pharmdb.rankTriplets(triplets, ranking);
	unsigned degrees[points.size()];
	memset(degrees, 0, sizeof(degrees));

	trips.push_back(triplets[mintrip]);
	incrementDegrees(trips.back(), degrees);
	unordered_set<unsigned> seen;
	seen.insert(mintrip);

	unsigned N = points.size();
	if (N > 3)
	{
		vector<PharmaPoint> unconnectedpoints; //point of each triplet not connected to successor
		while (trips.size() < N - 1)
		{
			QueryTriplet& prev = trips.back();
			mintrip = 0;
			unsigned bestj = 0;
			unsigned bestkindex = 0;
			double minval = HUGE_VAL;
			//for each side of the previous triangle
			for (unsigned i = 0; i < 3; i++)
			{
				unsigned iindex = prev.getPoints()[i].index;
				unsigned j = (i + 1) % 3;
				unsigned jindex = prev.getPoints()[j].index;
				//if one of the vertices has a degree less than 3
				if (degrees[iindex] < 3 || degrees[jindex] < 3)
				{
					//find a point with degree 0 to connect to, unless this
					//is the last triangle
					unsigned reqdeg = 0;
					if (trips.size() == N - 2)
						reqdeg = 2;
					for (unsigned kindex = 0; kindex < N; kindex++)
					{
						if (degrees[kindex] == reqdeg && kindex != iindex
								&& kindex != jindex)
						{
							unsigned tindex = tripIndex[iindex][jindex][kindex];
							if (ranking[tindex] < minval)
							{
								mintrip = tindex;
								minval = ranking[tindex];
								bestj = j;
								bestkindex = kindex;
							}
						}
					}
				}
			}
			assert(minval != HUGE_VAL);
			assert(seen.count(mintrip) == 0);
			seen.insert(mintrip);
			unsigned prevucon = (bestj + 1) % 3;
			trips.back().setNextUnconnected(prevucon);
			trips.push_back(triplets[mintrip]);
			unconnectedpoints.push_back(
					points[prev.getPoints()[prevucon].index]);
			trips.back().setPrevUnconnectedIndex(bestkindex, unconnectedpoints);
			incrementDegrees(trips.back(), degrees);
		}
	}
	expandedtrips.resize(trips.size());
	for (unsigned i = 0, n = trips.size(); i < n; i++)
	{
		trips[i].expand(expandedtrips[i], *this);
	}

	if (trips.size() == 1)
	{
		//don't bother with fingerprints
		BOOST_FOREACH(QueryTriplet& trp, expandedtrips[0])
		{
			trp.setSkipFingers(true);
		}
	}
}

//check to see if query threads are done, and if so deallocate them
void PharmerQuery::checkThreads()
{
	if (tripletMatchThread != NULL)
	{
		if (tripletMatchThread->timed_join(posix_time::time_duration(0, 0, 0)))
		{
			delete tripletMatchThread;
			tripletMatchThread = NULL; // thread is done
			if (!Quiet)
			{
				size_t mem = 0;
				MallocExtension::instance()->GetNumericProperty(
						"generic.current_allocated_bytes", &mem);
				double gb = round(10000.0 * mem / (1024.0 * 1024 * 1024))
						/ 10000.0;
				cout << "CORALLOC NUMCHUNKS " << coralloc.numChunks() << "\t"
						<< gb << "GB\n";
			}
		}
	}

	if (shapeMatchThread != NULL)
	{
		if (shapeMatchThread->timed_join(posix_time::time_duration(0, 0, 0)))
		{
			delete shapeMatchThread;
			shapeMatchThread = NULL; // thread is done
		}
	}

}

//launch Nthreads to match triplets
void PharmerQuery::thread_tripletMatches(PharmerQuery *query)
{
	//generate a thread for each set of databases
	thread_group tmthreads;

	for (unsigned t = 0; t < query->nthreads; t++)
	{
		tmthreads.add_thread(new boost::thread(thread_tripletMatch, query));
	}
	tmthreads.join_all();
}

//launch Nthreads to match shapes
void PharmerQuery::thread_shapeMatches(PharmerQuery *query)
{
	//generate a thread for each set of databases
	thread_group shthreads;

	for (unsigned t = 0; t < query->nthreads; t++)
	{
		shthreads.add_thread(new boost::thread(thread_shapeMatch, query));
	}
	shthreads.join_all();
}

//match all the triplets in a database
void PharmerQuery::thread_tripletMatch(PharmerQuery *query)
{
	unsigned db = 0;
	while (query->dbSearchQ.pop(db))
	{
		if (query->stopQuery)
			break;

		PharmerDatabaseSearcher& pharmdb = *query->databases[db];
		vector<vector<QueryTriplet> > trips;
		query->generateQueryTriplets(pharmdb, trips);
		TripletMatchAllocator tmalloc(trips.size());
		TripletMatches matches(tmalloc, query->params, trips.size(), 1);

		Timer t;
		pharmdb.generateTripletMatches(trips, matches,
				query->stopQuery);
		if (!Quiet)
			cout << "PMTime " << t.elapsed() << "\n";
		if (query->stopQuery)
			break;

		if (!Quiet)
		{
			cout << db << " TripletMatches ";
			matches.dumpCnts();
		}
		query->corrsQs[db].addProducer();

		Timer ct;
		Corresponder sponder(query->databases[db],
				db, query->databases.size(),
				query->points, trips, matches, query->coralloc, 0,
				query->corrsQs[db], query->params, query->excluder,
				query->stopQuery);
		sponder();
		if (!Quiet)
			cout << "CTime " << ct.elapsed() << "\n";

	}
}


//perform shape matching
void PharmerQuery::thread_shapeMatch(PharmerQuery *query)
{
	unsigned db = 0;
	while (query->dbSearchQ.pop(db))
	{
		if (query->stopQuery)
			break;

		PharmerDatabaseSearcher& pharmdb = *query->databases[db];
		MTQueue<CorrespondenceResult*>& corrQ =	query->corrsQs[db];
		corrQ.addProducer();

		ShapeResults shapes(query->databases[db], query->points, corrQ, query->coralloc, query->params, query->excluder, db, query->databases.size(), query->stopQuery);
		pharmdb.generateShapeMatches(query->excluder, shapes);
		corrQ.removeProducer();
	}
}


//execute the query, if block is true then perform synchronously
void PharmerQuery::execute(bool block)
{
	for (unsigned d = 0, nd = databases.size(); d < nd; d++)
		dbSearchQ.push(d);

	if(params.isshape)
	{
		coralloc.setSize(0); //no actuall correspondances
		assert(shapeMatchThread == NULL);
		shapeMatchThread = new boost::thread(thread_shapeMatches, this);

		if (block) //wait for completion
		{
			shapeMatchThread->join();
			delete shapeMatchThread;
			shapeMatchThread = NULL;
		}
	}
	else
	{
		assert(tripletMatchThread == NULL);
		tripletMatchThread = new boost::thread(thread_tripletMatches, this);

		if (block) //wait for completion
		{
			tripletMatchThread->join();
			delete tripletMatchThread;
			tripletMatchThread = NULL;
		}
	}
}

PharmerQuery::~PharmerQuery()
{
	checkThreads();
	if (tripletMatchThread != NULL)
		abort(); //threds must be done before we can destruct the results array
	if (shapeMatchThread != NULL)
		abort(); //threds must be done before we can destruct the results array
}

void PharmerQuery::cancel()
{
	stopQuery = true;
	cancelSmina();
}

void PharmerQuery::cancelSmina()
{
	if (sminaid > 0 && sminaServer.length() > 0)
	{
		boost::asio::ip::tcp::iostream strm(sminaServer, sminaPort);
		strm << "cancel\n" << sminaid << "\n";
	}
}

bool PharmerQuery::threadsDone()
{
	checkThreads();
	return tripletMatchThread == NULL && shapeMatchThread == NULL;
}

bool PharmerQuery::finished() //okay to deallocate
{
	if (!threadsDone())
		return false;
	if (inUseCnt > 0)
		return false;
	return true;
}

typedef bool (*QRCompare)(const QueryResult* lhs, const QueryResult* rhs);

static bool rmsdCompare(const QueryResult* lhs, const QueryResult* rhs)
{
	return lhs->c->val < rhs->c->val;
}

static bool mwCompare(const QueryResult* lhs, const QueryResult* rhs)
{
	return lhs->c->weight < rhs->c->weight;
}

static bool rbndCompare(const QueryResult* lhs, const QueryResult* rhs)
{
	return lhs->c->nRBnds < rhs->c->nRBnds;
}

class ReverseSort
{
	QRCompare cmp;
	public:
	ReverseSort(QRCompare c) :
			cmp(c)
	{
	}

	bool operator()(const QueryResult* lhs, const QueryResult* rhs)
	{
		return cmp(rhs, lhs);
	}
};

//sort according to srt, do not re-sort
//do stable sorts, which require an explicit reverse sort
void PharmerQuery::sortResults(SortTyp srt, bool reverse)
{
	if (currsort == srt && currrev == reverse)
		return;
	QRCompare cmp = NULL;
	switch (srt)
	{
		case SortType::RMSD:
			cmp = rmsdCompare;
			break;
		case SortType::MolWeight:
			cmp = mwCompare;
			break;
		case SortType::NRBnds:
			cmp = rbndCompare;
			break;
		default:
			break;
	}
	currsort = srt;
	currrev = reverse;
	if (cmp != NULL)
	{
		if (reverse)
			stable_sort(results.begin(), results.end(), ReverseSort(cmp));
		else
			stable_sort(results.begin(), results.end(), cmp);
	}
}

//reduce according to parameters, ie, fewer conformers
void PharmerQuery::reduceResults()
{
	if (params.reduceConfs != 0 && params.reduceConfs != UINT_MAX)
	{
		unordered_map<unsigned, unsigned> molCnts;
		vector<QueryResult*> newr;
		newr.reserve(results.size());
		for (unsigned i = 0, n = results.size(); i < n; i++)
		{
			unsigned molid = results[i]->c->molid;
			if (molCnts.count(molid) == 0)
				molCnts[molid] = 0;

			if (molCnts[molid]++ < params.reduceConfs)
			{
				newr.push_back(results[i]);
			}
		}
		swap(results, newr);
	}
}

//copy current set of results in to vector
//filter the vector based on query params
//return true if there are more results expected
bool PharmerQuery::loadResults()
{
	access();
	SpinLock lock(mutex);
	checkThreads();
	bool moretoread = !threadsDone();
	unsigned newcnt = 0; //number added to results
	for (unsigned i = 0, n = corrsQs.size(); i < n; i++)
	{
		vector<CorrespondenceResult*> corrs;
		moretoread |= corrsQs[i].popAll(corrs);

		for (unsigned j = 0, nc = corrs.size(); j < nc; j++)
		{
			CorrespondenceResult *c = corrs[j];
			currsort = SortType::Undefined;
			results.push_back(
					new (resalloc.alloc(sizeof(QueryResult))) QueryResult(c));
			newcnt++;
		}
	}

	//filter and sort if we've seen something new
	sortResults(params.sort, params.reverseSort);
	reduceResults();

	if (params.maxHits != 0 && results.size() > params.maxHits)
	{
		results.resize(params.maxHits);
	}

	return moretoread;
}

//return all current results subject to data parameters
//return true if results are still being produced
bool PharmerQuery::getResults(const DataParameters& dp,
		vector<QueryResult*>& out)
{
	bool notdone = loadResults();
	out.clear();

	SpinLock lock(mutex);
	if (dp.num > 0)
		out.reserve(dp.num);
	else
		out.reserve(results.size());

	sortResults(dp.sort, dp.reverseSort);

	unsigned end = dp.num == 0 ? results.size() : dp.num + dp.start;
	for (unsigned i = dp.start, n = results.size(); i < n && i < end; i++)
	{
		if (dp.extraInfo)
			setExtraInfo(*results[i]);
		out.push_back(results[i]);
	}

	return notdone;
}

//if necessary, load extra info into r
void PharmerQuery::setExtraInfo(QueryResult& r)
{
	if (!r.name[0])
	{
		std::shared_ptr<PharmerDatabaseSearcher> db;
		unsigned long loc = getLocation(&r, db);

		//TODO: make this more efficient (don't need to unpack full mol)
		MolData mdata;
		PMolReaderSingleAlloc pread;
		db->getMolData(loc, mdata, pread);
		r.name = mdata.mol->getTitle();
	}
}

//output text result of query (correspondences)
void PharmerQuery::outputData(const DataParameters& dp, ostream& out,
		bool jsonHeader)
{
	vector<QueryResult*> r;
	bool notdone = getResults(dp, r);

//json status line
	if (jsonHeader)
	{
		out << "{\"status\" : 1, \"done\" : ";
		if (!notdone)
			out << "true, ";
		else
			out << "false, ";

		out << "\"total\" : " << results.size() << " }\n";
	}

	for (unsigned i = 0, n = r.size(); i < n; i++)
	{
		out << i + dp.start << "," << r[i]->c->val << "," << r[i]->c->weight
				<< "," << r[i]->c->nRBnds << "," << r[i]->name << ","
				<< r[i]->c->molid << "," << r[i]->c->location << "\n";
	}
}

//including nan is invalid json, so make sure we wrap all divisions
static double fixNaN(double d)
{
	return isfinite(d) ? d : 0;
}
void PharmerQuery::computeBenchmarkStats(const vector<QueryResult*>& r, Json::Value& stat)
{
	//we only care about molecules, not conformations - count number of unique
	//molecules with active in the name
	unordered_set<unsigned> seenactives;
	unordered_set<unsigned> seendecoys;

	for(unsigned i = 0, n = r.size(); i < n; i++)
	{
		unsigned mid = r[i]->c->molid;
		if(!seenactives.count(mid) && !seendecoys.count(mid))
		{
			const string& name = r[i]->name;
			if(name.find("active") == string::npos)
			{
				//not active
				seendecoys.insert(mid);
			}
			else
			{
				//active
				seenactives.insert(mid);
			}
		}
	}

	//statistics
	int T = numactives;
	int F = totalmols-numactives;
	int TP = seenactives.size();
	int FP = seendecoys.size();
	int TN = F -seendecoys.size();
	int FN = T-seenactives.size();
	stat["TP"] = TP;
	stat["FP"] = FP;
	stat["TN"] = TN;
	stat["FN"] = FN;

	stat["TPR"] = fixNaN((double)TP/(double)T);
	stat["FPR"] = fixNaN((double)FP/(double)F);

	stat["ACC"] = fixNaN((TP+TN)/(double)(T+F));
	stat["F1"] = fixNaN(2*TP/(double)(2*TP+FP+FN));

	stat["EF"] = fixNaN(( TP/(double)(TP+FP) ) / (T/(double)totalmols));
}


void PharmerQuery::setDataJSON(const DataParameters& dp, Json::Value& data)
{
	vector<QueryResult*> r;
	bool notdone = getResults(dp, r);

	data["draw"] = dp.drawCode;
	data["recordsTotal"] = (unsigned) results.size();
	data["recordsFiltered"] = (unsigned) results.size();
	data["finished"] = !notdone;
	data["data"].resize(0); //make empty array

	if(!notdone && numactives > 0)
	{
		//when all done, if we are a benchmark library, compute stats
		//for all results
		vector<QueryResult*> all;
		DataParameters allparam;
		//need to retain original sort
		allparam.sort = dp.sort;
		allparam.reverseSort = dp.reverseSort;
		allparam.extraInfo = true; //need names
		getResults(allparam, all);
		computeBenchmarkStats(all, data["benchmark"]);
		getResults(dp, r); //resort
	}

	for (unsigned i = 0, n = r.size(); i < n; i++)
	{
		//column ordre is Name,RMSD, Mass,RBnds,location
		data["data"][i][0u] = r[i]->name;
		data["data"][i][1] = r[i]->c->val;
		data["data"][i][2] = r[i]->c->weight;
		data["data"][i][3] = r[i]->c->nRBnds;
		data["data"][i][4] = i + dp.start;
	}
}

static bool locationCompare(const QueryResult* lhs, const QueryResult* rhs)
{
	return lhs->c->location < rhs->c->location;
}

unsigned long PharmerQuery::getLocation(const QueryResult* r,
		std::shared_ptr<PharmerDatabaseSearcher>& db)
{
	unsigned dbid = r->c->location % databases.size();
	unsigned long loc = r->c->location / databases.size();
	db = databases[dbid];
	return loc;
}

//write out all results in sdf format - NOT sorted
void PharmerQuery::outputMols(ostream& out)
{
	loadResults();
	SpinLock lock(mutex);
//copy so we can sort weirdly and release access to results
	vector<QueryResult*> myres = results;
	lock.release();
	sort(myres.begin(), myres.end(), locationCompare);

	PMolReaderSingleAlloc pread;
	MolData mdata;
	vector<ASDDataItem> sddata;
	const char* dataname = "rmsd";
	if(params.isshape) dataname = "sim";

	for (unsigned i = 0, n = myres.size(); i < n && out; i++)
	{
		access();
		std::shared_ptr<PharmerDatabaseSearcher> db;
		unsigned long loc = getLocation(myres[i], db);

		sddata.clear();
		sddata.push_back(
				ASDDataItem(dataname, lexical_cast<string>(myres[i]->c->val)));

		db->getMolData(loc, mdata, pread);
		mdata.mol->writeSDF(out, sddata, myres[i]->c->rmsd);
	}
}

//output single mol in sdf format
void PharmerQuery::outputMol(const QueryResult* mol, ostream& out,
		bool minimize)
{
	access();
	PMolReaderSingleAlloc pread;
	MolData mdata;
	vector<ASDDataItem> sddata;
	const char* dataname = "rmsd";
	if(params.isshape) dataname = "sim";
	sddata.push_back(ASDDataItem(dataname, lexical_cast<string>(mol->c->val)));

	std::shared_ptr<PharmerDatabaseSearcher> db;
	unsigned long loc = getLocation(mol, db);

	db->getMolData(loc, mdata, pread);

//TODO: minimization if requested - openbabel isn't quite where I want it yet..
	mdata.mol->writeSDF(out, sddata, mol->c->rmsd);
}

//output single mol in sdf format, use index into results
void PharmerQuery::outputMol(unsigned index, ostream& out, bool jsonHeader,
		bool minimize)
{
	access();
	SpinLock lock(mutex);
	if (index >= results.size())
		return;
	QueryResult *m = results[index];
	lock.release();

	if (jsonHeader)
		out << "{\"status\": 1}\n";
	outputMol(m, out, minimize);
}

void PharmerQuery::print(ostream& out) const
		{
	Json::Value root;
	convertPharmaJson(root, points);
	out << root;
}

//use the ordering of sides from qt to set box box
static void setBoxMinMax(BoundingBox& box, double min, double max, unsigned i,
		unsigned j, const QueryTriplet& qt)
{
	if (qt.getPoints()[0].index == i)
	{
		if (qt.getPoints()[1].index == j)
		{
			box.setX(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));
		}
		else
		{
			box.setZ(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));
		}
	}
	else if (qt.getPoints()[1].index == i)
	{
		if (qt.getPoints()[0].index == j)
		{
			box.setX(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));
		}
		else
		{
			box.setY(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));
		}
	}
	else
	{
		if (qt.getPoints()[0].index == j)
		{
			box.setZ(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));
		}
		else
			box.setY(ThreePointData::reduceFloat(min),
					ThreePointData::reduceFloat(max));

	}
}

//return true if it is possible for there to be a triangle between points i,j and p where the sides ip and jp are
//constrained by the passed range
bool PharmerQuery::inRange(unsigned i, unsigned j, unsigned p, double minip,
		double maxip, double minjp, double maxjp) const
		{
	const QueryTriplet& trip = triplets[tripIndex[i][j][p]];
	BoundingBox box;

	double min = PharmaPoint::pharmaDist(points[i], points[j])
			- points[i].radius - points[j].radius;
	if (min < 0)
		min = 0;
	double max = PharmaPoint::pharmaDist(points[i], points[j])
			+ points[i].radius + points[j].radius;
	setBoxMinMax(box, min, max, i, j, trip);

	setBoxMinMax(box, minip, maxip, i, p, trip);
	setBoxMinMax(box, minjp, maxjp, j, p, trip);

	return trip.inRange(box);
}

//send smina data, for now do this single threaded - if its a bottleneck
//we can multithread
void PharmerQuery::thread_sendSmina(PharmerQuery *query, stream_ptr out,
		unsigned max)
{
	query->incrementUseCnt();
	SpinLock lock(query->mutex);
	vector<QueryResult*> rescopy(query->results);
	lock.release();

	try
	{

		//limit number of results
		if (max > 0 && rescopy.size() > max)
			rescopy.resize(max);

		//sort by location for sequential access
		sort(rescopy.begin(), rescopy.end(), locationCompare);

		PMolReaderSingleAlloc pread;
		MolData mdata;

		for (unsigned i = 0, n = rescopy.size(); i < n && *out; i++)
		{
			query->access();
			QueryResult *r = rescopy[i];
			std::shared_ptr<PharmerDatabaseSearcher> db;
			unsigned long loc = query->getLocation(r, db);

			//extract and zip rmsd transform
			stringstream str; //need to gzip rotation info
			iostreams::filtering_stream<iostreams::output> filter;
			filter.push(iostreams::gzip_compressor());
			filter.push(str);
			Matrix3d rotmat = r->c->rmsd.rotationMatrix().cast<double>();
			Vector3d trans = r->c->rmsd.translationVector().cast<double>();

			//smina expects things in row major order..
			rotmat.transposeInPlace();

			filter.write((char*) rotmat.data(), sizeof(double) * 9);
			filter.write((char*) trans.data(), sizeof(double) * 3);
			filter.flush();
			filter.pop();
			filter.pop();

			out->write(str.str().c_str(), str.str().length());
			db->getSminaData(loc, *out);

		}
	} catch (...) //don't let exceptions mess up usecnt or crash server
	{
	  std::cerr << "Exception in sendsmina\n";
	}
	query->decrementUseCnt();
}
