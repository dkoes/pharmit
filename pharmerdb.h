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
 * pharmerdb.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 *
 *     	Creation and search routines for the specially formatted pharmer database.
 */

#ifndef PHARMERDB_H_
#define PHARMERDB_H_


#include <iostream>
#include <string>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <sstream>
#include "boost/filesystem.hpp"
#include "boost/thread.hpp"

#include "pharmarec.h"
#include "ThreePointData.h"
#include "FloatCoord.h"
#include "PMol.h"
#include "MMappedRegion.h"
#include "ThreadCounter.h"
#include "BoundingBox.h"
#include "MTQueue.h"
#include "TripleIndexer.h"
#include "Triplet.h"
#include <openbabel/mol.h>
#include "CommandLine2/CommandLine.h"
#include "SPSCQueue.h"
#include "SimpleFingers.h"
#include "MolFilter.h"

using namespace std;
using namespace boost;
using namespace OpenBabel;

class TripletMatches;
extern cl::opt<bool> Quiet;

struct PointDataFile
{
	string name;
	FILE *file;
	long length;
	unsigned phTrip;
	PointDataFile(): file(NULL), length(0),  phTrip(0) {}
	PointDataFile(const string& nm, unsigned p): name(nm), file(NULL), length(0), phTrip(p)
	{
	}

	bool isValid() const { return file != NULL; }

	//only open file if actually written to
	void write(const void *ptr, size_t size, size_t nmemb)
	{
		if(nmemb == 0) return;
		if(size == 0) return;
		if(file == NULL)
		{
			file = fopen(name.c_str(), "w+");
			assert(file);
#ifdef POSIX_FADV_SEQUENTIAL
			posix_fadvise(fileno(file), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
		}

		fwrite(ptr, size, nmemb, file);
		length += size*nmemb;
	}

	void close()
	{
		if(file) fclose(file);
		file = NULL;
	}


};

struct MolDataHeader
{
		unsigned molID; /* unique identifier for mol */
		float molWeight;
		unsigned short confidx; /* which conformation*/

		MolDataHeader(): molID(0), molWeight(0),confidx(0)
		{
		}

		MolDataHeader(unsigned mid, float mw, unsigned short cid): molID(mid), molWeight(mw),  confidx(cid)
		{
		}
} __attribute__((__packed__));


/* Assembles all the data for a molecule and formats into a buffer for writing. */
class MolDataCreator
{
	struct ConfCreator
	{
		MolDataHeader header;
		string molData; //mol representation

		unsigned byteSize() //return number of bytes this takes up
		{
			return sizeof(header) + molData.size();
		}

		ConfCreator(unsigned mid, float mw, unsigned cid, OBMol& m);

		unsigned write(unsigned char *buf);
		bool isValid() const { return molData.length() > 0; }
	};

	const Pharmas& pharmas;
	const TripleIndexer& tindex;
	vector<ConfCreator> confs;
	vector<unsigned> confOffsets; //offset from start of mol for each conformer
	vector<vector<ThreePointData> > pdatas; //indexed by triplet type

	unsigned mid;
	unsigned char *buffer; //raw data for writing
	unsigned bufferSize;

	unsigned numConfs;

	void processMol(OBMol& mol, double mWeight, unsigned mid);
	void createBuffer();

	static unsigned maxIndex;


public:

	MolDataCreator(const Pharmas& ps, const TripleIndexer& t, OBMol& mol, double mw, unsigned m):
		pharmas(ps), tindex(t),
		pdatas(t.size()),
		mid(m), buffer(NULL), bufferSize(0), numConfs(0)
	{
		processMol(mol, mw, mid);
		createBuffer();
	}
	virtual ~MolDataCreator() { if(buffer) delete [] buffer; }
	unsigned write(FILE *molData, vector<PointDataFile>& pointDataFiles);
	unsigned NumPoints() const
	{
		unsigned ret = 0;
		for(unsigned i = 0, n = pdatas.size(); i < n; i++)
			ret += pdatas[i].size();
		return ret;
	}
	unsigned NumConfs() const { return numConfs; }
	static unsigned MaxIndex() { return maxIndex; }

	const vector<unsigned>& ConfOffsets() const { return confOffsets;}
};

//data for a single conformation, include full mol since this turns out to
//be faster than more compressed layouts
struct MolData
{
	PMol* mol; //just the conformation
	unsigned mid; //mol identifier
	float molWeight;
	unsigned short confidx; //which conformation

	MolData(): confidx(0) {}

	bool read(FILE *molData, PMolReader& areader);
	bool read(FILE *molData, unsigned long location, PMolReader& areader);
	bool read(unsigned char *molData, unsigned long location, PMolReader& areader);
	void readDataOnly(unsigned char *molData, unsigned long location);
	void clear();
};




struct GeoKDPage;


//functions to use when sorting point keys
extern bool comparePointDataX(const ThreePointData& lhs, const ThreePointData& rhs);
extern bool comparePointDataY(const ThreePointData& lhs, const ThreePointData& rhs);
extern bool comparePointDataZ(const ThreePointData& lhs, const ThreePointData& rhs);

typedef bool (*pointDataCompare)(const ThreePointData& lhs, const ThreePointData& rhs);

struct SplitInfo
{
	SplitType type;
	pointDataCompare func;
	unsigned short min;
	unsigned short max;

	SplitInfo() :
		type(NoSplit), func(NULL), min(0), max(0)
	{

	}

	SplitInfo(SplitType t, unsigned short n, unsigned short x) :
		type(t), min(n), max(x)
	{
		switch (type)
		{
		case SplitXFixed:
			func = comparePointDataX;
			break;
		case SplitYFixed:
			func = comparePointDataY;
			break;
		case SplitZFixed:
			func = comparePointDataZ;
			break;
		default:
			abort();
		}
	}

	//return the value from p that we're splitting on
	unsigned short splitVal(const ThreePointData& p) const
	{
		switch (type)
		{
		case SplitXFixed:
			return p.l1;
		case SplitYFixed:
			return p.l2;
		case SplitZFixed:
			return p.l3;
		default:
			abort();
		}
	}

	unsigned short getMin(const BoundingBox& box) const
	{
		switch (type)
		{
		case SplitXFixed:
			return box.minx;
		case SplitYFixed:
			return box.miny;
		case SplitZFixed:
			return box.minz;
		default:
			abort();
		}
	}

	unsigned short getMax(const BoundingBox& box) const
	{
		switch (type)
		{
		case SplitXFixed:
			return box.maxx;
		case SplitYFixed:
			return box.maxy;
		case SplitZFixed:
			return box.maxz;
		default:
			abort();
		}
	}
};

//a node in the internal kd tree of a page
struct GeoKDPageNode
{
	BoundingBox box; //bounding box of this entire node
	unsigned short splitVal; // value we are splitting on
	SplitType splitType:8; //how we are splitting, nosplit means we are a leaf
	//some extra space here if we can think of a use for it
	unsigned long splitData; //position within pointData where the split is

	GeoKDPageNode(): splitVal(0), splitType(NoSplit), splitData(0)
	{
	}
};

#define SPLITS_PER_GEOPAGE (256)
#define MAX_POINTDATAS_PER_LEAF (10)
#define MIN_POINTDATAS_PER_PAGE (2*SPLITS_PER_GEOPAGE)
#define MAX_UNIQUE_POINTS_PER_LEAF (2)

struct GeoKDPage
{
	GeoKDPageNode nodes[SPLITS_PER_GEOPAGE]; //index 0 is unused and can be overloaded
	unsigned long nextPages[SPLITS_PER_GEOPAGE];

	GeoKDPage()
	{
		memset(nodes, 0, sizeof(nodes));
		memset(nextPages, 0, sizeof(nextPages));
	}
};


class QueryPoint;

enum Stats
{
	NumMols,
	NumConfs,
	NumDbPoints,
	NumUniquePoints,
	NumInternalPages,
	LastStat
};


typedef vector<vector<PharmaPoint> > MCPoints;

#define MOLQ_CHUNK_SIZE (32)
#define LENGTH_BINS (32)
#define LENGTHDIV (500)

//interface to an anchor oriented database
class PharmerDatabaseCreator
{
private:
	filesystem::path dbpath; //path to database directory

	FILE *info; //small amount of metadata
	FILE *lookup; //human readable pharmacophore classes for my own benefit
	FILE *molData; //flat file of molecular data
	FILE *binData; //binned cnts of lengths
	FILE *midList; //index of mids corresponding to sequential index

	FILE *sminaIndex; //map from mol location to location in sminadata
	FILE *sminaData; //smina formated molecule

	vector<PointDataFile> pointDataFiles; //just pointdata objects; separate library for every pharma combo; indexed by triplet index
	MMappedRegion<ThreePointData> *pointDataArrays;

	vector<FILE*> geoDataFiles; //spatial indexes; indexed  by trip index

	//coarsely track length distributions using simple binning
	vector<array<array<array<unsigned, LENGTH_BINS>, LENGTH_BINS>, LENGTH_BINS> > binnedCnts;

	void initPointDataArrays();
	void initializeDatabases();

	void doSplitInPage(unsigned pharma, FILE *geoFile, GeoKDPage& page, unsigned pos, ThreePointData *start,
			ThreePointData *end, ThreePointData *begin, unsigned depth);

	unsigned long doSplitNewPage(unsigned pharma, FILE *geoFile, ThreePointData *start, ThreePointData *end, ThreePointData *begin, unsigned depth);

	void incrementBinCnt(const ThreePointData& t, unsigned pclass);
	void createIJKSpatialIndex(int p);

	void generateAtomData();

	static void thread_start_createPharmaSpatialIndex(PharmerDatabaseCreator *db);
	static void thread_start_doSplitInPage(PharmerDatabaseCreator *db, unsigned pharma, FILE *geoFile, GeoKDPage& page, unsigned pos,
			ThreePointData *start, ThreePointData *end, ThreePointData *begin, unsigned depth);

	void writeMIDs();

	unsigned long stats[LastStat];

	const Pharmas& pharmas;
	TripleIndexer tindex;

	shared_mutex fileAccessLock;
	long pdatasFitInMemory;
	ThreadCounter tcnt;

	MTQueue< vector<MolDataCreator*> > molDataWorkQ;

	vector<unsigned> mids;
	unsigned nThreads;

public:
	PharmerDatabaseCreator(const Pharmas& ps, const filesystem::path& dbp,
			 unsigned nt) :
				dbpath(dbp), info(NULL), molData(NULL), midList(NULL), pointDataArrays(NULL),
				  pharmas(ps), tindex(ps.size()),
				tcnt(nt), molDataWorkQ(32),
				nThreads(nt)
	{
		memset(&stats, 0, sizeof(stats));

		array<array<array<unsigned, LENGTH_BINS>, LENGTH_BINS>, LENGTH_BINS> zero;
		memset(zero.c_array(), 0, LENGTH_BINS*LENGTH_BINS*LENGTH_BINS*sizeof(unsigned));
		binnedCnts.resize(tindex.size(), zero);
		//create databases

		initializeDatabases();

		//amount of physical mem
#if defined(_SC_PHYS_PAGES)
		unsigned long memsz = sysconf (_SC_PHYS_PAGES) * sysconf (_SC_PAGESIZE);
		//allow 2/3 to be used to cache pointdatas
		pdatasFitInMemory = 2*memsz/sizeof(ThreePointData)/3;
#else
		pdatasFitInMemory = 1000000; //eh, os x
#endif
	//	cout << "PointData size " << sizeof(ThreePointData) << "\n";
		//cout << "GeoKDPage size " << sizeof(GeoKDPage) << "\n";
		//		cout << "pagesize " << sysconf(_SC_PAGESIZE) << "\n";

	}

	//ensure that all data is flushed
	~PharmerDatabaseCreator()
	{

		if (info)
			fclose(info);
		if (molData)
			fclose(molData);
		if(midList)
			fclose(midList);
		for(unsigned i = 0, n = pointDataFiles.size(); i < n; i++)
		{
			pointDataFiles[i].close();
		}
		for(unsigned i = 0, n = geoDataFiles.size(); i < n; i++)
		{
			if(geoDataFiles[i]) fclose(geoDataFiles[i]);
		}

		if(pointDataArrays != NULL) delete [] pointDataArrays;
	}

	void addMolToDatabase(OBMol& mol, double weight);
	//add all the conformations from infile into the database
	void addMolsToDatabase(istream& infile,
			OBFormat *format, unsigned start, unsigned stride,
			unsigned long readBytes, unsigned long lastByte);
	//create the spatial index
	//requirs pointdata to have been filled out
	void createSpatialIndex();

	void writeStats();

	//return number in database
	unsigned numMolecules() const
	{
		return stats[NumMols];
	}
	unsigned numConformations() const
	{
		return stats[NumConfs];
	}
};


//stuff that doesn't change through the recursive search of the
//spatial index
struct QueryInfo
{
	const QueryTriplet& triplet;
	const GeoKDPage *pages;
	const ThreePointData *data;
	const unsigned index;
	unsigned which; //which triplet ordering
	TripletMatches& M;
	volatile bool& stopEarly;
	QueryInfo(const QueryTriplet& trp, unsigned w, const GeoKDPage *p, const ThreePointData *d, unsigned i,
			TripletMatches& m, bool& stop):
				triplet(trp),  pages(p), data(d), index(i), which(w), M(m), stopEarly(stop)
	{

	}
};

//interface to an anchor oriented database - search only
class PharmerDatabaseSearcher
{
	typedef pair<unsigned long, unsigned long> StartEnd;

	Pharmas pharmas;
	TripleIndexer tindex;

	filesystem::path dbpath; //path to database directory

	FILE *info;
	MMappedRegion<unsigned char> molData;
	MMappedRegion<ThreePointData>  * tripletDataArrays;
	MMappedRegion<GeoKDPage> * geoDataArrays;
	MMappedRegion<unsigned> midList;

	MMappedRegion< pair<unsigned long, unsigned long> > sminaIndex; //maps moldata location to sminadata
	MMappedRegion<char> sminaData;

	unsigned goodChunkSize;

	void initializeDatabases();

	MMappedRegion<unsigned> binnedCnts;

	unsigned long stats[LastStat];
	bool valid;

	void queryProcessPoints(QueryInfo& t,
				unsigned long startLoc, unsigned long endLoc);
	void queryIndex(QueryInfo& t, const GeoKDPage *page,
			unsigned pos, unsigned long startLoc,
			unsigned long endLoc);

	unsigned getBinCnt(unsigned pclass, unsigned i, unsigned j, unsigned k);

	unsigned emptyCnt;
	unsigned processCnt;
	unsigned matchedCnt;
public:
	PharmerDatabaseSearcher(const filesystem::path& dbp) :
		dbpath(dbp), info(NULL),  valid(false), emptyCnt(0), processCnt(0), matchedCnt(0)
	{
		goodChunkSize = 100000; //what's life without a little magic (numbers)? - should probably be related to cache size
		memset(&stats, 0, sizeof(stats));

		//create databases
		initializeDatabases();
	}

	~PharmerDatabaseSearcher()
	{
		if(tripletDataArrays != NULL) delete [] tripletDataArrays;
		if(geoDataArrays != NULL) delete [] geoDataArrays;
	}

	//return number in database
	unsigned numMolecules() const
	{
		return stats[NumMols];
	}
	unsigned long numConformations() const
	{
		return stats[NumConfs];
	}

	//translate file-offset mid to true (lowest of cmpd) mid
	unsigned getBaseMID(unsigned lmid) const
	{
		unsigned len = midList.length();
		if(len == 0)
			return lmid;
		if(lmid >= len)
			return midList[len-1];
		return midList[lmid];
	}

	//generate ranking and return index of best triplet
	unsigned rankTriplets(const vector<QueryTriplet>& triplets, vector<double>& ranking);

	//put all matching triplets into Q
	void generateTripletMatches(const vector<vector<QueryTriplet> >& triplets, TripletMatches& Q, bool& stopEarly);

	//get mol data, a single conformation, at location
	bool getMolData(unsigned long location, MolData& mdata, PMolReader& reader)
	{
		return mdata.read(molData.begin(), location, reader);
	}

	//get mol info, don't parse in mol
	void getMolData(unsigned long location, MolData& mdata)
	{
		mdata.readDataOnly(molData.begin(), location);
	}

	void getSminaData(unsigned long location, ostream& out);


	const Pharmas& getPharmas() const { return pharmas; }

	string getName() { return dbpath.string(); }
	bool isValid() { return valid;}

};


#endif /* PHARMERDB_H_ */
