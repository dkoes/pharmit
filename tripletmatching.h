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
 * triplematching.h
 *
 *  Created on: Aug 6, 2010
 *      Author: dkoes
 *
 *      Classes and routines for performing triplet matching of queries.
 */

#ifndef TRIPLEMATCHING_H_
#define TRIPLEMATCHING_H_


#include "BumpAllocator.h"
#include "ThreePointData.h"
#include <boost/thread.hpp>
#include <vector>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/pool/object_pool.hpp>
#include "Triplet.h"
#include "params.h"
using namespace std;
using namespace boost;



struct PointCoords
{
	short x;
	short y;
	short z;

	PointCoords() :
		x(0), y(0), z(0)
	{
	}
	PointCoords(short a, short b, short c) :
		x(a), y(b), z(c)
	{
	}

	bool operator==(const PointCoords& rhs) const
	{
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}

	bool operator!=(const PointCoords& rhs) const
	{
		return !(*this == rhs);
	}

	bool operator<(const PointCoords& rhs) const
	{
		if(x != rhs.x) return x < rhs.x;
		if(y != rhs.y) return y < rhs.y;
		return z < rhs.z;
	}

	//set to largest coords
	void maximize()
	{
		x = y = z = SHRT_MAX;
	}

	unsigned distSQ(const PointCoords& r) const
	{
		int X = r.x-x;
		int Y = r.y-y;
		int Z = r.z-z;
		return X*X+Y*Y+Z*Z;
	}
};

class TripletMatch;
class TripletMatchInfoArray;

class TripletMatchAllocator
{
	const unsigned qsize; //number of triplets
	const unsigned PMsize; //size of a point match
	BumpAllocator<1024*1024, false> allocator; //1MB chunks, NOT THREAD SAFE
public:
	TripletMatchAllocator(unsigned qsz);
	TripletMatch* newTripletMatch(unsigned mid, const ThreePointData& tdata);
	TripletMatchInfoArray* newTripletMatchInfoArray(unsigned numEl);

	void clear()
	{
		allocator.clear();
	}

	unsigned numChunks() const { return allocator.numChunks(); }
};


//info of a single corresponding triplet
struct TripletMatchInfo
{
	array<PointCoords, 3> coords;
	array<unsigned char, 3> indices;
	unsigned char whichTripOrder:4;
	unsigned char unconnectedNextIndex:2; //0,1,2; which coord is not connected to the next x
	unsigned char unconnectedPrevIndex:2; //0,1,2: which oord is not connected to prev trip

	TripletMatchInfo(const ThreePointData& tdata, unsigned which, unsigned uncon, unsigned puncon): whichTripOrder(which), unconnectedNextIndex(uncon), unconnectedPrevIndex(puncon)
	{
		coords[0] = PointCoords(tdata.x1(), tdata.y1(), tdata.z1());
		coords[1] = PointCoords(tdata.x2(), tdata.y2(), tdata.z2());
		coords[2] = PointCoords(tdata.x3(), tdata.y3(), tdata.z3());
		indices[0] = tdata.i1;
		indices[1] = tdata.i2;
		indices[2] = tdata.i3;
	}

	TripletMatchInfo() {}

	void dump() const
	{
		cout << "( ";
		cout << (int)whichTripOrder << "," << (int)unconnectedNextIndex << " - ";
		for(unsigned i = 0; i < 3; i++)
			cout << (int)indices[i] << " : " << (int)coords[i].x << "," << (int)coords[i].y << "," << (int)coords[i].z << (i == 2 ? "" : "  ; ");
		cout << ")\n";
	}

	//return true if this and x share an edge
	bool isNext(const TripletMatchInfo& x) const
	{
		//return true;
		unsigned i = (unconnectedNextIndex+1)%3;
		unsigned j = (unconnectedNextIndex+2)%3;

		unsigned a = (x.unconnectedPrevIndex+1)%3;
		unsigned b = (x.unconnectedPrevIndex+2)%3;

		if(indices[i] == x.indices[a] && indices[j] == x.indices[b])
			return true;
		if(indices[i] == x.indices[b] && indices[j] == x.indices[a])
			return true;
		return false;
	}
};


#define TMINFO_ARRAY_STARTSIZE (2)
class TripletMatchInfoArray
{
	friend class TripletMatchAllocator;
	TripletMatchInfoArray *next; //for overflow
	unsigned num; //total number of elements, included attached chains
	TripletMatchInfo data[TMINFO_ARRAY_STARTSIZE]; //this size increases exponentially with each level of the chain

	static unsigned allocSize(unsigned numEl);

	//add to array
	void add(const TripletMatchInfo& info, TripletMatchAllocator& alloc, unsigned level)
	{
		unsigned pos = num++;
		unsigned thisSize = TMINFO_ARRAY_STARTSIZE<<level;
		if(pos < thisSize) //space here - most common case
		{
			data[pos] = info;
		}
		else if(pos == thisSize) //must create new array
		{
			TripletMatchInfoArray *newarray = alloc.newTripletMatchInfoArray(thisSize<<1);
			newarray->data[0] = info;
			newarray->num = 1;
			next = newarray;
		}
		else
		{
			//recurse down linkedlist
			next->add(info, alloc,level+1);
		}

	}

	//get element i, recursing through the exponentially increasing levels
	const TripletMatchInfo& get(unsigned i, unsigned level) const
	{
		unsigned thisSize = TMINFO_ARRAY_STARTSIZE<<level;
		if(i < thisSize)
			return data[i];
		return next->get(i-thisSize, level+1); //no bounds checking
	}

	void dump(unsigned level)
	{
		cerr << num << ":";
		unsigned thisSize = TMINFO_ARRAY_STARTSIZE<<level;
		for(unsigned i = 0; i < num && i < thisSize; i++)
		{
			cerr << " ";
			data[i].dump();
		}
		if(next) { cerr << "\n\t"; next->dump(); }
		else cerr << "\n";
	}
public:
	TripletMatchInfoArray():next(NULL), num(0)
	{

	}

	//add to array
	void add(const TripletMatchInfo& info, TripletMatchAllocator& alloc)
	{
		add(info, alloc, 0);
	}

	int size() const
	{
		return num;
	}

	const TripletMatchInfo& operator[](unsigned i) const
	{
		return get(i, 0);
	}

	void dump()
	{
		dump(0);
	}
};


typedef long int128_t __attribute__((mode(TI))) __attribute__((aligned(8)));

//stores all the matching triples for a unique conformer
struct TripletMatch
{
	//these can all probably be reduced to 64bits..
	int128_t	 mustMatch; //indices the current triple must contain two of
	int128_t	nextMustMatch; //indices for the next triple

	unsigned long id:TPD_MOLDATA_BITS; //mol we are matching
	unsigned molid:TPD_MOLID_BITS;

	float weight;
	short currIndex; //index of last matched triplet - must increase sequentially (for required)
	unsigned char numTriplets; //size of query
	unsigned char nRBnds;

	TripletMatchInfoArray matches[]; //indexed by query size
	//space must be reserved by allocator

	TripletMatch(unsigned mid, const ThreePointData& tdata, unsigned numtrip):
		mustMatch(0), nextMustMatch(-1),
		id(tdata.molPos), molid(mid), weight(ThreePointData::unreduceWeight(tdata.weight)),
		currIndex(-1), numTriplets(numtrip), nRBnds(tdata.nrot)
	{
		for (unsigned i = 0; i < numTriplets; i++)
		{
			new (&matches[i]) TripletMatchInfoArray(); //initialize
		}
	}


	bool valid() { return currIndex == numTriplets-1; }

	bool hasValidConnections(const ThreePointData& tdata, const QueryTriplet& trip, int tripIndex)
	{
		static int128_t one = 1;
		if(tripIndex == currIndex+1)
		{
			mustMatch = nextMustMatch;
			nextMustMatch = 0;
			currIndex++;
		}

		if(tripIndex == currIndex)
		{
			unsigned cnt = 0;
			cnt += !!((one << tdata.i1)&mustMatch);
			cnt += !!((one << tdata.i2)&mustMatch);
			cnt += !!((one << tdata.i3)&mustMatch);

			if(cnt < 2)
				return false; //have to correspond to two vertices of previously matched
		}
		else
			return false; //missed a match

		return true;
	}

	bool add(const TripletMatchInfo& tminfo, const QueryTriplet& trip, int tripIndex, TripletMatchAllocator& alloc)
	{
		static int128_t one = 1;
		if(tripIndex == currIndex+1)
		{
			mustMatch = nextMustMatch;
			nextMustMatch = 0;
			currIndex++;
		}

		//check distance between new point and unmatched point of previous triplet
		//the amount of pruning that results compensates for the extra overhead
		if(currIndex > 0)
		{
			for(int distback=0; distback < currIndex; distback++)
			{
				unsigned prev = currIndex - 1 - distback;
				unsigned puncon = trip.getPrevUnconnected();
				unsigned i = 0;
				unsigned n = matches[prev].size();
				for (; i < n; i++)
				{
					const TripletMatchInfo& pmatch = matches[prev][i];

					//immediate successors must line up points
					if(distback > 0 || pmatch.isNext(tminfo))
					{
						unsigned d =
								pmatch.coords[pmatch.unconnectedNextIndex].distSQ(
									tminfo.coords[puncon]);
						if (trip.goodKDistance(d, distback))
							break;
					}
				}
				if (i == n) //no good distances
				{
					return false;
				}
			}
		}
		unsigned badi = tminfo.unconnectedNextIndex;
		unsigned i1 = (badi+1)%3;
		unsigned i2 = (badi+2)%3;

		nextMustMatch |= (one << tminfo.indices[i1]) | (one << tminfo.indices[i2]);

		matches[tripIndex].add(tminfo,alloc);
		return true;
	}

	void dump(int qsz)
	{
		cerr << id << "\n";
		for(int i = 0; i < qsz; i++)
			matches[i].dump();
	}
};


//an open address hash table
//intentionally not thread-safe (parrallel matching doesn't work so well, so ditch the overhead)
class TripletMatchHash
{
	TripletMatch **table;
	unsigned prime_index;
	unsigned long table_size;
	unsigned long num_elements;

	TripletMatchAllocator& alloc;

	const static unsigned long  prime_list[];
	static unsigned const endPrimeIndex;

	void grow();
	unsigned long hash(unsigned long id)
	{
		//hash int using murmerhash2
		const uint64_t m = 0xc6a4a7935bd1e995;
		const int r = 47;

		uint64_t h = 8 * m;
		uint64_t k = id;
		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;

		h ^= h >> r;
		h *= m;
		h ^= h >> r;

		return h % table_size;
	}

	unsigned long getPos(unsigned long id)
	{
		unsigned long pos = hash(id);
		while (table[pos] != NULL)
		{
			if (table[pos]->id == id)
				return pos;
			pos++;
			if (pos >= table_size)
				pos = 0;
		}
		return pos;
	}

public:
	TripletMatchHash(TripletMatchAllocator& a): prime_index(12), table_size(prime_list[prime_index]), num_elements(0), alloc(a)
	{
		table = (TripletMatch**)malloc(sizeof(TripletMatch*)*table_size);
		memset(table, 0, sizeof(TripletMatch*)*table_size);
	}

	~TripletMatchHash()
	{
		free(table);
	}

	//return non-null if pm already exists
	TripletMatch* exists(unsigned long id)
	{
		return table[getPos(id)];
	}

	//allocates a triplet match and returns it
	//alternatively, if returns an already created triplet match
	TripletMatch* create(unsigned mid, const ThreePointData& tdata)
	{
		unsigned long id = tdata.molPos;

		unsigned long pos = hash(id);
		while (table[pos] != NULL)
		{
			if (table[pos]->id == id)
			{
				return table[pos];
			}
			pos++;
			if (pos >= table_size)
				pos = 0;
		}
		TripletMatch *match = alloc.newTripletMatch(mid, tdata);
		table[pos] = match;
		num_elements++;

		if(num_elements >= table_size/8)
		{
			grow();
		}
		return match;
	}

	unsigned long size() const { return table_size; };

	TripletMatch* operator[](unsigned long pos) { return table[pos]; }
};

//keep track of all our  matches in a thread safe manner; however
// all threads must finish matching one query point before moving onto the next
// one; also I find that splitting the database up and having separate threads
// matching different databases makes better use of parallelism
class TripletMatches
{
	const QueryParameters& params;
	TripletMatchAllocator& alloc; //separating allocator from matcher allows results to outlive the matcher
	TripletMatchHash seenMatches;

	typedef union {
		unsigned long head;
		long padding[16]; //eliminate false sharing
	} HeadType;

	vector<HeadType> heads;

	vector<unsigned long> counts;
	unsigned curIndex;
	unsigned qsize;

public:
	TripletMatches(TripletMatchAllocator& a, const QueryParameters& p, unsigned sz, unsigned nthreads) : params(p), alloc(a), seenMatches(a),
	 heads(nthreads), counts(sz, 0), curIndex(0), qsize(sz)
	{
		for(unsigned i = 0; i < nthreads; i++)
			heads[i].head = i;
	}

	virtual ~TripletMatches()
	{
	}

	//register that we are now processing the next triplet in the query
	//return true if there is still a possibility of matching something
	bool nextIndex()
	{
		bool stillgood = counts[curIndex];
		curIndex++;
		assert(curIndex <= qsize); /* equal at end */
		return stillgood;
	}

	//add point, return true if triplet is actual valid
	bool add(unsigned mid, const ThreePointData& tdata, const QueryTriplet& trip, unsigned which)
	{
		//see if we've seen this match already
		unsigned long key = tdata.molPos;
		TripletMatch *match = NULL;

		//check against query params
		if(tdata.nrot < params.minRot)
			return false;
		if(tdata.nrot > params.maxRot)
			return false;

		if(tdata.weight < params.reducedMinWeight)
			return false;
		if(tdata.weight > params.reducedMaxWeight)
			return false;

		//do fast checks before calling isMatch
		if(curIndex > 0)
		{
			//must already exist
			match = seenMatches.exists(key);
			if(match == NULL)
				return false;
			if(!match->hasValidConnections(tdata, trip, curIndex))
				return false;
			if(!trip.isMatch(tdata))
				return false;
		}
		else
		{
			if(!trip.isMatch(tdata))
				return false;
			//create the match
			match = seenMatches.create(mid, tdata);
		}
		//create match info
		assert(match != NULL);
		TripletMatchInfo info(tdata, which, trip.getNextUnconnected(), trip.getPrevUnconnected());
		if(match->add(info, trip, curIndex, alloc))
		{
			counts[curIndex]++;
			return true;
		}
		return false;
	}

	//pop from the t'th queue
	bool pop(TripletMatch*& match, unsigned t)
	{
		unsigned nthreads = heads.size();
		assert(t < nthreads);
		while(heads[t].head < seenMatches.size())
		{
			unsigned long pos = heads[t].head;
			heads[t].head += nthreads;

			if(seenMatches[pos] != NULL)
			{
				match = seenMatches[pos];
				if(match->valid())
				{
					return true;
				}
			}
		}
		return false;
	}


	//number of matches in all
	size_t total() { return counts.back(); }
	size_t empty(unsigned t) {return seenMatches.size() <= heads[t].head; }

	void dumpCnts()
	{
		for(unsigned i = 0, n = counts.size(); i < n; i++)
			cout << counts[i] << " ";
		cout << "\n";
	}

	unsigned getQSize() const { return qsize; }

	unsigned numChunks() const { return alloc.numChunks(); }
};


#endif /* TRIPLEMATCHING_H_ */
