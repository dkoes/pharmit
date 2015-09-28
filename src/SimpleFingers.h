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
 * SimpleFingers.h
 *
 *  Created on: Jan 14, 2011
 *      Author: dkoes
 *
 *      Classes and routines for computing and manipulating triplet relative fingerprints
 *      which is what I call a very coarse finger print of the entire molecule defined
 *      relative to a specific triplet.  The coarseness makes it more appropriate binning
 *      molecule, although of course they will not be evenly distributed.
 */

#ifndef PHARMITSERVER_SIMPLEFINGERS_H_
#define PHARMITSERVER_SIMPLEFINGERS_H_

#if 0
This didn't have anywhere close to the hoped for effect. In most cases doesn't
even offset the increased fragementation.

#include "pharmarec.h"
#include <bm/bm.h>

typedef bm::bvector<bm::standard_allocator> bvect;

//for each bit in a simplefinger, precompute all the possible bitsets
class SimpleFingerCollection
{
	vector<unsigned> pharmaStart;
	vector<bvect> checkIfSet; //indexed by simplefingerprint, all collections overlapping those bits
	unsigned N; //number of bits in SF

public:

	SimpleFingerCollection(): N(0) {}

	SimpleFingerCollection(const Pharmas& pharmas): pharmaStart(pharmas.size())
	{
		setPharmas(pharmas);
	}

	void setPharmas(const Pharmas& pharmas)
	{
		pharmaStart.resize(pharmas.size());
		N = 0;
		for(unsigned i = 0, n = pharmas.size(); i < n; i++)
		{
			pharmaStart[i] = N;
			N += pharmas[i]->simpleFingerBits;
		}

		//store what numbers have each bit set for easy manipulation
		unsigned max = 1<<N;
		checkIfSet.assign(max, bvect(max));

		for(unsigned i = 0; i < max; i++)
		{
			for(unsigned v = 0; v < max; v++)
			{
				if((v&i)) //i overlaps v
				{
					checkIfSet[i].set(v); //v is a collection we must consider for i
				}
			}
		}
	}

	const bvect& getCollections(unsigned b) const
	{
		assert(b < checkIfSet.size());
		return checkIfSet[b];
	}

	unsigned numBits() const { return N; }
	unsigned numCollections() const { return 1 << N; }

	unsigned getPharmaStart(const Pharma *p) const { return pharmaStart[p->index]; }
};

//compute the fingerprint bit of a single point
static unsigned computeSimpleFingerBit(const SimpleFingerCollection& SFC, const Pharma* p, bool underplane, double dist)
{
	unsigned index = SFC.getPharmaStart(p);
	unsigned ret = 0;
	switch(p->simpleFingerBits)
	{
	case 0:
		return ret;
	case 1: //present/not present
		ret |= (1<<index);
		break;
	case 4: //chiral/not chiral and in/out of radius
		//bit 0: underplane and short dist
		//bit 1: underplane and long dist
		//bit 2: overplane and short dist
		//bit 3: overplane and long dist
		if(underplane)
		{
			if(dist <= p->simpleFingerRadius)
				ret |= (1<<(index));
			else
				ret |= (1<<(index+1));
		}
		else
		{
			if(dist <= p->simpleFingerRadius)
				ret |= (1<<(index+2));
			else
				ret |= (1<<(index+3));
		}
		break;
	case 3: //chiral/not chiral and inside radius
		//bit 0: underplane and long dist
		//bit 1: overplane and long dist
		//bit 3: short dist
		if(dist <= p->simpleFingerRadius)
			ret |= 1<<(index+2);
		else if(underplane)
			ret |= 1<<index;
		else
			ret |= 1<<(index+1);
		break;
	case 2: //chiral/not chiral
		if(underplane)
			ret |= 1<<(index);
		else
			ret |= 1<<(index+1);
		break;
	default:
		abort();
	}
	return ret;
}

class SimpleFinger
{
	const SimpleFingerCollection& SFC;

	bvect finalMask; //for whole mol
	unsigned currMask; //for current unfinished point
public:

	SimpleFinger(const SimpleFingerCollection& sfc): SFC(sfc), finalMask(sfc.numCollections()), currMask(0)
	{
		assert(SFC.numBits()< sizeof(currMask)*8);
		finalMask.set(); //start all set since intersect between points
	}


	//add one piece of data for a query point; these all get or'ed together until
	//finishpoint which will find an intersection
	void addQueryPointData(const Pharma* p, bool underplane, double mindist, double maxdist)
	{
		currMask |= computeSimpleFingerBit(SFC, p, underplane, mindist);
		currMask |= computeSimpleFingerBit(SFC, p, underplane, maxdist);
	}


	//interesect all collections required for this point with finalMask, hopefully
	//reducing the set
	void finishPoint()
	{
		finalMask &= SFC.getCollections(currMask);
		currMask = 0;
	}

	//return a bitset indicated what collections are required for a search
	const bvect& getCollections() const
	{
		return finalMask;
	}

};

#endif
#endif /* PHARMITSERVER_SIMPLEFINGERS_H_ */
