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
 * TripletFingerprint.h
 *
 *  Created on: Nov 2, 2010
 *      Author: dkoes
 *
 *  Compute, store, and manipulate triplet fingerprints, which represent
 *  all the other pharmacophore points of the molecule relative to a given
 *  triplet.
 */

#ifndef PHARMITSERVER_TRIPLETFINGERPRINT_H_
#define PHARMITSERVER_TRIPLETFINGERPRINT_H_
#include <vector>

#include "pharmarec.h"
#include "SimpleFingers.h"
#include "CommandLine2/CommandLine.h"

using namespace std;

extern cl::opt<unsigned> BloomBitsLarge;
extern cl::opt<unsigned> BloomBitsSmall;


typedef unsigned long uint128_t __attribute__((mode(TI)));
typedef long int128_t __attribute__((mode(TI)));

struct TripletFingerprint
{
	uint128_t f1;
	uint128_t f2; //fingerprint bits - 256

	TripletFingerprint() :
		f1(0), f2(0)
	{
	}
	~TripletFingerprint()
	{
	}

	//set all the bits
	//return simplefinger
	void set(unsigned i, unsigned j, unsigned k, const vector<PharmaPoint>& points);

	//all bits set in rhs are set in this
	bool contains(const TripletFingerprint& rhs) const
	{
		return ((f1 & rhs.f1) == rhs.f1) && ((f2 & rhs.f2) == rhs.f2);
	}

	bool overlaps(const TripletFingerprint& rhs) const
	{
		return (f1 & rhs.f1) || (f2 & rhs.f2);
	}

	unsigned firstNBits(unsigned n) const
	{
		unsigned mask = (1 << n) - 1;
		return mask & f1;
	}

	void setAll()
	{
		int128_t ones = -1;
		f1 = ones;
		f2 = ones;
	}

	bool isZero() const
	{
		return !f1 && !f2;
	}

	bool isFull() const
	{
		return !(~f1) && !(~f2);
	}

	bool getBit(unsigned pos) const
	{
		uint128_t one = 1;
		if (pos < 128)
			return f1 & (one << pos);
		else
			return f2 & (one << (pos - 128));
	}

	static void bloom(unsigned long val, uint128_t& f1, uint128_t& f2, unsigned bits);

	TripletFingerprint& operator|=(const TripletFingerprint& rhs)
	{
		f1 |= rhs.f1;
		f2 |= rhs.f2;
		return *this;
	}

	TripletFingerprint& operator&=(const TripletFingerprint& rhs)
	{
		f1 &= rhs.f1;
		f2 &= rhs.f2;
		return *this;
	}

	TripletFingerprint& operator^=(const TripletFingerprint& rhs)
	{
		f1 ^= rhs.f1;
		f2 ^= rhs.f2;
		return *this;
	}

	TripletFingerprint operator|(const TripletFingerprint& rhs) const
	{
		TripletFingerprint ret;
		ret.f1 = f1 | rhs.f1;
		ret.f2 = f2 | rhs.f2;
		return ret;
	}

	TripletFingerprint operator&(const TripletFingerprint& rhs) const
	{
		TripletFingerprint ret;
		ret.f1 = f1 & rhs.f1;
		ret.f2 = f2 & rhs.f2;
		return ret;
	}

	bool operator==(const TripletFingerprint& rhs) const
	{
		return f1 == rhs.f1 && f2 == rhs.f2;
	}

	unsigned bitcnt() const
	{
		//kernigans alg
		unsigned int c; // c accumulates the total bits set in v
		uint128_t v = f1;
		for (c = 0; v; c++)
		{
			v &= v - 1; // clear the least significant bit set
		}
		v = f2;
		for (; v; c++)
		{
			v &= v - 1;
		}

		return c;
	}

	void printBin() const
	{
		uint128_t v = f1;
		for (unsigned i = 0; i < 128; i++)
		{
			cout << (int) (v & 1);
			v >>= 1;
		}
		v = f2;
		for (unsigned i = 0; i < 128; i++)
		{
			cout << (int) (v & 1);
			v >>= 1;
		}
	}
	static const double MAXDISTS[];
	static const double DISTSPACES[];
	static const unsigned NUMDISTSPACES;
	static const unsigned MAXPHARMA;
};

//tally up counts for various choices of threshold for each pharmacophore feature
class ThresholdComputer
{
	unsigned steps;
	double incr;

	vector<vector<vector<unsigned long> > > counts; //indexed by pharma, then by threshold choice, then by equivalence class

	vector<vector<unsigned> > curr; //bit mask for current, indexed by pharma then threshold choice

public:

	ThresholdComputer(unsigned nump, unsigned s, double inc) :
		steps(s), incr(inc), counts(nump, vector<vector<unsigned long> > (s,
				vector<unsigned long> (4, 0))), curr(nump, vector<unsigned> (s,
				0))
	{

	}

	//add bits of a single point
	void addVal(unsigned p, bool chiral, double distance)
	{
		assert(p < curr.size());
		for (unsigned i = 0; i < steps; i++)
		{
			double threshold = i * incr;
			//only set a bit if greater than threshold
			if (distance > threshold)
			{
				curr[p][i] |= 1 << chiral;
			}
		}
	}

	//count classes within accumulated fingerprint
	void finishFinger()
	{
		for (unsigned i = 0, n = curr.size(); i < n; i++)
		{
			for (unsigned j = 0, nt = curr[i].size(); j < nt; j++)
			{
				unsigned val = curr[i][j];
				assert(val < 4);
				counts[i][j][val]++;

				curr[i][j] = 0;
			}
		}
	}

	//for each pharma, find the threshold where the minimum of |00|,|01|+|11|, and |10|+|11| is maximized
	//since these are the collections that will be searched
	void printBestThresholds()
	{
		for (unsigned p = 0, np = counts.size(); p < np; p++)
		{
			unsigned long maxmin = 0;
			unsigned bestt = 0;
			for (unsigned t = 0, nt = counts[p].size(); t < nt; t++)
			{
				unsigned long min = std::min(counts[p][t][0], counts[p][t][1]
						+ counts[p][t][3]);
				min = std::min(min, counts[p][t][2] + counts[p][t][3]);
				if (min > maxmin)
				{
					maxmin = min;
					bestt = t;
				}
			}

			cout << "THRESH " << p << " " << bestt * incr << " "
					<< counts[p][bestt][0] << " " << counts[p][bestt][1] << " "
					<< counts[p][bestt][2] << " " << counts[p][bestt][3] << "\n";
		}
	}
};

extern ThresholdComputer thresholdComputer;

#endif /* PHARMITSERVER_TRIPLETFINGERPRINT_H_ */
