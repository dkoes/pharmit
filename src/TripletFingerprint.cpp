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
 * TripletFingerprint.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: dkoes
 */

#include "TripletFingerprint.h"
#include "openbabel/math/vector3.h"
#include "CommandLine2/CommandLine.h"

using namespace OpenBabel;

const double TripletFingerprint::MAXDISTS[2] = {15, 32.0};
const double TripletFingerprint::DISTSPACES[2] = {4.0, 1.0};
const unsigned TripletFingerprint::NUMDISTSPACES = 2;
const unsigned TripletFingerprint::MAXPHARMA = 32;


cl::opt<bool> ComputeThresholds("compute-thresholds", cl::desc("compute good thresholds for binning fingerprints"), cl::Hidden);
cl::opt<unsigned> BloomBitsLarge("bloom-large", cl::desc("Number of bits for coarse discretization fingerprint."), cl::init(6),cl::Hidden);
cl::opt<unsigned> BloomBitsSmall("bloom-small", cl::desc("Number of bits for fine discretization fingerprint."), cl::init(8),cl::Hidden);


ThresholdComputer thresholdComputer(6, 32, .5);

static uint64_t bloomhash(unsigned long val)
{
	//hash int using murmerhash2
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = 8 * m;
	uint64_t k = val;
	k *= m;
	k ^= k >> r;
	k *= m;

	h ^= k;
	h *= m;

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}
//take an integer and set an 8-hash 256-bit bloom filter
//does not clear existing bits of f1 and f2
void TripletFingerprint::bloom(unsigned long val, uint128_t& f1, uint128_t& f2, unsigned bits)
{
	uint64_t h = bloomhash(val);
	uint64_t horig = h;
	assert(bits <= 16);
	//each byte sets a bit in the bloom filter
	uint128_t one = 1;
	for(unsigned i = 0; i < 8 && i < bits; i++)
	{
		unsigned pos = h&0xff;
		h >>= 8;
		if(pos < 128)
			f1 |= (one << pos);
		else
			f2 |= (one << (pos-128));
	}

	if(bits > 8)
	{
		bits -= 8;
		uint64_t h2 = bloomhash(horig);
		for(unsigned i = 0; i < 8 && i < bits; i++)
		{
			unsigned pos = h2&0xff;
			h2 >>= 8;
			if(pos < 128)
				f1 |= (one << pos);
			else
				f2 |= (one << (pos-128));
		}
	}
}

//add points to bloom filter, include i,j,k distance, point type, and chirality
void TripletFingerprint::set(unsigned i, unsigned j, unsigned k, const vector<PharmaPoint>& points)
{
	f1 = 0;
	f2 = 0;
	//compute cross product for chirality
	//k to j and k to i
	vector3 vk(points[k].x, points[k].y, points[k].z);
	vector3 vi(points[i].x, points[i].y, points[i].z);
	vector3 vj(points[j].x, points[j].y, points[j].z);

	vector3 ki = vi-vk;
	vector3 kj = vj-vk;
	vector3 cr = cross(ki, kj);
	bool iszero = cr.length_2() == 0;

	for (unsigned p = 0, n = points.size(); p < n; p++)
	{
		if (p == i)
			continue;
		if (p == j)
			continue;
		if (p == k)
			continue;

		for (unsigned d = 0; d < NUMDISTSPACES; d++)
		{
			//compute a val that is a combination of distances, type
			unsigned long val = 0;

			double dist = trunc(PharmaPoint::pharmaDist(points[i], points[p])
					/ DISTSPACES[d]);
			if (dist > MAXDISTS[d])
				dist = MAXDISTS[d];
			val += dist;
			val *= (MAXDISTS[d] + 1);
			dist = trunc(PharmaPoint::pharmaDist(points[j], points[p])
					/ DISTSPACES[d]);
			if (dist > MAXDISTS[d])
				dist = MAXDISTS[d];
			val += dist;
			val *= (MAXDISTS[d] + 1);
			dist = trunc(PharmaPoint::pharmaDist(points[k], points[p])
					/ DISTSPACES[d]);
			if (dist > MAXDISTS[d])
				dist = MAXDISTS[d];
			val += dist;
			val *= MAXPHARMA;
			val += points[p].pharma->index;

			val *= 2;

			//chirality - pointing in direction of cross product
			vector3 pv(points[p].x, points[p].y, points[p].z);
			pv -= vk;
			bool underplane = false;
			if (!iszero && vectorAngle(cr, pv) > 90)
			{
				underplane = true;
				val++;
			}

//			printf("FINGERC %d,%d,%d %d:%d  %ld  (%f,%f,%f)\n",i,j,k,p,d,val,l1,l2,l3);
			if (ComputeThresholds)
				thresholdComputer.addVal(points[p].pharma->index, underplane,
						pv.length());
			bloom(val, f1, f2, (d == 0 ? BloomBitsLarge : BloomBitsSmall));
		}
	}

	if(ComputeThresholds)
		thresholdComputer.finishFinger();
}
