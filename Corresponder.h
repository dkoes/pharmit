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
 * Corresponder.h
 *
 *  Created on: Aug 10, 2010
 *      Author: dkoes
 */

#ifndef CORRESPONDER_H_
#define CORRESPONDER_H_

#include "pharmerdb.h"
#include "Triplet.h"
#include "cors.h"
#include "MTQueue.h"
#include "RMSD.h"
#include "params.h"
#include "CommandLine2/CommandLine.h"
#include "Excluder.h"

extern cl::opt<bool> UnWeightedRMSD;
typedef long int128_t __attribute__((mode(TI)));

class Corresponder
{
	shared_ptr<PharmerDatabaseSearcher> dbptr;
	unsigned dbid;
	unsigned numdbids;
	const vector<PharmaPoint>& points;
	const vector<vector<QueryTriplet> >& triplets;
	TripletMatches& inQ;
	CorAllocator& alloc;
	unsigned threadQ;
	MTQueue<CorrespondenceResult*>& resultQ;
	const QueryParameters& qparams;
	const Excluder& excluder;

	CorrespondenceResult *tmpresult; //used for building up result
	TripletMatch *tm; //current match being processed
	vector<double> pointCoords;
	vector<double> weights;
	vector<double> molCoords; //coordinates of current match
	vector<double> tmpCoords;

	unsigned thisConfCnt;
	unsigned processedCnt;
	unsigned matchedCnt;

	bool& stopEarly;


	//return true if result will fall in exclusion zone - this can be expensive
	//since it must load the molecular data
	bool isExcluded(CorrespondenceResult *result)
	{
		MolData mdata;
		PMolReaderSingleAlloc pread;
		dbptr->getMolData(tm->id, mdata, pread);

		vector<FloatCoord> coords;
		mdata.mol->getCoords(coords, result->rmsd);
		for (unsigned i = 0, n = coords.size(); i < n; i++)
		{
			if (excluder.isExcluded(coords[i]))
				return true;
		}

		return false;
	}

	//recursively generate all possible legal correspondences
	//return true if should keep going
	bool generate(int pos, int128_t alreadyMatched)
	{
		if (pos < 0) //base case
		{
			//comptue rmsd
			assert(
					pointCoords.size() == molCoords.size() && pointCoords.size() % 3 == 0);
			unsigned n = weights.size();
			assert(n == pointCoords.size()/3);

			if (UnWeightedRMSD)
				tmpresult->rmsd = calculateRMSD(&pointCoords[0], &molCoords[0],
						n);
			else
				tmpresult->rmsd = calculateRMSD(&pointCoords[0], &molCoords[0],
						&weights[0], n);
			if (UnWeightedRMSD || tmpresult->rmsd.value() <= 1.0)
			{
				//rmsd is weighted so the max legal value is 1.0
				//check distance of each point
				tmpCoords = molCoords;
				tmpresult->rmsd.reorient(n, &tmpCoords[0]);
				bool isgood = true;
				double rms = 0;
				for (unsigned i = 0; i < n && isgood; i++)
				{
					double val = 0;
					for (unsigned j = 0; j < 3; j++)
					{
						double v = tmpCoords[3 * i + j]
								- pointCoords[3 * i + j];
						v *= v;
						val += v;
					}
					//weight is 1/r^2; val must be <= r^2
					rms += val;
					if (weights[i] * val > 1.0)
					{
						isgood = false;
					}
				}
				if (isgood)
				{
					tmpresult->val = sqrt(rms / n);

					if (tmpresult->val <= qparams.maxRMSD)
					{
						if (!excluder.isDefined() || !isExcluded(tmpresult))
						{
							resultQ.push(alloc.newCorResult(*tmpresult));
							thisConfCnt++;
							if (thisConfCnt >= qparams.orientationsPerConf)
							{
								return false; //terminate early
							}
						}
					}
				}
			}
			else if (false)
			{
				cout << "BADLENS " << tmpresult->location;
				;
				for (unsigned i = 0; i < 3; i++)
				{
					double d = 0;
					for (unsigned j = 0; j < 3; j++)
					{
						double df = molCoords[3 * i + j]
								- molCoords[3 * ((i + 1) % 3) + j];
						df *= df;
						d += df;
					}
					cout << " " << sqrt(d);
				}
				cout << "\n";
			}

			return true;
		}
		else
		{
			const int128_t one = 1;
			const vector<QueryTriplet>& trips = triplets[pos];
			//for each 3 points that match triplet pos
			for (unsigned i = 0, n = tm->matches[pos].size(); i < n; i++)
			{
				const TripletMatchInfo& info = tm->matches[pos][i];
				const Triplet& trip = trips[info.whichTripOrder];
				//figure out what points of this triple are new
				int newqpoints[3] =
				{ -1, -1, -1 };
				int newmpoints[3] =
				{ -1, -1, -1 };
				bool valid = true;
				//points must be unmatched or identically matched
				for (unsigned p = 0; p < 3; p++)
				{
					unsigned qpoint = trip.getPoints()[p].index;
					unsigned mpoint = info.indices[p];

					int curmatch = tmpresult->cor.matchIndex(qpoint);
					if (curmatch < 0) //query point unmatched
					{
						if (alreadyMatched & (one << mpoint)) //already assigned to a different query point
						{
							valid = false;
							break;
						}
						newqpoints[p] = qpoint;
						newmpoints[p] = mpoint;
					}
					else if ((unsigned) curmatch != mpoint) //not identical
					{
						valid = false;
						break;
					}
				}

				if (valid) //recursively descend
				{
					//set cor
					int128_t newmatches = 0;
					for (unsigned p = 0; p < 3; p++)
					{
						if (newqpoints[p] >= 0)
						{
							tmpresult->cor.setIndex(newqpoints[p],
									newmpoints[p]);
							molCoords.push_back(
									ThreePointData::unreduceFloat(
											info.coords[p].x));
							molCoords.push_back(
									ThreePointData::unreduceFloat(
											info.coords[p].y));
							molCoords.push_back(
									ThreePointData::unreduceFloat(
											info.coords[p].z));

							pointCoords.push_back(trip.getPoints()[p].point->x);
							pointCoords.push_back(trip.getPoints()[p].point->y);
							pointCoords.push_back(trip.getPoints()[p].point->z);

							weights.push_back(
									trip.getPoints()[p].point->radiusWeight());
							newmatches |= (one << newmpoints[p]);
						}
					}
					if (!generate(pos - 1, alreadyMatched | newmatches))
					{
						return false;
					}

					//undo chagnes
					for (unsigned p = 0; p < 3; p++)
					{
						if (newqpoints[p] >= 0)
						{
							tmpresult->cor.setIndex(newqpoints[p], -1);

							molCoords.resize(molCoords.size() - 3);
							pointCoords.resize(pointCoords.size() - 3);
							weights.pop_back();
						}
					}
				}

			}
		}
		return true;
	}

public:
	Corresponder(shared_ptr<PharmerDatabaseSearcher>& dptr, unsigned dbid_,
			unsigned ndbids, const vector<PharmaPoint>& pts,
			const vector<vector<QueryTriplet> >& trips, TripletMatches& m,
			CorAllocator& ca, unsigned t, MTQueue<CorrespondenceResult*> & Q,
			const QueryParameters& qp, const Excluder& ex, bool& stop) :
			dbptr(dptr), dbid(dbid_), numdbids(ndbids), points(pts), triplets(
					trips), inQ(m), alloc(ca), threadQ(t), resultQ(Q), qparams(
					qp), excluder(ex), tmpresult(NULL), tm(NULL), thisConfCnt(
					0), processedCnt(0), matchedCnt(0), stopEarly(stop)
	{
		tmpresult = alloc.newCorResult();
		pointCoords.reserve(points.size() * 3);
		molCoords.reserve(points.size() * 3);
		weights.reserve(points.size());
	}

	virtual ~Corresponder()
	{
	}

	//for every t'th match compute all possible, legal correspondences
	void operator()()
	{
		tm = NULL;

		while (inQ.pop(tm, threadQ))
		{
			if (stopEarly)
				break;
			//now recursively greedily enumerate correspondences
			//require a one-to-one correspondence
			tmpresult->reinitialize(*tm, dbid, numdbids);
			thisConfCnt = 0;
			processedCnt++;
			if (!generate(triplets.size() - 1, 0))
			{
				//early termination, clear bookkeepping arrays
				pointCoords.clear();
				molCoords.clear();
				weights.clear();
			}
			if (thisConfCnt > 0)
				matchedCnt++;

		}
		resultQ.removeProducer();

		if (!Quiet)
			cout << "Correspondence ProccessedCnt " << processedCnt
					<< " MatchedCnt " << matchedCnt << " ("
					<< 100 * (double) matchedCnt / (double) processedCnt
					<< "%)\n";
	}
};

#endif /* CORRESPONDER_H_ */
