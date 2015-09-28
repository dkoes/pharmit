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
 * pharminfo.cpp
 *
 *  Created on: Aug 11, 2015
 *      Author: dkoes
 */
#include "pharminfo.h"

using namespace OpenBabel;

struct ReducedFeatureVec
{
	float x;
	float y;
	float z;

	ReducedFeatureVec(): x(0), y(0), z(0) {}
	ReducedFeatureVec(const vector3& v): x(v.x()), y(v.y()), z(v.z()) {}
};


struct ReducedFeature
{
	float x;
	float y;
	float z;
	short vecstart;
	char veclen;
	char size;

	ReducedFeature() :
			x(0), y(0), z(0), vecstart(0), veclen(0), size(0)
	{
	}

	ReducedFeature(const PharmaPoint& pt, unsigned vecst) :
			x(pt.x), y(pt.y), z(pt.z), vecstart(vecst), veclen(pt.vecs.size()), size(
					pt.size)
	{
	}

	bool matches(const PharmaPoint& pt, const ReducedFeatureVec* vecs) const
	{
		//distance short circuiting
		double r = pt.radius;
		double xdiff = fabs(pt.x - x);
		if(xdiff > r) return false;

		double ydiff = fabs(pt.y-y);
		if(ydiff > r) return false;

		double zdiff = fabs(pt.z-z);
		if(zdiff > r) return false;

		double dist2 = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff;
		double r2 = r*r;
		if(dist2 > r2) return false;

		if(pt.minSize != pt.maxSize)
		{
			if((unsigned)size < pt.minSize) return false;
			if(pt.maxSize != UINT_MAX && (unsigned)size > pt.maxSize) return false;
		}

		if(pt.vecs.size() > 0 && veclen > 0)
		{
			//only look at first vector and only if we actually have vectors
			vector3 qvec = pt.vecs[0];
			bool hasvecmatch = false;
			for(unsigned i = vecstart, n = vecstart+veclen; i < n; i++)
			{
				const ReducedFeatureVec& vec = vecs[i];
				vector3 v(vec.x, vec.y, vec.z);
				double val = dot(v, qvec);
				if(val >= 0) //same general direction
				{
					hasvecmatch = true;
					break;
				}
			}
			if(!hasvecmatch) return false;
		}
		return true;
	}
};

//write points to file and return the starting offset
// first write offsets into different kinds of features
//then write the features themselves, with offsets for vectors
//then write the vectors
unsigned long writePharmacophoreInfo(FILE *f, const vector<PharmaPoint>& points, const Pharmas& pharmas)
{
	unsigned long ret = ftell(f);
	//first the offsets indexed by pharma kind
	//offsets[i] is where pharma kind i-1 stops (ie, a kind >= i is there)
	vector<short> offsets(pharmas.size(),0);
	unsigned currindex = 0;
	vector<ReducedFeature> features; features.reserve(points.size());
	vector<ReducedFeatureVec> vectors;
	for(unsigned i = 0, n = points.size(); i < n; i++)
	{
		const PharmaPoint& pt = points[i];
		unsigned indx = pt.pharma->index;
		if(indx < currindex) abort(); //assume sorted
		while(indx > currindex)
		{
			offsets[currindex] = i;
			currindex++;
		}

		//fill out feature and vector arrays
		features.push_back(ReducedFeature(pt, vectors.size()));
		for(unsigned v = 0, nv = pt.vecs.size(); v < nv; v++)
		{
			vectors.push_back(ReducedFeatureVec(pt.vecs[v]));
		}
	}
	while(currindex < offsets.size()) {
		offsets[currindex] = points.size();
		currindex++;
	}

	//write out data
	fwrite(&offsets[0], sizeof(short), offsets.size(), f);
	fwrite(&features[0], sizeof(ReducedFeature), features.size(), f);
	fwrite(&vectors[0], sizeof(ReducedFeatureVec), vectors.size(), f);

	return ret;
}

//takes a readonly pointer to a pharmacophore written in the above format and checks for a match to the query
bool pharmacophoreMatchesQuery(const char *pharmacophore, const vector<PharmaPoint>& querypoints, const Pharmas& pharmas)
{
	short *offsets = (short*)pharmacophore;
	unsigned n = pharmas.size();
	unsigned numfeatures = offsets[n-1];
	const ReducedFeature *features = (ReducedFeature*)(offsets+n);
	const ReducedFeatureVec *vecs = (ReducedFeatureVec*)(features+numfeatures);

	//every querypoint must match some feature (same one can be matched by multiple querypoints if they overlap)
	for(unsigned i = 0, nq = querypoints.size(); i < nq; i++)
	{
		const PharmaPoint& qpt = querypoints[i];
		unsigned phindx = qpt.pharma->index;
		unsigned start = 0; //index of first feature of same type
		if(phindx > 0) start = offsets[phindx-1];
		unsigned end = offsets[phindx]; //end of range
		bool hasmatch = false;
		for(unsigned f = start; f < end; f++)
		{
			const ReducedFeature& feature = features[f];
			if(feature.matches(qpt, vecs))
			{
				hasmatch = true;
				break;
			}
		}

		if(!hasmatch) //this query point did not overlap anything
			return false;
	}
	return true;
}

