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
 * Triplet.h
 *
 *  Created on: Aug 5, 2010
 *      Author: dkoes
 *
 *   Three points make a triplet.  Has a canonical form.
 */

#ifndef TRIPLET_H_
#define TRIPLET_H_
#include "pharmarec.h"
#include "ThreePointData.h"
#include <boost/array.hpp>
#include "BoundingBox.h"
#include "RMSD.h"
#include "QueryTripletFingerprint.h"
#include "SimpleFingers.h"
#include "basis.h"
#include "SphereGrid.h"
using namespace boost;


//convenience class for dealing with indexed pharma points
struct PharmaIndex
{
	const PharmaPoint* point;
	unsigned index;
	unsigned pharmaIndex;

	PharmaIndex()
	{
	}
	PharmaIndex(const PharmaPoint& p, unsigned i) :
		point(&p), index(i), pharmaIndex(p.pharma->index)
	{

	}

	bool operator<(const PharmaIndex& rhs) const
	{
		return pharmaIndex < rhs.pharmaIndex;
	}
};

struct TripletRange
{
	unsigned short min;
	unsigned short max;
	unsigned short length;
};

class Triplet
{
protected:
	array<PharmaIndex, 3> PIs;
	array<TripletRange, 3> range;
	unsigned nextUnconnectIndex; //vertex not connect to very next triplet
	unsigned prevUnconnectIndex; //vertex not conneted to previous triplet

	enum DupKind
	{
		AllDifferent, FirstTwoSame, LastTwoSame, AllSame
	} kind;

	void computeLengths()
	{
		range[0].length = ThreePointData::reduceFloatUnsigned(
				PharmaPoint::pharmaDist(*PIs[0].point, *PIs[1].point));
		range[1].length = ThreePointData::reduceFloatUnsigned(
				PharmaPoint::pharmaDist(*PIs[1].point, *PIs[2].point));
		range[2].length = ThreePointData::reduceFloatUnsigned(
				PharmaPoint::pharmaDist(*PIs[2].point, *PIs[0].point));

		//mins
		range[0].min = 0;
		range[1].min = 0;
		range[2].min = 0;

		unsigned short dist = ThreePointData::reduceFloatUnsigned(
				PIs[0].point->radius + PIs[1].point->radius);
		if (dist < range[0].length)
			range[0].min = range[0].length - dist;

		dist = ThreePointData::reduceFloatUnsigned(PIs[1].point->radius
				+ PIs[2].point->radius);
		if (dist < range[1].length)
			range[1].min = range[1].length - dist;

		dist = ThreePointData::reduceFloatUnsigned(PIs[2].point->radius
				+ PIs[0].point->radius);
		if (dist < range[2].length)
			range[2].min = range[2].length - dist;

		//maxs
		range[0].max = range[0].length + ThreePointData::reduceFloatUnsigned(
				PIs[0].point->radius + PIs[1].point->radius);
		range[1].max = range[1].length + ThreePointData::reduceFloatUnsigned(
				PIs[1].point->radius + PIs[2].point->radius);
		range[2].max = range[2].length + ThreePointData::reduceFloatUnsigned(
				PIs[2].point->radius + PIs[0].point->radius);

	}

	//swap points i and j in p and update ranges appropriately
	void swap(unsigned i, unsigned j)
	{
		::swap(PIs[i], PIs[j]);

		if (nextUnconnectIndex == i)
			nextUnconnectIndex = j;
		else if (nextUnconnectIndex == j)
			nextUnconnectIndex = i;

		if (prevUnconnectIndex == i)
			prevUnconnectIndex = j;
		else if (prevUnconnectIndex == j)
			prevUnconnectIndex = i;

		switch (i)
		{
		case 0:
			if (j == 1)
			{
				::swap(range[1], range[2]);
			}
			else if (j == 2)
			{
				::swap(range[0], range[1]);
			}
			break;
		case 1:
			if (j == 0)
			{
				::swap(range[1], range[2]);
			}
			else if (j == 2)
			{
				::swap(range[0], range[2]);
			}
			break;
		case 2:
			if (j == 0)
			{
				::swap(range[0], range[1]);
			}
			else if (j == 1)
			{
				::swap(range[0], range[2]);
			}
			break;
		}
	}

	//put in canonical order
	void canonize()
	{
		//4 cases
		//all different kinds
		if (PIs[0].pharmaIndex != PIs[1].pharmaIndex && PIs[0].pharmaIndex
				!= PIs[2].pharmaIndex && PIs[1].pharmaIndex
				!= PIs[2].pharmaIndex)
		{
			//use the current order, it's fully constrained
			kind = AllDifferent;
		}
		//all the same kind
		else if (PIs[0].pharmaIndex == PIs[1].pharmaIndex && PIs[1].pharmaIndex
				== PIs[2].pharmaIndex)
		{
			kind = AllSame;
			//sort completely on length
			//do manual "sort"
			//the first point is between the smallest and largest length
			if (range[0].length <= range[1].length && range[0].length
					<= range[2].length)
			{
				if (range[1].length <= range[2].length)
				{
					//perfect, no change necessary
				}
				else
				{
					//need to swap first and second point
					swap(0, 1);
				}
			}
			else if (range[1].length <= range[0].length && range[1].length
					<= range[2].length)
			{
				if (range[0].length <= range[2].length)
				{
					//swap first and third
					swap(0, 2);
				}
				else
				{
					swap(0, 1);
					swap(1, 2);
				}
			}
			else //l2 is shortest
			{
				if (range[0].length <= range[1].length)
				{
					swap(0, 2);
					swap(1, 2);
				}
				else
				{
					swap(1, 2);
				}
			}

		}
		//first two are the same kind
		else if (PIs[0].pharmaIndex == PIs[1].pharmaIndex)
		{
			kind = FirstTwoSame;
			//may need to swap
			if (range[1].length > range[2].length)
			{
				swap(0, 1);
			}

		}
		//last two the same
		else if (PIs[1].pharmaIndex == PIs[2].pharmaIndex)
		{
			kind = LastTwoSame;
			//may need to swap
			if (range[0].length > range[2].length)
			{
				swap(1, 2);
			}
		}
		else
			abort(); //huh?
	}

	Triplet(const array<PharmaIndex, 3>& pis, const array<TripletRange, 3>& r,
			unsigned nc, unsigned pc) :
		PIs(pis), range(r), nextUnconnectIndex(nc), prevUnconnectIndex(pc)
	{
	}
public:

	Triplet()
	{
	}

	Triplet(const vector<PharmaPoint>& points, unsigned i, unsigned j,
			unsigned k) :
		nextUnconnectIndex(0), prevUnconnectIndex(0)
	{
		//points are always ordered by pharma kind first
		PIs[0] = PharmaIndex(points[i], i);
		PIs[1] = PharmaIndex(points[j], j);
		PIs[2] = PharmaIndex(points[k], k);
		sort(PIs.begin(), PIs.end());

		computeLengths();
		canonize();
	}
	virtual ~Triplet()
	{
	}

	//index k is not is not connected to the next triplet
	//and the next's triplet's unconnected point is kpoint
	void setNextUnconnected(unsigned k)
	{
		nextUnconnectIndex = k;

	}

	//volume
	double getV() const
	{
		double ret = 1.0;
		for (unsigned i = 0; i < 3; i++)
		{
			ret *= ThreePointData::unreduceFloat(range[i].max - range[i].min);
		}
		return ret;
	}

	//return which point of this triplet ordering is not connected to the next triplet
	unsigned getNextUnconnected() const
	{
		return nextUnconnectIndex;
	}

	unsigned getPrevUnconnected() const
	{
		return prevUnconnectIndex;
	}

	const array<PharmaIndex, 3>& getPoints() const
	{
		return PIs;
	}
	const array<TripletRange, 3>& getRanges() const
	{
		return range;
	}

	unsigned getPharma(unsigned p) const
	{
		return PIs[p].point->pharma->index;
	}

	void dump() const
	{
		for (unsigned i = 0; i < 3; i++)
		{
			cout << PIs[i].index << " " << PIs[i].pharmaIndex << " "
					<< range[i].length << " " << range[i].min << " "
					<< range[i].max << " | ";
		}
		cout << "N" << nextUnconnectIndex << " P" << prevUnconnectIndex << "\n";
	}

};

static unsigned sharedR[3][3] =
{
{ 0, 1, 0 },
{ 1, 0, 2 },
{ 0, 2, 0 } }; //radii index shared be lengths i/j
static unsigned unsharedR[3][3] =
{
{ 0, 0, 1 },
{ 2, 0, 1 },
{ 2, 0, 0 } }; //radii index of i not shared with j
static unsigned sharedA[3][3] =
{
{ 0, 2, 1 },
{ 2, 0, 0 },
{ 1, 0, 0 } }; //index of shared angle between lengths i/j
class QueryTriplet: public Triplet
{
protected:
	//used by search
	BoundingBox mybox;
	array<double, 9> ref;
	double a, b, c; //lengths, shortest to longest
	double aang, bang, cang; //opposing angle
	double r0, r1, r2;
	double smallAthresholdsmallB; // if a <=, b a bounded
	double smallAthresholdlargeB;
	double smallAthresholdlargeC;
	double largeCthresholdsmallA;

	double lengths[3];
	double angles[3];
	double radii[3];
	double tangd[3];
	double tangdM[3];
	double smallLowerThresholds[3][3]; //if first index is smaller than threshold then second has lower bound
	double smallUpperThresholds[3][3]; //if first index is smaller than threshold then second has upper bound
	double largeLowerThresholds[3][3]; //if first index is larger than threshold then second has lower bound
	double largeUpperThresholds[3][3]; //if first index is larger than threshold then second has upper bound

	QueryTripletFingerprint fingerprint;

	//min/max distances to previous triplet's unconnects to successor triplet points
	//indexed by number of hops back (so prev triplet is at index 0, it's prev is at 1)
	vector<unsigned> minkdistsq;
	vector<unsigned> maxkdistsq;

	bool skipfingers;

	//optional constraints, if zero ignore
	unsigned minSize[3];
	unsigned maxSize[3];
	unsigned long vectorMask[3];
	bool hasExtra;
	//set parameters for any extra constraints
	void setExtra(const PharmaPoint* p, unsigned& min, unsigned& max, unsigned long& mask)
	{
		if(p->maxSize > 0) //has a range
		{
			max = p->maxSize;
			min = p->minSize;
			hasExtra = true;
		}

		if(p->vecs.size() > 0)
		{
			//just take first!
			CoordinateBasis basis(*PIs[0].point, *PIs[1].point, *PIs[2].point, true);
			if(basis.hasValidBasis())
			{
				hasExtra = true;
				vector3 v = p->vecs[0];
				float x = 0;
				float y = 0;
				float z = 0;
				basis.setTranslate(vector3(0,0,0));
				basis.replot(v.x(), v.y(), v.z(), x, y, z);

				unsigned g = sphereGrid.pointToGrid(x ,y, z);

				//very approximate and poorly calibrated to query tolerances masks
				if(p->pharma->name == "Aromatic")
				{
					//can be in either direction, so use 45 degrees
					mask = sphereGrid.searchMask(g, 45);
					unsigned g2 = sphereGrid.pointToGrid(-x,-y,-z);
					mask |= sphereGrid.searchMask(g2, 45);
				}
				else
				{
					mask = sphereGrid.searchMask(g, 90);
				}
				//0 is reserved
				mask <<= 1;
			}
		}
	}

	//return true if extra isn't appropriate for this query
	bool extraFails(unsigned min, unsigned max, unsigned long mask, unsigned extra) const
	{
		if(max > 0)
		{
			if(extra >= min && extra <= max)
				return false;
			return true;
		}
		else if(mask != 0 && extra != 0)
		{
			if((1<<extra)&mask)
				return false;
			return true;
		}
		return false;
	}

	void setExtras()
	{
		hasExtra = false;
		memset(minSize, 0, sizeof(minSize));
		memset(maxSize, 0, sizeof(maxSize));
		memset(vectorMask, 0, sizeof(vectorMask));

		for(unsigned i = 0; i < 3; i++)
		{
			setExtra(PIs[i].point, minSize[i], maxSize[i], vectorMask[i]);
		}
	}

	//the distance between a center c and the closest tangent plane of the other two spheres
	double computeTangDMin(double ra, double rb, double ab, double theta, double ac)
	{
		//find distance between circle and plane (from maple);

		double t1, t3, t4, t8, t10, t11, t12, t15, t16, t19, t22, t25, t28;
		double x3 = cos(theta)*ac;
		double y3 = sin(theta)*ac;

		t1 = ra - rb;
		t3 = 1.0 / ab;
		t4 = t3 * t1 * ra;
		t8 = ab + t3 * t1 * rb - t4;
		t10 = (t1 * t1);
		t11 = ab * ab;
		t12 = 1.0 / t11;
		t15 = sqrt((double) (1 - t12 * t10));
		t16 = t15 * ra;
		t19 = t15 * rb - t16;
		t22 = (double) t12 * (t8 * (x3 - t4) + t19 * (y3 - t16));
		t25 = x3 - t4 - t8 * t22;
		t25 *= t25;
		t28 = y3 - t16 - t19 * t22;
		t28 *= t28;
		return sqrt(t25 + t28);
	}

	double computeTangDMax(double ra, double rb, double ab, double theta, double ac)
	{
		double t1, t3, t4, t8, t10, t11, t12, t15, t16, t19, t22, t25, t28;
		double x3 = cos(theta)*ac;
		double y3 = sin(theta)*ac;

		t1 = ra - rb;
		t3 = 1.0 / ab;
		t4 = t3 * t1 * ra;
		t8 = ab + t3 * t1 * rb - t4;
		t10 = (t1 * t1);
		t11 = ab * ab;
		t12 = 1.0 / t11;
		t15 = sqrt((double) (1 - t12 * t10));
		t16 = t15 * ra;
		t19 = -t15 * rb + t16;
		t22 = (t8 * (x3 - t4) + (y3 + t16) * t19) * (double) t12;
		t25 = x3 - t4 - t8 * t22;
		t25 *= t25;
		t28 = y3 + t16 - t19 * t22;
		t28 *= t28;
		return sqrt(t25 + t28);
	}
	//cache data used by search, but not creation
	void computeSearchData()
	{
		setExtras();

		mybox.minx = range[0].min;
		mybox.maxx = range[0].max;
		mybox.miny = range[1].min;
		mybox.maxy = range[1].max;
		mybox.minz = range[2].min;
		mybox.maxz = range[2].max;

		for (unsigned i = 0; i < 3; i++)
		{
			ref[3 * i] = PIs[i].point->x;
			ref[3 * i + 1] = PIs[i].point->y;
			ref[3 * i + 2] = PIs[i].point->z;
		}

		//compute the minimum and maximum perimeters of valid triangles
		//first get the angles (law of cos)
		a = range[0].length;
		b = range[1].length;
		c = range[2].length;

		r0 = ThreePointData::reduceFloatUnsigned(PIs[0].point->radius);
		r1 = ThreePointData::reduceFloatUnsigned(PIs[1].point->radius);
		r2 = ThreePointData::reduceFloatUnsigned(PIs[2].point->radius);

		radii[0] = r0;
		radii[1] = r1;
		radii[2] = r2;

		lengths[0] = a;
		lengths[1] = b;
		lengths[2] = c;

		cang = acos((c * c - a * a - b * b) / (-2 * a * b));
		bang = acos((b * b - a * a - c * c) / (-2 * a * c));
		aang = acos((a * a - b * b - c * c) / (-2 * b * c));

		angles[0] = aang;
		angles[1] = bang;
		angles[2] = cang;


		tangd[2] = computeTangDMin(r0,r1,a,bang,c)-r2;
		tangd[0] = computeTangDMin(r1,r2,b,cang,a)-r0;
		tangd[1] = computeTangDMin(r2,r0,c,aang,b)-r1;

		tangdM[2] = computeTangDMax(r0,r1,a,bang,c)+r2;
		tangdM[0] = computeTangDMax(r1,r2,b,cang,a)+r0;
		tangdM[1] = computeTangDMax(r2,r0,c,aang,b)+r1;

		//if a is short enough, it provides a lower bound on the other distances
		//must be shorter than connecting to other distances shortest point
		smallAthresholdsmallB = floor(sqrt(a * a + r1 * r1 - 2 * a * r1 * cos(
				cang)) - r0);
		smallAthresholdlargeB = floor(sqrt(a * a + r1 * r1 - 2 * a * r1 * cos(
				M_PI - cang)) - r0);
		smallAthresholdlargeC = sqrt(a * a + r0 * r0 - 2 * a * r0 * cos(M_PI
				- bang)) - r1;

		//if c is large enough provides bounds on other distances
		largeCthresholdsmallA = sqrt(c * c + r0 * r0 - 2 * c * r0 * cos(bang))
				+ r2;

		memset(smallLowerThresholds, 0, sizeof(smallLowerThresholds));
		memset(smallUpperThresholds, 0, sizeof(smallUpperThresholds));
		memset(largeLowerThresholds, 0, sizeof(largeLowerThresholds));
		memset(largeUpperThresholds, 0, sizeof(largeUpperThresholds));

		//compute 2D thresholds
		//in theory you can further narrow by considering all three lengths, but that's really complicated
		for (unsigned i = 0; i < 3; i++)
		{
			for (unsigned j = 0; j < 3; j++)
			{
				if (i == j)
					continue;
				//if i's length is less than a certain threshold, j has a lower bound
				double sR = radii[sharedR[i][j]];
				double l = lengths[i];
				double sA = angles[sharedA[i][j]];
				double iR = radii[unsharedR[i][j]];
				smallLowerThresholds[i][j] = floor(sqrt(l * l + sR * sR - 2 * l
						* sR * cos(sA)) - iR);
				//and an upper bound
				smallUpperThresholds[i][j] = floor(sqrt(l * l + sR * sR - 2 * l
						* sR * cos(M_PI - sA)) - iR);

				//if i is larger than threshold, j has lower bound
				largeLowerThresholds[i][j] = ceil(sqrt(l * l + sR * sR - 2 * l
						* sR * cos(sA)) + iR);
				//and an upper bound
				largeUpperThresholds[i][j] = ceil(sqrt(l * l + sR * sR - 2 * l
						* sR * cos(M_PI-sA)) + iR);
			}
		}

	}

	QueryTriplet(const array<PharmaIndex, 3>& pis,
			const array<TripletRange, 3>& r, unsigned nc, unsigned pc) :
		Triplet(pis, r, nc, pc), aang(0), bang(0), cang(0),minkdistsq(0), maxkdistsq(0), skipfingers(false)
	{
		computeSearchData();
	}

public:
	void setFingerPrints(const PharmerQuery& query)
	{
		fingerprint.set(PIs[0].index, PIs[1].index, PIs[2].index, query);
	}
protected:

	//do not use in a ranged context
	bool distancekTooLargeForij(unsigned i, double Li, unsigned j, double Lj, unsigned k, double Lk) const
	{
		//consider the plane tangent to two circles, and the line segment perpindicular to
		//this plane that intersects the third; use lengths i and j to find a maximum bound on k
		double d = tangd[sharedR[i][j]];
		if(d <= 0)
			return false; //no gap

		double maxk = sqrt(Li*Li-d*d)+sqrt(Lj*Lj-d*d);
		if(Lk > maxk)
		{
			return true;
		}

		//additionally consider large i/j
		double rc = radii[sharedR[i][j]];
		if(Li > lengths[i]+rc && Lj > lengths[j]+rc)
		{
			//the mutually constrain each other - don't attempt a perfect calculation
			//though, just consider each individual constraint and add
			double ac = lengths[i];
			double ra = radii[sharedR[i][k]];
			double t1 = ac * ac;
			double t2 = ra * ra;
			double t4 = Li - rc;
			t4 *= t4;
			double t11 = acos((t1 + t2 - t4) / ac / ra / 0.2e1);
			double aa = M_PI - t11;

			double A = angles[sharedA[k][k]];
			if(aa < A)
			{
				double bc = lengths[j];
				double rb = radii[sharedR[j][k]];

				t1 = bc * bc;
				t2 = rb * rb;
				t4 = Lj - rc;
				t4 *= t4;
				t11 = acos((t1 + t2 - t4) / bc / rb / 0.2e1);
				double bb = M_PI - t11;
				double B = angles[sharedA[j][k]];
				if(bb < B)
				{
					double ab = lengths[k];
					double Aang = A-aa;
					double Bang = B-bb;
					t1 = cos(Bang);
					double t3 = cos(Aang);
					double t6 = ab + rb * t1 + ra * t3;
					t6 *= t6;
					double t7 = sin(Bang);
					double t9 = sin(Aang);
					double t12 = -rb * t7 + ra * t9;
					t12*=t12;
					double m = sqrt(t6 + t12);

					if(Lk >= m)
						return true;
				}
			}

		}
		return false;
	}

	//don't use this in a ranged context
	bool distancekTooSmallForij(unsigned i, double Li, unsigned j, double Lj, unsigned k, double Lk) const
	{
		//consider the plane tangent to bottom two circles, and the line segment perpindicular to
		//this plane that intersects the third; use lengths i and j to find a minimum bound on k
		double d = tangdM[sharedR[i][j]];
		if(Li < d || Lj < d)
			return false; //not long enough to be constrained

		double v1 = sqrt(Li*Li-d*d);
		double v2 = sqrt(Lj*Lj-d*d);
		double mink = v1+v2;
		if(Lk < mink)
		{
			//actually two possible traingles, one of which is unconstraining
			double limit = lengths[k] - radii[sharedR[k][j]] - radii[sharedR[k][i]];
			//alternate triangle must be useless
			if(fabs(v1-v2) < limit)
			{
				return true;
			}
		}
		return false;
	}

	//return true if bounds are not met
	bool distancejTooSmallForSmalli(unsigned i, double Li, unsigned j,
			double Lj) const
	{
		if (Li <= smallLowerThresholds[i][j])
		{
			double A = angles[sharedA[i][j]];
			if (isfinite(A))
			{
				double iR = radii[unsharedR[i][j]];
				double jR = radii[unsharedR[j][i]];
				double sR = radii[sharedR[i][j]];
				double l = Li + iR;
				double l2 = Lj + jR;
				double iL = lengths[i];
				double jL = lengths[j];
				double min = -2 * cos(-acos(.5 * (-l * l + iL * iL + sR * sR)
						/ (iL * sR)) + A) * jL * sR + jL * jL + sR * sR;
				if (l2 * l2 < min)
					return true;
			}
		}
		return false;
	}

	//return true if bounds are not met
	bool distancejTooSmallForLargei(unsigned i, double Li, unsigned j,
			double Lj) const
	{
		abort(); //this function is redudant in most cases to distancejTooLargeForSmalli
		if (Li >= largeLowerThresholds[i][j])
		{
			double A = angles[sharedA[i][j]];
			if (isfinite(A))
			{
				double iR = radii[unsharedR[i][j]];
				double jR = radii[unsharedR[j][i]];
				double sR = radii[sharedR[i][j]];
				double l = Li - iR;
				double l2 = Lj + jR;
				double iL = lengths[i];
				double jL = lengths[j];
				double min = -2 * cos(acos(.5 * (-l * l + iL * iL + sR * sR)
						/ (iL * sR)) - A) * jL * sR + jL * jL + sR * sR;
				if (l2 * l2 < min)
					return true;
			}
		}
		return false;
	}

	//return true if j is too large for i
	bool distancejTooLargeForSmalli(unsigned i, double Li, unsigned j,
			double Lj) const
	{
		if (Li <= smallUpperThresholds[i][j])
		{
			double A = angles[sharedA[i][j]];
			if (isfinite(A))
			{
				double iR = radii[unsharedR[i][j]];
				double jR = radii[unsharedR[j][i]];
				double sR = radii[sharedR[i][j]];
				double l = Li + iR;
				double l2 = Lj - jR;
				double iL = lengths[i];
				double jL = lengths[j];

				double max = jL * jL + sR * sR - 2.0 * jL * sR * cos(A + acos(
						(-l * l + iL * iL + sR * sR) / iL / sR / 2.0));
				if (l2 * l2 > max)
					return true;
			}
		}
		return false;
	}

	bool distancejTooLargeForLargei(unsigned i, double Li, unsigned j,
			double Lj) const
	{
		if (Li >= largeUpperThresholds[i][j])
		{
			double A = angles[sharedA[i][j]];
			double iR = radii[unsharedR[i][j]];
			double jR = radii[unsharedR[j][i]];
			double l = Li - iR;
			double l2 = Lj - jR;
			if (l > 0 && l2 > 0 && isfinite(A))
			{
				double sR = radii[sharedR[i][j]];
				double iL = lengths[i];
				double jL = lengths[j];
				double max = -2 * cos(2*M_PI - acos(.5 * (-l * l + iL * iL + sR * sR)
						/ (iL * sR)) - A) * jL * sR + jL * jL + sR * sR;
				if (l2 * l2 > max)
					return true;
			}
		}
		return false;
	}

	//return true if any distances don't work
	//note the relations are not symetric
	//distances must be in order
	bool distancesAreBad(double Li, double Lj, double Lk) const
	{
		double d[3] =
		{ Li, Lj, Lk };

		for (unsigned i = 0; i < 3; i++)
		{
			for (unsigned j = 0; j < 3; j++)
			{
				if (i != j)
				{
					if (distancejTooLargeForSmalli(i, d[i], j, d[j]))
						return true;
					if (distancejTooSmallForSmalli(i, d[i], j, d[j]))
						return true;
					if (distancejTooLargeForLargei(i, d[i], j, d[j]))
						return true;
				}
			}
		}

		if(distancekTooLargeForij(0, Li, 1, Lj, 2, Lk))
			return true;
		if(distancekTooLargeForij(0, Li, 2, Lk, 1, Lj))
			return true;
		if(distancekTooLargeForij(1, Lj, 2, Lk, 0, Li))
			return true;

		if(distancekTooSmallForij(0, Li, 1, Lj, 2, Lk))
			return true;
		if(distancekTooSmallForij(0, Li, 2, Lk, 1, Lj))
			return true;
		if(distancekTooSmallForij(1, Lj, 2, Lk, 0, Li))
			return true;

		return false;
	}

	//return true if l1/l2 are valid based on precomputed thresholds for sides AB
	bool BtooSmallForSmallA(unsigned short L1, unsigned short L2) const
	{
		if (L1 < smallAthresholdsmallB)
		{
			double l = L1 + r0;
			double l2 = L2 + r2;
			double min = -2 * cos(-acos((.5) * (-l * l + a * a + r1 * r1) / (a
					* r1)) + cang) * b * r1 + b * b + r1 * r1;

			if (l2 * l2 < min)
			{
				return true;
			}
		}
		return false;
	}

	bool BtooBigForSmallA(unsigned short L1, unsigned short L2) const
	{
		if (L1 < smallAthresholdlargeB)
		{
			double l = L1 + r0;
			double l2 = L2 - r2;
			double max = b * b + r1 * r1 - 2.0 * b * r1 * cos(cang + acos((-l
					* l + a * a + r1 * r1) / a / r1 / 2.0));

			if (l2 * l2 > max)
			{
				return true;
			}
		}
		return false;
	}

	bool CtooBigForSmallA(unsigned short L1, unsigned short L3) const
	{
		if (L1 < smallAthresholdlargeC)
		{
			float l = L1 + r1;
			float l2 = L3 - r2;
			float max = c * c + r0 * r0 - 2.0 * c * r0 * cos(bang + acos((-l
					* l + a * a + r0 * r0) / a / r0 / 2.0));

			if (l2 * l2 > max)
			{
				return true;
			}
		}
		return false;
	}

public:

	QueryTriplet()
	{
	}

	QueryTriplet(const vector<PharmaPoint>& points, unsigned i, unsigned j,
			unsigned k) :
		Triplet(points, i, j, k), skipfingers(false)
	{
		computeSearchData();
	}
	virtual ~QueryTriplet()
	{
	}

	const QueryTripletFingerprint& getFingerPrint() const { return fingerprint; }

	//return true if order is at least partially canonical (ranges overlap appropriately)
	bool validOverlap() const
	{
		switch (kind)
		{
		case AllDifferent:
			//only one ordering
			return true;
		case AllSame:
			if (range[0].min < range[1].max && range[1].min < range[2].max)
				return true;
			return false;
		case FirstTwoSame:
			if (range[0].min < range[1].max)
				return true;
			return false;

		case LastTwoSame:
			if (range[1].min < range[2].max)
				return true;
			return false;
		}
		return false;
	}

	//return number of useful orderings
	unsigned numOrderings() const
	{
		switch (kind)
		{
		case AllDifferent:
			//only one ordering
			return 1;
		case AllSame:
			return 6;
		case FirstTwoSame:
		case LastTwoSame:
			return 2;
		}
		return 1;
	}

	//generate all useful orderings
	//an ordering if it can be a valid canonical ordering for some values
	//within the range the pharma points can take
	void expand(vector<QueryTriplet>& result, const PharmerQuery& query) const
	{
		QueryTriplet tmp(*this);
		result.clear();
		result.reserve(6);
		//always include canonical version
		tmp.canonize();
		result.push_back(tmp);

		//todo - make this more efficient
		switch (kind)
		{
		case AllDifferent:
			//only one ordering
			break;
		case AllSame:

			tmp.swap(1, 2);//132
			if (tmp.validOverlap())
				result.push_back(tmp);

			tmp.swap(0, 2); //231
			if (tmp.validOverlap())
				result.push_back(tmp);

			tmp.swap(1, 2); //213
			if (tmp.validOverlap())
				result.push_back(tmp);

			tmp.swap(0, 2); //312
			if (tmp.validOverlap())
				result.push_back(tmp);

			tmp.swap(1, 2);//321
			if (tmp.validOverlap())
				result.push_back(tmp);
			break;

		case FirstTwoSame:
			//try swapping
			tmp.swap(0, 1);
			if (tmp.validOverlap())
				result.push_back(tmp);
			break;
		case LastTwoSame:
			tmp.swap(1, 2);
			if (tmp.validOverlap())
				result.push_back(tmp);
			break;
		}

		for(unsigned i = 0, n = result.size(); i < n; i++)
		{
			result[i].setFingerPrints(query);
			result[i].setExtras();
		}
	}

	//return true if there is any possibility of overlap
	//between the triplet 3-space and the box
	bool inRange(const BoundingBox& box) const
	{
		if (mybox.hasOverlap(box))
		{
			//the two most likely culprits
			if (distancejTooSmallForSmalli(0, box.maxx, 2, box.maxz))
				return false;
			if(distancejTooLargeForLargei(2, box.minz, 1, box.miny))
				return false;

			//every other combination
			if (distancejTooSmallForSmalli(0, box.maxx, 1, box.maxy))
				return false;

			if (distancejTooLargeForSmalli(0, box.maxx, 1, box.miny))
				return false;
			if (distancejTooLargeForSmalli(0, box.maxx, 2, box.minz))
				return false;

			if (distancejTooSmallForSmalli(1, box.maxy, 0, box.maxx))
				return false;
			if (distancejTooSmallForSmalli(1, box.maxy, 2, box.maxz))
				return false;

			if (distancejTooLargeForSmalli(1, box.maxy, 0, box.minx))
				return false;
			if (distancejTooLargeForSmalli(1, box.maxy, 2, box.minz))
				return false;

			if (distancejTooSmallForSmalli(2, box.maxz, 0, box.maxx))
				return false;
			if (distancejTooSmallForSmalli(2, box.maxz, 1, box.maxy))
				return false;

			if (distancejTooLargeForSmalli(2, box.maxz, 0, box.minx))
				return false;
			if (distancejTooLargeForSmalli(2, box.maxz, 1, box.miny))
				return false;

			if(distancejTooLargeForLargei(0, box.minx, 1, box.miny))
				return false;
			if(distancejTooLargeForLargei(0, box.minx, 2, box.minz))
				return false;

			if(distancejTooLargeForLargei(1, box.miny, 0, box.minx))
				return false;
			if(distancejTooLargeForLargei(1, box.miny, 2, box.minz))
				return false;

			if(distancejTooLargeForLargei(2, box.minz, 0, box.minx))
				return false;

			return true;
		}
		return false;
	}

	//return true if the bounding box is fully within the 3space
	//defined by the triplet
	bool allInRange(const BoundingBox& box) const
	{
		if(box.containedIn(mybox))
		{
			//additional checks, all coords in the box should be valid
			if(distancejTooSmallForSmalli(0, box.maxx, 1, box.miny))
				return false;
			if(distancejTooSmallForSmalli(0, box.maxx, 2, box.minz))
				return false;

			if(distancejTooSmallForSmalli(1, box.maxy, 0, box.minx))
				return false;
			if(distancejTooSmallForSmalli(1, box.maxy, 2, box.minz))
				return false;

			if(distancejTooSmallForSmalli(2, box.maxz, 0, box.minx))
				return false;
			if(distancejTooSmallForSmalli(2, box.maxz, 1, box.miny))
				return false;

			if(distancejTooLargeForSmalli(0, box.maxx, 1, box.maxy))
				return false;
			if(distancejTooLargeForSmalli(0, box.maxx, 2, box.maxz))
				return false;

			if(distancejTooLargeForSmalli(1, box.maxy, 0, box.maxx))
				return false;
			if(distancejTooLargeForSmalli(1, box.maxy, 2, box.maxz))
				return false;

			if(distancejTooLargeForSmalli(2, box.maxz, 1, box.maxy))
				return false;
			if(distancejTooLargeForSmalli(2, box.maxz, 0, box.maxx))
				return false;

			if(distancejTooLargeForLargei(0, box.minx, 1, box.miny))
				return false;
			if(distancejTooLargeForLargei(0, box.minx, 2, box.minz))
				return false;

			if(distancejTooLargeForLargei(1, box.miny, 0, box.minx))
				return false;
			if(distancejTooLargeForLargei(1, box.miny, 2, box.minz))
				return false;

			if(distancejTooLargeForLargei(2, box.minz, 0, box.minx))
				return false;
			if(distancejTooLargeForLargei(2, box.minz, 1, box.miny))
				return false;

			return true;
		}
		return false;


	}

	//check to see if minimal rmsd fits within query - slow and inaccurate, used for debugging
	bool rmsdMatch(const ThreePointData& tdata) const
	{
		double a[9] = {0,};
		double b[9] = {0,};

		a[0] = ThreePointData::unreduceFloat(tdata.x1()); a[1] = ThreePointData::unreduceFloat(tdata.y1()); a[2] = ThreePointData::unreduceFloat(tdata.z1());
		a[3] = ThreePointData::unreduceFloat(tdata.x2()); a[4] = ThreePointData::unreduceFloat(tdata.y2()); a[5] = ThreePointData::unreduceFloat(tdata.z2());
		a[6] = ThreePointData::unreduceFloat(tdata.x3()); a[7] = ThreePointData::unreduceFloat(tdata.y3()); a[8] = ThreePointData::unreduceFloat(tdata.z3());

		b[0] = PIs[0].point->x; b[1] = PIs[0].point->y; b[2] = PIs[0].point->z;
		b[3] = PIs[1].point->x; b[4] = PIs[1].point->y; b[5] = PIs[1].point->z;
		b[6] = PIs[2].point->x; b[7] = PIs[2].point->y; b[8] = PIs[2].point->z;

		RMSDResult r = calculateRMSD(b, a, 3);

		r.reorient(3, a);

		for(unsigned i = 0; i < 3; i++)
		{
			double d = sqrt(SQR(a[3*i]-b[3*i]) + SQR(a[3*i+1]-b[3*i+1]) + SQR(a[3*i+2]-b[3*i+2]));
			if(d > PIs[i].point->radius)
				return false;
		}

		return true;
	}

	//return true if tdata matches this triplet
	//all checks need to be computational simple
	bool isMatch(const ThreePointData& tdata) const
	{

		if (tdata.l1 >= range[0].min && tdata.l1 <= range[0].max && tdata.l2
				>= range[1].min && tdata.l2 <= range[1].max && tdata.l3
				>= range[2].min && tdata.l3 <= range[2].max)
		{
			if(hasExtra)
			{
				if(extraFails(minSize[0], maxSize[0], vectorMask[0], tdata.extra1))
					return false;
				if(extraFails(minSize[1], maxSize[1], vectorMask[1], tdata.extra2))
					return false;
				if(extraFails(minSize[2], maxSize[2], vectorMask[2], tdata.extra3))
					return false;
			}

			if (distancesAreBad(tdata.l1, tdata.l2, tdata.l3))
				return false;

			if(!skipfingers && !fingerprint.isValid(tdata.fingerprint))
			{
				return false;
			}

	/*		if(!rmsdMatch(tdata))
				return false; */
			return true;
		}
		return false;
	}

	//index k is the point index of the point that is not connected to the previous triplet
	//kpoints are the unconnected points of the previous triplets with the farthest back first
	void setPrevUnconnectedIndex(unsigned k, const vector<PharmaPoint>& kpoints)
	{
		for (unsigned i = 0; i < 3; i++)
		{
			if (PIs[i].index == k)
			{
				prevUnconnectIndex = i;

				for (int j = kpoints.size()-1; j >= 0; j--)
				{
					//store distance ranges between unconnected point of this trip and
					//unconnected point of prev trip
					double dist =
							PharmaPoint::pharmaDist(*PIs[i].point, kpoints[j]);
					double rrange = PIs[i].point->radius + kpoints[j].radius;
					double mind = dist - rrange;
					double maxd = dist + rrange;
					if (mind < 0)
						mind = 0;
					//+/- 1 to provide some slack
					double min = ThreePointData::reduceFloat(mind) - 1;
					if(min < 0) min = 0;
					min *= min;
					minkdistsq.push_back(min);
					double max = ThreePointData::reduceFloat(maxd) + 1;
					max *= max;
					maxkdistsq.push_back(max);
				}
				break;
			}
		}
	}

	bool goodKDistance(unsigned kdist, unsigned distback) const
	{
		if(minkdistsq.size() <= distback) return true;
		return kdist >= minkdistsq[distback] && kdist <= maxkdistsq[distback];
	}

	//set whether or not to skip fingerprint checking; ie for a single triplet query
	void setSkipFingers(bool val)
	{
		skipfingers = val;
	}
};

#endif /* TRIPLET_H_ */
